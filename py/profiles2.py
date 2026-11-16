"""Simple robust vertical design: smoothed terrain + grade-limited envelopes over clearance constraints."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
from profiles import *
from scipy.ndimage import gaussian_filter1d
PRIO=['b0','b6','b10','b8','b9','b23','b1','b25','b3','b7','b5','b2']
def gs(a,sig): return gaussian_filter1d(a,sig/DS,mode='nearest')
def envelope(m,g):
    e=m.copy(); st=g*DS
    for i in range(1,len(e)): e[i]=max(e[i],e[i-1]-st)
    for i in range(len(e)-2,-1,-1): e[i]=max(e[i],e[i+1]-st)
    return e
def relations(chains,plan,cls,verbose=False):
    decks=plan['decks']; ids=list(chains.keys())
    lines={c:sg.LineString(cls[c]['P']) for c in ids}; tang={}
    for c in ids:
        P=cls[c]['P']; t=np.gradient(P,axis=0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9; tang[c]=t
    trees={c:cKDTree(cls[c]['P']) for c in ids}
    over=[]; pins=[]; par=[]
    for a_ in range(len(ids)):
        for b_ in range(a_+1,len(ids)):
            A,Bc=ids[a_],ids[b_]
            if lines[A].distance(lines[Bc])>14: continue
            PA,PB=cls[A]['P'],cls[Bc]['P']; sa,sb=cls[A]['st'],cls[Bc]['st']
            X=lines[A].intersection(lines[Bc])
            pts=[g for g in (X.geoms if hasattr(X,'geoms') else [X]) if g.geom_type=='Point']
            stA=[]
            for p in pts:
                i=min(int(round(lines[A].project(p)/DS)),len(PA)-1); j=min(int(round(lines[Bc].project(p)/DS)),len(PB)-1)
                if abs(tang[A][i]@tang[Bc][j])>0.985: continue
                ka,kb=sa[i],sb[j]
                if ka==kb: pins.append((A,i,Bc,j)); continue
                if ka>=0 and kb<0: up,lo=(A,i),(Bc,j)
                elif kb>=0 and ka<0: up,lo=(Bc,j),(A,i)
                else:
                    la=1 if decks[ka].get('iso') else 0; lb=1 if decks[kb].get('iso') else 0
                    if la==lb:
                        if verbose: print('ambiguous',A,Bc,(round(p.x),round(p.y)),ka,kb)
                        continue
                    up,lo=((A,i),(Bc,j)) if la>lb else ((Bc,j),(A,i))
                over.append((up[0],up[1],lo[0],lo[1])); stA.append(i)
                if verbose: print('crossing',up[0],'over',lo[0],(round(p.x),round(p.y)))
            dd,jj=trees[Bc].query(PA)
            for i in np.where(dd<14)[0]:
                j=jj[i]
                if sa[i]!=sb[j]: continue
                if any(abs(i-q)*DS<150 for q in stA): continue
                if abs(tang[A][i]@tang[Bc][j])>0.94: par.append((A,i,Bc,j,dd[i]))
    return over,pins,par
def design2(chains,plan,cls,clear=6.8,grade=0.045,passes=4,verbose=False):
    over,pins,par=relations(chains,plan,cls,verbose)
    ids=list(chains.keys()); pr={c:(PRIO.index(c) if c in PRIO else 50+k) for k,c in enumerate(ids)}
    Gex=plan.get('Gexist'); Gb=Gex.buffer(1.0) if Gex is not None else None
    base={}; need={}; ktrees={c:cKDTree(cls[c]['P']) for c in ids}
    for c in ids:
        d=cls[c]; P=d['P']; z0=terr(P[:,0],P[:,1]); d['zt']=z0; m=len(z0); ch=chains[c]
        sig=np.full(m,ch.get('sig',70.0))
        ex=shapely.contains(Gb,shapely.points(P)) if Gb is not None else np.zeros(m,bool)
        b_new=gs(z0,ch.get('sig',70.0)); b_ex=gs(z0,ch.get('sig_ex',22.0))
        wex=gs(ex.astype(float),40.0)
        base[c]=b_ex*wex+b_new*(1-wex)
        dk=d['st']>=0; idx=np.arange(m); nd=np.full(m,1e9)
        ends=[i for i in range(m) if dk[i] and (i==0 or not dk[i-1] or i==m-1 or not dk[i+1])]
        for e in ends: nd=np.minimum(nd,np.abs(idx-e)*DS)
        zs=gs(z0,40.0)
        need[c]=np.where(dk,np.minimum(z0,zs+1.0)+np.clip(0.06*(nd-10.0),0,ch.get('deckclear',5.5)),-1e9)
        if ch.get('deckclear',0)>20:   # existing bridge: straight chord between abutments
            k=np.where(dk)[0]; a,b=k[0],k[-1]; need[c][dk]=np.interp(k,[a,b],[base[c][a],base[c][b]])
        d['z']=base[c].copy()
    # where a bridge crosses one of today's roads in the 2018-design area, let the lower road dip under it
    for (u,i,l,j) in over:
        if not (l.startswith('x') or l.startswith('e_south')): continue
        zu=max(base[u][i],need[u][i]); lack=base[l][j]+clear-zu
        if lack>0:
            m=len(base[l]); t=np.clip(1-np.abs(np.arange(m)-j)*DS/110.0,0,1); base[l]=base[l]-min(lack,5.0)*t*t*(3-2*t)
    for it in range(passes):
        for c in sorted(ids,key=lambda c:pr[c]):
            d=cls[c]; m=len(d['z']); mm=need[c].copy()
            for (u,i,l,j) in over:
                if u!=c: continue
                for di in range(-4,5):
                    q=i+di
                    if 0<=q<m: mm[q]=max(mm[q],cls[l]['z'][j]+clear)
            z=np.maximum(base[c],envelope(mm,grade))
            z=gs(z,35.0); z=np.maximum(z,mm); z=gs(z,12.0)
            # pins at same-level crossings and parallel ties toward higher-priority chains
            delta=np.zeros(m); wt=np.zeros(m)
            for (A,i,B,j,dist) in par:
                if A==c and pr[B]<pr[A]: k,o,oj=i,B,j
                elif B==c and pr[A]<pr[B]: k,o,oj=j,A,i
                else: continue
                dz_=cls[o]['z'][oj]-z[k]
                w_=min(1.0,(14-dist)/7.0)*float(np.clip(1-(abs(dz_)-3)/3,0,1))
                if w_>wt[k]: wt[k]=w_; delta[k]=dz_
            if wt.any():
                # fill delta between tied stations smoothly, decay outside
                k=np.where(wt>0)[0]; dfull=np.interp(np.arange(m),k,delta[k]); 
                dist=np.abs(np.arange(m)[:,None]-k[None,:]).min(1)*DS
                wfull=np.interp(np.arange(m),k,wt[k])*np.clip(1-dist/120.0,0,1)
                z=z+gs(dfull*wfull,15.0)
            for (A,i,B,j) in pins:
                if A==c and pr[B]<pr[A]: k,o,oj=i,B,j
                elif B==c and pr[A]<pr[B]: k,o,oj=j,A,i
                else: continue
                dz=cls[o]['z'][oj]-z[k]; z=z+dz*np.clip(1-np.abs(np.arange(m)-k)*DS/90.0,0,1)**2*(3-2*np.clip(1-np.abs(np.arange(m)-k)*DS/90.0,0,1))
            # ends that merge into / diverge from a higher-priority road take its elevation
            for e in (0,m-1):
                p=d['P'][e]; te=d['P'][min(e+1,m-1)]-d['P'][max(e-1,0)]; te=te/(np.hypot(*te)+1e-9); best=None
                for o in ids:
                    if pr[o]>=pr[c]: continue
                    dd_,j=ktrees[o].query(p)
                    if dd_>26: continue
                    Po=cls[o]['P']; to=Po[min(j+1,len(Po)-1)]-Po[max(j-1,0)]; to=to/(np.hypot(*to)+1e-9)
                    if abs(te@to)<0.9 and dd_>16: continue
                    if (cls[o]['st'][j]>=0)!=(d['st'][e]>=0): continue
                    if best is None or dd_<best[0]: best=(dd_,cls[o]['z'][j])
                if best is not None:
                    t=np.clip(1-np.abs(np.arange(m)-e)*DS/min(260.0,0.6*m*DS),0,1); z=z+(best[1]-z[e])*t*t*(3-2*t)
            if chains[c].get('linear'):
                line=np.linspace(z[0],z[-1],m); z=line+(1-chains[c]['linear'])*(z-line)
            for (u,i,l,j) in over:   # hard clearance repair
                if u!=c: continue
                lack=cls[l]['z'][j]+clear-z[max(0,i-3):i+4].min()
                if lack>0: 
                    t=np.clip(1-np.abs(np.arange(m)-i)*DS/(40+lack/grade),0,1); z=z+lack*t*t*(3-2*t)
            d['z']=gs(z,12.0)
    return cls
HWD={'b0':22.0,'b6':14.0,'b10':12.0,'b23':12.0,'b8':9.0,'b9':9.0,'b1':8.0,'b3':7.0,'b7':7.0,'b2':6.0,'b5':6.0}
def hw(c,chains=None):
    if c in HWD: return HWD[c]
    if c.startswith('e_'): return 9.0
    if c.startswith('x'): return 8.0 if chains and chains[c].get('cls')=='motorway' else 5.0
    if c.startswith('y'): return 6.0
    return 4.5
def patch_plan(chains,plan,cls):
    """fill holes in pavement along alignments (under crossing bridges, at plan seams, between bridge pieces)"""
    import shapely.ops as so
    G=plan['G']; decks=plan['decks']; DU=shapely.unary_union([d['poly'] for d in decks]); Gb=G.buffer(0.6); shapely.prepare(Gb)
    addG=[]; nd=0
    for c,d in cls.items():
        xr=c.startswith('x') and chains[c].get('ride')
        if (c.startswith('x') and not xr) or c.startswith('y'): continue
        P=d['P']; st=d['st']; pts=shapely.points(P); h=hw(c,chains); ls=sg.LineString(P)
        inG=shapely.contains(Gb,pts); inD=shapely.contains(DU.buffer(0.6),pts)
        m=(st<0)&~inG
        if xr: m&=shapely.distance(shapely.unary_union([decks[k]['poly'] for k in chains[c]['ride']]),pts)<45
        idx=np.where(m)[0]
        for run in (np.split(idx,np.where(np.diff(idx)>1)[0]+1) if len(idx) else []):
            if len(run)<2: continue
            a,b=d['s'][run[0]]-7,d['s'][run[-1]]+7; seg=so.substring(ls,max(a,0),min(b,ls.length))
            under=inD[run].mean()>0.5
            if under: addG.append(seg.buffer(h,cap_style=2).intersection(DU.buffer(2.0)))
            elif b-a<110: addG.append(seg.buffer(min(h,8.0),cap_style=2))
        for k in np.unique(st[st>=0]):
            own=decks[k]['poly'].buffer(1.0); mm=(st==k)&~shapely.contains(own,pts); idx=np.where(mm)[0]
            for run in (np.split(idx,np.where(np.diff(idx)>1)[0]+1) if len(idx) else []):
                if len(run)<2: continue
                a,b=d['s'][run[0]]-7,d['s'][run[-1]]+7; seg=so.substring(ls,max(a,0),min(b,ls.length))
                dp=seg.buffer(h,cap_style=2).difference(shapely.unary_union([dd['poly'] for dd in decks if not dd.get('iso')]).buffer(0.05))
                dp=dp.buffer(-0.5).buffer(0.5)
                if dp.is_empty: continue
                if dp.geom_type!='Polygon': dp=max(dp.geoms,key=lambda g:g.area)
                if dp.area<40: continue
                decks.append({'poly':dp,'iso':False,'src':'patch','abut':sg.Polygon(),'width':2*h,'own':c}); st[run]=len(decks)-1; nd+=1
    if addG: G=shapely.unary_union([G]+[g for g in addG if not g.is_empty]).buffer(0)
    plan['G']=G; print('patched: G pieces',len(addG),'deck pieces',nd)
    return plan
def axis_of(p):
    r=p.minimum_rotated_rectangle; c=np.array(r.exterior.coords)[:4]; e=[c[1]-c[0],c[2]-c[1]]
    a=e[0] if np.hypot(*e[0])>=np.hypot(*e[1]) else e[1]; L=np.hypot(*a); return a/L,L,c
def end_zones(p):
    a,L,c=axis_of(p); cen=np.array(p.centroid.coords[0]); n=np.array([-a[1],a[0]])
    return shapely.unary_union([sg.LineString([cen+a*s*L/2-n*40,cen+a*s*L/2+n*40]).buffer(5.0) for s in (-1,1)])
# (bridge centre, point the bridge must reach): Brevard Road also spans the southern collector roadway
BRIDGE_REACH=[((283864,207096),(283867,207036))]
def add_axis_chains(chains,plan,cls):
    """a bridge nobody rides: bend the nearest parallel local road onto its long axis so it crosses on the bridge"""
    riders=set()
    for d in cls.values(): riders|=set(int(k) for k in np.unique(d['st']) if k>=0)
    n=0
    for k,d in enumerate(plan['decks']):
        if k in riders or d.get('src')!='kmz2': continue
        a,L,c=axis_of(d['poly']); cen=np.array(d['poly'].centroid.coords[0]); nrm=np.array([-a[1],a[0]]); best=None
        for cid,ch in chains.items():
            if not cid.startswith('x') or ch.get('cls')=='motorway': continue
            P=cls[cid]['P']; dd=np.hypot(*(P-cen).T); i=dd.argmin()
            if dd[i]>60: continue
            t=P[min(i+2,len(P)-1)]-P[max(i-2,0)]; t=t/max(np.hypot(*t),1e-9)
            if abs(t@a)<0.85: continue
            if best is None or dd[i]<best[0]: best=(dd[i],cid)
        if best is None: continue
        cid=best[1]; P=cls[cid]['P'].copy(); al=(P-cen)@a; lat=(P-cen)@nrm
        u=np.clip((np.abs(al)-(L/2+15))/110.0,0,1); w=1-u*u*(3-2*u)
        P=P-np.outer(lat*w,nrm)
        # lengthen the bridge over any roadway that crosses just beyond its drawn ends
        E=(np.array(d['poly'].exterior.coords)-cen); wd=(E@nrm).max()-(E@nrm).min(); lo,hi=-L/2,L/2
        from scipy.spatial import cKDTree
        mw=np.vstack([cls[c2]['P'] for c2 in cls if c2.startswith('e_') or chains[c2].get('cls')=='motorway'])
        mt=cKDTree(mw); Gp=plan['G']
        for sgn in (-1,1):
            for t in np.arange(1.0,41.0,2.0):
                q=[cen+a*sgn*(L/2+t)+nrm*s2*(wd/2+10) for s2 in (-1,1)]
                if all(Gp.contains(sg.Point(*x)) for x in q) and max(mt.query(q)[0])<26:
                    if sgn>0: hi=L/2+t-1
                    else: lo=-(L/2+t-1)
        for c0,pt in BRIDGE_REACH:      # manual: roadways known to pass under this bridge
            if np.hypot(*(cen-np.array(c0)))<40: t=(np.array(pt)-cen)@a; lo=min(lo,t); hi=max(hi,t)
        if hi>L/2 or lo<-L/2:
            ext=sg.Polygon([cen+a*lo-nrm*wd/2,cen+a*hi-nrm*wd/2,cen+a*hi+nrm*wd/2,cen+a*lo+nrm*wd/2])
            d['poly']=shapely.unary_union([d['poly'],ext]).buffer(0.5).buffer(-0.5); d['abut']=end_zones(d['poly']); print('bridge',k,'lengthened to',round(hi-lo),'m from',round(L))
        al=(P-cen)@a
        if al[0]>al[-1]: P=P[::-1]; al=al[::-1]
        if al[0]>lo-40: P=np.vstack([cen+a*(lo-60),P[al>lo+5]])
        al=(P-cen)@a
        if al[-1]<hi+40: P=np.vstack([P[al<hi-5],cen+a*(hi+60)])
        chains[cid]=dict(chains[cid],xy=P,ride=[k],decks=[k],sig=45.0); n+=1; print('bridge',k,'carried by',cid)
    print('axis chains',n); return chains
def load_all():
    import ext
    plan=pickle.load(open(W+'/sheets_plan.pkl','rb'))
    chains=load_chains_B()
    for c,pts in RIDE.items(): chains[c]['ride']=[deck_at(plan,p) for p in pts]
    chains['b22']['linear']=0.85
    G=plan['G']
    for c in ('b8','b9'):
        s_,P=resample(chains[c]['xy']); z0=terr(P[:,0],P[:,1])
        low=np.where(z0<612)[0]; a=max(0,low[0]-12); b=min(len(s_)-1,low[-1]+12)
        seg=sg.LineString(P[a:b+1]); dp=seg.buffer(12.5,cap_style=2).intersection(G)
        if dp.geom_type!='Polygon': dp=max(dp.geoms,key=lambda g:g.area)
        ab=shapely.unary_union([sg.Point(*P[a]).buffer(16),sg.Point(*P[b]).buffer(16)])
        plan['decks'].append({'poly':dp,'iso':False,'src':'exist','abut':ab,'width':20,'own':c})
        G=G.difference(dp.buffer(0.05)); chains[c]['ride']=[len(plan['decks'])-1]; chains[c]['deckclear']=30.0
    # extension roads (Overture centrelines) with synthesized pavement
    E=ext.build_ext2(); synth=[]
    hw={'south':0,'i240e':7.6,'north':7.6,'pate':9.2,'patw':9.2}
    for k,v in E.items():
        for dname in ('out','in'):
            if dname=='in' and not v['dual'] and hw[k]>8: continue
            chains[f'e_{k}_{dname}']={'xy':v[dname],'sig':45.0,'sig_ex':45.0,'ext':True}
            if hw[k]>0: synth.append(sg.LineString(v[dname]).buffer(hw[k] if (v['dual'] or hw[k]>8) else 12.0,cap_style=2))
    # where the plan shows lanes next to an extension, size its pavement from the plan's lane layout instead (lanes_ext.py)
    import lanes_ext as LE
    LT=LE.build(chains=chains,E=E,hwf=globals()['hw']); pickle.dump(LT,open(W+'/lanes_ext.pkl','wb'))
    if LT:
        keep=[g for g in synth if not any(g.distance(sg.LineString(E[k]['out']))<1.0 for k in LT)]
        synth=keep+[LE.pavement(LT)]
    allp=shapely.unary_union([G]+[d['poly'] for d in plan['decks']])
    S=shapely.unary_union(synth).difference(allp.buffer(3.0)).buffer(3.2).difference(shapely.unary_union([d['poly'] for d in plan['decks']]))
    plan['synth']=shapely.unary_union(synth)
    G=shapely.unary_union([G,S]).buffer(0)
    # bridges drawn in the 2018 design that the deck extraction missed (e.g. Brevard Road over I-26)
    K=pickle.load(open(W+'/kmz_sp.pkl','rb')); DU0=shapely.unary_union([d['poly'] for d in plan['decks']])
    miss=[p.buffer(0) for p in K['S'] if p.centroid.y<208250 and p.area>40 and p.buffer(0).intersection(DU0).area<0.3*p.area]
    MU=shapely.unary_union(miss).buffer(0.3).buffer(-0.3) if miss else sg.Polygon(); nmiss=0
    for g in (MU.geoms if hasattr(MU,'geoms') else [MU]):
        if g.is_empty or g.area<300 or g.geom_type!='Polygon': continue
        g=sg.Polygon(g.exterior)
        plan['decks'].append({'poly':g,'iso':False,'src':'kmz2','abut':end_zones(g),'width':20,'own':None}); nmiss+=1
    print('added missing bridges',nmiss)
    # Section A (2018 design, no alignments): today's road centrelines give every roadway, ramp and cross street a profile
    R=ext.load(('motorway','trunk','primary','secondary','tertiary','residential','unclassified'))
    south=sg.box(270000,200000,300000,208200.0)
    GA=G.intersection(south); dA=[d['poly'] for d in plan['decks'] if d['poly'].centroid.y<208300]
    zmain=GA.buffer(25); zdeck=shapely.unary_union(dA).buffer(60); clipz=shapely.unary_union([GA.buffer(140),zdeck.buffer(350)]).intersection(south)
    rank={'motorway':0,'trunk':1,'primary':2,'secondary':3,'tertiary':4,'residential':5,'unclassified':5}
    R.sort(key=lambda r:rank[r['cls']]); n=0
    for r in R:
        ls=sg.LineString(r['xy'])
        if not (ls.intersects(zdeck) or (rank[r['cls']]<=4 and ls.intersects(zmain))): continue
        c=ls.intersection(clipz)
        for g in (c.geoms if hasattr(c,'geoms') else [c]):
            if g.geom_type=='LineString' and g.length>60: chains[f'x{n}']={'xy':np.array(g.coords),'sig':45.0,'sig_ex':30.0,'cls':r['cls']}; n+=1
    print('x chains',n)
    plan['G']=G
    pickle.dump(plan,open(W+'/sheets_plan2.pkl','wb'))
    return plan,chains
if __name__=='__main__':
    plan,chains=load_all()
    cls=classify(chains,plan,noride={'b1','b10'},straddle=hw)
    chains=add_axis_chains(chains,plan,cls); cls=classify(chains,plan,noride={'b1','b10'},straddle=hw)
    plan=patch_plan(chains,plan,cls); pickle.dump(plan,open(W+'/sheets_plan2.pkl','wb'))
    cls=design2(chains,plan,cls,verbose=True)
    pickle.dump({'chains':chains,'cls':cls},open(W+'/profiles.pkl','wb'))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    show=sys.argv[1].split(',')
    fg,axs=plt.subplots(len(show),1,figsize=(16,2.6*len(show)),squeeze=False)
    for ax,c in zip(axs[:,0],show):
        d=cls[c]; ax.plot(d['s'],d['zt'],'g-',lw=0.8); ax.plot(d['s'],d['z'],'k-',lw=1.2); dk=d['st']>=0; ax.plot(d['s'][dk],d['z'][dk],'r.',ms=3)
        g=np.gradient(d['z'],d['s'])*100; ax.set_title(f'{c}: max grade {np.abs(g).max():.1f}%'); ax.grid(True)
    os.makedirs(W+'/img',exist_ok=True)
    fg.tight_layout(); fg.savefig(W+'/img/prof3.png',dpi=62)
    for c in cls:
        g=np.abs(np.gradient(cls[c]['z'],cls[c]['s'])).max()*100
        if g>7: print(c,'grade %.1f'%g,'len %.0f'%cls[c]['s'][-1])
