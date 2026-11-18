"""Simple robust vertical design: smoothed terrain + grade-limited envelopes over clearance constraints."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
from profiles import *
from scipy.ndimage import gaussian_filter1d
PRIO=['b0','b6','b10','b8','b9','b23','b1','b25','b3','b7','b5','b2','b22']   # b22 last: it sets the level of the ramp it shares with b18 north of Patton Avenue
def gs(a,sig): return gaussian_filter1d(a,sig/DS,mode='nearest')
def envelope(m,g):
    e=m.copy(); st=g*DS
    for i in range(1,len(e)): e[i]=max(e[i],e[i-1]-st)
    for i in range(len(e)-2,-1,-1): e[i]=max(e[i],e[i+1]-st)
    return e
def relations(chains,plan,cls,verbose=False):
    decks=plan['decks']; ids=list(chains.keys()); stack=set(plan.get('stack',()))
    lines={c:sg.LineString(cls[c]['P']) for c in ids}; tang={}
    for c in ids:
        P=cls[c]['P']; t=np.gradient(P,axis=0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9; tang[c]=t
    trees={c:cKDTree(cls[c]['P']) for c in ids}
    over=[]; pins=[]; par=[]
    for a_ in range(len(ids)):
        for b_ in range(a_+1,len(ids)):
            A,Bc=ids[a_],ids[b_]
            # roads count as running side by side out to 14 m, or as far as their carriageways overlap if that is more
            # (a ramp's lane along the edge of a wide mainline is 20 m and more from the mainline's centreline)
            lim=14.0 if A.startswith(('x','e_')) or Bc.startswith(('x','e_')) else max(14.0,hw(A,chains)+hw(Bc,chains))
            if lines[A].distance(lines[Bc])>lim: continue
            PA,PB=cls[A]['P'],cls[Bc]['P']; sa,sb=cls[A]['st'],cls[Bc]['st']
            X=lines[A].intersection(lines[Bc])
            pts=[g for g in (X.geoms if hasattr(X,'geoms') else [X]) if g.geom_type=='Point']
            stA=[]
            for p in pts:
                i=min(int(round(lines[A].project(p)/DS)),len(PA)-1); j=min(int(round(lines[Bc].project(p)/DS)),len(PB)-1)
                if abs(tang[A][i]@tang[Bc][j])>0.985: continue
                ka,kb=sa[i],sb[j]
                if (A,Bc) in stack or (Bc,A) in stack:      # stacked crossing on one merged bridge outline (split_stacked)
                    up,lo=((A,i),(Bc,j)) if (A,Bc) in stack else ((Bc,j),(A,i))
                    over.append((up[0],up[1],lo[0],lo[1])); stA.append(i)
                    if verbose: print('crossing',up[0],'over',lo[0],(round(p.x),round(p.y)),'(stacked)')
                    # the upper bridge must clear the lower deck over the whole patch where the two overlap, not only
                    # where the centrelines meet: a wide lower deck crossed at a skew is under the upper one for a long
                    # way, and its far side sits at the height of the lower road some way up or down its grade
                    U_,iu,L_,jl=up[0],up[1],lo[0],lo[1]
                    ku=set(int(k) for k in cls[U_]['st'][max(0,iu-60):iu+61] if k>=0)
                    kl=set(int(k) for k in cls[L_]['st'][max(0,jl-60):jl+61] if k>=0)-ku
                    if ku and kl:
                        ov=shapely.unary_union([decks[k]['poly'] for k in ku]).intersection(shapely.unary_union([decks[k]['poly'] for k in kl]))
                        if not ov.is_empty:
                            x0,y0,x1,y1=ov.bounds; gx,gy=np.meshgrid(np.arange(x0,x1,3.0),np.arange(y0,y1,3.0)); Q=np.c_[gx.ravel(),gy.ravel()]
                            Q=Q[shapely.contains(ov,shapely.points(Q))]
                            if len(Q):
                                prs=set(zip(trees[U_].query(Q)[1].tolist(),trees[L_].query(Q)[1].tolist()))
                                for (qu,ql) in sorted(prs): over.append((U_,int(qu),L_,int(ql)))
                                if verbose: print('  clearance held over %.0f m2 where the decks overlap (%d station pairs)'%(ov.area,len(prs)))
                    continue
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
                # a skewed crossing is long: the upper road must clear the lower one wherever the lower one is under
                # the bridge and within the upper carriageway, not only where the centrelines meet
                U_,iu,L_,jl=up[0],up[1],lo[0],lo[1]; su=cls[U_]['st']
                if su[iu]>=0 and not U_.startswith(('x','e_')):
                    ks=set(int(k) for k in su[max(0,iu-60):iu+61] if k>=0)
                    Ud=shapely.unary_union([decks[k]['poly'] for k in ks]).buffer(0.5)
                    a0,a1=max(0,jl-80),min(len(cls[L_]['P']),jl+81); Pl=cls[L_]['P'][a0:a1]
                    du,ii=trees[U_].query(Pl); ins=shapely.contains(Ud,shapely.points(Pl))&(du<=hw(U_,chains)+2.0)&(su[ii]>=0)
                    for q in np.where(ins)[0][::2]:
                        if abs(int(ii[q])-iu)>3: over.append((U_,int(ii[q]),L_,a0+int(q)))
            dd,jj=trees[Bc].query(PA)
            lim2=14.0 if stA else lim       # roads that cross at different levels are not one roadway nearby
            for i in np.where(dd<lim2)[0]:
                j=jj[i]
                if sa[i]!=sb[j]: continue
                if any(abs(i-q)*DS<150 for q in stA): continue
                if abs(tang[A][i]@tang[Bc][j])>0.94: par.append((A,i,Bc,j,dd[i]))
    return over,pins,par
def design2(chains,plan,cls,clear=6.8,grade=0.045,passes=4,verbose=False):
    over,pins,par=relations(chains,plan,cls,verbose)
    crossed=set((min(u,l),max(u,l)) for (u,i,l,j) in over)
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
        # a road that lands on the ground for a short way between two bridges (profiles.lands) rests on that ground:
        # the crest it lands on, not the valleys either side that the smoothed terrain would average in
        ag=np.where(~dk)[0]
        for r in (np.split(ag,np.where(np.diff(ag)>1)[0]+1) if len(ag) else []):
            if r[0]>0 and r[-1]<m-1 and len(r)*DS<=60: need[c][r]=z0[r]
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
            bs=base[c]
            for (o,s0,s1,bl) in chains[c].get('follow',()):      # hand hint: run at a neighbour's level, not the hill's
                s_=d['s']; w_=np.clip((s_-s0)/bl,0,1)*np.clip((s1-s_)/bl,0,1); w_=w_*w_*(3-2*w_)
                jo=ktrees[o].query(d['P'])[1]; bs=bs*(1-w_)+cls[o]['z'][jo]*w_
            z=np.maximum(bs,envelope(mm,grade))
            z=gs(z,35.0); z=np.maximum(z,mm); z=gs(z,12.0)
            # pins at same-level crossings and parallel ties toward higher-priority chains
            delta=np.zeros(m); wt=np.zeros(m); hard=np.zeros(m,bool)
            for (A,i,B,j,dist) in par:
                if A==c and pr[B]<pr[A]: k,o,oj=i,B,j
                elif B==c and pr[A]<pr[B]: k,o,oj=j,A,i
                else: continue
                dz_=cls[o]['z'][oj]-z[k]
                w_=max(0.0,min(1.0,(14-dist)/7.0))*float(np.clip(1-(abs(dz_)-3)/3,0,1))
                # carriageways that overlap are one roadway (the shared trunk of a ramp before it splits): same level,
                # whatever the two profiles would do on their own
                if dist<=hw(c,chains)+hw(o,chains) and not c.startswith(('x','e_')) and not o.startswith(('x','e_')) and (dist<14 or (min(c,o),max(c,o)) not in crossed): w_=1.0; hard[k]=True
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
                # straight grade from end to end; a ramp that starts on a trunk shared with another ramp keeps the
                # trunk's level and runs straight from the split, and one that ends merging alongside another road
                # (the Patton Avenue ramp joining the I-26 to I-240 East ramp on one bridge) keeps that road's level
                # from where their lanes meet: two halves of one deck at different heights leave holes between them
                hk=np.where(hard[:m//2])[0]; a0=int(hk[-1]) if len(hk) and hk[0]*DS<80 else 0
                he=np.where(hard[m//2:])[0]+m//2; b1=m-1
                if len(he) and (m-1-he[-1])*DS<80:
                    br=np.where(np.diff(he)>1)[0]; b1=int(he[br[-1]+1] if len(br) else he[0])
                line=np.linspace(z[a0],z[b1],b1-a0+1); z[a0:b1+1]=line+(1-chains[c]['linear'])*(z[a0:b1+1]-line)
            for (u,i,l,j) in over:   # hard clearance repair
                if u!=c: continue
                lack=cls[l]['z'][j]+clear-z[max(0,i-3):i+4].min()
                if lack>0: 
                    t=np.clip(1-np.abs(np.arange(m)-i)*DS/(40+lack/grade),0,1); z=z+lack*t*t*(3-2*t)
            d['z']=gs(z,12.0)
    return cls
HWD={'b0':22.0,'b6':14.0,'b10':12.0,'b23':12.0,'b8':9.0,'b9':9.0,'b25':9.0,'b1':8.0,'b3':7.0,'b7':7.0,'b2':6.0,'b5':6.0}
def hw(c,chains=None):
    if c in HWD: return HWD[c]
    if c.startswith('e_'): return 9.0
    if c.startswith('x'): return 8.0 if chains and chains[c].get('cls')=='motorway' else 5.0
    if c.startswith('y'): return 6.0
    return 4.5
def split_stacked(chains,plan,cls,verbose=True):
    """Two roads that cross at an angle while both sit on the same bridge outline are stacked, not level: the plan's
    structure fills merge where one bridge passes over another. The road with the longer run on that outline is the
    bridge itself and goes over; the other road's stretch is cut out into its own deck, below."""
    decks=plan['decks']; ids=[c for c in cls if not c.startswith(('x','e_'))]
    lines={c:sg.LineString(cls[c]['P']) for c in ids}; stack=set(plan.get('stack',()))
    def run(c,i):
        st=cls[c]['st']; k=st[i]; a=i; b=i
        while a>0 and st[a-1]==k: a-=1
        while b<len(st)-1 and st[b+1]==k: b+=1
        return a,b
    n=0
    for x in range(len(ids)):
        for y in range(x+1,len(ids)):
            A,B=ids[x],ids[y]
            if lines[A].distance(lines[B])>1: continue
            X=lines[A].intersection(lines[B])
            for p in [g for g in (X.geoms if hasattr(X,'geoms') else [X]) if g.geom_type=='Point']:
                i=min(int(round(lines[A].project(p)/DS)),len(cls[A]['P'])-1); j=min(int(round(lines[B].project(p)/DS)),len(cls[B]['P'])-1)
                k=cls[A]['st'][i]
                if k<0 or cls[B]['st'][j]!=k: continue
                tA=np.gradient(cls[A]['P'],axis=0)[i]; tB=np.gradient(cls[B]['P'],axis=0)[j]
                if abs(tA@tB)/(np.hypot(*tA)*np.hypot(*tB)+1e-9)>0.985: continue      # merging or diverging, not crossing
                ra,rb=run(A,i),run(B,j)
                (up,ru),(lo,rl)=((A,ra),(B,rb)) if ra[1]-ra[0]>=rb[1]-rb[0] else ((B,rb),(A,ra))
                stl=cls[lo]['st']; nbk=[stl[q] for q in (rl[0]-1,rl[1]+1) if 0<=q<len(stl) and stl[q]>=0 and stl[q]!=k]
                if (rl[1]-rl[0]+1)*DS<60 and nbk:
                    # a short clip of the lower road, carried on from its own bridge next door: keep it on that bridge
                    stl[rl[0]:rl[1]+1]=nbk[0]; stack.add((up,lo)); n+=1
                    if verbose: print('stacked',up,'over',lo,'on bridge',k,'at',(round(p.x),round(p.y)),'-> back onto its own bridge',nbk[0])
                    continue
                Plo=cls[lo]['P'][rl[0]:rl[1]+1]
                if len(Plo)<2: continue
                hl=hw(lo,chains); hu=hw(up,chains)
                dk=decks[k]['poly']
                L=polyonly(dk.intersection(sg.LineString(Plo).buffer(hl+1.0,cap_style=2)))
                band=sg.LineString(cls[up]['P'][ru[0]:ru[1]+1]).buffer(hu+1.0,cap_style=2)
                rest=polyonly(dk.difference(L.difference(band)).buffer(-0.3).buffer(0.3))
                if L.area<40 or rest.area<40: continue
                # the two roads are separate structures here: drop what is left of the merged outline beside the lower
                # road (the gore between it and its neighbours, and the slivers along its edges) wherever no other road runs
                bands=[lines[c].buffer(hw(c,chains)+2.0,cap_style=2) for c in ids if c!=lo and lines[c].distance(L)<20.0]
                gore=polyonly(rest.intersection(L.buffer(14.0)).difference(shapely.unary_union(bands)))
                gore=[g for g in getattr(gore,'geoms',[gore]) if not g.is_empty and g.distance(L)<0.5]
                if gore:
                    rest2=polyonly(rest.difference(shapely.unary_union(gore).buffer(0.2)).buffer(-0.75).buffer(0.75))
                    if rest2.area>40:
                        if verbose: print('  removed %.0f m2 of deck'%(rest.area-rest2.area),k,'beside',lo)
                        rest=rest2
                decks[k]=dict(decks[k],poly=rest)
                decks.append({'poly':L,'iso':decks[k].get('iso',False),'src':'split','abut':polyonly(decks[k]['abut'].intersection(L.buffer(6))) if not decks[k]['abut'].is_empty else sg.Polygon(),
                              'width':2*hl,'own':lo})
                cls[lo]['st'][rl[0]:rl[1]+1]=len(decks)-1; stack.add((up,lo)); n+=1
                if verbose: print('stacked',up,'over',lo,'on bridge',k,'at',(round(p.x),round(p.y)),'-> new deck',len(decks)-1)
    plan['stack']=sorted(stack)
    return plan,cls
def polyonly(g):
    if g.geom_type in ('Polygon','MultiPolygon'): return g
    q=[x for x in getattr(g,'geoms',[]) if x.geom_type in ('Polygon','MultiPolygon')]
    return shapely.unary_union(q) if q else sg.Polygon()
def trim_spurs(plan,cls,r=6.0,verbose=True):
    """Bridge outlines traced from the plans sometimes carry short side tongues (wing walls, hatching, a drawn stub of
    the structure being replaced). Cut off the parts narrower than 2r that stick out of a wide deck, unless a road that
    rides the deck runs out along them."""
    for k,d in enumerate(plan['decks']):
        p=d['poly']
        if d.get('src')!='pdf' or p.geom_type!='Polygon': continue
        # notches and slots under 5 m wide are drawing artefacts (joint lines, gaps between traced fills), not openings
        q=p.buffer(2.5,join_style=2).buffer(-2.5,join_style=2)
        if q.geom_type=='Polygon' and q.area-p.area>3:
            if verbose: print('deck',k,'closed %.0f m2 of narrow notches'%(q.area-p.area))
            p=d['poly']=q
        runs=[]
        for c,v in cls.items():
            idx=np.where(v['st']==k)[0]
            for r_ in (np.split(idx,np.where(np.diff(idx)>1)[0]+1) if len(idx) else []):
                if len(r_)>1: runs.append(sg.LineString(v['P'][r_]))
        if not runs: continue
        core=p.buffer(-r,join_style=2).buffer(r,join_style=2).intersection(p)
        if core.area<0.6*p.area: continue
        rl=shapely.unary_union(runs).buffer(1.0); spurs=polyonly(p.difference(core.buffer(0.05)))
        drop=[g for g in getattr(spurs,'geoms',[spurs]) if g.area>15 and not g.intersects(rl)
              and g.boundary.intersection(core.buffer(0.3)).length<0.35*g.boundary.length]     # it sticks out: joined to the deck only at its root
        if not drop: continue
        new=polyonly(p.difference(shapely.unary_union(drop).buffer(0.1)).buffer(-0.3).buffer(0.3))
        if new.geom_type!='Polygon': new=max(new.geoms,key=lambda g:g.area)
        if verbose: print('deck',k,'trimmed %d side tongue(s), %.0f m2'%(len(drop),p.area-new.area),'at',[round(v) for v in p.centroid.coords[0]])
        d['poly']=new
        # the ground-level pavement has a hole where the tongues were drawn over it: pave the carriageways of the
        # roads passing under at grade, on and right around the tongues
        tg=shapely.unary_union(drop).buffer(12.0); bands=[]
        for c,v in cls.items():
            idx=np.where(v['st']<0)[0]
            for r_ in (np.split(idx,np.where(np.diff(idx)>1)[0]+1) if len(idx) else []):
                if len(r_)>1 and sg.LineString(v['P'][r_]).distance(tg)<1.0: bands.append(sg.LineString(v['P'][r_]).buffer(hw(c),cap_style=2))
        if bands:
            patch=shapely.unary_union(bands).intersection(tg).difference(new.buffer(-0.5))
            plan['G']=polyonly(shapely.unary_union([plan['G'],patch]).buffer(0))
    # small gaps fully enclosed by bridge deck (where one traced outline stops short of the next) are deck too:
    # give each to the deck it shares the longest edge with
    decks=plan['decks']; DU=shapely.unary_union([d['poly'] for d in decks])
    for g in getattr(DU,'geoms',[DU]):
        for h in (g.interiors if g.geom_type=='Polygon' else []):
            hp=sg.Polygon(h)
            if not 3<hp.area<400: continue
            # ...unless a road passes through it at ground level: then it is the opening that road runs under
            if any(((v['st']<0)&shapely.contains(hp.buffer(2.0),shapely.points(v['P']))).any() for v in cls.values()): continue
            k=max(range(len(decks)),key=lambda k:hp.boundary.intersection(decks[k]['poly'].buffer(0.3)).length)
            new=polyonly(shapely.unary_union([decks[k]['poly'],hp.buffer(0.3)]).buffer(0))
            if new.geom_type=='Polygon':
                decks[k]['poly']=new
                if verbose: print('deck',k,'filled an enclosed gap of %.0f m2 at'%hp.area,[round(v) for v in hp.centroid.coords[0]])
    return plan
def join_decks(chains,plan,cls,g=4.0,ge=7.0,verbose=True):
    """Close the gaps the traced outlines leave along one road's bridge:
      * between consecutive (or side-by-side) deck pieces the road rides, up to 2g apart;
      * between the end of a deck and the ground-level pavement the road continues on, up to 2ge apart.
    Pavement of roads running at another level (2 m or more above or below this road there) is left out, so the deck
    does not grow sideways toward a road it passes over. Needs the road heights, so it runs after the design."""
    decks=plan['decks']; G=plan['G']; ids=[c for c in cls if not c.startswith(('x','e_','y'))]
    lines={c:sg.LineString(cls[c]['P']) for c in ids}
    def close(U,r): return polyonly(U.buffer(r,join_style=2).buffer(-r,join_style=2).difference(U.buffer(0.01)))
    def give(k,pieces):
        new=polyonly(shapely.unary_union([decks[k]['poly']]+[q.buffer(0.05) for q in pieces]).buffer(0))
        if new.geom_type=='Polygon' or decks[k]['poly'].geom_type!='Polygon': decks[k]['poly']=new; return True
        return False
    done=set(); n1=n2=0; a1=a2=0.0
    for c in ids:
        st=cls[c]['st']; P=cls[c]['P']; h=hw(c,chains); ctree=cKDTree(P)
        for i in range(1,len(st)):
            a,b=int(st[i-1]),int(st[i])
            if a<0 or b<0 or a==b or (min(a,b),max(a,b)) in done: continue
            done.add((min(a,b),max(a,b)))
            if 'split' in (decks[a].get('src'),decks[b].get('src')) or decks[a].get('drawn') or decks[b].get('drawn'): continue
            A,B=decks[a]['poly'],decks[b]['poly']
            add=[q for q in getattr(close(shapely.unary_union([A,B]),g),'geoms',[]) if q.area>1 and q.distance(A)<0.5 and q.distance(B)<0.5]
            if add and give(a,add): n1+=len(add); a1+=sum(q.area for q in add)
        for k in np.unique(st[st>=0]):
            k=int(k)
            if decks[k].get('src')=='split' or decks[k].get('drawn'): continue
            idx=np.where(st==k)[0]; lo,hi=max(0,idx[0]-12),min(len(st)-1,idx[-1]+12)
            if hi-lo<2: continue
            run=sg.LineString(P[lo:hi+1]); band=run.buffer(h+1.0,cap_style=2); own=sg.LineString(P[idx[0]:idx[-1]+1]) if len(idx)>1 else None
            under=[]
            for o in ids:
                if o==c or lines[o].distance(run)>h+hw(o,chains)+6.0: continue
                Po=cls[o]['P']; jc=ctree.query(Po)[1]; off=(cls[o]['st']<0)&(np.abs(cls[o]['z']-cls[c]['z'][jc])>2.0)
                io=np.where(off)[0]
                for r_ in (np.split(io,np.where(np.diff(io)>1)[0]+1) if len(io) else []):
                    if len(r_)>1: under.append(sg.LineString(Po[r_]).buffer(hw(o,chains)+4.0))
            near=polyonly(G.intersection(band).difference(shapely.unary_union(under)) if under else G.intersection(band))
            if near.is_empty: continue
            K=decks[k]['poly']
            add=[q for q in getattr(close(shapely.unary_union([K,near]),ge),'geoms',[]) if 2<q.area<1200 and q.distance(K)<0.5 and q.distance(near)<0.5]
            add=[polyonly(q.intersection(band)) for q in add]; add=[q for q in add if not q.is_empty and q.area>2]
            if add and give(k,add):
                n2+=len(add); a2+=sum(q.area for q in add)
                for q in add:
                    if verbose and q.area>150: print('  deck',k,'(%s) extended %.0f m2 to the pavement at'%(c,q.area),[round(v) for v in q.centroid.coords[0]])
    if verbose: print('joined decks: %d gaps between deck pieces (%.0f m2), %d gaps to the pavement at deck ends (%.0f m2)'%(n1,a1,n2,a2))
    return plan
def lateral(D,p,nv,W,step=0.5):
    """extent of deck D across a road at one station: offsets (lo,hi) along the normal nv of the stretch of the
    cross-section inside D that holds the centreline (or comes nearest to it); None if the section misses D"""
    o=np.arange(-W,W+step/2,step); ins=shapely.contains(D,shapely.points(p+nv*o[:,None]))
    if not ins.any(): return None
    idx=np.where(ins)[0]; rs=np.split(idx,np.where(np.diff(idx)>1)[0]+1)
    r=min(rs,key=lambda r:0.0 if o[r[0]]<=0<=o[r[-1]] else min(abs(o[r[0]]),abs(o[r[-1]])))
    return float(o[r[0]]),float(o[r[-1]])
def span_gap(c,P,st,run,k,decks,h,verbose=True):
    """A road's bridge drawn in two pieces with a gap between them: the plans draw only the top structure where another
    bridge passes over (the I-240 West to I-26 South flyover over the I-26 river bridges), and traced fills sometimes
    stop short of each other. Carry the deck across the gap at the width it has on either side, following the road's
    curve: from the last full-width section before the gap to the first one after it, edges interpolated between the
    two. The fill joins the deck the road is assigned to there. Returns False when the road is not on a bridge on
    both sides of the gap (the caller then patches the gap along the carriageway instead)."""
    r0,r1=int(run[0]),int(run[-1]); n=len(P)
    if r0<1 or r1>n-2 or st[r0-1]<0 or st[r1+1]<0: return False
    ka,kb=int(st[r0-1]),int(st[r1+1]); A,B=decks[ka]['poly'],decks[kb]['poly']
    t=np.gradient(P,axis=0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9; nv=np.c_[-t[:,1],t[:,0]]; W=h+15.0
    def full(D,steps):
        # walk away from the gap until the deck's section stops changing: nearer stations are cut by the skewed end
        prev=None
        for i in steps:
            if not 0<=i<n: break
            e=lateral(D,P[i],nv[i],W)
            if e is None: prev=None; continue
            if prev is not None and abs(e[0]-prev[1][0])<0.6 and abs(e[1]-prev[1][1])<0.6: return prev
            prev=(i,e)
        return None
    fa=full(A,range(r0-1,r0-25,-1)); fb=full(B,range(r1+1,r1+25))
    if fa is None or fb is None: return False
    (ia,ea),(ib,eb)=fa,fb
    q=np.arange(ia,ib+1); u=(q-ia)/max(ib-ia,1)
    lo=ea[0]+(eb[0]-ea[0])*u; hi=ea[1]+(eb[1]-ea[1])*u
    fill=sg.Polygon(np.vstack([P[q]+nv[q]*lo[:,None],(P[q]+nv[q]*hi[:,None])[::-1]])).buffer(0)
    if fill.is_empty: return False
    if k==ka: new=shapely.unary_union([A,fill.difference(B.buffer(0.05)) if kb!=ka else fill])
    elif k==kb: new=shapely.unary_union([B,fill.difference(A.buffer(0.05))])
    else: new=None
    if new is not None:
        new=polyonly(new.buffer(0.05,join_style=2).buffer(-0.05,join_style=2))
        if new.geom_type!='Polygon': new=None
    if new is not None:
        added=new.area-decks[k]['poly'].area; decks[k]=dict(decks[k],poly=new)
    else:      # the road rides some other deck here: the fill is a deck of its own
        piece=polyonly(fill.difference(shapely.unary_union([A,B]).buffer(0.05)).buffer(-0.3).buffer(0.3))
        if piece.geom_type!='Polygon':
            if piece.is_empty: return False
            piece=max(piece.geoms,key=lambda g:g.area)
        decks.append({'poly':piece,'iso':decks[ka].get('iso',False),'src':'patch','abut':sg.Polygon(),'width':2*h,'own':c})
        st[r0:r1+1]=len(decks)-1; added=piece.area; k=len(decks)-1
    if verbose: print('bridge',c,'carried across a %.0f m gap between decks %d and %d (deck %d, +%.0f m2) at'%((r1-r0+2)*DS,ka,kb,k,added),P[(r0+r1)//2].round())
    return True
def fill_notches(chains,plan,cls,win=5,verbose=True):
    """Bridge outlines traced from the plans are notched where the surface roads underneath are drawn across the
    structure fill (the I-240 West to I-26 South flyover over the US 19/23 and Riverside Drive lanes): along a road's
    run on a deck, the deck's edge dips in and comes back out. Close those dips: at each station the deck reaches at
    least as far out as the lower of the widest sections within win stations (25 m) before and after it, but no
    further than the road's own carriageway (+1.5 m), and never into another deck. Tapers and deck ends, where the
    edge moves one way only, are left as drawn."""
    decks=plan['decks']; n=0; tot=0.0
    def close(a):
        m=len(a); l=np.array([np.nanmax(a[max(0,i-win):i+1]) for i in range(m)]); r=np.array([np.nanmax(a[i:i+win+1]) for i in range(m)])
        return np.minimum(l,r)
    for c,v in cls.items():
        if c.startswith(('x','e_','y')): continue
        st=v['st']; P=v['P']; h=hw(c,chains)+1.5
        t=np.gradient(P,axis=0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9; nv=np.c_[-t[:,1],t[:,0]]
        for k in np.unique(st[st>=0]):
            k=int(k); idx=np.where(st==k)[0]
            if decks[k].get('drawn'): continue
            for run in np.split(idx,np.where(np.diff(idx)>1)[0]+1):
                if len(run)<win: continue
                D=decks[k]['poly']; E=[lateral(D,P[i],nv[i],h+10.0) for i in run]
                if sum(e is None for e in E)>len(E)//3: continue
                lo=np.array([np.nan if e is None else max(-e[0],0.0) for e in E]); hi=np.array([np.nan if e is None else max(e[1],0.0) for e in E])
                if np.isnan(lo).all(): continue
                lo2=np.fmax(lo,np.minimum(close(lo),h)); hi2=np.fmax(hi,np.minimum(close(hi),h))
                ok=np.isfinite(lo2)&np.isfinite(hi2)
                dip=ok&(np.nan_to_num(np.maximum(lo2-lo,hi2-hi),nan=0.0)>0.4)
                if not dip.any(): continue
                others=shapely.unary_union([d['poly'] for j,d in enumerate(decks) if j!=k]).buffer(0.3)
                add=[]
                dq=np.where(dip)[0]
                for g_ in np.split(dq,np.where(np.diff(dq)>2)[0]+1):      # each notch: its stations and one either side
                    q=np.arange(max(0,g_[0]-1),min(len(run),g_[-1]+2)); q=q[ok[q]]
                    if len(q)<2: continue
                    i_=run[q]
                    band=sg.Polygon(np.vstack([P[i_]-nv[i_]*lo2[q][:,None],(P[i_]+nv[i_]*hi2[q][:,None])[::-1]])).buffer(0)
                    a_=polyonly(band.difference(D.buffer(0.05)).difference(others))
                    add+=[g for g in getattr(a_,'geoms',[a_]) if not g.is_empty and g.area>2.0 and g.distance(D)<0.2]
                if not add: continue
                new=polyonly(shapely.unary_union([D]+[g.buffer(0.05) for g in add]).buffer(0.05,join_style=2).buffer(-0.05,join_style=2))
                if new.geom_type!='Polygon': continue
                a=new.area-D.area; decks[k]=dict(decks[k],poly=new); n+=len(add); tot+=a
                if verbose: print('  deck',k,'(%s) %d notch(es) filled, +%.0f m2 at'%(c,len(add),a),P[run[dip]].mean(0).round())
    if verbose: print('notches filled: %d (%.0f m2)'%(n,tot))
    return plan
def patch_plan(chains,plan,cls):
    """fill holes in pavement along alignments (under crossing bridges, at plan seams, between bridge pieces)"""
    import shapely.ops as so
    G=plan['G']; decks=plan['decks']; DU=shapely.unary_union([d['poly'] for d in decks]); Gb=G.buffer(0.6); shapely.prepare(Gb)
    addG=[]; nd=0
    for c,d in cls.items():
        xr=c.startswith('x') and chains[c].get('ride')
        # (an Overture freeway riding a bridge over a street runs down the median here: pave_between_edges paves the
        # carriageways from the plan's edge lines instead)
        if xr and all(decks[k].get('over_street') for k in chains[c]['ride']): xr=False
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
            own=decks[k]['poly'].buffer(1.0); mm=(st==k)&~shapely.contains(own,pts)
            if mm.any():      # ...unless the deck still carries part of the carriageway there (centreline along its edge)
                tg=np.gradient(P,axis=0); tg/=np.hypot(tg[:,0],tg[:,1])[:,None]+1e-9; nn=np.c_[-tg[:,1],tg[:,0]]; offs=np.linspace(-h,h,23)
                for i_ in np.where(mm)[0]:
                    if shapely.contains(own,shapely.points(P[i_]+nn[i_]*offs[:,None])).mean()>=0.2: mm[i_]=False
            idx=np.where(mm)[0]
            for run in (np.split(idx,np.where(np.diff(idx)>1)[0]+1) if len(idx) else []):
                if len(run)<2: continue
                if span_gap(c,P,st,run,k,decks,h): continue
                a,b=d['s'][run[0]]-7,d['s'][run[-1]]+7; seg=so.substring(ls,max(a,0),min(b,ls.length))
                dp=seg.buffer(h,cap_style=2).difference(shapely.unary_union([dd['poly'] for dd in decks if not dd.get('iso')]).buffer(0.05))
                dp=dp.buffer(-0.5).buffer(0.5)
                if dp.is_empty: continue
                if dp.geom_type!='Polygon': dp=max(dp.geoms,key=lambda g:g.area)
                if dp.area<40: continue
                decks.append({'poly':dp,'iso':False,'src':'patch','abut':sg.Polygon(),'width':2*h,'own':c}); st[run]=len(decks)-1; nd+=1
    if addG: G=shapely.unary_union([G]+[g for g in addG if not g.is_empty]).buffer(0)
    # the plans draw a road's pavement fill straight through under its own bridges; no road is down there. Keep
    # ground-level pavement beneath a deck only along roads that really pass under it at grade.
    DU=shapely.unary_union([d['poly'] for d in decks]); keep=[]
    for c,d in cls.items():
        h=hw(c,chains)+3.0; idx=np.where(d['st']<0)[0]
        for run in (np.split(idx,np.where(np.diff(idx)>1)[0]+1) if len(idx) else []):
            if len(run)>1 and sg.LineString(d['P'][run]).distance(DU)<h: keep.append(sg.LineString(d['P'][run]).buffer(h,cap_style=2))
    gone=polyonly(G.intersection(DU.buffer(-0.5)).difference(shapely.unary_union(keep)).buffer(-0.5).buffer(0.7))
    if not gone.is_empty:
        G=polyonly(G.difference(gone).buffer(-0.3).buffer(0.3)); print('removed %.0f m2 of ground-level pavement under bridges'%gone.area)
        plan['Gunder']=gone
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
# Section A bridges that carry the freeway over a street. The 2018 design draws each as one outline that the street
# crosses on the skew, so the street (the longer diagonal run) was being put on the deck, above the freeway.
FREEWAY_OVER=[(284865,208127)]      # I-26 / I-240 over State Street (Overture flags the freeway as the bridge)
def freeway_over(plan,chains):
    """the freeway roadways and ramps crossing a FREEWAY_OVER bridge ride it; every other road passes under it"""
    for pt in FREEWAY_OVER:
        k=deck_at(plan,pt); D=plan['decks'][k]['poly']; nd=len(plan['decks'])
        plan['decks'][k]['drawn']=True     # its 2018 outline is complete: join_decks and fill_notches leave it as drawn
        plan['decks'][k]['over_street']=True   # (meshz.py: its clearance cap shapes the ground under it, not the pavement)
        for c,ch in chains.items():
            if not sg.LineString(ch['xy']).intersects(D): continue
            # (any stretch on it counts: a plan alignment that starts on the bridge rides it from its first station)
            if c.startswith(('b','e_')) or ch.get('cls')=='motorway': ch['ride']=sorted(set(ch.get('ride',[]))|{k}); ch['ride_min']=0
            else: ch['decks']=[j for j in (ch['decks'] if ch.get('decks') is not None else range(nd)) if j!=k]
            print('bridge',k,'over a street at',pt,':',c,'rides it' if k in ch.get('ride',()) else 'passes under it')
# (median alignment, box, clear below y): stretches where the 2018 design draws both carriageways' edge-of-travel lines
# but its pavement fill is missing (north of the State Street bridge, up to where the Section B sheets pave the northbound
# side again); south of 'clear below y' the KMZ's odd pavement pieces there, drawn 3 m east of the plan sheets' lines,
# are cleared from beside the new carriageways (north of it the Haywood Road on-ramp joins the southbound side)
EOT_FILL=[('b0',(284840,208140,285000,208400),208270)]
EOT_TAPER=40.0      # m over which the carriageways move across from the KMZ's edge lines to the plan sheets'
def eot_sections(lines,P,nv,box):
    """per station of the median alignment: (outer, inner) edge offsets of the right carriageway and (inner, outer) of
    the left one, cast across to the edge lines; None where they are not two clean pairs"""
    EU=shapely.unary_union([sg.LineString(l) for l in lines if len(l)>1 and sg.LineString(l).intersects(box.buffer(40))])
    sec=[]
    for i in range(len(P)):
        if not box.contains(sg.Point(*P[i])): sec.append(None); continue
        X=sg.LineString([P[i]-nv[i]*30,P[i]+nv[i]*30]).intersection(EU)
        o=sorted(float((np.array(g.coords[0])-P[i])@nv[i]) for g in getattr(X,'geoms',[X]) if g.geom_type=='Point')
        o=[v for j,v in enumerate(o) if j==0 or v-o[j-1]>0.3]      # one line drawn twice (overlapping plan sheets)
        lo=[v for v in o if v<0][::-1]; hi=[v for v in o if v>0]
        ok=len(lo)>=2 and len(hi)>=2 and -lo[0]<10 and hi[0]<10 and 6<lo[0]-lo[1]<22 and 6<hi[1]-hi[0]<22
        sec.append(np.array((lo[1],lo[0],hi[0],hi[1])) if ok else None)
    return sec
def pave_between_edges(plan,chains):
    """pave each carriageway between its two edge-of-travel lines, cast across from the median alignment: the nearest
    line on each side is the carriageway's inner edge, the next one out its outer edge.
    The plan sheets' lines are used where they exist: the lane lines and the median alignment come from the same
    sheets. The 2018 KMZ draws this stretch about 3 m further east, as it does the bridge outline, so the carriageways
    leave the bridge on the KMZ's lines and move across to the plan sheets' over the first EOT_TAPER m."""
    K=pickle.load(open(W+'/kmz_sp.pkl','rb')); DU=shapely.unary_union([d['poly'] for d in plan['decks']])
    LY=pickle.load(open(W+'/pdf_layers.pkl','rb'))
    PE=[np.asarray(it[1] if isinstance(it,tuple) else it,float) for it in LY.get(('Prop EOT','s'),[])]
    for c,bx,yclr in EOT_FILL:
        box=sg.box(*bx)
        s_,P=resample(chains[c]['xy'],2.0); t=np.gradient(P,axis=0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9; nv=np.c_[-t[:,1],t[:,0]]
        sk=eot_sections(K['EOT'],P,nv,box); sp=eot_sections(PE,P,nv,box)
        iP=[i for i in range(len(P)) if sp[i] is not None]; s0=s_[iP[0]] if iP else np.inf
        sec=[]; wt=[]
        for i in range(len(P)):
            w=float(np.clip((s_[i]-s0)/EOT_TAPER,0,1)) if sp[i] is not None else 0.0
            if sk[i] is None: w=1.0 if sp[i] is not None else 0.0
            sec.append(None if sk[i] is None and sp[i] is None else sp[i] if w==1.0 else sk[i] if w==0.0 else (1-w)*sk[i]+w*sp[i])
            wt.append(w)
        quads=[]; qp=[]
        for i in range(len(P)-1):
            a,b=sec[i],sec[i+1]
            if a is None or b is None: continue
            for u,v in ((0,1),(2,3)):
                quads.append(sg.Polygon([P[i]+nv[i]*a[u],P[i+1]+nv[i+1]*b[u],P[i+1]+nv[i+1]*b[v],P[i]+nv[i]*a[v]]).buffer(0.05))
                if min(wt[i],wt[i+1])>=0.5: qp.append(quads[-1])
        if not quads: continue
        zone=shapely.unary_union(quads)
        # clear the KMZ's pieces from the verges and the median beside the new carriageways (only the pieces there change)
        RM=zone.buffer(6.0).intersection(sg.box(bx[0],bx[1],bx[2],yclr)).difference(zone.buffer(0.05)).difference(DU.buffer(2.0))
        parts=[]
        for q in getattr(plan['G'],'geoms',[plan['G']]):
            if not q.intersects(RM): parts.append(q); continue
            r=q.difference(RM); parts+=[g for g in getattr(r,'geoms',[r]) if g.geom_type=='Polygon' and g.area>1.0]
            print('cleared %.0f m2 of pavement beside the new carriageways at'%q.intersection(RM).area,np.round(q.representative_point().coords[0]))
        fill=polyonly(zone.difference(DU.buffer(0.1)).buffer(-0.3).buffer(0.3))
        fill=polyonly(fill.difference(sg.MultiPolygon(parts)))
        merge_fill(plan,parts,fill,box)
        # for marks.py: lines inside the zone are kept to its pavement, and where the plan sheets' edges govern, the
        # KMZ's own edge lines (3 m away) are left out
        plan.setdefault('eot_fill',[]).append({'zone':zone.difference(DU),'pdf':shapely.unary_union(qp).difference(DU) if qp else sg.Polygon()})
        print('paved %.0f m2 between edge lines along'%fill.area,c,'in',bx,'(plan sheets\' lines from %.0f m along it)'%s0)
def merge_fill(plan,parts,fill,box):
    """merge new pavement into the pieces of G it touches only; the rest of G stays exactly as it was (a full re-union
    reorders vertices, which shifts later simplified outlines, such as the fine-terrain tile set, far from here)"""
    near=[q.intersects(fill.buffer(0.5)) for q in parts]
    m=shapely.unary_union([q for q,t in zip(parts,near) if t]+[fill]); m=[q for q in getattr(m,'geoms',[m]) if q.geom_type=='Polygon']
    # (slivers left where a cleared piece met the new pavement)
    m=[sg.Polygon(q.exterior,[h for h in q.interiors if not (sg.Polygon(h).area<20 and box.contains(sg.Polygon(h)))]) for q in m]
    i0=near.index(True) if any(near) else len(parts)
    plan['G']=sg.MultiPolygon([q for q,t in zip(parts[:i0],near[:i0]) if not t]+m+[q for q,t in zip(parts[i0:],near[i0:]) if not t])
# (ramp's plan alignment, box): ramps whose first stretch the 2018 KMZ leaves unpaved though the plan sheets draw its
# edge-of-travel lines; it is paved between them, the shoulders widening to meet the KMZ pavement it runs into
RAMP_FILL=[('b12',(284940,208350,284990,208500))]     # exit from I-26 West / I-240 East to Haywood Road: first 130 m
def pave_ramp_gaps(plan,chains):
    LY=pickle.load(open(W+'/pdf_layers.pkl','rb'))
    PE=[np.asarray(it[1] if isinstance(it,tuple) else it,float) for it in LY.get(('Prop EOT','s'),[])]
    for c,bx in RAMP_FILL:
        box=sg.box(*bx); G=plan['G']; parts=list(getattr(G,'geoms',[G]))
        EU=shapely.unary_union([sg.LineString(l) for l in PE if len(l)>1 and sg.LineString(l).intersects(box.buffer(20))])
        s_,P=resample(chains[c]['xy'],2.0); t=np.gradient(P,axis=0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9; nv=np.c_[-t[:,1],t[:,0]]
        sec={}; strip={}
        for i in range(len(P)):
            if not box.contains(sg.Point(*P[i])): continue
            X=sg.LineString([P[i]-nv[i]*12,P[i]+nv[i]*12]).intersection(EU)
            o=sorted(float((np.array(g.coords[0])-P[i])@nv[i]) for g in getattr(X,'geoms',[X]) if g.geom_type=='Point')
            a=[v for v in o if -1.5<v<4]                 # the ramp's left edge, next to its alignment
            if not a: continue
            a=min(a,key=abs); b=[v for v in o if a-10<v<a-4]
            if not b: continue
            g_=[v for v in o if a<v<=a+4]                # the gore, paved up to where the through road's edge is 4 m away
            # past the gore's nose: where two plan sheets draw the through road's edge up to 2.5 m apart, pave between
            th=[v for v in o if v>a+4][:2]
            strip[i]=(th[0],th[1]) if not g_ and len(th)==2 and th[1]-th[0]<2.5 else None
            sec[i]=(max(g_) if g_ else a,a,max(b))
        if len(sec)<2: continue
        if min(sec)>0: sec[min(sec)-1]=sec[min(sec)]          # (back over the station where the edge lines start)
        for i in sorted(k for k in strip if strip[k] and not strip.get(k-1)): strip[i-1]=strip[i-2]=strip[i]   # (and its start)
        ii=sorted(sec); i1=ii[-1]
        # shoulders of the pavement it runs into, measured across its last station
        X=sg.LineString([P[i1]-nv[i1]*15,P[i1]+nv[i1]*15]).intersection(G); mid=P[i1]+nv[i1]*(sec[i1][1]+sec[i1][2])/2
        sl=sr=0.0
        for g in getattr(X,'geoms',[X]):
            if g.geom_type=='LineString' and g.distance(sg.Point(*mid))<0.1:
                e=sorted(float((np.array(q)-P[i1])@nv[i1]) for q in g.coords)
                sl=float(np.clip(e[-1]-sec[i1][1],0,4)); sr=float(np.clip(sec[i1][2]-e[0],0,4))
        quads=[]
        for i,j in zip(ii[:-1],ii[1:]):
            if j!=i+1: continue
            w=[(s_[k]-s_[ii[0]])/max(s_[i1]-s_[ii[0]],1e-9) for k in (i,j)]
            L=[max(sec[k][0],sec[k][1]+wk*sl) for k,wk in zip((i,j),w)]; R_=[sec[k][2]-wk*sr for k,wk in zip((i,j),w)]
            quads.append(sg.Polygon([P[i]+nv[i]*L[0],P[j]+nv[j]*L[1],P[j]+nv[j]*R_[1],P[i]+nv[i]*R_[0]]).buffer(0.05))
            if strip.get(i) and strip.get(j):
                u,v=strip[i],strip[j]; quads.append(sg.Polygon([P[i]+nv[i]*u[0],P[j]+nv[j]*v[0],P[j]+nv[j]*v[1],P[i]+nv[i]*u[1]]).buffer(0.05))
        fill=polyonly(shapely.unary_union(quads).buffer(-0.3).buffer(0.3))
        fill=polyonly(fill.difference(sg.MultiPolygon(parts)))
        merge_fill(plan,parts,fill,box)
        print('paved %.0f m2 of ramp between edge lines along'%fill.area,c,'in',bx,'(shoulders to %.1f / %.1f m)'%(sl,sr))
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
    chains['b18']['linear']=0.85     # on structure from the split to the I-240 flyover it joins: no dip into the valley
    # north of the river bridges I-26 comes down to the level of the US 19/23 lanes beside it (cut into the Montford
    # hillside), instead of riding over the hill: stations along b0 from the bridge end to where those lanes join it
    chains['b0']['follow']=[('b2',3530.0,4000.0,150.0)]
    # and the ramp from I-240 West climbs straight to I-26's level at the bridge end, without a hump over that hill
    chains['b5']['follow']=[('b0',900.0,1300.0,80.0)]
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
    freeway_over(plan,chains)
    plan['G']=G; pave_between_edges(plan,chains); pave_ramp_gaps(plan,chains)
    pickle.dump(plan,open(W+'/sheets_plan2.pkl','wb'))
    return plan,chains
if __name__=='__main__':
    plan,chains=load_all()
    cls=classify(chains,plan,noride={'b1','b10'},straddle=hw)
    plan=trim_spurs(plan,cls)
    chains=add_axis_chains(chains,plan,cls); cls=classify(chains,plan,noride={'b1','b10'},straddle=hw)
    plan,cls=split_stacked(chains,plan,cls)
    plan=patch_plan(chains,plan,cls); pickle.dump(plan,open(W+'/sheets_plan2.pkl','wb'))
    cls=design2(chains,plan,cls,verbose=True)
    plan=join_decks(chains,plan,cls); plan=fill_notches(chains,plan,cls); pickle.dump(plan,open(W+'/sheets_plan2.pkl','wb'))
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
