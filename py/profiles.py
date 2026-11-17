"""1-D vertical design along alignments (beam on elastic foundation with clearance constraints)."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely, shapely.geometry as sg
import scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.spatial import cKDTree
from scipy.ndimage import map_coordinates, gaussian_filter
DEM=np.load(W+'/dem_sp5.npy'); DX0,DY1,DRES=277500.0,215200.0,5.0
DEMs=gaussian_filter(DEM,1.2)
def terr(x,y):
    c=(np.asarray(x)-DX0)/DRES-0.5; r=(DY1-np.asarray(y))/DRES-0.5
    return map_coordinates(DEMs,[r,c],order=1,mode='nearest')
DS=5.0
LAM_NEW=260.0
def resample(c,ds=DS):
    ls=sg.LineString(c); n=max(2,int(round(ls.length/ds))); s=np.linspace(0,ls.length,n+1)
    return s,np.array([ls.interpolate(t).coords[0] for t in s])
def longside(p):
    r=p.minimum_rotated_rectangle; xs,ys=r.exterior.xy
    return max(np.hypot(xs[k+1]-xs[k],ys[k+1]-ys[k]) for k in range(4))
def deck_at(plan,pt):
    p=sg.Point(*pt); return min(range(len(plan['decks'])),key=lambda k:plan['decks'][k]['poly'].distance(p))
def continue_runs(cid,s,P,st,decks,dpoly,h,ext=150.0):
    """A road that is on a bridge stays on the structure it runs into:
      * a stretch inside a deck that touches the deck the road is already on counts as on the bridge, however short
        (a divided road whose directions split onto two parallel structures runs along the seam between them);
      * past the end of its deck run, the road stays up while that same deck still covers part of its carriageway
        (half-width h), until the structure ends.
    Roads that merely pass under or beside a bridge never start a deck run, so they are not affected."""
    t=np.gradient(P,axis=0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9; n=np.c_[-t[:,1],t[:,0]]
    offs=np.linspace(-h,h,23); pts=shapely.points(P); before=st.copy()
    def runs(mask):
        idx=np.where(mask)[0]
        return [r for r in np.split(idx,np.where(np.diff(idx)>1)[0]+1) if len(r)] if len(idx) else []
    changed=True
    while changed:
        changed=False
        for k,dp in enumerate(dpoly):
            ins=shapely.contains(dp,pts)&(st<0)
            for r in runs(ins):
                for e in (r[0]-1,r[0]-2,r[0]-3,r[-1]+1,r[-1]+2,r[-1]+3):
                    if 0<=e<len(s) and st[e]>=0 and st[e]!=k and dpoly[st[e]].distance(dp)<1.0:
                        st[r]=k; changed=True; break
    def cov(k,i): return shapely.contains(dpoly[k],shapely.points(P[i]+n[i]*offs[:,None])).mean()
    for r in runs(st>=0):
        for i0,step in ((r[0],-1),(r[-1],1)):
            k=st[i0]; i=i0+step
            while 0<=i<len(s) and st[i]<0 and abs(s[i]-s[i0])<ext and cov(k,i)>=0.25:
                st[i]=k; i+=step
    for r in runs(st!=before):
        print('continued',cid,'on deck',sorted(set(int(k) for k in st[r])),'%.0f m'%(s[r[-1]]-s[r[0]]+DS),'at',P[r[0]].round(),P[r[-1]].round())
    return st
def lands(decks,ka,pa,kb,pb,r=8.0):
    """The road comes off one bridge and onto another a short way on, and the plans draw an abutment at both of those
    bridge ends: the road lands on the ground in between (I-26 and both I-240 ramps on the hill between Emma Road and
    the French Broad). A gap with no abutment at one side is a break in the traced outline, and the bridge goes on."""
    if ka==kb: return False
    for k,p in ((ka,pa),(kb,pb)):
        a=decks[k].get('abut')
        if a is None or a.is_empty or a.distance(sg.Point(*p))>r: return False
    return True
def classify(chains,plan,noride=(),straddle=None):
    """chains: dict id -> dict(xy=Nx2, kind) ; returns per chain s, P, status(-1 at grade / deck index)"""
    decks=plan['decks']; out={}
    dpoly=[d['poly'].buffer(1.0) for d in decks]; dlong=[longside(d['poly']) for d in decks]
    for cid,ch in chains.items():
        s,P=resample(ch['xy']); st=np.full(len(s),-1)
        pts=shapely.points(P)
        allowed=ch.get('decks')    # None = any
        if cid not in noride:
            for k,dp in enumerate(dpoly):
                if allowed is not None and k not in allowed: continue
                ins=shapely.contains(dp,pts)
                if not ins.any(): continue
                # contiguous runs
                i=0
                while i<len(s):
                    if ins[i]:
                        j=i
                        while j+1<len(s) and ins[j+1]: j+=1
                        run=s[j]-s[i]
                        if run>=150 or run>=0.6*dlong[k] or (k in ch.get('ride',()) and run>ch.get('ride_min',30)):
                            st[i:j+1]=k
                        i=j+1
                    else: i+=1
        # stay on the structure a road is already on (see continue_runs)
        if cid not in noride and straddle is not None and not cid.startswith(('x','e_')):
            st=continue_runs(cid,s,P,st,decks,dpoly,straddle(cid,chains))
        # fill short at-grade gaps between deck runs, except where the plans end both bridges at the gap (see lands)
        i=0
        while i<len(s):
            if st[i]<0:
                j=i
                while j+1<len(s) and st[j+1]<0: j+=1
                if i>0 and j<len(s)-1 and (s[j]-s[i])<60:
                    if lands(decks,st[i-1],P[i-1],st[j+1],P[j+1]):
                        print('lands',cid,'at grade for %.0f m between decks %d and %d at'%(s[j]-s[i]+DS,st[i-1],st[j+1]),P[(i+j)//2].round())
                    else: st[i:j+1]=st[i-1]
                i=j+1
            else: i+=1
        out[cid]={'s':s,'P':P,'st':st}
    return out
def design(chains,plan,cls,clear=6.8,verbose=True):
    decks=plan['decks']; Gex=plan.get('Gexist')
    ids=list(chains.keys()); off={}; n=0
    for c in ids: off[c]=n; n+=len(cls[c]['s'])
    N=n
    zt=np.zeros(N); w=np.zeros(N); isdeck=np.zeros(N,bool); need=np.zeros(N)
    rows=[];cols=[];vals=[]
    def add(i,j,v): rows.append(i);cols.append(j);vals.append(v)
    for c in ids:
        d=cls[c]; o=off[c]; m=len(d['s']); P=d['P']; ch=chains[c]
        z0=terr(P[:,0],P[:,1]); zt[o:o+m]=z0
        lam=np.full(m,ch.get('lam',LAM_NEW))
        for s0,s1,lv in ch.get('lamseg',()): lam[(d['s']>=s0)&(d['s']<=s1)]=lv
        if Gex is not None and not Gex.is_empty and not ch.get('new'):
            ex=shapely.contains(Gex.buffer(1.0),shapely.points(P)); lam[ex]=np.minimum(lam[ex],ch.get('lam_ex',35.0))
        ww=1.0/lam**4; dk=d['st']>=0; ww[dk]=0.0
        w[o:o+m]=ww; isdeck[o:o+m]=dk
        # distance (along chain) to nearest deck end -> terrain clearance requirement
        idx=np.arange(m); nd=np.full(m,1e9)
        ends=[i for i in range(m) if dk[i] and (i==0 or not dk[i-1] or i==m-1 or not dk[i+1])]
        for e in ends: nd=np.minimum(nd,np.abs(idx-e)*DS)
        need[o:o+m]=np.where(dk,np.clip(0.10*(nd-10.0),0,ch.get('deckclear',6.0)),0)
        B=ch.get('B',1.0); T=ch.get('T',1.0/90.0**2)
        for i in range(1,m-1):   # bending
            a=[(o+i-1,1.0),(o+i,-2.0),(o+i+1,1.0)]; k=B/DS**3
            for p,vp in a:
                for q,vq in a: add(p,q,k*vp*vq)
        for i in range(m-1):     # tension
            k=T/DS; add(o+i,o+i,k);add(o+i+1,o+i+1,k);add(o+i,o+i+1,-k);add(o+i+1,o+i,-k)
    K=sp.coo_matrix((vals,(rows,cols)),shape=(N,N)).tocsr()
    # ---- relations between chains
    ties=[]   # (i,j,w)
    over=[]   # (upper i, lower j)
    lines={c:sg.LineString(cls[c]['P']) for c in ids}
    tang={}
    for c in ids:
        P=cls[c]['P']; t=np.gradient(P,axis=0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9; tang[c]=t
    trees={c:cKDTree(cls[c]['P']) for c in ids}
    for a_ in range(len(ids)):
        for b_ in range(a_+1,len(ids)):
            A,Bc=ids[a_],ids[b_]
            if lines[A].distance(lines[Bc])>14: continue
            PA,PB=cls[A]['P'],cls[Bc]['P']; sa,sb=cls[A]['st'],cls[Bc]['st']
            dd,jj=trees[Bc].query(PA)
            # crossings first
            X=lines[A].intersection(lines[Bc])
            pts=[g for g in (X.geoms if hasattr(X,'geoms') else [X]) if g.geom_type=='Point']
            stackedA=[]
            for p in pts:
                i=int(round(lines[A].project(p)/DS)); j=int(round(lines[Bc].project(p)/DS))
                i=min(i,len(PA)-1); j=min(j,len(PB)-1)
                if abs(tang[A][i]@tang[Bc][j])>0.985: continue      # merging, handled by ties
                ka,kb=sa[i],sb[j]
                if ka==kb:
                    for di in range(-1,2):
                        ties.append((off[A]+min(max(i+di,0),len(PA)-1),off[Bc]+j,2e-3*DS))
                    continue
                if ka>=0 and kb<0: up,lo=(A,i),(Bc,j)
                elif kb>=0 and ka<0: up,lo=(Bc,j),(A,i)
                else:
                    ia,ib=decks[ka].get('iso',False),decks[kb].get('iso',False)
                    la,lb=decks[ka].get('level',1 if ia else 0),decks[kb].get('level',1 if ib else 0)
                    if la==lb:
                        if verbose: print('  ambiguous deck-deck crossing',A,Bc,'at',(round(p.x),round(p.y)),'decks',ka,kb)
                        continue
                    up,lo=((A,i),(Bc,j)) if la>lb else ((Bc,j),(A,i))
                cu,iu=up; cl,il=lo
                stackedA.append(i)
                span=5   # +-25 m along the upper road
                for di in range(-span,span+1):
                    q=iu+di
                    if 0<=q<len(cls[cu]['s']): over.append((off[cu]+q,off[cl]+il))
                if verbose: print('  crossing: chain',cu,'over',cl,'at',(round(p.x),round(p.y)))
            # parallel ties on same sheet (not near a stacked crossing of the same pair)
            for i in np.where(dd<14)[0]:
                j=jj[i]
                if sa[i]!=sb[j]: continue
                if any(abs(i-q)*DS<300 for q in stackedA): continue
                par=abs(tang[A][i]@tang[Bc][j])
                if par>0.94: ties.append((off[A]+i,off[Bc]+j,1e-4*(1-dd[i]/14)**4*DS))
    ties=np.array(ties) if ties else np.zeros((0,3)); over=np.array(over,int) if over else np.zeros((0,2),int)
    ti=ties[:,0].astype(int); tj=ties[:,1].astype(int); tw=ties[:,2]
    Tm=sp.coo_matrix((np.concatenate([tw,tw,-tw,-tw]),(np.concatenate([ti,tj,ti,tj]),np.concatenate([ti,tj,tj,ti]))),shape=(N,N)).tocsr()
    A0=K+Tm
    act_o=np.zeros(len(over),bool); act_t=np.zeros(N,bool); z=None
    PEN=0.5
    for it in range(80):
        d=w*DS; rhs=w*DS*zt
        ti_=np.where(act_t)[0]; d[ti_]+=PEN; rhs[ti_]+=PEN*(zt[ti_]+need[ti_])
        if act_o.any():
            u=over[act_o,0]; l=over[act_o,1]; tgt=z[l]+clear
            np.add.at(d,u,PEN); np.add.at(rhs,u,PEN*tgt)
        zn=spl.spsolve((A0+sp.diags(d)).tocsc(),rhs)
        dz=np.abs(zn-z).max() if z is not None else 99
        z=zn
        if len(over):
            gap=z[over[:,0]]-z[over[:,1]]-clear
            v1=(gap<-0.03)&~act_o; r1=(gap>0.002)&act_o
        else: v1=r1=np.zeros(0,bool)
        gt=z-(zt+need)
        v2=isdeck&(gt<-0.03)&~act_t; r2=(gt>0.002)&act_t
        if verbose and it%5==0: print('iter',it,'add',int(v1.sum()),int(v2.sum()),'release',int(r1.sum()),int(r2.sum()),'dz %.3f'%dz)
        if v1.sum()==0 and v2.sum()==0 and r1.sum()==0 and r2.sum()==0 and dz<0.02: break
        if len(over): act_o=(act_o|v1)&~r1
        act_t=(act_t|v2)&~r2
    if verbose: print('done iter',it,'dz %.3f'%dz)
    for c in ids:
        o=off[c]; m=len(cls[c]['s']); cls[c]['z']=z[o:o+m]; cls[c]['zt']=zt[o:o+m]
    return cls
def load_chains_B():
    raw=pickle.load(open(W+'/align_chains.pkl','rb'))
    EXCL={4,16,17,21,35,36,20,27}
    ch={}
    for i,c in enumerate(raw):
        if i in EXCL or sg.LineString(c).length<45: continue
        L=sg.LineString(c).length
        ch[f'b{i}']={'xy':c,'lam':LAM_NEW if L>450 else 45.0}
    ch['b0']['lamseg']=[(0,1300,60.0),(3400,9999,60.0)]
    for k,v in {'b6':60,'b1':40,'b5':90,'b2':90,'b11':80,'b12':80,'b13':90,'b24':90,'b33':90,'b14':110,'b10':40,'b8':40,'b9':40,'b23':40,'b15':45,'b25':40}.items():
        if k in ch: ch[k]['lam']=float(v)
    return ch
RIDE={}   # (an earlier hand override put b2, the US 19/23 southbound lanes, on the I-26 bridge; it passes under it)
if __name__=="__main__":
    plan=pickle.load(open(W+'/sheets_plan.pkl','rb'))
    chains=load_chains_B()
    for c,pts in RIDE.items(): chains[c]['ride']=[deck_at(plan,p) for p in pts]
    cls=classify(chains,plan,noride={'b1','b10'})
    cls=design(chains,plan,cls)
    pickle.dump({'chains':chains,'cls':cls},open(W+'/profiles.pkl','wb'))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    show=sys.argv[1].split(',') if len(sys.argv)>1 else ['b0','b3','b7','b2','b5','b10']
    fg,axs=plt.subplots(len(show),1,figsize=(16,3.0*len(show)))
    for ax,c in zip(axs,show):
        d=cls[c]; ax.plot(d['s'],d['zt'],'g-',lw=0.8); ax.plot(d['s'],d['z'],'k-',lw=1.2); dk=d['st']>=0; ax.plot(d['s'][dk],d['z'][dk],'r.',ms=3)
        g=np.gradient(d['z'],d['s'])*100; ax.set_title(f'{c}: max grade {np.abs(g).max():.1f}%  decks {sorted(set(d["st"][dk]))}'); ax.grid(True)
    fg.tight_layout(); fg.savefig(W+'/img/prof2.png',dpi=70)
