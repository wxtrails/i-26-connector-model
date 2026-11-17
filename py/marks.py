"""Lane lines, edge lines, barriers and walls draped on the road surface."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import os
from surf import Surf
import shapely, shapely.geometry as sg
from scipy.spatial import cKDTree
OX,OY=285000.0,210000.0
class Picker:
    """choose the pavement sheet (ground or a given deck) for linework using the nearest parallel alignment"""
    def __init__(self,surf):
        self.surf=surf; pr=pickle.load(open(W+'/profiles.pkl','rb'))['cls']
        P=[];T=[];S=[];Z=[]
        for c,d in pr.items():
            t=np.gradient(d['P'],axis=0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9
            P.append(d['P']);T.append(t);S.append(d['st']);Z.append(d['z'])
        self.P=np.vstack(P); self.T=np.vstack(T); self.S=np.concatenate(S); self.Z=np.concatenate(Z); self.tree=cKDTree(self.P)
        self.gs=[i for i,s in enumerate(surf.sheets) if s['kind']=='G']
        m=surf.m; stv=np.concatenate([np.full(len(s['v']),-1 if s['kind']=='G' else s['deck']) for s in surf.sheets])
        self.vt={}
        for k in np.unique(stv):
            mm=stv==k; self.vt[int(k)]=(cKDTree(m['XY'][mm]),m['z'][mm])
        self.vall=(cKDTree(m['XY']),m['z'])
        self.ds={}
        for i,s in enumerate(surf.sheets):
            if s['kind']=='D': self.ds.setdefault(s['deck'],[]).append(i)
    def near(self,p,t):
        """(status, station index) of the alignment a line at p running along t belongs to"""
        idx=self.tree.query_ball_point(p,32.0)
        if not idx: return None
        idx=np.array(idx); d=np.hypot(*(self.P[idx]-p).T); c=np.clip(np.abs(self.T[idx]@t),0,1)
        # direction counts for more than distance: lines of a road that passes under a bridge at a shallow angle lie
        # close to the bridge's alignment but are never parallel to it (2.5 m per degree)
        k=idx[np.argmin(d+2.5*np.degrees(np.arccos(c)))]
        return int(self.S[k]),int(k)
    def status(self,p,t):
        r=self.near(p,t); return None if r is None else r[0]
    def z(self,P):
        """heights for a polyline, choosing ground or deck per stretch from the nearest parallel alignment"""
        P=np.asarray(P,float)
        if len(P)<2: return self.surf.drape(P)[0]
        t=np.gradient(P,axis=0); n=np.hypot(t[:,0],t[:,1]); t=t/(n[:,None]+1e-9)
        res=[self.near(P[i],t[i]) if n[i]>1e-6 else None for i in range(len(P))]
        st=np.array([-2 if r_ is None else r_[0] for r_ in res]); zr=np.array([np.nan if r_ is None else self.Z[r_[1]] for r_ in res])
        # no flicker: stretches shorter than 3 vertices take the status of the stretch before them
        for _ in range(2):
            b=np.r_[0,np.where(np.diff(st)!=0)[0]+1,len(st)]
            for a_,b_ in zip(b[:-1],b[1:]):
                if b_-a_<3 and a_>0: st[a_:b_]=st[a_-1]
        z=np.full(len(P),np.nan); b=np.r_[0,np.where(np.diff(st)!=0)[0]+1,len(st)]
        for a_,b_ in zip(b[:-1],b[1:]):
            lo,hi=max(0,a_-1),min(len(P),b_+1)     # overlap one vertex so the pieces join
            zz=self._z(P[lo:hi],None if st[a_]==-2 else int(st[a_]),zr[lo:hi])
            seg=z[lo:hi]; m=np.isnan(seg); seg[m]=zz[m]; z[lo:hi]=seg
        return z
    def _z(self,P,st,zr=None):
        if len(P)<2: return self.surf.drape(P)[0]
        cand=self.surf.cands(P); order=[]
        if st is not None:
            pref=self.gs if st<0 else self.ds.get(st,[])
            order=[i for i in cand if i in pref]
        z=np.full(len(P),np.nan)
        for i in order:
            zi=self.surf.z_on(i,P); mm=np.isnan(z)&np.isfinite(zi); z[mm]=zi[mm]
        if np.isnan(z).any() and st is not None and st in self.vt:
            mm=np.isnan(z); d,j=self.vt[st][0].query(P[mm],distance_upper_bound=7.0); ok=np.isfinite(d)
            idx=np.where(mm)[0][ok]; z[idx]=self.vt[st][1][j[ok]]
        if np.isnan(z).any():
            z2,_=self.surf.drape(P); mm=np.isnan(z)
            if st is not None and st<0 and zr is not None:
                # a ground-level road's line with no ground pavement under it (a barrier along the edge of a road that
                # passes under a bridge, where the plan's pavement stops short of it): when the only surface there is
                # the deck overhead, it stays at the height of the road it runs along
                off=mm&np.isfinite(z2)&np.isfinite(zr)&(z2>np.nan_to_num(zr)+2.5); z[off]=zr[off]; mm&=~off
            z[mm]=z2[mm]
        if np.isnan(z).any():
            mm=np.isnan(z); d,j=self.vall[0].query(P[mm],distance_upper_bound=7.0); ok=np.isfinite(d)
            idx=np.where(mm)[0][ok]; z[idx]=self.vall[1][j[ok]]
        return z
def densify(P,step=6.0):
    P=np.asarray(P,float); out=[P[0]]
    for a,b in zip(P[:-1],P[1:]):
        n=max(1,int(np.ceil(np.hypot(*(b-a))/step)))
        for k in range(1,n+1): out.append(a+(b-a)*k/n)
    return np.array(out)
class Strips:
    def __init__(self): self.P=[];self.I=[];self.n=0
    def add(self,P,z,w,lift=0.05):
        P=np.asarray(P,float); ok=np.isfinite(z)
        if ok.sum()<2: return
        # split at nan
        idx=np.where(ok)[0]; runs=np.split(idx,np.where(np.diff(idx)>1)[0]+1)
        for r in runs:
            if len(r)<2: continue
            p=P[r]; zz=z[r]+lift; t=np.gradient(p,axis=0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9; nr=np.c_[-t[:,1],t[:,0]]*w/2
            a=p+nr; b=p-nr; m=len(p)
            v=np.empty((2*m,3),np.float32); v[0::2]=np.c_[a[:,0]-OX,zz,-(a[:,1]-OY)]; v[1::2]=np.c_[b[:,0]-OX,zz,-(b[:,1]-OY)]
            k=np.arange(m-1)*2; tri=np.c_[k,k+1,k+3,k,k+3,k+2].reshape(-1,3)+self.n
            self.P.append(v); self.I.append(tri); self.n+=2*m
    def wall(self,P,z,h,base=0.0):
        P=np.asarray(P,float); ok=np.isfinite(z)
        idx=np.where(ok)[0]
        if len(idx)<2: return
        for r in np.split(idx,np.where(np.diff(idx)>1)[0]+1):
            if len(r)<2: continue
            p=P[r]; zz=z[r]; m=len(p); v=np.empty((2*m,3),np.float32)
            v[0::2]=np.c_[p[:,0]-OX,zz+base,-(p[:,1]-OY)]; v[1::2]=np.c_[p[:,0]-OX,zz+h,-(p[:,1]-OY)]
            k=np.arange(m-1)*2; tri=np.c_[k,k+1,k+3,k,k+3,k+2].reshape(-1,3)+self.n
            self.P.append(v); self.I.append(tri); self.n+=2*m
    def get(self): return (np.vstack(self.P),np.vstack(self.I).astype(np.uint32)) if self.P else (np.zeros((0,3),np.float32),np.zeros((0,3),np.uint32))
def lines_of(L,key):
    out=[]
    for it in L.get(key,[]):
        a=np.asarray(it[1] if isinstance(it,tuple) else it,float)
        if len(a)>=2: out.append(a)
    return out
def dedupe(lines,tol=0.6):
    """drop duplicates from overlapping sheets (same midpoint & length)"""
    seen=set(); out=[]
    for a in lines:
        m=a[len(a)//2]; L=np.hypot(*np.diff(a,axis=0).T).sum(); k=(round(m[0]/tol),round(m[1]/tol),round(L/2))
        if k in seen: continue
        seen.add(k); out.append(a)
    return out
def simp(a,tol=0.06):
    if len(a)<=2: return a
    return np.array(sg.LineString(a).simplify(tol).coords)
def build(surf=None):
    surf=surf or Surf(); pk=Picker(surf)
    L=pickle.load(open(W+'/pdf_layers.pkl','rb')); K=pickle.load(open(W+'/kmz_sp.pkl','rb'))
    white=Strips(); conc=Strips()
    ll=dedupe(lines_of(L,('Prop Lane Lines','s'))); print('lane lines',len(ll))
    for a in ll:
        a=densify(simp(a),8.0); white.add(a,pk.z(a),0.15)
    eot=dedupe(lines_of(L,('Prop EOT','s')))+[a for a in K['EOT'] if a[:,1].mean()<208250]
    # chain EOT pieces coarsely: just draw each piece
    # an edge line with pavement on both sides is a plan match line drawn across a road that continues: skip those parts
    Gf=pickle.load(open(W+'/sheets_plan2.pkl','rb'))['G']; shapely.prepare(Gf)
    def open_edge(a):
        t=np.gradient(a,axis=0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9; nr=np.c_[-t[:,1],t[:,0]]
        both=shapely.contains(Gf,shapely.points(a+nr*1.5))&shapely.contains(Gf,shapely.points(a-nr*1.5))
        _,j=pk.tree.query(a); across=np.abs(np.sum(pk.T[j]*t,axis=1))<0.5      # runs across the nearest alignment
        return ~(both&across)
    n=0; ncut=0
    for a in eot:
        if np.hypot(*np.diff(a,axis=0).T).sum()<0.5: continue
        a=densify(simp(a,0.08),2.0); keep=open_edge(a); z=pk.z(a)
        if not keep.all(): ncut+=1
        idx=np.where(keep)[0]
        for r in (np.split(idx,np.where(np.diff(idx)>1)[0]+1) if len(idx) else []):
            if len(r)<2: continue
            b=densify(simp(a[r],0.08),8.0); white.add(b,pk.z(b),0.15); n+=1
    print('eot',n,'edge lines trimmed where the road continues',ncut)
    for a in dedupe(lines_of(L,('Prop Paint Stripe','s'))):
        a=densify(simp(a),8.0); white.add(a,pk.z(a),0.2)
    nb=0
    for a in dedupe(lines_of(L,('Prop Conc Barrier','s'))):
        if np.hypot(*np.diff(a,axis=0).T).sum()<0.8: continue
        a=densify(simp(a,0.08),8.0); conc.wall(a,pk.z(a),0.85); nb+=1
    for a in dedupe(lines_of(L,('Prop Noise Wall','s'))):
        a=densify(simp(a,0.08),8.0); z=pk.z(a)
        if np.isnan(z).all():
            from terrain import sample; T=np.load(W+'/terrain_conf.npy'); z=sample(T,a[:,0],a[:,1])
        conc.wall(a,z,4.5)
    print('barrier lines',nb)
    # lanes carried past the plan limits (lanes_ext.py)
    import lanes_ext as LE
    yellow=Strips(); T=pickle.load(open(W+'/lanes_ext.pkl','rb')) if os.path.exists(W+'/lanes_ext.pkl') else {}
    def dashes(a,z,period=12.0,on=3.0):
        s=np.r_[0,np.cumsum(np.hypot(*np.diff(a,axis=0).T))]
        for s0 in np.arange(0,s[-1]-on,period):
            ss=np.linspace(s0,s0+on,3); q=np.c_[np.interp(ss,s,a[:,0]),np.interp(ss,s,a[:,1])]; white.add(q,np.interp(ss,s,z),0.15)
    nt=0
    for a,kd in LE.paint(T):
        a=densify(simp(a,0.05),4.0); z=pk.z(a); nt+=1
        if kd=='dash': dashes(a,z)
        elif kd=='white': white.add(a,z,0.15)
        else: yellow.add(a,z,0.15)
    Graw=pickle.load(open(W+'/sheets_plan.pkl','rb'))['G']; Gf=pickle.load(open(W+'/sheets_plan2.pkl','rb'))['G']
    stubs=LE.stub_extensions(Graw,Gf,[sg.LineString(a) for a in ll],exclude=LE.pavement(T) if T else None)
    for a in stubs:
        a=densify(a,4.0); dashes(a,pk.z(a))
    print('carried lane lines',nt,'stub continuations',len(stubs))
    return white.get(),conc.get(),yellow.get()
if __name__=='__main__':
    (wp,wi),(cp,ci),(yp,yi)=build()
    print(wp.shape,wi.shape,cp.shape,ci.shape,yp.shape)
    pickle.dump({'white':(wp,wi),'conc':(cp,ci),'yellow':(yp,yi)},open(W+'/marks.pkl','wb'))
