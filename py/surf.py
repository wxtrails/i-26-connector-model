"""Sheet-aware surface elevation queries."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely, shapely.geometry as sg
import matplotlib.tri as mtri
class Surf:
    def __init__(self,mesh=None):
        m=mesh or pickle.load(open(W+'/mesh.pkl','rb'))
        self.m=m; self.sheets=m['sheets']; self.off=m['off']; self.z=m['z']
        self.polys=[s['poly'] for s in self.sheets]
        self.tree=shapely.STRtree(self.polys)
        self.tri=[None]*len(self.sheets); self.finder=[None]*len(self.sheets)
    def _prep(self,i):
        if self.tri[i] is None:
            s=self.sheets[i]; t=mtri.Triangulation(s['v'][:,0],s['v'][:,1],s['t'])
            self.tri[i]=t; self.finder[i]=t.get_trifinder()
    def z_on(self,i,P):
        """z of sheet i at points P (n,2); nan outside"""
        self._prep(i); P=np.asarray(P,float); s=self.sheets[i]
        ti=self.finder[i](P[:,0],P[:,1]); out=np.full(len(P),np.nan)
        ok=ti>=0
        if ok.any():
            tr_=s['t'][ti[ok]]; v=s['v']; zz=self.z[self.off[i]:self.off[i+1]]
            a=v[tr_[:,0]];b=v[tr_[:,1]];c=v[tr_[:,2]]; p=P[ok]
            d=(b[:,1]-c[:,1])*(a[:,0]-c[:,0])+(c[:,0]-b[:,0])*(a[:,1]-c[:,1])
            w0=((b[:,1]-c[:,1])*(p[:,0]-c[:,0])+(c[:,0]-b[:,0])*(p[:,1]-c[:,1]))/d
            w1=((c[:,1]-a[:,1])*(p[:,0]-c[:,0])+(a[:,0]-c[:,0])*(p[:,1]-c[:,1]))/d
            out[ok]=w0*zz[tr_[:,0]]+w1*zz[tr_[:,1]]+(1-w0-w1)*zz[tr_[:,2]]
        return out
    def cands(self,P):
        P=np.asarray(P,float)
        box=sg.box(P[:,0].min()-1,P[:,1].min()-1,P[:,0].max()+1,P[:,1].max()+1)
        return [int(i) for i in self.tree.query(box)]
    def drape(self,P,prefer=None):
        """polyline P -> (z, sheet index per vertex). Chooses for the whole polyline the sheet containing most vertices; other vertices fall back to other sheets."""
        P=np.asarray(P,float); best=None; res={}
        for i in self.cands(P):
            z=self.z_on(i,P); n=np.isfinite(z).sum()
            if n: res[i]=z
        if not res: return np.full(len(P),np.nan),np.full(len(P),-1)
        order=sorted(res,key=lambda i:(-(np.isfinite(res[i]).sum()+(1000 if i==prefer else 0)),0 if self.sheets[i]['kind']=='D' else 1))
        z=np.full(len(P),np.nan); sh=np.full(len(P),-1)
        for i in order:
            m=np.isnan(z)&np.isfinite(res[i]); z[m]=res[i][m]; sh[m]=i
        return z,sh
    def track(self,P,start=None):
        """path following: keep to the current sheet while it contains the point; else switch to a sheet containing it (prefer closest z)."""
        P=np.asarray(P,float); cs=self.cands(P); Z={i:self.z_on(i,P) for i in cs}
        z=np.full(len(P),np.nan); sh=np.full(len(P),-1); cur=start; lastz=None; lastk=0
        for k in range(len(P)):
            if cur is not None and cur in Z and np.isfinite(Z[cur][k]): pass
            else:
                opts=[i for i in cs if np.isfinite(Z[i][k])]
                if opts:
                    if lastz is None: cur=min(opts,key=lambda i:(0 if self.sheets[i]['kind']=='G' else 1))
                    else:
                        cur=min(opts,key=lambda i:abs(Z[i][k]-lastz))
                        if abs(Z[cur][k]-lastz)>3.0+0.09*np.hypot(*(P[k]-P[lastk])): cur=None
                else: cur=None
            if cur is not None: z[k]=Z[cur][k]; sh[k]=cur; lastz=z[k]; lastk=k
        return z,sh
