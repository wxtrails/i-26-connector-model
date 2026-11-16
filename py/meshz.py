"""Assign pavement mesh elevations from 1-D profiles (inverse-distance across alignments)."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import mesh as M, shapely
from scipy.spatial import cKDTree
def run(planfile='/sheets_plan2.pkl'):
    plan=pickle.load(open(W+planfile,'rb')); pr=pickle.load(open(W+'/profiles.pkl','rb')); cls=pr['cls']
    sheets=M.build(plan)
    off=np.cumsum([0]+[len(s['v']) for s in sheets]); XY=np.vstack([s['v'] for s in sheets])
    sid=np.concatenate([np.full(len(s['v']),i) for i,s in enumerate(sheets)])
    zt=M.terr(XY[:,0],XY[:,1]); z=zt.copy()
    groups={}
    for c,d in cls.items():
        for k in np.unique(d['st']):
            m=d['st']==k; groups.setdefault(int(k),[]).append((cKDTree(d['P'][m]),d['z'][m],c))
    from profiles2 import hw
    chains=pr['chains']
    allst=[(cKDTree(d['P']),d['z'],c) for c,d in cls.items()]
    allup=[(cKDTree(d['P'][d['st']>=0]),d['z'][d['st']>=0],c) for c,d in cls.items() if (d['st']>=0).any()]
    DU=shapely.unary_union([d['poly'] for d in plan['decks']]).buffer(1.0); shapely.prepare(DU)
    import caps as CP
    CAPS=CP.build(plan,cls); print('clearance caps',len(CAPS)); nlow=0
    ndrop=0
    for i,s in enumerate(sheets):
        P=s['v']; key=-1 if s['kind']=='G' else s['deck']
        gl=groups.get(key)
        fallback=False
        if not gl:
            fallback=True; best=None      # deck without a riding alignment: use the alignment with the longest run across it
            for c,d in cls.items():
                n_in=int(shapely.contains(s['poly'].buffer(1.0),shapely.points(d['P'])).sum())
                if n_in and (best is None or n_in>best[0]): best=(n_in,c)
            gl=[(cKDTree(cls[best[1]]['P']),cls[best[1]]['z'],best[1])] if best else allst
        under=shapely.contains(DU,shapely.points(P)) if key<0 else None
        cand=gl+(allup if key<0 else []); D=np.full((len(cand),len(P)),np.inf); Zc=np.zeros((len(cand),len(P))); Mr=np.full((len(cand),len(P)),np.inf); Ed=np.full((len(cand),len(P)),np.inf)
        for gi,(tree,zz,c) in enumerate(cand):
            d,j=tree.query(P,distance_upper_bound=400.0 if key>=0 else 60.0)
            ok=np.isfinite(d)
            if key<0 and gi>=len(gl): ok&=~under
            D[gi,ok]=d[ok]; Zc[gi,ok]=zz[j[ok]]; Mr[gi,ok]=d[ok]/hw(c,chains); Ed[gi,ok]=np.maximum(d[ok]-hw(c,chains),0.0)
            # ground-level pavement next to a bridge is a ground-level road's unless none is within 8 m of it
            if key<0 and gi>=len(gl): Ed[gi,ok]+=8.0
        if len(cand):
            # a point belongs to the road whose carriageway it lies in, else to the one whose edge is nearest (a wide
            # road must not claim the aprons of a narrow ramp beside it just because its half-width is large)
            n0=np.argmin(np.where(np.isfinite(Mr),Ed+0.01*Mr,np.inf),axis=0); dmin=D.min(axis=0); z0=Zc[n0,np.arange(len(P))]
            Wt=np.where(np.isfinite(D)&(np.abs(Zc-z0[None,:])<4.0),1.0/(np.where(np.isfinite(Mr),Mr,0)**2+0.05)**2,0.0)
            num=(Wt*Zc).sum(0); den=Wt.sum(0)
        else:
            num=np.zeros(len(P)); den=np.zeros(len(P)); dmin=np.full(len(P),1e9)
        zi=np.where(den>0,num/np.maximum(den,1e-30),np.nan)
        if key<0:
            wa=np.clip((50.0-dmin)/30.0,0,1); wa=wa*wa*(3-2*wa)
            zz_=np.where(den>0,zi*wa+zt[off[i]:off[i+1]]*(1-wa),zt[off[i]:off[i+1]])
            cf=CP.field(CAPS,P[:,0],P[:,1]); nlow+=int((zz_>cf+0.05).sum()); zz_=np.minimum(zz_,cf)
        else:
            zz_=np.where(den>0,zi,zt[off[i]:off[i+1]]+6)
            if fallback or (den==0).any(): print('deck',key,'fallback' if fallback else '','unassigned',int((den==0).sum()),'of',len(P),s['poly'].representative_point().coords[0])
        z[off[i]:off[i+1]]=zz_
        # drop impossible (very steep) triangles: they join pavements that are really at different levels
        v=s['v']; t=s['t']; pa=np.c_[v[t[:,0]],zz_[t[:,0]]]; pb=np.c_[v[t[:,1]],zz_[t[:,1]]]; pc=np.c_[v[t[:,2]],zz_[t[:,2]]]
        nrm=np.cross(pb-pa,pc-pa); slope=np.hypot(nrm[:,0],nrm[:,1])/(np.abs(nrm[:,2])+1e-9)
        keep=slope<0.40
        if not keep.any(): keep[:]=True
        ndrop+=int((~keep).sum()); s['t']=t[keep]; s['tw']=t[~keep]
    pickle.dump({'sheets':[{k:v for k,v in s.items() if k!='meta'} for s in sheets],'XY':XY,'z':z,'sid':sid,'off':off,'zt':zt,'caps':CAPS},open(W+'/mesh.pkl','wb'))
    print('lowered',nlow); print('verts',len(XY),'dropped steep triangles',ndrop)
if __name__=='__main__': run()
