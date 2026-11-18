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
    groups={}; runs={}
    for c,d in cls.items():
        st=d['st']; on=st>=0
        t=np.gradient(d['P'],axis=0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9
        # where each ground-level stretch starts and ends, to tell a point beside a road's ground stretch from one past
        # its end (I-26 is 44 m wide: by distance alone, pavement under its bridge 20 m before it touches down is "inside")
        g=st==-1; runs[c]=(np.where(g)[0],d['P'],t,np.r_[True,~g[:-1]],np.r_[~g[1:],True])
        for k in np.unique(st):
            m=st==k
            if k>=0:
                # a road's run on this deck, plus the stations either side where it carries on along the next deck
                # (up to 60 m): where two pieces of one bridge meet at a skewed joint, both sides then take the height
                # of the station they are actually beside, not the end of their own run, and the joint stays level
                idx=np.where(m)[0]; ext=m.copy()
                for i0,step in ((idx[0],-1),(idx[-1],1)):
                    i=i0+step
                    while 0<=i<len(st) and on[i] and abs(i-i0)<=12: ext[i]=True; i+=step
                m=ext
            groups.setdefault(int(k),[]).append((cKDTree(d['P'][m]),d['z'][m],c))
    from profiles2 import hw
    chains=pr['chains']; ICC={c for c in chains if chains[c].get('ic')}; ICZ=plan.get('ic_zone'); IC_FADE=15.0
    if ICZ is not None: shapely.prepare(ICZ)
    allst=[(cKDTree(d['P']),d['z'],c) for c,d in cls.items()]
    allup=[(cKDTree(d['P'][d['st']>=0]),d['z'][d['st']>=0],c) for c,d in cls.items() if (d['st']>=0).any()]
    DU=shapely.unary_union([d['poly'] for d in plan['decks']]).buffer(1.0); shapely.prepare(DU)
    import caps as CP
    CAPS=CP.build(plan,cls); print('clearance caps',len(CAPS)); nlow=0
    # a freeway bridge over a street keeps its cap for the ground under it (terrain.py cuts the street's corridor), but
    # not for pavement: past both ends of the deck the cap would pull the freeway's own approach down to street level
    CAPG=[c for c in CAPS if not plan['decks'][c['k']].get('over_street')]
    # behind the end walls of a bridge standing on walls (profiles2.BR_DECK) the pavement is its own road's, at its level:
    # the roads it spans (at grade under it, or along its ends within 6 m) do not count there, however near their
    # centrelines run (the lanes of the road below end at the foot of the wall)
    import gapfill as GF
    WB=GF.wall_backs(plan); QB=shapely.unary_union([q for _,_,_,q in WB]) if WB else None; SPAN=set()
    if QB is not None: shapely.prepare(QB)
    for k in sorted({k for k,_,_,_ in WB}):
        near=plan['decks'][k]['poly'].buffer(6.0)
        for c,d in cls.items():
            if (d['st']==k).any(): continue
            if (shapely.contains(near,shapely.points(d['P']))&(d['st']<0)).any(): SPAN.add(c)
    if WB: print('walled bridge ends: pavement behind them leaves out',', '.join(sorted(SPAN)))
    ndrop=0
    for i,s in enumerate(sheets):
        P=s['v']; key=-1 if s['kind']=='G' else s['deck']
        gl=groups.get(key)
        fallback=False
        if not gl:
            fallback=True; best=None      # deck without a riding alignment: use the alignment with the longest run across it
            for c,d in cls.items():
                if c in ICC: continue
                n_in=int(shapely.contains(s['poly'].buffer(1.0),shapely.points(d['P'])).sum())
                if n_in and (best is None or n_in>best[0]): best=(n_in,c)
            gl=[(cKDTree(cls[best[1]]['P']),cls[best[1]]['z'],best[1])] if best else allst
        under=shapely.contains(DU,shapely.points(P)) if key<0 else None
        behind=GF.behind(WB,P) if key<0 and QB is not None and s['poly'].intersects(QB.buffer(1.0)) else None
        cand=gl+(allup if key<0 else []); D=np.full((len(cand),len(P)),np.inf); Zc=np.zeros((len(cand),len(P))); Mr=np.full((len(cand),len(P)),np.inf); Ed=np.full((len(cand),len(P)),np.inf); Er=np.full((len(cand),len(P)),np.inf)
        for gi,(tree,zz,c) in enumerate(cand):
            d,j=tree.query(P,distance_upper_bound=400.0 if key>=0 else 60.0)
            ok=np.isfinite(d)
            if key<0 and gi>=len(gl): ok&=~under
            if behind is not None and c in SPAN: ok&=~behind
            D[gi,ok]=d[ok]; Zc[gi,ok]=zz[j[ok]]; Mr[gi,ok]=d[ok]/hw(c,chains); Ed[gi,ok]=np.maximum(d[ok]-hw(c,chains),0.0)
            if key<0 and gi<len(gl):
                # ground pavement well past the end of a road's ground stretch (where it goes onto a bridge, or ends) is not inside
                # that road's carriageway, however wide the road: it counts as outside by however far it lies beyond
                # 10 m (a skewed joint leaves the road's own pavement up to about that far past its centreline's
                # transition), so a road passing under the bridge whose carriageway holds it takes it
                ix,PP,tt,first,last=runs[c]; i_=ix[j[ok]]; a_=np.einsum('ij,ij->i',P[ok]-PP[i_],tt[i_])
                past=np.where(first[i_]&(a_<0),-a_,0.0)+np.where(last[i_]&(a_>0),a_,0.0)
                Ed[gi,ok]=np.maximum(Ed[gi,ok],past-10.0)
            Er[gi,ok]=Ed[gi,ok]
            # ground-level pavement next to a bridge is a ground-level road's unless none is within 8 m of it
            if key<0 and gi>=len(gl): Ed[gi,ok]+=8.0
        def mix(sel):
            if not sel.any(): return np.full(len(P),np.nan),np.zeros(len(P)),np.full(len(P),1e9)
            D_,Zc_,Mr_,Ed_,Er_=D[sel],Zc[sel],Mr[sel],Ed[sel],Er[sel]
            # a point belongs to the road whose carriageway it lies in, else to the one whose edge is nearest (a wide
            # road must not claim the aprons of a narrow ramp beside it just because its half-width is large)
            n0=np.argmin(np.where(np.isfinite(Mr_),Ed_+0.01*Mr_,np.inf),axis=0); dmin=D_.min(axis=0); z0=Zc_[n0,np.arange(len(P))]
            # blend only with roads whose edge is about as near as that road's: pavement inside one carriageway keeps its
            # level even with another road a few metres higher beyond the verge (I-26 and the I-26 to I-240 East ramp
            # landing side by side on the hill west of the river)
            # (fading back to the full blend over the first 6 m outside the carriageway, so there is no crease at its edge)
            e0=Er_[n0,np.arange(len(P))][None,:]
            far=np.clip((np.where(np.isfinite(Er_),Er_,1e9)-e0-2.0)/3.0,0,1); keepw=1-far*(1-np.clip(e0/6.0,0,1)) if key<0 else 1.0
            Wt=np.where(np.isfinite(D_)&(np.abs(Zc_-z0[None,:])<4.0),keepw/(np.where(np.isfinite(Mr_),Mr_,0)**2+0.05)**2,0.0)
            num=(Wt*Zc_).sum(0); den=Wt.sum(0)
            return np.where(den>0,num/np.maximum(den,1e-30),np.nan),den,dmin
        icm=np.array([c in ICC for (_,_,c) in cand],bool)
        if key>=0: zi,den,dmin=mix(np.ones(len(cand),bool))
        else:
            zi,den,dmin=mix(~icm)
            # the interchange south of State Street (profiles2.ic_build): its own roads' levels, today's roads' outside
            # it, blended over the zone's outer IC_FADE m
            wz=np.zeros(len(P))
            if ICZ is not None and icm.any() and s['poly'].intersects(ICZ):
                pts=shapely.points(P); inz=shapely.contains(ICZ,pts)
                if inz.any(): u=np.clip(shapely.distance(ICZ.boundary,pts[inz])/IC_FADE,0,1); wz[inz]=u*u*(3-2*u)
            if (wz>0).any():
                zb,denb,dminb=mix(icm); m=(wz>0)&(denb>0); both=m&(den>0)
                with np.errstate(invalid='ignore'):     # np.where works out both branches: inf*0 in cells it then discards
                    zi=np.where(both,wz*zb+(1-wz)*np.nan_to_num(zi),np.where(m,zb,zi))
                    dmin=np.where(both,wz*dminb+(1-wz)*dmin,np.where(m,dminb,dmin)); den=np.where(m,np.maximum(den,denb),den)
        if key<0:
            wa=np.clip((50.0-dmin)/30.0,0,1); wa=wa*wa*(3-2*wa)
            zz_=np.where(den>0,zi*wa+zt[off[i]:off[i+1]]*(1-wa),zt[off[i]:off[i+1]])
            cf=CP.field(CAPG,P[:,0],P[:,1]); nlow+=int((zz_>cf+0.05).sum()); zz_=np.minimum(zz_,cf)
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
