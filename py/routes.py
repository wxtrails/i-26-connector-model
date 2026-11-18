"""Drive-through routes: lane paths along alignments."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
from surf import Surf
import shapely, shapely.geometry as sg, shapely.ops as so
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter1d
OX,OY=285000.0,210000.0
FW,RP,ST=27.0,17.0,15.5      # speeds m/s: freeway, ramp, street
# leg: (chain, dir, offset right of travel or 'auto', twoway, vmax)
def L(c,d=1,off=0.0,two=False,v=RP): return (c,d,off,two,v)
I26N=[L('e_south_in',1,0,False,FW),L('b0',1,9.0,True,FW)]
NORTH_OUT=[L('e_north_out',1,'auto',True,FW)]; NORTH_IN=[L('e_north_in',1,'auto',True,FW)]
SOUTH_OUT=[L('e_south_out',1,0,False,FW)]
E240_IN=[L('e_i240e_in',1,'auto',False,FW),L('b6',-1,7.0,True,24.0)]
E240_OUT=[L('b6',1,7.0,True,24.0),L('e_i240e_out',1,'auto',False,FW)]
PATW_IN=[L('e_patw_in',1,'auto',True,ST),L('b10',1,5.2,True,ST)]
PATW_OUT=[L('b10',-1,5.2,True,ST),L('e_patw_out',1,'auto',True,ST)]
PATE_IN=[L('e_pate_in',1,'auto',False,ST),L('b23',-1,5.2,True,ST),L('b9',-1,0,False,ST)]
PATE_OUT=[L('b8',1,0,False,ST),L('b23',1,5.2,True,ST),L('e_pate_out',1,'auto',False,ST)]
ROUTES=[
 ('i26n','I-26','I-26 West / northbound: I-40 side to Weaverville side',I26N+NORTH_OUT),
 ('i26s','I-26','I-26 East / southbound: Weaverville side to I-40 side',NORTH_IN+[L('b0',-1,9.0,True,FW)]+SOUTH_OUT),
 ('i240w_i26n','I-240 and I-26','I-240 West from downtown to I-26 north',E240_IN+[L('b5',1,0,False,20.0),L('b0',1,9.0,True,FW)]+NORTH_OUT),
 ('i26s_i240e','I-240 and I-26','I-26 from the north to I-240 East into downtown',NORTH_IN+[L('b0',-1,9.0,True,FW),L('b2',-1,0,False,20.0)]+E240_OUT),
 ('i240w_i26s','I-240 and I-26','I-240 West from downtown to I-26 south',E240_IN+[L('b7',-1,0,False,22.0),L('b0',-1,9.0,True,FW)]+SOUTH_OUT),
 ('i26n_i240e','I-240 and I-26','I-26 from the south to I-240 East into downtown',I26N[:1]+[L('b0',1,9.0,True,FW),L('b3',1,0,False,22.0)]+E240_OUT),
 ('pat_we','Patton Avenue','Patton Avenue: West Asheville to downtown',PATW_IN+PATE_OUT),
 ('pat_ew','Patton Avenue','Patton Avenue: downtown to West Asheville',PATE_IN+PATW_OUT),
 ('patw_i26n','Patton Avenue to I-26','Patton Avenue (from West Asheville) onto I-26 north',PATW_IN+[L('b8',1,0,False,ST),L('b22',1,0,False,RP),L('b0',1,9.0,True,FW)]+NORTH_OUT),
 ('pate_i26n','Patton Avenue to I-26','Patton Avenue (from downtown) onto I-26 north',PATE_IN+[L('b22',1,0,False,RP),L('b0',1,9.0,True,FW)]+NORTH_OUT),
 ('patw_i26s','Patton Avenue to I-26','Patton Avenue (from West Asheville) onto I-26 south',PATW_IN+[L('b24',-1,0,False,RP),L('b0',-1,9.0,True,FW)]+SOUTH_OUT),
 ('pate_i26s','Patton Avenue to I-26','Patton Avenue (from downtown) onto I-26 south',PATE_IN+[L('b10',-1,5.2,True,ST),L('b24',-1,0,False,RP),L('b0',-1,9.0,True,FW)]+SOUTH_OUT),
]
# I-26 south of State Street: the south-bound drives follow today's south-bound centreline (e_south_out), which runs on
# the new north-bound carriageway from just north of the interchange 650 m south of State Street to the bridge about
# 850 m south-west of it (over a side street). There they follow the middle of the new south-bound carriageway instead:
# these points (traced from its pavement, south to north), then the interchange's own centreline for it (profiles2 ic_sb)
SB_GUIDE=[(283694,207045),(283743,207062),(283792,207079),(283843,207099),(283890,207120),(283940,207144),(283990,207170)]
def tangents(P):
    t=np.gradient(P,axis=0); return t/(np.hypot(t[:,0],t[:,1])[:,None]+1e-9)
def resamp(P,ds):
    ls=sg.LineString(P); n=max(2,int(round(ls.length/ds))); return np.array([ls.interpolate(t).coords[0] for t in np.linspace(0,ls.length,n+1)])
class Builder:
    def __init__(self):
        self.cls=pickle.load(open(W+'/profiles.pkl','rb'))['cls']; self.surf=Surf()
        for c in list(self.cls):
            if c.endswith('_out') and c[:-4]+'_in' not in self.cls: self.cls[c[:-4]+'_in']={'P':self.cls[c]['P'][::-1].copy(),'same':True}
        self.guide_sb()
        self.plan=pickle.load(open(W+'/sheets_plan2.pkl','rb'))
        U=shapely.unary_union([self.plan['G']]+[d['poly'] for d in self.plan['decks']]).buffer(0.5); self.U=U; shapely.prepare(self.U)
        A=pickle.load(open(W+'/arrows.pkl','rb')); a=np.vstack([np.array(A['pdf'])[:,:4],np.array(A['kmz'])[:,:4]])
        n=np.hypot(a[:,2],a[:,3]); self.AP=a[:,:2]; self.AD=a[:,2:4]/n[:,None]; self.AT=cKDTree(self.AP)
    def guide_sb(self,fade=80.0):
        """move today's south-bound I-26 centreline onto the new south-bound carriageway's right-hand lane along SB_GUIDE
        (fading in and out over the guide's first and last fade m)"""
        if 'ic_sb' not in self.cls or 'e_south_out' not in self.cls: return
        G=np.vstack([np.array(SB_GUIDE,float),self.cls['ic_sb']['P']]); t=tangents(G)
        G=sg.LineString(G+np.c_[-t[:,1],t[:,0]]*1.8)      # in the right-hand lane (south-bound runs against the guide's direction)
        P=self.cls['e_south_out']['P']; pts=shapely.points(P)
        s=shapely.line_locate_point(G,pts); d=shapely.distance(G,pts); Q=shapely.get_coordinates(shapely.line_interpolate_point(G,s))
        ss=lambda u:np.clip(u,0,1)**2*(3-2*np.clip(u,0,1))
        w=ss(np.minimum(s,G.length-s)/fade)*ss((60.0-d)/10.0)
        self.cls['e_south_out']=dict(self.cls['e_south_out'],P=P+(Q-P)*w[:,None])
        print('south-bound I-26 guided onto its new carriageway over %.0f m (moved up to %.1f m)'%((w>0.5).sum()*np.median(np.hypot(*np.diff(P,axis=0).T)),(np.hypot(*(Q-P).T)*w).max()))
    def ctree(self,c):
        if not hasattr(self,'_ct'): self._ct={}
        if c not in self._ct: self._ct[c]=cKDTree(self.cls[c]['P'])
        return self._ct[c]
    def leg(self,c,d,off):
        if c not in self.cls and c.endswith('_in'): c=c[:-3]+'_out'; d=-d
        P=self.cls[c]['P']; P=P[::-1] if d<0 else P
        t=tangents(P); nr=np.c_[t[:,1],-t[:,0]]
        if off=='auto':
            opp=c.replace('_in','_X').replace('_out','_in').replace('_X','_out')
            if opp in self.cls and opp!=c and not self.cls[c].get('same') and not self.cls[opp].get('same'): dd=shapely.distance(sg.LineString(self.cls[opp]['P']),shapely.points(P))
            else: dd=np.zeros(len(P))
            o=5.4*np.clip(1-dd/9.0,0,1); o=gaussian_filter1d(o,6,mode='nearest'); return P+nr*o[:,None]
        return P+nr*off
    def join(self,A,B):
        lA,lB=sg.LineString(A),sg.LineString(B); X=lA.intersection(lB)
        pts=[g for g in (X.geoms if hasattr(X,'geoms') else [X]) if g.geom_type=='Point']
        if pts:
            p=max(pts,key=lambda g:lA.project(g)); sa,sb=lA.project(p),lB.project(p)
        else:
            dA=lA.distance(sg.Point(*B[0])); dB=lB.distance(sg.Point(*A[-1]))
            if dA<=dB: sa,sb=lA.project(sg.Point(*B[0])),0.0
            else: sa,sb=lA.length,lB.project(sg.Point(*A[-1]))
        a=np.array(so.substring(lA,0,max(sa,1.0)).coords); b=np.array(so.substring(lB,min(sb,lB.length-1.0),lB.length).coords)
        return a,b,float(np.hypot(*(a[-1]-b[0])))
    def snap(self,P,tol,sig):
        t=tangents(P); nr=np.c_[t[:,1],-t[:,0]]; sh=np.full(len(P),np.nan)
        for i,idx in enumerate(self.AT.query_ball_point(P,14.0)):
            if not idx: continue
            idx=np.array(idx); al=self.AD[idx]@t[i]; idx=idx[al>0.93]
            if not len(idx): continue
            lat=(self.AP[idx]-P[i])@nr[i]; k=np.argmin(np.abs(lat))
            if abs(lat[k])<=tol: sh[i]=lat[k]
        ok=np.isfinite(sh)
        if ok.sum()<3: return P
        # stations far from any arrow decay toward 0 shift only beyond 400 m
        full=np.interp(np.arange(len(P)),np.where(ok)[0],sh[ok])
        full=gaussian_filter1d(full,sig/4.0,mode='nearest')
        return P+nr*full[:,None]
    def pull(self,P):
        """keep the path on pavement: dynamic programming over lateral offsets (prefers staying put, then right)"""
        if not hasattr(self,'U2'):
            self.U2=self.U.buffer(-2.3); shapely.prepare(self.U2); self.U3=self.U.buffer(-3.9); shapely.prepare(self.U3)
        t=tangents(P); nr=np.c_[t[:,1],-t[:,0]]; offs=np.arange(-33.0,33.5,1.0); n=len(P); K=len(offs)
        C=np.empty((n,K))
        for k,o in enumerate(offs):
            pts=shapely.points(P+nr*o); i3=shapely.contains(self.U3,pts); i2=shapely.contains(self.U2,pts)
            C[:,k]=np.where(i3,0.0,np.where(i2,6.0,40.0))+(0.15 if o>=0 else 0.3)*abs(o)
        trans=6.0*(offs[:,None]-offs[None,:])**2
        cost=C[0].copy(); back=np.zeros((n,K),int)
        for i in range(1,n):
            m=cost[:,None]+trans; back[i]=np.argmin(m,axis=0); cost=m[back[i],np.arange(K)]+C[i]
        k=int(np.argmin(cost)); path=np.zeros(n,int)
        for i in range(n-1,-1,-1): path[i]=k; k=back[i][k]
        best=gaussian_filter1d(offs[path],4,mode='nearest')
        return P+nr*best[:,None]
    def build(self,legs):
        paths=[self.leg(c,d,off) for (c,d,off,two,v) in legs]
        segs=[paths[0]]; joints=[]
        for k in range(1,len(paths)):
            a,b,gap=self.join(segs[-1],paths[k]); segs[-1]=a; segs.append(b)
        P=[];leg_id=[]
        for k,s in enumerate(segs):
            s=resamp(s,4.0) if len(s)>1 else s
            if k>0: joints.append(len(np.vstack(P)))
            P.append(s); leg_id+= [k]*len(s)
        P=np.vstack(P); leg_id=np.array(leg_id)
        # resample whole to 4 m (keep leg ids by nearest)
        Q=resamp(P,4.0); tr=cKDTree(P); _,j=tr.query(Q); lid=leg_id[j]
        jn=np.where(np.diff(lid)!=0)[0]
        w=np.zeros(len(Q))
        for q in jn: w=np.maximum(w,np.clip(1-np.abs(np.arange(len(Q))-q)*4.0/140.0,0,1))
        w=gaussian_filter1d(w,5)
        sm=lambda X,s: np.c_[gaussian_filter1d(X[:,0],s/4.0,mode='nearest'),gaussian_filter1d(X[:,1],s/4.0,mode='nearest')]
        Q=sm(Q,24.0)*w[:,None]+sm(Q,7.0)*(1-w[:,None])
        Q=self.snap(Q,14.0,120.0); Q=self.pull(Q); Q=self.snap(Q,2.3,45.0)
        Q=sm(Q,6.0); Q=resamp(Q,4.0); _,j=tr.query(Q); lid=leg_id[j]
        zm,sh=self.surf.track(Q)
        zp=np.zeros(len(Q)); wsum=np.zeros(len(Q))     # profile elevation of the legs' alignments
        for k,(c,d,off,two,vv) in enumerate(legs):
            cc=c if 'z' in self.cls[c] else c[:-3]+'_out'
            dd,j=self.ctree(cc).query(Q); w_=np.clip(1-np.abs(np.arange(len(Q))-np.clip(np.arange(len(Q)),*( (np.where(lid==k)[0][[0,-1]]) if (lid==k).any() else (0,0))))*4.0/60.0,0,1)*(dd<40)
            if not (lid==k).any(): w_[:]=0
            zp+=w_*self.cls[cc]['z'][j]; wsum+=w_
        zp=np.where(wsum>0,zp/np.maximum(wsum,1e-9),np.nan); okp=np.isfinite(zp); zp=np.interp(np.arange(len(Q)),np.where(okp)[0],zp[okp])
        isG=np.array([s>=0 and self.surf.sheets[s]['kind']=='G' for s in sh])
        use=np.isfinite(zm)&((np.abs(zm-zp)<1.5)|(isG&(np.abs(zm-zp)<6.0))); nanfrac=1-use.mean()
        z=np.where(use,zm,zp); z=gaussian_filter1d(z,4.0,mode='nearest')
        # speed
        t=tangents(Q); hd=np.unwrap(np.arctan2(t[:,1],t[:,0])); kap=np.abs(gaussian_filter1d(np.gradient(hd)/4.0,3))
        vmax=np.array([legs[k][4] for k in lid],float)
        v=np.minimum(vmax,np.sqrt(2.6/np.maximum(kap,1e-5))); v=np.maximum(v,7.0)
        # forward/backward accel limits
        for i in range(1,len(v)): v[i]=min(v[i],np.sqrt(v[i-1]**2+2*1.6*4.0))
        for i in range(len(v)-2,-1,-1): v[i]=min(v[i],np.sqrt(v[i+1]**2+2*2.0*4.0))
        v=gaussian_filter1d(v,4,mode='nearest')
        # adjacent lane availability
        nr=np.c_[t[:,1],-t[:,0]]; mask=np.zeros(len(Q),np.uint8)
        for bit,sgn in ((1,-1.0),(2,1.0)):
            a=shapely.contains(self.U,shapely.points(Q+nr*sgn*3.6))&shapely.contains(self.U,shapely.points(Q+nr*sgn*6.0))
            if sgn<0:
                for k,(c,d,off,two,vv) in enumerate(legs):
                    if not two: continue
                    m=lid==k
                    if not m.any(): continue
                    cl=sg.LineString(self.cls[c]['P'] if d>0 else self.cls[c]['P'][::-1]); pts=Q[m]+nr[m]*sgn*5.4
                    # signed right distance from centreline
                    sgnd=[]
                    for p in pts:
                        s_=cl.project(sg.Point(*p)); q0=np.array(cl.interpolate(s_).coords[0]); q1=np.array(cl.interpolate(min(s_+1,cl.length)).coords[0]); tt=q1-q0
                        if np.hypot(*tt)<1e-6: q1=q0; q0=np.array(cl.interpolate(max(s_-1,0)).coords[0]); tt=q1-q0
                        sgnd.append((tt[1]*(p[0]-q0[0])-tt[0]*(p[1]-q0[1]))/(np.hypot(*tt)+1e-9))
                    a[m]&=np.array(sgnd)>0.8
            # require runs of at least 200 m
            a=a.astype(float); a=gaussian_filter1d(a,6)>0.98
            mask|=(a.astype(np.uint8)*bit)
        inp=shapely.contains(self.U,shapely.points(Q)).mean()
        return {'P':Q,'z':z,'v':v,'mask':mask,'lid':lid,'nan':nanfrac,'inpave':inp}
if __name__=='__main__':
    B=Builder(); out=[]
    for rid,grp,name,legs in ROUTES:
        r=B.build(legs); r.update(id=rid,group=grp,name=name); out.append(r)
        g=np.abs(np.gradient(r['z'])/4.0).max()*100
        print(rid,'len %.0f m'%(4*len(r['P'])),'time %.0f s'%(4.0/r['v']).sum(),'offpave %.1f%%'%(100*(1-r['inpave'])),'nan %.1f%%'%(100*r['nan']),'maxgrade %.1f%%'%g,'left %.0f%% right %.0f%%'%(100*(r['mask']&1).mean(),100*((r['mask']&2)>0).mean()))
    pickle.dump(out,open(W+'/routes.pkl','wb'))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    views=[(283700,207000,284300,207500),(284300,207300,284900,208200),(282850,206700,283700,207100),(285000,209300,285500,209900)]
    fg,axs=plt.subplots(2,2,figsize=(22,22))
    for ax,(x0,y0,x1,y1) in zip(axs.ravel(),views):
        box=sg.box(x0,y0,x1,y1)
        for g in B.plan['G'].geoms:
            if g.intersects(box):
                ax.fill(*g.exterior.xy,color='#ddd')
                for h in g.interiors: ax.fill(*h.xy,color='w')
        for d in B.plan['decks']:
            if d['poly'].intersects(box):
                for g_ in getattr(d['poly'],'geoms',[d['poly']]): ax.fill(*g_.exterior.xy,color='#f5c8c8',alpha=0.7)
        m=(B.AP[:,0]>x0)&(B.AP[:,0]<x1)&(B.AP[:,1]>y0)&(B.AP[:,1]<y1)
        ax.quiver(B.AP[m,0],B.AP[m,1],B.AD[m,0],B.AD[m,1],color='#888',scale=90,width=0.0015)
        for k,r in enumerate(out): ax.plot(r['P'][:,0],r['P'][:,1],lw=0.9,color=plt.cm.tab20(k))
        ax.set_xlim(x0,x1); ax.set_ylim(y0,y1); ax.set_aspect('equal')
    fg.tight_layout(); fg.savefig(W+'/img/routes.png',dpi=55)
