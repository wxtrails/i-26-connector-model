"""Clearance caps: everything in the corridor a short overpass spans must sit low enough to pass under it."""
import numpy as np
def build(plan,cls,ymax=208300.0,maxlen=130.0,clear=6.8):
    caps=[]
    for k,d in enumerate(plan['decks']):
        if d['poly'].centroid.y>ymax: continue
        best=None
        for c,v in cls.items():
            m=v['st']==k
            if m.sum()>=3 and (best is None or m.sum()>best[0]): best=(m.sum(),v['P'][m],v['z'][m])
        if best is None: continue
        P=best[1]; a=P[-1]-P[0]; L0=np.hypot(*a)
        if L0<10: continue
        a=a/L0; n=np.array([-a[1],a[0]]); cen=np.array(d['poly'].centroid.coords[0])
        # direction of the roads passing underneath (bridges are often skewed to them)
        import shapely
        pb=d['poly'].buffer(1.0); ub=None
        for c,v in cls.items():
            ins=shapely.contains(pb,shapely.points(v['P']))&(v['st']!=k)
            if ins.sum()>=3 and (ub is None or ins.sum()>ub[0]): Q=v['P'][ins]; ub=(ins.sum(),Q[-1]-Q[0])
        u=n
        if ub is not None and np.hypot(*ub[1])>5:
            u=ub[1]/np.hypot(*ub[1])
            if abs(u@a)>0.8: u=n
        Mi=np.linalg.inv(np.c_[a,u])
        E=(np.array(d['poly'].exterior.coords)-cen)@Mi.T; al=E[:,0]; la=E[:,1]
        L=al.max()-al.min(); wd=la.max()-la.min()
        if L>maxlen: continue
        cen=cen+a*(al.max()+al.min())/2+u*(la.max()+la.min())/2
        caps.append({'k':k,'cen':cen,'Mi':Mi,'L':L,'wd':wd,'z':float(best[2].min())-clear})
    return caps
def field(caps,X,Y,grade=0.045):
    out=np.full(np.shape(X),np.inf)
    for c in caps:
        dx=X-c['cen'][0]; dy=Y-c['cen'][1]
        al=np.abs(dx*c['Mi'][0,0]+dy*c['Mi'][0,1]); la=np.abs(dx*c['Mi'][1,0]+dy*c['Mi'][1,1])
        v=c['z']+grade*np.maximum(la-c['wd']/2,0)+1.0*np.maximum(al-(c['L']/2-2.0),0)
        v=np.where(la<c['wd']/2+300,v,np.inf)
        out=np.minimum(out,v)
    return out
