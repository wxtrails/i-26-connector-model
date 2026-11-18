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
        if d.get('walled') and d.get('ends'):     # along its end walls (today's road under it can run a few degrees off them)
            q=np.array(d['ends'][0].coords[-1])-np.array(d['ends'][0].coords[0]); u=q/np.hypot(*q)
        Mi=np.linalg.inv(np.c_[a,u])
        E=(np.array(d['poly'].exterior.coords)-cen)@Mi.T; al=E[:,0]; la=E[:,1]
        L=al.max()-al.min(); wd=la.max()-la.min()
        if L>maxlen: continue
        cen=cen+a*(al.max()+al.min())/2+u*(la.max()+la.min())/2
        c={'k':k,'cen':cen,'Mi':Mi,'L':L,'wd':wd,'z':float(best[2].min())-clear,'walled':bool(d.get('walled'))}
        if c['walled']:      # each end wall's line and its outward normal (the cap stops at the walls)
            dc=np.array(d['poly'].centroid.coords[0]); c['ends']=[]
            for E in d.get('ends',[]):
                p0=np.array(E.coords[0]); q=np.array(E.coords[-1])-p0; nO=np.array([-q[1],q[0]])/np.hypot(*q)
                c['ends'].append((p0,nO if (dc-p0)@nO<0 else -nO))
        caps.append(c)
    return caps
def field(caps,X,Y,grade=0.045):
    out=np.full(np.shape(X),np.inf)
    for c in caps:
        dx=X-c['cen'][0]; dy=Y-c['cen'][1]
        al=np.abs(dx*c['Mi'][0,0]+dy*c['Mi'][0,1]); la=np.abs(dx*c['Mi'][1,0]+dy*c['Mi'][1,1])
        v=c['z']+grade*np.maximum(la-c['wd']/2,0)+1.0*np.maximum(al-(c['L']/2-2.0),0)
        # a bridge standing on walls at its ends (profiles2.BR_DECK) spans just the corridor between them: behind the
        # walls the road carries on at its own level on fill, so the cap stops at the ends instead of easing out past them
        if c.get('walled'):
            past=np.zeros(np.shape(X),bool)
            for p0,nO in c.get('ends',()): past|=(X-p0[0])*nO[0]+(Y-p0[1])*nO[1]>0
            v=np.where((al<=c['L']/2+3.0)&~past,c['z']+grade*np.maximum(la-c['wd']/2,0),np.inf)
        v=np.where(la<c['wd']/2+300,v,np.inf)
        out=np.minimum(out,v)
    return out
