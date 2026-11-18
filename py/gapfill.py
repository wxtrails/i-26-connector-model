"""Ground between closely adjacent roads at different levels: fill up to the higher road, retaining wall at the lower one.

The terrain is a height grid (2.5 m raster, drawn on a 6 m grid), so it cannot stand a wall next to a road: beside a
low road it stays low for a grid step or two, and a ramp a few metres away and a few metres up is left hanging over a
pit. Here the strips between such roads get their own surface: level with the higher road for 6 m, then a 2:1 slope,
and where that still ends above the lower road, a vertical wall along the lower road's edge."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely, shapely.geometry as sg
from scipy.ndimage import distance_transform_edt, map_coordinates, binary_opening, binary_dilation
from rasterio.features import rasterize, shapes
from rasterio.transform import from_origin
import mesh as M
def walls(g,lev,ground,WV,WT,wn):
    """wall panels round a fill polygon, from its top down to the ground just outside, where that drop is 1 m or more
    (the exterior runs counter-clockwise and holes clockwise, so "to the right of travel" points out of the fill)"""
    for ring in [g.exterior]+list(g.interiors):
        ls=sg.LineString(ring.coords); m=max(3,int(ls.length/3.0)); P=np.array([ls.interpolate(s).coords[0] for s in np.linspace(0,ls.length,m+1)])
        tg=np.gradient(P,axis=0); tg/=np.hypot(tg[:,0],tg[:,1])[:,None]+1e-9; out=np.c_[tg[:,1],-tg[:,0]]
        top=lev(P[:,0],P[:,1]); Q=P+out*1.5; bot=np.minimum(top,ground(Q[:,0],Q[:,1]))-0.6
        for k in range(len(P)-1):
            if max(top[k]-bot[k],top[k+1]-bot[k+1])<1.0: continue
            WV.append(np.array([[P[k,0],P[k,1],top[k]],[P[k+1,0],P[k+1,1],top[k+1]],[P[k+1,0],P[k+1,1],bot[k+1]],[P[k,0],P[k,1],bot[k]]]))
            WT.append(np.array([[0,1,2],[0,2,3]])+wn); wn+=4
    return wn
def wall_backs(plan,depth=22.0,wing=12.0):
    """for each end of a bridge standing on walls (profiles2.BR_DECK): (deck index, end line, outward normal, the strip
    behind the wall: the end line, 'wing' m longer each way, swept 'depth' m back from the bridge)"""
    out=[]
    for k,d in enumerate(plan['decks']):
        if not d.get('walled'): continue
        cen=np.array(d['poly'].centroid.coords[0])
        for E in d.get('ends',[]):
            a,b=np.array(E.coords[0]),np.array(E.coords[-1]); e=(b-a)/E.length; nO=np.array([-e[1],e[0]])
            if (cen-a)@nO>0: nO=-nO           # away from the bridge
            q=sg.Polygon([a-e*wing,b+e*wing,b+e*wing+nO*depth,a-e*wing+nO*depth]).difference(d['poly'].buffer(0.05))
            out.append((k,E,nO,max(getattr(q,'geoms',[q]),key=lambda g:g.area)))
    return out
def behind(WB,P,depth=22.0,wing=12.0):
    """points past a walled bridge's end line (by any amount, up to 'depth' m), within 'wing' m beyond its ends"""
    m=np.zeros(len(P),bool)
    for k,E,nO,q in WB:
        a=np.array(E.coords[0]); e=(np.array(E.coords[-1])-a)/E.length; v=(P-a)@nO; u=(P-a)@e
        m|=(v>0)&(v<=depth)&(u>=-wing)&(u<=E.length+wing)
    return m
def abutments(zr,plan,SX0,SY1,RES,depth=22.0,wing=12.0,drop=0.3):
    """Fill behind the end walls of a bridge that stands on walls (profiles2.BR_DECK): beyond each end its road carries
    on at its own level, while the road it spans runs along the foot of the wall. The height grid cannot stand a wall
    (beside the low road it stays low for a grid step or two, a pit under the high road's pavement), so the strip behind
    each end gets its own surface, level with the pavement behind the wall (pavement more than 4 m below that, the road
    under the bridge, does not count), and a wall face down to the ground in front. Returns (polygon, level function)."""
    h,w=zr.shape; out=[]
    for k,E,nO,quad in wall_backs(plan,depth,wing):
        a=np.array(E.coords[0]); e=(np.array(E.coords[-1])-a)/E.length
        x0,y0,x1,y1=quad.buffer(30.0).bounds
        c0=max(0,int((x0-SX0)/RES)); c1=min(w,int((x1-SX0)/RES)+1); r0=max(0,int((SY1-y1)/RES)); r1=min(h,int((SY1-y0)/RES)+1)
        X,Y=np.meshgrid(SX0+(np.arange(c0,c1)+0.5)*RES,SY1-(np.arange(r0,r1)+0.5)*RES)
        v=(X-a[0])*nO[0]+(Y-a[1])*nO[1]; u=(X-a[0])*e[0]+(Y-a[1])*e[1]; zz=zr[r0:r1,c0:c1]
        ok=np.isfinite(zz)&(v>=1.5)
        near=ok&(v<10.0)&(u>=0)&(u<=E.length)
        if not near.any(): continue
        ok&=zz>=np.median(zz[near])-4.0
        _,(ir,ic)=distance_transform_edt(~ok,return_indices=True); L=(zz[ir,ic]-drop).astype(np.float32)
        def lev(x,y,L=L,c0=c0,r0=r0): return map_coordinates(L,[(SY1-np.asarray(y))/RES-0.5-r0,(np.asarray(x)-SX0)/RES-0.5-c0],order=1,mode='nearest')
        out.append((quad,lev))
        cx,cy=quad.centroid.coords[0]
        print('bridge end wall: %.0f m long, fill behind it %.0f m2 at %.1f m'%(E.length,quad.area,float(lev(np.array([cx]),np.array([cy]))[0])))
    return out
def build(H,zr,plan,SX0,SY1,RES,reach=16.0,span=34.0,mindz=1.2,maxdz=7.5,flat=6.0,drop=0.3):
    h,w=zr.shape; mask=np.isfinite(zr); tf=from_origin(SX0,SY1,RES,RES)
    AB=abutments(zr,plan,SX0,SY1,RES,drop=drop)
    DU=shapely.unary_union([d['poly'] for d in plan['decks']])
    deck=rasterize([(DU.buffer(1.0),1)],out_shape=(h,w),transform=tf,fill=0).astype(bool)
    d1,(ir,ic)=distance_transform_edt(~mask,return_indices=True); d1=d1*RES
    cand=np.where((~mask)&(~deck)&(d1<=reach)&(d1>0))
    r=cand[0].astype(float); c=cand[1].astype(float); nr=ir[cand]; nc=ic[cand]
    ur=r-nr; uc=c-nc; n=np.hypot(ur,uc)+1e-9; ur/=n; uc/=n
    z1=zr[nr,nc]; z2=np.full(len(r),np.nan); d2=np.full(len(r),np.nan); blocked=np.zeros(len(r),bool)
    for t in np.arange(1.0,span/RES+1.0,0.5):          # march on, away from the nearest road, to the road on the other side
        rr=np.rint(r+ur*t).astype(int); cc=np.rint(c+uc*t).astype(int); ok=(rr>=0)&(rr<h)&(cc>=0)&(cc<w)&np.isnan(z2)&~blocked
        k=np.where(ok)[0]
        blocked[k[deck[rr[k],cc[k]]]]=True              # not across a bridge
        k=k[mask[rr[k],cc[k]]&~deck[rr[k],cc[k]]]
        z2[k]=zr[rr[k],cc[k]]; d2[k]=t*RES
    d1c=d1[cand]; hit=np.isfinite(z2)&(d1c+d2<=span)&(np.abs(z1-z2)>=mindz)&(np.abs(z1-z2)<=maxdz)     # taller breaks stay graded slopes
    zhi=np.maximum(z1,z2); zlo=np.minimum(z1,z2); dhi=np.where(z1>=z2,d1c,d2)
    L=np.maximum(zlo-drop,zhi-drop-0.5*np.maximum(dhi-flat,0.0))
    hit&=L>H[cand]+0.3                                  # only where the height grid falls short of it
    C=np.zeros((h,w),bool); C[cand[0][hit],cand[1][hit]]=True
    C=binary_dilation(binary_opening(C,iterations=1),iterations=1)&~mask&~deck
    Lr=np.full((h,w),np.nan,np.float32); Lr[cand[0][hit],cand[1][hit]]=L[hit]
    _,(jr,jc)=distance_transform_edt(np.isnan(Lr),return_indices=True); Lf=Lr[jr,jc]     # nearest fill level everywhere
    def lev(x,y): return map_coordinates(Lf,[(SY1-np.asarray(y))/RES-0.5,(np.asarray(x)-SX0)/RES-0.5],order=1,mode='nearest')
    def ground(x,y): return map_coordinates(H,[(SY1-np.asarray(y))/RES-0.5,(np.asarray(x)-SX0)/RES-0.5],order=1,mode='nearest')
    e3=np.zeros((0,3)); empty={'v':e3,'t':e3.astype(int),'wv':e3,'wt':e3.astype(int)}
    polys=[sg.shape(g) for g,v in shapes(C.astype(np.uint8),mask=C,transform=tf) if v==1]
    if not polys and not AB: return empty
    road=shapely.unary_union([plan['G'].buffer(0.3),DU.buffer(0.5)])
    U=shapely.unary_union(polys).buffer(2.0).buffer(-1.0).difference(road).buffer(-0.6).buffer(0.6).difference(road) if polys else sg.Polygon()
    if AB:      # the fill behind a bridge's end walls has its own level (only the pieces it overlaps: see mesh.build)
        K=shapely.unary_union([q for q,_ in AB]); parts=[]
        for g in (U.geoms if hasattr(U,'geoms') else [U]):
            if not g.intersects(K): parts.append(g); continue
            r=g.difference(K); parts+=[q for q in getattr(r,'geoms',[r]) if q.geom_type=='Polygon']
        U=sg.MultiPolygon([q for q in parts if q.geom_type=='Polygon'])
    V=[];T=[];n=0; WV=[];WT=[];wn=0; area=0.0
    for g,lv in [(g,lev) for g in (U.geoms if hasattr(U,'geoms') else [U])]+AB:
        if g.geom_type!='Polygon' or g.area<25: continue
        g=sg.polygon.orient(g,1.0); res=M.tri_poly(g,25.0)
        if res is None: continue
        v,t=res; z=lv(v[:,0],v[:,1]); V.append(np.c_[v,z]); T.append(t+n); n+=len(v); area+=g.area
        wn=walls(g,lv,ground,WV,WT,wn)
    print('gap fill: %d pieces, %.0f m2, %d wall panels'%(len(V),area,len(WT)))
    if not V: return empty
    return {'v':np.vstack(V),'t':np.vstack(T),'wv':np.vstack(WV) if WV else e3,'wt':np.vstack(WT) if WT else e3.astype(int)}
