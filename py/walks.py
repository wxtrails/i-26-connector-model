"""Sidewalks and multi-use paths from the plans.

The plans colour new concrete red ('Prop PHM Concrete Structures Shape': bridges, sidewalks, curbs, median islands) and
the multi-use paths magenta ('Prop PHM Misc Impact Area Hatching'). A red shape is a sidewalk when it is a strip a
sidewalk's width (1.2-4.5 m) edged by the plan's sidewalk lines ('Prop Conc 4in Sidewalk') with pavement on at most one
side (a median island has road on both). Each walk is laid at the top of the curb beside its road (0.15 m above the
pavement edge, rising 2% away from it) and follows the smoothed ground where it leaves the road; a short skirt drops from
its edges so the coarser terrain grid never shows a gap beneath it. terrain.py grades the ground to the walks."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely, shapely.geometry as sg
from shapely.validation import make_valid
from scipy.spatial import cKDTree
import scipy.sparse as sp, scipy.sparse.linalg as spl
import mesh as M
from surf import Surf
CURB=0.15; SLOPE=0.02; SKIRT=0.8
SIDEWALK_W=(1.2,4.5)

def _parts(g):
    if g.is_empty: return []
    if g.geom_type=='Polygon': return [g]
    if hasattr(g,'geoms'): return [p for h in g.geoms for p in _parts(h)]
    return []

def plan_polys():
    """(red, magenta) filled shapes from all three plan sheets, in state plane"""
    G=georef(); red=[]; mag=[]
    for n in (1,2,3):
        for d in load_pdf(n):
            if d['type'] not in ('f','fs'): continue
            col=tuple(round(x,2) for x in (d.get('fill') or ()))
            if d['layer']=='Prop PHM Concrete Structures Shape' and col==(1.0,0.0,0.0): out=red
            elif d['layer']=='Prop PHM Misc Impact Area Hatching': out=mag
            else: continue
            for pl in path_polys(d):
                if len(pl)<4: continue
                out+=[h for h in _parts(make_valid(sg.Polygon(pdf2sp(pl,G[n])))) if h.area>0.5]
    return red,mag

def select(red,mag,G,decks):
    import marks as MK
    L=pickle.load(open(W+'/pdf_layers.pkl','rb'))
    SW=shapely.unary_union([sg.LineString(a) for a in MK.lines_of(L,('Prop Conc 4in Sidewalk','s'))]).buffer(0.6)
    shapely.prepare(SW); shapely.prepare(G)
    side=[]
    for h in red:
        wd=2*h.area/max(h.length,1e-9)
        if not SIDEWALK_W[0]<=wd<=SIDEWALK_W[1]: continue
        b=h.exterior
        if b.intersection(SW).length<0.4*b.length: continue
        ring=h.buffer(3.0).difference(h.buffer(0.3))
        if ring.intersection(G).area>0.6*ring.area: continue          # road on both sides: a median island
        side.append(h)
    close=lambda P: shapely.unary_union([p.buffer(0.4,join_style='mitre') for p in P]).buffer(-0.4,join_style='mitre')
    U=close(mag).difference(G)                    # the path pieces are drawn between its centre-line dashes
    S=close(side).difference(G).difference(U)
    # under a bridge the plans' red bridge fill is cut along the path's edge lines: those strips are the path
    DU=shapely.unary_union([d['poly'] for d in decks]); shapely.prepare(DU); Ub=U.buffer(1.0)
    sl=[p for p in _parts(S) if p.intersection(DU).area>0.5*p.area and p.intersects(Ub)]
    if sl: U=shapely.unary_union([U]+sl); S=S.difference(shapely.unary_union(sl).buffer(0.01))
    keep=lambda g,wmin: [p for p in _parts(g) if p.area>=8 and 2*p.area/p.length>=wmin]
    print('red strips under bridges joined to the path',len(sl))
    return keep(S,1.0),keep(U,1.5)

def road_edges(surf,kind,step=1.0):
    """points along the edges of the ground pavement (kind 'G') or the bridge decks ('D'), with height and direction"""
    P=[];Tg=[];Z=[]
    for i,s in enumerate(surf.sheets):
        if s['kind']!=kind: continue
        zz=surf.z[surf.off[i]:surf.off[i+1]]; vt=cKDTree(s['v'])
        for r in [s['poly'].exterior]+list(s['poly'].interiors):
            n=max(int(r.length/step),4); q=np.array([r.interpolate(t).coords[0] for t in np.linspace(0,r.length,n,endpoint=False)])
            t=np.roll(q,-1,0)-np.roll(q,1,0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9
            z=surf.z_on(i,q); bad=~np.isfinite(z)
            if bad.any(): _,j=vt.query(q[bad]); z[bad]=zz[j]
            P.append(q);Tg.append(t);Z.append(z)
    P,Tg,Z=np.vstack(P),np.vstack(Tg),np.concatenate(Z)
    return P,Tg,Z,cKDTree(P)

def lay(poly,EG,ED,maxarea=6.0,near=6.0,reach=20.0):
    """triangulate a walk and give it heights: the curb top beside a road running alongside it (a road it meets end on,
    or one it passes under, does not count), the smoothed ground where no road runs alongside"""
    r=M.tri_poly(poly,maxarea)
    if r is None: return None
    V,T=r; nv=len(V)
    # the walk's direction at each vertex: main axis of its outline nearby
    b=poly.exterior; n=max(int(b.length/0.7),8); B=np.array([b.interpolate(t).coords[0] for t in np.linspace(0,b.length,n,endpoint=False)])
    bt=cKDTree(B); u=np.zeros((nv,2))
    for k,p in enumerate(V):
        Q=B[bt.query_ball_point(p,7.0)]
        if len(Q)<3: u[k]=(1,0); continue
        Q=Q-Q.mean(0); w,e=np.linalg.eigh(Q.T@Q); u[k]=e[:,1]
    def cands(E):
        Ep,Et,Ez,Etree=E; out=[]
        for p,uu in zip(V,u):
            idx=np.array(Etree.query_ball_point(p,reach),int)
            if len(idx): idx=idx[np.abs(Et[idx]@uu)>0.7]
            out.append((np.hypot(*(Ep[idx]-p).T),Ez[idx]+CURB) if len(idx) else (np.zeros(0),np.zeros(0)))
        return out
    CG=cands(EG); CD=cands(ED)
    zg=M.terr(V[:,0],V[:,1])
    e=np.vstack([T[:,[0,1]],T[:,[1,2]],T[:,[2,0]]]); e=np.unique(np.sort(e,1),axis=0)
    A=sp.coo_matrix((np.ones(len(e)),(e[:,0],e[:,1])),shape=(nv,nv)); A=(A+A.T).tocsr()
    Lp=sp.diags(np.asarray(A.sum(1)).ravel())-A
    def fit(ref,tol,smooth,drop=np.inf):
        tgt=np.zeros(nv); wt=np.zeros(nv); free=np.zeros(nv,bool)
        for k in range(nv):
            d=np.r_[CG[k][0],CD[k][0]]; zc=np.r_[CG[k][1],CD[k][1]]
            if not len(d): continue
            sc=d+4.0*np.maximum(np.abs(zc-ref[k])-tol,0); j=np.argmin(sc)
            if abs(zc[j]-ref[k])>drop: free[k]=True; continue     # no road alongside at this level: follow the walk itself
            tgt[k]=zc[j]+SLOPE*min(d[j],4.0); wt[k]=1.0 if d[j]<=near else 0.25+0.75*(reach-d[j])/(reach-near)
        wg=np.where(free,0.0,0.03*(1-wt)+0.002)     # soft fit: curb targets strongly, the smoothed ground weakly,
        wg[free]=1e-6                                # neighbours tied together
        return spl.spsolve((sp.diags(wt+wg)+smooth*Lp).tocsc(),wt*tgt+wg*zg),tgt,wt
    # the nearest road alongside, unless it runs far above or below the ground here while another is near ground level
    # (a bridge overhead); then again against the walk's own overall level, so a sidewalk reaching a bridge end takes the
    # bridge, not the road passing under it
    z1,t1,w1=fit(zg,4.0,0.5)
    # the walk's own level nearby: median of those targets over 30 m, each 3 m cell counted once (the triangulation
    # is denser at the narrow ends, where a walk meets a bridge or a road passing under it)
    ref=z1.copy(); ok=w1>0
    if ok.any():
        cell={}
        for c,t in zip(map(tuple,np.floor(V[ok]/3.0).astype(int)),t1[ok]): cell.setdefault(c,[]).append(t)
        C=np.array(list(cell.keys()),float)*3.0+1.5; Cz=np.array([np.mean(v) for v in cell.values()]); ct=cKDTree(C)
        for k,nb in enumerate(ct.query_ball_point(V,30.0)):
            if nb: ref[k]=np.median(Cz[nb])
    z=fit(ref,1.0,0.5,drop=2.5)[0]
    # a lone vertex left far off its neighbours (a tip where a walk meets a bridge end) takes their median
    nb=[[] for _ in range(nv)]
    for a_,b_ in e: nb[a_].append(b_); nb[b_].append(a_)
    for _ in range(3):
        m=np.array([np.median(z[n]) if n else z[k] for k,n in enumerate(nb)]); bad=np.abs(z-m)>0.8; z[bad]=m[bad]
    return V,T,z

def skirt(V,T,z):
    """vertical faces down from the walk's outline"""
    e=np.vstack([T[:,[0,1]],T[:,[1,2]],T[:,[2,0]]]); key=np.sort(e,1)
    u,inv,cnt=np.unique(key,axis=0,return_inverse=True,return_counts=True); be=e[cnt[inv.ravel()]==1]
    P=[];I=[]
    for k,(a,b) in enumerate(be):
        P.append([[V[a,0],V[a,1],z[a]],[V[b,0],V[b,1],z[b]],[V[b,0],V[b,1],z[b]-SKIRT],[V[a,0],V[a,1],z[a]-SKIRT]])
        I.append(np.array([[0,1,2],[0,2,3]])+4*k)
    return (np.array(P).reshape(-1,3),np.vstack(I)) if P else (np.zeros((0,3)),np.zeros((0,3),int))

def build():
    surf=Surf(); plan=pickle.load(open(W+'/sheets_plan2.pkl','rb')); G=plan['G']
    red,mag=plan_polys(); S,U=select(red,mag,G,plan['decks'])
    EG=road_edges(surf,'G'); ED=road_edges(surf,'D')
    out={}
    for name,polys in (('walk',S),('mup',U)):
        pieces=[]
        for p in polys:
            r=lay(p,EG,ED)
            if r is None: continue
            V,T,z=r; sv,si=skirt(V,T,z); pieces.append({'poly':p,'v':V,'t':T,'z':z,'sv':sv,'st':si})
        out[name]=pieces
        print(name,len(pieces),'pieces',round(sum(p['poly'].area for p in pieces)),'m2',round(sum(p['poly'].length/2 for p in pieces)),'m long')
    pickle.dump(out,open(W+'/walks.pkl','wb'))

if __name__=='__main__':
    build()
