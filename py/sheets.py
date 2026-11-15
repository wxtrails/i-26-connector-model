"""Build pavement 'sheets': at-grade polygons (G) and bridge decks (D) with abutment zones."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely, shapely.geometry as sg, shapely.ops as so
from shapely.validation import make_valid
LY=pickle.load(open(W+'/pdf_layers.pkl','rb'))
KM=pickle.load(open(W+'/kmz_sp.pkl','rb'))
chains=pickle.load(open(W+'/align_chains.pkl','rb'))
CUT_Y=208250.0
def polyonly(g):
    g=make_valid(g)
    if g.geom_type in('Polygon','MultiPolygon'): return g
    if hasattr(g,'geoms'):
        q=[x for x in g.geoms if x.geom_type in('Polygon','MultiPolygon')]
        return shapely.unary_union(q) if q else sg.Polygon()
    return sg.Polygon()
def pdfpolys(layer,minarea=1.0):
    out=[]
    for n,p in LY.get((layer,'f'),[]):
        if len(p)>=4:
            g=polyonly(sg.Polygon(p))
            if g.area>minarea: out.append(g)
    return out
southbox=sg.box(270000,200000,300000,CUT_Y)
# ---------- at grade
Gp=pdfpolys('Prop PHM Roadway Shape',2)+pdfpolys('Prop PHM Road to be Resurfaced Shape',2)
Gk=[polyonly(g.intersection(southbox)) for g in KM['R']]
Gk=[g for g in Gk if g.area>5]
G_plan=shapely.unary_union([g.buffer(0.3) for g in Gp+Gk]).buffer(-0.3)
print('G plan area',round(G_plan.area))
# ---------- decks (PDF): pieces -> groups
M=pickle.load(open(W+'/struct_merged.pkl','rb'))
P=[polyonly(m) for m in M]
P=[polyonly(p.buffer(-3).buffer(3.2).intersection(p)) for p in P]
P=[p for p in P if p.area>60]
chl=[sg.LineString(c) for c in chains]
EXCL={4,16,17,19,21,35,36,20,27}   # greenways / non-road alignments
def owners(p):
    best=(0,None)
    r=p.minimum_rotated_rectangle; xs,ys=r.exterior.xy
    e=[np.hypot(xs[i+1]-xs[i],ys[i+1]-ys[i]) for i in range(4)]; longside=max(e)
    for i,l in enumerate(chl):
        if i in EXCL: continue
        if l.intersects(p):
            ln=l.intersection(p).length
            if ln>best[0]: best=(ln,i)
    return {best[1]} if best[1] is not None and best[0]>=0.4*longside else set()
# KMZ deck pieces (south only)
Pk=[polyonly(g.intersection(southbox)) for g in KM['S']]
Pk=[p for p in Pk if p.area>40 and p.buffer(-1.5).buffer(1.5).area>40]
src=['pdf']*len(P)+['kmz']*len(Pk)
P=P+Pk
own=[owners(p) if s=='pdf' else set() for p,s in zip(P,src)]
def pdir(i,c):
    if own[i]:
        l=chl[list(own[i])[0]]; s_=l.project(c); a=l.interpolate(max(0,s_-3)); b=l.interpolate(min(l.length,s_+3))
        return np.degrees(np.arctan2(b.y-a.y,b.x-a.x))%180
    r=P[i].minimum_rotated_rectangle; xs,ys=r.exterior.xy
    e=[(np.hypot(xs[k+1]-xs[k],ys[k+1]-ys[k]),k) for k in range(4)]; k=max(e)[1]
    return np.degrees(np.arctan2(ys[k+1]-ys[k],xs[k+1]-xs[k]))%180
tree=shapely.STRtree(P)
par=list(range(len(P)))
def f(i):
    while par[i]!=i: par[i]=par[par[i]]; i=par[i]
    return i
ISO={7,2}       # owner chains whose decks cross over other decks: keep isolated
def rel(i,j):
    g,h=P[i],P[j]
    if g.distance(h)>3.0: return None
    ov=g.intersection(h).area
    if ov<3.0:
        return 'touch' if g.buffer(1.6).intersection(h.buffer(1.6)).area>12 else None
    if ov>0.6*min(g.area,h.area): return 'dup'
    c=g.intersection(h).centroid
    dang=abs((pdir(i,c)-pdir(j,c)+90)%180-90)
    return 'par' if dang<25 else 'cross'
pairs=[]
for i,g in enumerate(P):
    for j in tree.query(g.buffer(3.0)):
        if j>i:
            r=rel(i,j)
            if r: pairs.append((i,j,r))
# pass 1: isolated family via strict relations
iso=set(i for i in range(len(P)) if own[i]&ISO)
def longside(p):
    r=p.minimum_rotated_rectangle; xs,ys=r.exterior.xy
    return max(np.hypot(xs[k+1]-xs[k],ys[k+1]-ys[k]) for k in range(4))
for i in range(len(P)):
    if i in iso: continue
    for c7 in ISO:
        l=chl[c7]; c=P[i].centroid
        if l.distance(c)<7 and P[i].intersection(l).length>0.7*longside(P[i]):
            iso.add(i)
# raw pieces (unfiltered) for robust unions
RAW=[polyonly(m) for m in M]+[polyonly(g.intersection(southbox)) for g in KM['S']]
RAW=[r for r in RAW if r.area>20]
isoU=shapely.unary_union([P[i] for i in iso]) if iso else sg.Polygon()
def comps(g):
    return [q for q in (g.geoms if hasattr(g,'geoms') else [g]) if q.geom_type=='Polygon' and not q.is_empty]
decks=[]
# isolated (flyover) family
isoC=polyonly(isoU.buffer(2.5,join_style=2).buffer(-2.5,join_style=2))
for q in comps(isoC):
    if q.area>250: decks.append({'poly':q,'iso':True,'src':'pdf'})
# everything else: union -> close gaps -> remove thin slivers
rest=[r for r in RAW if not (r.intersection(isoU).area>0.5*r.area)]
U=shapely.unary_union([r.buffer(0.2) for r in rest]).buffer(-0.2)
U=polyonly(U.buffer(2.0,join_style=2).buffer(-2.0,join_style=2))
U=polyonly(U.buffer(-3.0).buffer(3.1).intersection(U))
for q in comps(U):
    if q.area>300:
        decks.append({'poly':q,'iso':False,'src':'kmz' if q.centroid.y<CUT_Y else 'pdf'})
for d in decks:
    p=d['poly']; d['width']=2*p.area/p.length
    d['own']=[i for i,l in enumerate(chl) if i not in EXCL and l.intersects(p) and l.intersection(p).length>30]
decks.sort(key=lambda d:-d['poly'].area)
# ---------- approach slabs -> abutment zones
AS=[sg.LineString(p) for n,p in LY.get(('Prop Approach Slab','s'),[]) if len(p)>=2]+[sg.LineString(p) for p in KM['AS'] if len(p)>=2]
AS=[a for a in AS if a.length>0.3]
asu=shapely.unary_union([a.buffer(1.0) for a in AS])
slabs=[polyonly(g.convex_hull) for g in (asu.geoms if hasattr(asu,'geoms') else [asu])]
slabs=[s for s in slabs if s.area>20 and s.centroid.distance(sg.Point(285381,209895))>25]
print('slabs',len(slabs))
stree=shapely.STRtree(slabs)
for k,d in enumerate(decks):
    z=[slabs[j] for j in stree.query(d['poly'].buffer(4)) if slabs[j].distance(d['poly'])<4]
    d['abut']=shapely.unary_union([s.buffer(5) for s in z]) if z else sg.Polygon()
    c=d['poly'].representative_point()
    print(k,d['src'],'ISO' if d['iso'] else '','area',round(d['poly'].area),'w',round(d['width'],1),'own',d['own'],'nslabs',len(z),'at',(round(c.x),round(c.y)))
for d in decks:
    if d['abut'].is_empty: continue
    cut=polyonly(G_plan.intersection(d['poly']).intersection(d['abut'].buffer(18)))
    if cut.area>1: G_plan=polyonly(G_plan.difference(cut.buffer(0.3)))
G_plan=polyonly(G_plan.buffer(-0.4).buffer(0.4))
decks=[d for d in decks if d['iso'] or not d['abut'].is_empty or (d['poly'].area>1500 and d['width']>8)]
Gex=shapely.unary_union([g.buffer(0.3) for g in pdfpolys('Prop PHM Road to be Resurfaced Shape',2)]).buffer(-0.3)
print('decks kept',len(decks))
pickle.dump({'G':G_plan,'decks':decks,'slabs':slabs,'Gexist':Gex},open(W+'/sheets_plan.pkl','wb'))
