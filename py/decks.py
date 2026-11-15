import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely, shapely.geometry as sg
from shapely.validation import make_valid
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
LY=pickle.load(open(W+'/pdf_layers.pkl','rb'))
M=pickle.load(open(W+'/struct_merged.pkl','rb'))
def clean(g):
    g=make_valid(g)
    if g.geom_type=='GeometryCollection': g=shapely.unary_union([q for q in g.geoms if q.geom_type in('Polygon','MultiPolygon')])
    return g
M=[clean(m) for m in M]
# drop slivers: opening
P=[]
for m in M:
    o=m.buffer(-1.5).buffer(1.5)
    if o.area>40: P.append(m)
print('pieces',len(P))
tree=shapely.STRtree(P)
par=list(range(len(P)))
def f(i):
    while par[i]!=i: par[i]=par[par[i]]; i=par[i]
    return i
for i,g in enumerate(P):
    for j in tree.query(g.buffer(0.4)):
        if j<=i: continue
        h=P[j]; ov=g.intersection(h).area
        if ov>0.25*min(g.area,h.area): continue     # stacked (crossing)
        # shared boundary length
        sh=g.buffer(0.4).intersection(h.buffer(0.4))
        if sh.area>0.8*4: # at least ~4m long contact
            par[f(i)]=f(j)
groups={}
for i in range(len(P)): groups.setdefault(f(i),[]).append(i)
D=[]
for v in groups.values():
    u=shapely.unary_union([P[i].buffer(0.25) for i in v]).buffer(-0.25)
    if u.area>250: D.append(u)
D.sort(key=lambda d:-d.area)
print('decks',len(D),[round(d.area) for d in D])
pickle.dump(D,open(W+'/decks_pdf.pkl','wb'))
AS=[p for n,p in LY.get(('Prop Approach Slab','s'),[])]
print('approach slab paths',len(AS))
def fig(name,bbox,size):
    fg,ax=plt.subplots(figsize=(size,size*(bbox[3]-bbox[1])/(bbox[2]-bbox[0])))
    cm=plt.get_cmap('tab20')
    for i,m in enumerate(D):
        for g in (m.geoms if hasattr(m,'geoms') else [m]):
            x,y=g.exterior.xy; ax.fill(x,y,color=cm(i%20),alpha=0.5,lw=0.4,ec='k')
        c=m.representative_point()
        if bbox[0]<c.x<bbox[2] and bbox[1]<c.y<bbox[3]: ax.text(c.x,c.y,str(i),fontsize=10,weight='bold')
    for b in AS: ax.plot(b[:,0],b[:,1],'b-',lw=0.8)
    for n,p in LY[('Prop EOT','s')]: ax.plot(p[:,0],p[:,1],'k-',lw=0.2)
    ax.set_xlim(bbox[0],bbox[2]); ax.set_ylim(bbox[1],bbox[3]); ax.set_aspect(1); ax.grid(True,lw=0.3)
    ax.xaxis.set_major_locator(plt.MultipleLocator(100)); ax.yaxis.set_major_locator(plt.MultipleLocator(100)); ax.tick_params(labelsize=6); plt.setp(ax.get_xticklabels(),rotation=90)
    fg.savefig(W+'/img/'+name,dpi=75,bbox_inches='tight'); plt.close(fg)
fig('dk.png',(284800,208400,286500,211600),20)
