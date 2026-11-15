import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely, shapely.geometry as sg
from shapely.validation import make_valid
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
LY=pickle.load(open(W+'/pdf_layers.pkl','rb'))
def polys(layer,minarea=1.0):
    out=[]
    for n,p in LY.get((layer,'f'),[]):
        if len(p)>=4:
            g=make_valid(sg.Polygon(p))
            if g.area>minarea: out.append((n,g))
    return out
S=polys('Prop PHM Concrete Structures Shape',20)
print(len(S))
# merge duplicates: union-find on high mutual overlap
par=list(range(len(S)))
def f(i):
    while par[i]!=i: par[i]=par[par[i]]; i=par[i]
    return i
tree=shapely.STRtree([g for n,g in S])
pairs=[]
for i,(n,g) in enumerate(S):
    for j in tree.query(g):
        if j<=i: continue
        h=S[j][1]; a=g.intersection(h).area
        if a<1: continue
        fr=a/min(g.area,h.area)
        if fr>0.7: par[f(i)]=f(j)
        else: pairs.append((i,j,round(a),round(fr,2)))
groups={}
for i in range(len(S)): groups.setdefault(f(i),[]).append(i)
M=[shapely.unary_union([S[i][1] for i in v]) for v in groups.values()]
M=[m for m in M if m.area>20]
print('merged',len(M),'areas',sorted([round(m.area) for m in M],reverse=True)[:40])
pickle.dump(M,open(W+'/struct_merged.pkl','wb'))
# bridge outlines
BL=[p for n,p in LY[('Prop Roadway Bridge','s')]]
def fig(name,bbox,size):
    fg,ax=plt.subplots(figsize=(size,size*(bbox[3]-bbox[1])/(bbox[2]-bbox[0])))
    cm=plt.get_cmap('tab20')
    for i,m in enumerate(M):
        for g in (m.geoms if hasattr(m,'geoms') else [m]):
            if g.geom_type!='Polygon': continue
            x,y=g.exterior.xy; ax.fill(x,y,color=cm(i%20),alpha=0.55,lw=0.3,ec='k')
        c=m.representative_point()
        if bbox[0]<c.x<bbox[2] and bbox[1]<c.y<bbox[3] and m.area>60: ax.text(c.x,c.y,str(i),fontsize=8,weight='bold')
    for b in BL: ax.plot(b[:,0],b[:,1],'r-',lw=0.4)
    ax.set_xlim(bbox[0],bbox[2]); ax.set_ylim(bbox[1],bbox[3]); ax.set_aspect(1); ax.grid(True,lw=0.3)
    ax.xaxis.set_major_locator(plt.MultipleLocator(100)); ax.yaxis.set_major_locator(plt.MultipleLocator(100)); ax.tick_params(labelsize=6); plt.setp(ax.get_xticklabels(),rotation=90)
    fg.savefig(W+'/img/'+name,dpi=75,bbox_inches='tight'); plt.close(fg)
fig('st_m.png',(284900,209600,286500,211500),22)
