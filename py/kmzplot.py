import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely, shapely.geometry as sg
from shapely.validation import make_valid
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
k=json.load(open(W+'/kml.json'))
def kp(name,minarea):
    out=[]
    for rings in k[name]['polys']:
        a=np.array(rings[0]); x,y=ll2sp(a[:,0],a[:,1])
        if len(x)>=4:
            g=sg.Polygon(np.c_[x,y])
            if g.area>minarea: out.append(make_valid(g))
    return out
def kl(name):
    out=[]
    for l in k[name]['lines']:
        a=np.array(l); x,y=ll2sp(a[:,0],a[:,1]); out.append(np.c_[x,y])
    return out
R=kp('Prop PHM Roadway Shape',20); S=kp('Prop PHM Concrete Structures Shape',20); A=kp('Prop Lane Arrow CELL',0.5)
print(len(R),sum(g.area for g in R),len(S),sum(g.area for g in S),len(A))
print('roadway lines',len(k['Prop PHM Roadway Shape']['lines']),'struct lines',len(k['Prop PHM Concrete Structures Shape']['lines']))
pickle.dump({'R':R,'S':S,'A':A,'EOT':kl('Prop EOT'),'BR':kl('Prop Roadway Bridge'),'RL':kl('Prop PHM Roadway Shape'),'SL':kl('Prop PHM Concrete Structures Shape'),'LL':kl('Prop Lane Lines'),'AS':kl('Prop Approach Slab')},open(W+'/kmz_sp.pkl','wb'))
def fig(name,bbox,size):
    fg,ax=plt.subplots(figsize=(size,size*(bbox[3]-bbox[1])/(bbox[2]-bbox[0])))
    for g in R:
        for q in (g.geoms if hasattr(g,'geoms') else [g]):
            if q.geom_type=='Polygon': x,y=q.exterior.xy; ax.fill(x,y,color='#eeee88',lw=0.2,ec='#999900')
    for g in S:
        for q in (g.geoms if hasattr(g,'geoms') else [g]):
            if q.geom_type=='Polygon': x,y=q.exterior.xy; ax.fill(x,y,color='#ff9999',alpha=0.7,lw=0.2,ec='r')
    ax.add_collection(LineCollection(kl('Prop EOT'),colors='k',linewidths=0.4))
    ax.add_collection(LineCollection(kl('Prop Roadway Bridge'),colors='r',linewidths=0.6))
    for g in A:
        if g.geom_type=='Polygon': x,y=g.exterior.xy; ax.fill(x,y,color='b')
    ax.set_xlim(bbox[0],bbox[2]); ax.set_ylim(bbox[1],bbox[3]); ax.set_aspect(1); ax.grid(True,lw=0.3)
    ax.xaxis.set_major_locator(plt.MultipleLocator(200)); ax.yaxis.set_major_locator(plt.MultipleLocator(200)); ax.tick_params(labelsize=6); plt.setp(ax.get_xticklabels(),rotation=90)
    fg.savefig(W+'/img/'+name,dpi=75,bbox_inches='tight'); plt.close(fg)
fig('ka.png',(279300,205000,285200,208600),24)
