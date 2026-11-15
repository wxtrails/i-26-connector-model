import json, sys, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
W=sys.argv[1]
d=json.load(open(W+'/kml.json'))
r=d['Prop PHM Roadway Shape']
import collections
print('poly sizes',collections.Counter(len(p[0]) for p in r['polys']).most_common(8))
print('rings>1',sum(1 for p in r['polys'] if len(p)>1))
print('line sizes',collections.Counter(len(l) for l in r['lines']).most_common(8))
# area stats
def area(ring):
    a=np.array(ring)[:,:2]; x=(a[:,0]+82.6)*90600; y=(a[:,1]-35.58)*110900
    return 0.5*abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
ar=np.array([area(p[0]) for p in r['polys']])
print('area pct',np.percentile(ar,[5,25,50,75,95,99,100]), ar.sum())
def draw(ax,name,color,lw=0.4,fill=None):
    L=d[name]
    if L['polys']:
        pc=PolyCollection([np.array(p[0])[:,:2] for p in L['polys']],facecolors=fill or color,edgecolors='none' if fill else color,linewidths=lw,alpha=0.8)
        ax.add_collection(pc)
    if L['lines']:
        ax.add_collection(LineCollection([np.array(l)[:,:2] for l in L['lines']],colors=color,linewidths=lw))
def fig(name,bbox,size=(20,20)):
    f,ax=plt.subplots(figsize=size)
    draw(ax,'Prop PHM Roadway Shape','#bbbbbb',fill='#cccccc')
    draw(ax,'Prop PHM Concrete Structures Shape','#ff8888',fill='#ffaaaa')
    draw(ax,'Prop EOT','k',0.5)
    draw(ax,'Prop Roadway Bridge','r',0.8)
    draw(ax,'Prop Lane Arrow CELL','b',0.3,fill='b')
    ax.set_xlim(bbox[0],bbox[2]); ax.set_ylim(bbox[1],bbox[3]); ax.set_aspect(1/np.cos(np.radians(35.58)))
    ax.grid(True,lw=0.3); ax.xaxis.set_major_locator(plt.MultipleLocator(0.005)); ax.yaxis.set_major_locator(plt.MultipleLocator(0.005))
    ax.tick_params(labelsize=7); plt.setp(ax.get_xticklabels(),rotation=90)
    f.savefig(W+'/img/'+name,dpi=90,bbox_inches='tight'); plt.close(f)
fig('all.png',(-82.645,35.545,-82.56,35.618),(18,18))
