import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import shapely.geometry as sg
k=json.load(open(W+'/kml.json'))
def klines(name):
    out=[]
    for l in k[name]['lines']:
        a=np.array(l); x,y=ll2sp(a[:,0],a[:,1]); out.append(np.c_[x,y])
    return out
seg=json.load(open(W+'/ov/segment.geojson'))
roads=[]
for f in seg['features']:
    if f['properties'].get('subtype')!='road': continue
    if f['geometry']['type']!='LineString': continue
    a=np.array(f['geometry']['coordinates']); x,y=ll2sp(a[:,0],a[:,1]); roads.append(np.c_[x,y])
def fig(name,bbox,size=18):
    f,ax=plt.subplots(figsize=(size,size*(bbox[3]-bbox[1])/(bbox[2]-bbox[0])))
    ax.add_collection(LineCollection(roads,colors='#88aaff',linewidths=0.6))
    ax.add_collection(LineCollection(klines('Prop EOT'),colors='#00aa00',linewidths=0.5))
    for n,c in ((1,'r'),(2,'m'),(3,'#ff8800')):
        ax.add_collection(LineCollection(layer_lines(n,'Prop EOT',('s',)),colors='k',linewidths=0.4))
        ax.add_collection(LineCollection(layer_lines(n,'Prop Horizontal Alignment',('s',)),colors=c,linewidths=0.9))
    ax.set_xlim(bbox[0],bbox[2]); ax.set_ylim(bbox[1],bbox[3]); ax.set_aspect(1); ax.grid(True,lw=0.3)
    ax.tick_params(labelsize=7); plt.setp(ax.get_xticklabels(),rotation=90)
    f.savefig(W+'/img/'+name,dpi=80,bbox_inches='tight'); plt.close(f)
args=sys.argv[1:]
fig(args[0],tuple(float(v) for v in args[1:5]),float(args[5]) if len(args)>5 else 18)
