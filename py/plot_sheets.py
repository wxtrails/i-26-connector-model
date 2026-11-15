import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
S=pickle.load(open(W+'/sheets_plan.pkl','rb'))
def fig(name,bbox,size):
    fg,ax=plt.subplots(figsize=(size,size*(bbox[3]-bbox[1])/(bbox[2]-bbox[0])))
    for g in S['G'].geoms:
        x,y=g.exterior.xy; ax.fill(x,y,color='#dddddd',lw=0.3,ec='#888888')
        for h in g.interiors: x,y=h.xy; ax.fill(x,y,color='w',lw=0.3,ec='#888888')
    cm=plt.get_cmap('tab20')
    for i,d in enumerate(S['decks']):
        m=d['poly']
        for q in (m.geoms if hasattr(m,'geoms') else [m]):
            x,y=q.exterior.xy; ax.fill(x,y,color=cm(i%20),alpha=0.55,lw=0.5,ec='k')
        c=m.representative_point()
        if bbox[0]<c.x<bbox[2] and bbox[1]<c.y<bbox[3]: ax.text(c.x,c.y,str(i),fontsize=10,weight='bold')
        a=d['abut']
        for q in (a.geoms if hasattr(a,'geoms') else [a]):
            if not q.is_empty: x,y=q.exterior.xy; ax.plot(x,y,'b-',lw=0.8)
    ax.set_xlim(bbox[0],bbox[2]); ax.set_ylim(bbox[1],bbox[3]); ax.set_aspect(1); ax.grid(True,lw=0.3)
    g_=100; ax.xaxis.set_major_locator(plt.MultipleLocator(g_)); ax.yaxis.set_major_locator(plt.MultipleLocator(g_)); ax.tick_params(labelsize=6); plt.setp(ax.get_xticklabels(),rotation=90)
    fg.savefig(W+'/img/'+name,dpi=75,bbox_inches='tight'); plt.close(fg)
a=sys.argv
fig(a[1],tuple(float(v) for v in a[2:6]),float(a[6]))
