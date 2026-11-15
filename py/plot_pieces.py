import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
D=pickle.load(open(W+'/pieces.pkl','rb')); P=D['P']; own=D['own']
S=pickle.load(open(W+'/sheets_plan.pkl','rb'))
chains=pickle.load(open(W+'/align_chains.pkl','rb'))
AR=pickle.load(open(W+'/arrows.pkl','rb'))['pdf']
def fig(name,bbox,size):
    fg,ax=plt.subplots(figsize=(size,size*(bbox[3]-bbox[1])/(bbox[2]-bbox[0])))
    for g in S['G'].geoms:
        x,y=g.exterior.xy; ax.fill(x,y,color='#e4e4e4',lw=0.3,ec='#888888')
        for h in g.interiors: x,y=h.xy; ax.fill(x,y,color='w',lw=0.3,ec='#888888')
    cm=plt.get_cmap('tab20')
    for i,m in enumerate(P):
        for q in (m.geoms if hasattr(m,'geoms') else [m]):
            x,y=q.exterior.xy; ax.fill(x,y,color=cm(i%20),alpha=0.5,lw=0.5,ec='k')
        c=m.representative_point()
        if bbox[0]<c.x<bbox[2] and bbox[1]<c.y<bbox[3]: ax.text(c.x,c.y,f"{i}:{','.join(map(str,own[i])) or '-'}",fontsize=7,weight='bold')
    for i,c in enumerate(chains):
        ax.plot(c[:,0],c[:,1],color='k',lw=0.8,ls='--')
        inb=(c[:,0]>bbox[0])&(c[:,0]<bbox[2])&(c[:,1]>bbox[1])&(c[:,1]<bbox[3])
        idx=np.where(inb)[0]
        for t in (0.2,0.5,0.8):
            if len(idx): m=c[idx[int(len(idx)*t)]]; ax.text(m[0],m[1],str(i),fontsize=8,color='r')
    ax.quiver(AR[:,0],AR[:,1],AR[:,2],AR[:,3],color='b',scale_units='xy',scale=0.15,width=0.0015,pivot='mid',zorder=5)
    for a in S['slabs']: x,y=a.exterior.xy; ax.plot(x,y,'b-',lw=1)
    ax.set_xlim(bbox[0],bbox[2]); ax.set_ylim(bbox[1],bbox[3]); ax.set_aspect(1); ax.grid(True,lw=0.3)
    g_=50; ax.xaxis.set_major_locator(plt.MultipleLocator(g_)); ax.yaxis.set_major_locator(plt.MultipleLocator(g_)); ax.tick_params(labelsize=6); plt.setp(ax.get_xticklabels(),rotation=90)
    fg.savefig(W+'/img/'+name,dpi=75,bbox_inches='tight'); plt.close(fg)
a=sys.argv
fig(a[1],tuple(float(v) for v in a[2:6]),float(a[6]))
