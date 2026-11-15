import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
chains=pickle.load(open(W+'/'+(sys.argv[7] if len(sys.argv)>7 else 'align_chains.pkl'),'rb'))
cache=W+'/pdf_layers.pkl'
import os
if os.path.exists(cache): LY=pickle.load(open(cache,'rb'))
else:
    LY={}
    for n in (1,2,3):
        for d in load_pdf(n):
            if d['layer'] in ('Prop EOT','Prop Roadway Bridge','Prop PHM Concrete Structures Shape','Prop PHM Roadway Shape','Prop PHM Road to be Resurfaced Shape','Prop Lane Lines','Prop Conc Barrier','Prop Lane Arrow CELL','Prop Retaining Wall','Prop Paved Shoulder','Exist PHM Roadway Shape','Prop Paint Stripe','Prop Conc 2ft 6in Curb and Gutter','Prop Conc 4in Sidewalk','Prop Noise Wall','Prop Guardrail','Exist PHM Water Shape','Prop PHM Preferred Alternate 9 Shape','Prop Approach Slab','Prop Conc Expressway Gutter','Prop PHM Black Shape','Prop PHM White Shape','Prop Conc 5in Monolithic Keyed In Island','Prop Conc 5in Monolithic Surface Mount Island','Prop Conc 1ft 6in Curb and Gutter'):
                for pl in path_polys(d):
                    LY.setdefault((d['layer'],d['type'][0]),[]).append((n,pdf2sp(pl,georef()[n])))
    pickle.dump(LY,open(cache,'wb'))
def L(layer,t='s'): return [p for n,p in LY.get((layer,t),[])]
def fig(name,bbox,size=18,labels=True):
    f,ax=plt.subplots(figsize=(size,size*(bbox[3]-bbox[1])/(bbox[2]-bbox[0])))
    ax.add_collection(PolyCollection(L('Prop PHM Roadway Shape','f'),facecolors='#ffffaa',edgecolors='none'))
    ax.add_collection(PolyCollection(L('Prop PHM Road to be Resurfaced Shape','f'),facecolors='#ffd9a0',edgecolors='none'))
    ax.add_collection(PolyCollection(L('Prop PHM Concrete Structures Shape','f'),facecolors='#ffb0b0',edgecolors='none'))
    ax.add_collection(PolyCollection(L('Exist PHM Water Shape','f'),facecolors='#b0e0ff',edgecolors='none'))
    ax.add_collection(LineCollection(L('Prop EOT'),colors='k',linewidths=0.5))
    ax.add_collection(LineCollection(L('Prop Lane Lines'),colors='#888888',linewidths=0.3))
    ax.add_collection(LineCollection(L('Prop Roadway Bridge'),colors='r',linewidths=0.7))
    AR=pickle.load(open(W+'/arrows.pkl','rb'))['pdf']
    ax.quiver(AR[:,0],AR[:,1],AR[:,2],AR[:,3],color='b',scale_units='xy',scale=0.12,width=0.0018,headwidth=4,pivot='mid',zorder=5)
    ax.add_collection(LineCollection(L('Prop Conc Barrier'),colors='#aa00aa',linewidths=0.5))
    cm=plt.get_cmap('tab20')
    for i,c in enumerate(chains):
        col=cm(i%20)
        ax.plot(c[:,0],c[:,1],color=col,lw=1.6)
        if labels:
            inb=(c[:,0]>bbox[0])&(c[:,0]<bbox[2])&(c[:,1]>bbox[1])&(c[:,1]<bbox[3])
            if inb.any():
                idx=np.where(inb)[0]; m=c[idx[len(idx)//2]]
                ax.text(m[0],m[1],str(i),fontsize=11,color='k',weight='bold',bbox=dict(fc=col,alpha=0.8,pad=1,lw=0))
                for e,t in ((0,'a'),(-1,'b')):
                    if inb[e]: ax.text(c[e,0],c[e,1],f'{i}{t}',fontsize=7,color='#0000aa')
    ax.set_xlim(bbox[0],bbox[2]); ax.set_ylim(bbox[1],bbox[3]); ax.set_aspect(1); ax.grid(True,lw=0.3)
    g_=50 if (bbox[2]-bbox[0])<900 else 100; ax.xaxis.set_major_locator(plt.MultipleLocator(g_)); ax.yaxis.set_major_locator(plt.MultipleLocator(g_))
    ax.tick_params(labelsize=6); plt.setp(ax.get_xticklabels(),rotation=90)
    f.savefig(W+'/img/'+name,dpi=float(sys.argv[8]) if len(sys.argv)>8 else 75,bbox_inches='tight'); plt.close(f)
a=sys.argv
fig(a[1],tuple(float(v) for v in a[2:6]),float(a[6]))
