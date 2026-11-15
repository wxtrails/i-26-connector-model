import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
from surf import Surf
from mesh import terr
import shapely.geometry as sg
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
surf=Surf(); chains=pickle.load(open(W+'/align_chains.pkl','rb'))
ids=[int(a) for a in sys.argv[2].split(',')]; offs=[float(a) for a in sys.argv[3].split(',')] if len(sys.argv)>3 else [0]*len(ids)
fg,axs=plt.subplots(len(ids),1,figsize=(16,3.2*len(ids)))
if len(ids)==1: axs=[axs]
for ax,ci,o in zip(axs,ids,offs):
    ls=sg.LineString(chains[ci]); ls=ls.parallel_offset(abs(o),'left' if o>0 else 'right') if o else ls
    if ls.geom_type!='LineString': ls=max(ls.geoms,key=lambda g:g.length)
    s=np.arange(0,ls.length,5.0); P=np.array([ls.interpolate(t).coords[0] for t in s])
    z,sh=surf.track(P); zt=terr(P[:,0],P[:,1],False)
    ax.plot(s,zt,'g-',lw=0.8,label='terrain'); ax.plot(s,z,'k-',lw=1.2,label='road')
    isd=np.array([surf.sheets[i]['kind']=='D' if i>=0 else False for i in sh])
    ax.plot(s[isd],z[isd],'r.',ms=3,label='bridge')
    g=np.gradient(z,s)*100; ax.set_title(f'chain {ci} off {o}: max grade {np.nanmax(np.abs(g)):.1f}%'); ax.grid(True); ax.legend(fontsize=7)
fg.tight_layout(); fg.savefig(W+'/img/'+sys.argv[1],dpi=70)
