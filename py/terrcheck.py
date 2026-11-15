"""How far does the rendered terrain rise above the pavement? Exact viewer interpolation (fine tiles, 5 cm quantisation)."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
T=np.load(W+'/terrain_tiles.npz'); X0,Y0,X1,Y1=T['bbox']; TS=240.0; FS=6.0
fh=np.floor(np.clip((T['fh']-500)*20,0,65535))/20+500; n1=fh.shape[1]
nx=int(round((X1-X0)/TS)); idx=-np.ones((int(round((Y1-Y0)/TS)),nx),int)
for k,(i,j) in enumerate(T['tiles']): idx[j,i]=k
def th(x,y):
    ti=np.floor((x-X0)/TS).astype(int); tj=np.floor((y-Y0)/TS).astype(int); k=idx[tj,ti]
    fx=(x-X0-ti*TS)/FS; fy=(y-Y0-tj*TS)/FS; i=np.minimum(n1-2,np.floor(fx).astype(int)); j=np.minimum(n1-2,np.floor(fy).astype(int)); u=fx-i; v=fy-j
    kk=np.maximum(k,0); a=fh[kk,j,i]; b=fh[kk,j,i+1]; c=fh[kk,j+1,i]; d=fh[kk,j+1,i+1]
    h=np.where(u>=v, a+u*(b-a)+v*(d-b), a+u*(d-c)+v*(c-a))
    return np.where(k>=0,h,np.nan)
if __name__=='__main__':
    m=pickle.load(open(W+'/mesh.pkl','rb')); bad=[]; tot=0
    for i,s in enumerate(m['sheets']):
        v=s['v']; z=m['z'][m['off'][i]:m['off'][i+1]]; t=s['t']
        if len(t)==0: continue
        P=np.vstack([v,(v[t[:,0]]+v[t[:,1]]+v[t[:,2]])/3]); Z=np.concatenate([z,(z[t[:,0]]+z[t[:,1]]+z[t[:,2]])/3])
        e=th(P[:,0],P[:,1])-Z; tot+=len(P)
        if s['kind']=='D': e=e+1.0       # decks: ground should stay a metre or more below
        ok=np.isfinite(e)&(e>0.0); bad.append(np.c_[P[ok],e[ok],np.full(ok.sum(),i)])
    B=np.vstack(bad)
    if len(B)==0: print('pavement points',tot,'terrain above road: none'); sys.exit(0)
    print('pavement points',tot,'terrain above road',len(B),'>0.5m',int((B[:,2]>0.5).sum()),'>2m',int((B[:,2]>2).sum()),'max',round(B[:,2].max(),1))
    np.save(W+'/terrbad.npy',B)
    o=np.argsort(-B[:,2])[:15]; print(np.round(B[o],1))
