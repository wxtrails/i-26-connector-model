import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
from scipy.spatial import cKDTree
k=json.load(open(W+'/kml.json'))
def kpts(names):
    P=[p[:2] for name in names for l in k[name]['lines'] for p in l]
    a=np.array(P); x,y=ll2sp(a[:,0],a[:,1]); return np.c_[x,y]
K=np.unique(np.round(kpts(['Exist Right of Way Line','Exist Right of Way Control of Access Line','Exist Property Line']),2),axis=0)
tree=cKDTree(K)
S=0.42325
g=json.load(open(W+'/georef.json'))
for n in (1,3):
    P=[]
    for d in load_pdf(n):
        if d['layer'] in('Exist Right of Way Line',) and d['type']=='s':
            for it in d['items']:
                if it[0]=='l': P+= [it[1],it[2]]
    P=np.unique(np.round(np.array(P),2),axis=0)
    x=P[:,0]; y=-P[:,1]
    cands=[]
    for th in np.radians(np.arange(-125,-55,0.25)):
        c,s_=S*np.cos(th),S*np.sin(th)
        X=c*x+s_*y; Y=-s_*x+c*y
        sel=slice(None,None,max(1,len(P)//1500))
        dx=(K[:,0][None,:]-X[sel,None]).ravel(); dy=(K[:,1][None,:]-Y[sel,None]).ravel()
        key=np.floor(dx/3).astype(np.int64)*1000003+np.floor(dy/3).astype(np.int64)
        u,cn=np.unique(key,return_counts=True)
        i=cn.argmax(); m=key==u[i]
        cands.append((cn[i],th,dx[m].mean(),dy[m].mean()))
    cands.sort(key=lambda c:-c[0])
    print(n,[(c[0],round(np.degrees(c[1]),2),round(c[2]),round(c[3])) for c in cands[:6]])
    cnt,th,tx,ty=cands[0]
    for it in range(15):
        c,s_=S*np.cos(th),S*np.sin(th)
        Q=np.c_[c*x+s_*y+tx,-s_*x+c*y+ty]
        dd,ii=tree.query(Q); m=dd<(3.0 if it<5 else 0.6)
        T=K[ii[m]]; A=np.c_[x[m],y[m]]
        Ac=A-A.mean(0); Tc=T-T.mean(0)
        num=(Ac[:,1]*Tc[:,0]-Ac[:,0]*Tc[:,1]).sum(); den=(Ac[:,0]*Tc[:,0]+Ac[:,1]*Tc[:,1]).sum()
        th=np.arctan2(num,den)
        c,s_=S*np.cos(th),S*np.sin(th)
        tx=T[:,0].mean()-(c*A[:,0].mean()+s_*A[:,1].mean()); ty=T[:,1].mean()-(-s_*A[:,0].mean()+c*A[:,1].mean())
    print(n,'inliers',m.sum(),'of',len(P),'rms',np.sqrt((dd[m]**2).mean()),'rot',np.degrees(th),'extent',A[:,0].min(),A[:,0].max(),A[:,1].min(),A[:,1].max())
    g[str(n)]=[S*np.cos(th),S*np.sin(th),tx,ty]
json.dump(g,open(W+'/georef.json','w'))
