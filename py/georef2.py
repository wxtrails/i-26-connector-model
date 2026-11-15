import pickle, json, sys, numpy as np, pyproj
from scipy.spatial import cKDTree
W=sys.argv[1]
tr=pyproj.Transformer.from_crs(4326,32119,always_xy=True)
k=json.load(open(W+'/kml.json'))
def kpts(names):
    P=[p[:2] for name in names for l in k[name]['lines'] for p in l]
    a=np.array(P); x,y=tr.transform(a[:,0],a[:,1]); return np.c_[x,y]
K=np.unique(np.round(kpts(['Exist Right of Way Line','Exist Right of Way Control of Access Line','Exist Property Line']),2),axis=0)
tree=cKDTree(K)
g=json.load(open(W+'/georef.json'))
S=0.42325
res={}
for n in (1,2,3):
    s=pickle.load(open(f'{W}/pdf{n}.pkl','rb'))
    P=[]
    for d in s:
        if d['layer'] in('Exist Right of Way Line',) and d['type']=='s':
            for it in d['items']:
                if it[0]=='l': P+= [it[1],it[2]]
    P=np.unique(np.round(np.array(P),2),axis=0)
    a,b,tx,ty=g[str(n)]
    th=np.arctan2(b,a)
    x=P[:,0]; y=-P[:,1]
    for it in range(15):
        c,s_=S*np.cos(th),S*np.sin(th)
        Q=np.c_[c*x+s_*y+tx,-s_*x+c*y+ty]
        dd,ii=tree.query(Q); m=dd<(2.0 if it<5 else 0.6)
        T=K[ii[m]]; A=np.c_[x[m],y[m]]
        # procrustes rotation with fixed scale
        Ac=A-A.mean(0); Tc=T-T.mean(0)
        # rotation R(th): [c s; -s c]
        num=(Ac[:,1]*Tc[:,0]-Ac[:,0]*Tc[:,1]).sum(); den=(Ac[:,0]*Tc[:,0]+Ac[:,1]*Tc[:,1]).sum()
        th=np.arctan2(num,den)
        c,s_=S*np.cos(th),S*np.sin(th)
        tx=T[:,0].mean()-(c*A[:,0].mean()+s_*A[:,1].mean()); ty=T[:,1].mean()-(-s_*A[:,0].mean()+c*A[:,1].mean())
    print(n,'inliers',m.sum(),'rms',np.sqrt((dd[m]**2).mean()),'rot',np.degrees(th), 'spread of inliers', A[:,0].min(),A[:,0].max(),A[:,1].min(),A[:,1].max())
    res[n]=[S*np.cos(th),S*np.sin(th),tx,ty]
json.dump(res,open(W+'/georef.json','w'))
