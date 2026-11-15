import pickle, json, sys, numpy as np, pyproj
W=sys.argv[1]
tr=pyproj.Transformer.from_crs(4326,32119,always_xy=True)
k=json.load(open(W+'/kml.json'))
def kpts(name):
    P=[p[:2] for l in k[name]['lines'] for p in l]
    a=np.array(P); x,y=tr.transform(a[:,0],a[:,1]); return np.c_[x,y]
K=np.unique(np.round(kpts('Exist Right of Way Line'),2),axis=0)
print('kml pts',len(K), K.min(0),K.max(0))
from scipy.spatial import cKDTree
tree=cKDTree(K)
res={}
for n in (1,2,3):
    s=pickle.load(open(f'{W}/pdf{n}.pkl','rb'))
    P=[]
    for d in s:
        if d['layer']=='Exist Right of Way Line' and d['type']=='s':
            for it in d['items']:
                if it[0]=='l': P+= [it[1],it[2]]
    P=np.unique(np.round(np.array(P),2),axis=0)
    # keep within map frame
    print(n,'pdf pts',len(P))
    rng=np.random.default_rng(0)
    Ps=P[rng.choice(len(P),min(700,len(P)),replace=False)]
    S=0.3048*100/72  # m per pt assuming 1in=100ft
    best=(0,)
    for th in np.radians(np.arange(0,360,0.5)):
        c,s_=np.cos(th),np.sin(th)
        # pdf y down -> flip
        X=S*(c*Ps[:,0]+s_*(-Ps[:,1])); Y=S*(-s_*Ps[:,0]+c*(-Ps[:,1]))
        dx=(K[:,0][None,:]-X[:,None]).ravel(); dy=(K[:,1][None,:]-Y[:,None]).ravel()
        bx=np.floor(dx/4).astype(np.int64); by=np.floor(dy/4).astype(np.int64)
        key=bx*1000003+by
        u,cn=np.unique(key,return_counts=True)
        i=cn.argmax()
        if cn[i]>best[0]:
            m=key==u[i]; best=(cn[i],th,dx[m].mean(),dy[m].mean())
    print(n,'best',best[0],np.degrees(best[1]),best[2],best[3])
    # refine with least squares similarity using nearest neighbours, iterating
    cnt,th,tx,ty=best
    def apply(P,prm):
        a,b,tx,ty=prm; x=P[:,0]; y=-P[:,1]
        return np.c_[a*x+b*y+tx, -b*x+a*y+ty]
    prm=np.array([S*np.cos(th),S*np.sin(th),tx,ty])
    for it in range(8):
        Q=apply(P,prm); dd,ii=tree.query(Q)
        m=dd<(3 if it<3 else 1.0)
        x=P[m,0]; y=-P[m,1]; T=K[ii[m]]
        A=np.zeros((2*m.sum(),4)); A[0::2,0]=x; A[0::2,1]=y; A[0::2,2]=1; A[1::2,0]=y; A[1::2,1]=-x; A[1::2,3]=1
        prm=np.linalg.lstsq(A,T.ravel(),rcond=None)[0]
    Q=apply(P,prm); dd,ii=tree.query(Q); m=dd<1.0
    print(n,'inliers',m.sum(),'of',len(P),'rms',np.sqrt((dd[m]**2).mean()),'scale m/pt',np.hypot(prm[0],prm[1]),'rot',np.degrees(np.arctan2(prm[1],prm[0])))
    res[n]=prm.tolist()
json.dump(res,open(W+'/georef.json','w'))
