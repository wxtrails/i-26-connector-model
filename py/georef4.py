import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
from scipy.spatial import cKDTree
import shapely.geometry as sg
b=json.load(open(W+'/ov/building.geojson'))
C=[]
for f in b['features']:
    g=sg.shape(f['geometry']); c=g.centroid; C.append((c.x,c.y))
C=np.array(C); x,y=ll2sp(C[:,0],C[:,1]); C=np.c_[x,y]
kk=json.load(open(W+'/kml.json'))
E=np.array([p[:2] for l in kk['Prop EOT']['lines'] for p in l[::3]]); ex,ey=ll2sp(E[:,0],E[:,1]); et=cKDTree(np.c_[ex,ey])
dd_,_=et.query(C); m=(dd_<350)&(C[:,1]>207500); C=C[m]
print('ov bldgs',len(C))
ctree=cKDTree(C)
k=json.load(open(W+'/kml.json'))
def kpts(names):
    P=[p[:2] for name in names for l in k[name]['lines'] for p in l]
    a=np.array(P); x,y=ll2sp(a[:,0],a[:,1]); return np.c_[x,y]
K=np.unique(np.round(kpts(['Exist Right of Way Line','Exist Right of Way Control of Access Line','Exist Property Line']),2),axis=0)
ktree=cKDTree(K)
S=0.42325
g=json.load(open(W+'/georef.json'))
def fit(A,T,th):
    Ac=A-A.mean(0); Tc=T-T.mean(0)
    num=(Ac[:,1]*Tc[:,0]-Ac[:,0]*Tc[:,1]).sum(); den=(Ac[:,0]*Tc[:,0]+Ac[:,1]*Tc[:,1]).sum()
    th=np.arctan2(num,den); c,s_=S*np.cos(th),S*np.sin(th)
    tx=T[:,0].mean()-(c*A[:,0].mean()+s_*A[:,1].mean()); ty=T[:,1].mean()-(-s_*A[:,0].mean()+c*A[:,1].mean())
    return th,tx,ty
def app(A,th,tx,ty):
    c,s_=S*np.cos(th),S*np.sin(th); return np.c_[c*A[:,0]+s_*A[:,1]+tx,-s_*A[:,0]+c*A[:,1]+ty]
for n in (1,2,3):
    B=[]
    for d in load_pdf(n):
        if d['layer']=='Exist PHM Building Shape' and d['fill'] is not None:
            for pl in path_polys(d):
                if len(pl)>=4:
                    try:
                        p=sg.Polygon(pl)
                        if p.area*S*S>30: B.append((p.centroid.x,-p.centroid.y))
                    except Exception: pass
    B=np.array(B); print(n,'pdf bldgs',len(B))
    cands=[]
    for th in np.radians(np.arange(-180,180,1.0)):
        Q=app(B,th,0,0)
        dx=(C[:,0][None,:]-Q[:,0][:,None]).ravel(); dy=(C[:,1][None,:]-Q[:,1][:,None]).ravel()
        ix=np.floor((dx-285500+12000)/6).astype(np.int64); iy=np.floor((dy-210500+12000)/6).astype(np.int64)
        ok=(ix>=0)&(ix<4000)&(iy>=0)&(iy<4000); key=ix[ok]*4000+iy[ok]; dx=dx[ok]; dy=dy[ok]
        cn=np.bincount(key,minlength=16000000); i=cn.argmax(); mm=key==i; cn={i:cn[i]}
        cands.append((cn[i],th,dx[mm].mean(),dy[mm].mean()))
    cands.sort(key=lambda c:-c[0])
    print(n,[(int(c[0]),round(float(np.degrees(c[1])),2),round(c[2]),round(c[3])) for c in cands[:4]])
    cnt,th,tx,ty=cands[0]
    for it in range(10):
        Q=app(B,th,tx,ty); dd,ii=ctree.query(Q); mm=dd<(8 if it<4 else 4)
        th,tx,ty=fit(B[mm],C[ii[mm]],th)
    print(n,'bldg inliers',mm.sum(),'rms',np.sqrt((dd[mm]**2).mean()),'rot',np.degrees(th),tx,ty)
    # fine tune with ROW points
    P=[]
    for d in load_pdf(n):
        if d['layer'] in('Exist Right of Way Line','Exist Property Line') and d['type']=='s':
            for it in d['items']:
                if it[0]=='l': P+= [it[1],it[2]]
    P=np.unique(np.round(np.array(P),2),axis=0); P[:,1]*=-1
    for it in range(12):
        Q=app(P,th,tx,ty); dd,ii=ktree.query(Q); mm=dd<(2.5 if it<4 else 0.8)
        if mm.sum()<30: print('too few'); break
        th,tx,ty=fit(P[mm],K[ii[mm]],th)
    print(n,'row inliers',mm.sum(),'rms',np.sqrt((dd[mm]**2).mean()),'rot',np.degrees(th),tx,ty,'extent',np.ptp(P[mm],axis=0))
    g[str(n)]=[S*np.cos(th),S*np.sin(th),tx,ty]
json.dump(g,open(W+'/georef.json','w'))
