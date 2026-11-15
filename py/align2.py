import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely.geometry as sg, shapely.ops as so, shapely
raw=pickle.load(open(W+'/align_raw.pkl','rb'))
def densify(l,step=5.0):
    ls=sg.LineString(l); n=max(1,int(np.ceil(ls.length/step)))
    return np.array([ls.interpolate(i/n,normalized=True).coords[0] for i in range(n+1)])
by={1:[],2:[],3:[]}
for n,l,ln in raw:
    if ln>0.5: by[n].append(sg.LineString(l))
base=shapely.unary_union(by[2]).buffer(1.2)
pieces=list(by[2])
for n in (1,3):
    for ls in by[n]:
        d=ls.difference(base)
        for g in (d.geoms if hasattr(d,'geoms') else [d]):
            if g.length>2.0: pieces.append(g)
pieces=[np.array(p.coords) for p in pieces]
print('pieces',len(pieces))
# chain
def enddir(p,end):
    if end==0: v=p[0]-p[min(len(p)-1,1)]
    else: v=p[-1]-p[-2]
    return v/ (np.hypot(*v)+1e-9)   # pointing outward
used=[False]*len(pieces); chains=[]
ends=[(i,e,(p[0] if e==0 else p[-1]),enddir(p,e)) for i,p in enumerate(pieces) for e in (0,1)]
def find(pt,dr,skip):
    best=None
    for i,e,q,d in ends:
        if used[i] or i in skip: continue
        dist=np.hypot(*(q-pt))
        if dist<3.0:
            # outward dirs should be opposite
            cosang=-(d@dr)
            if cosang>np.cos(np.radians(15)):
                if best is None or dist<best[0]: best=(dist,i,e)
    return best
order=sorted(range(len(pieces)),key=lambda i:-sg.LineString(pieces[i]).length)
for i in order:
    if used[i]: continue
    used[i]=True; ch=pieces[i].copy()
    for side in (1,0):
        while True:
            pt=ch[-1] if side==1 else ch[0]
            dr=(ch[-1]-ch[-2]) if side==1 else (ch[0]-ch[1]); dr=dr/(np.hypot(*dr)+1e-9)
            b=find(pt,dr,())
            if not b: break
            _,j,e=b; used[j]=True; p=pieces[j] if e==0 else pieces[j][::-1]
            if side==1: ch=np.vstack([ch,p])
            else: ch=np.vstack([p[::-1],ch])
    chains.append(ch)
chains=[c for c in chains if sg.LineString(c).length>25]
chains.sort(key=lambda c:-sg.LineString(c).length)
print(len(chains),[round(sg.LineString(c).length) for c in chains])
pickle.dump(chains,open(W+'/align_chains.pkl','wb'))
