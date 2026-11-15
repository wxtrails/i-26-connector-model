import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely, shapely.geometry as sg
LY=pickle.load(open(W+'/pdf_layers.pkl','rb'))
K=pickle.load(open(W+'/kmz_sp.pkl','rb'))
def arrow(p):
    p=np.asarray(p)
    if len(p)<4: return None
    g=sg.Polygon(p)
    if not g.is_valid: g=g.buffer(0)
    if g.is_empty or g.area<0.3 or g.area>30: return None
    c=np.array(g.centroid.coords[0])
    q=p-p.mean(0); u,s,vt=np.linalg.svd(q,full_matrices=False); ax=vt[0]
    t=(p-c)@ax; L=t.max()-t.min()
    if L<1.5 or L>12: return None
    mid=(t.max()+t.min())/2   # bbox centre relative to centroid
    # width profile: head end is wider
    w=np.abs((p-c)@np.array([-ax[1],ax[0]]))
    hi=w[t>mid+L*0.15].max(initial=0); lo=w[t<mid-L*0.15].max(initial=0)
    d=ax if hi>lo else -ax
    return (c[0],c[1],d[0],d[1],L,g.area)
A=[]
for n,p in LY[('Prop Lane Arrow CELL','f')]:
    a=arrow(p)
    if a: A.append(a+(n,))
print('pdf arrows',len(A))
B=[]
for g in K['A']:
    if g.geom_type=='Polygon':
        a=arrow(np.array(g.exterior.coords))
        if a: B.append(a+(0,))
print('kmz arrows',len(B))
pickle.dump({'pdf':np.array(A),'kmz':np.array(B)},open(W+'/arrows.pkl','wb'))
chains=pickle.load(open(W+'/align_chains.pkl','rb'))
A=np.array(A)
for i,c in enumerate(chains):
    ls=sg.LineString(c); res=[]
    for a in A:
        pt=sg.Point(a[0],a[1]); d=ls.distance(pt)
        if d<14:
            s=ls.project(pt); p1=np.array(ls.interpolate(max(0,s-1)).coords[0]); p2=np.array(ls.interpolate(min(ls.length,s+1)).coords[0]); t=(p2-p1); t/=np.hypot(*t)+1e-9
            side=np.sign(t[0]*(a[1]-p1[1])-t[1]*(a[0]-p1[0]))  # +1 left
            dot=t@a[2:4]
            if abs(dot)>0.8: res.append((round(s),int(side),int(np.sign(dot)),round(d,1)))
    res.sort()
    fw=sum(1 for r in res if r[2]>0); bw=len(res)-fw
    L_f=sum(1 for r in res if r[2]>0 and r[1]>0); R_f=sum(1 for r in res if r[2]>0 and r[1]<0); L_b=sum(1 for r in res if r[2]<0 and r[1]>0); R_b=sum(1 for r in res if r[2]<0 and r[1]<0)
    print(i,'len',round(ls.length),'a=(%d,%d) b=(%d,%d)'%(c[0,0],c[0,1],c[-1,0],c[-1,1]),'arrows fwd',fw,'(L%d R%d)'%(L_f,R_f),'back',bw,'(L%d R%d)'%(L_b,R_b))
