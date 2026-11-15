import os, pickle, json, numpy as np, pyproj
W=os.environ.get('I26_WORK') or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'work')
_tr=pyproj.Transformer.from_crs(4326,32119,always_xy=True)
_inv=pyproj.Transformer.from_crs(32119,4326,always_xy=True)
def ll2sp(lon,lat): return _tr.transform(lon,lat)
def sp2ll(x,y): return _inv.transform(x,y)
def georef(): return {int(k):v for k,v in json.load(open(W+'/georef.json')).items()}
def pdf2sp(P,prm):
    P=np.asarray(P,float); a,b,tx,ty=prm; x=P[...,0]; y=-P[...,1]
    return np.stack([a*x+b*y+tx,-b*x+a*y+ty],-1)
def bez(p0,p1,p2,p3,n=8):
    t=np.linspace(0,1,n+1)[:,None]; p0,p1,p2,p3=map(np.array,(p0,p1,p2,p3))
    return (1-t)**3*p0+3*(1-t)**2*t*p1+3*(1-t)*t**2*p2+t**3*p3
def path_polys(d,n=8):
    """return list of polylines (np arrays, pdf coords) for a drawing; consecutive items chained"""
    out=[];cur=None
    for it in d['items']:
        if it[0]=='l': seg=np.array([it[1],it[2]])
        elif it[0]=='c': seg=bez(it[1],it[2],it[3],it[4],n)
        elif it[0]=='re': (x0,y0),(x1,y1)=it[1],it[2]; seg=np.array([(x0,y0),(x1,y0),(x1,y1),(x0,y1),(x0,y0)])
        elif it[0]=='qu': seg=np.array([it[1],it[2],it[3],it[4],it[1]])
        else: continue
        if cur is not None and np.hypot(*(cur[-1]-seg[0]))<1e-3: cur=np.vstack([cur,seg[1:]])
        else:
            if cur is not None: out.append(cur)
            cur=seg
    if cur is not None: out.append(cur)
    return out
def load_pdf(n): return pickle.load(open(f'{W}/pdf{n}.pkl','rb'))
def layer_lines(n,layer,types=('s','f','fs'),prm=None):
    prm=prm or georef()[n]
    res=[]
    for d in load_pdf(n):
        if d['layer']==layer and d['type'] in types:
            for pl in path_polys(d): res.append(pdf2sp(pl,prm))
    return res
