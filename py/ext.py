"""Overture road centrelines -> extension chains and route graph."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely, shapely.geometry as sg, shapely.ops as so
def load(classes=('motorway','trunk','primary','secondary','tertiary')):
    d=json.load(open(W+'/ov/segment.geojson')); out=[]
    for f in d['features']:
        p=f['properties']
        if p.get('subtype')!='road' or p.get('class') not in classes: continue
        g=f['geometry']; cs=[g['coordinates']] if g['type']=='LineString' else g['coordinates']
        n=p.get('names'); n=json.loads(n) if isinstance(n,str) else n
        for c in cs:
            c=np.array(c); x,y=ll2sp(c[:,0],c[:,1]); out.append({'xy':np.c_[x,y],'cls':p['class'],'name':(n or {}).get('primary'),'sub':p.get('subclass'),'flags':p.get('road_flags')})
    return out
if __name__=='__main__':
    R=load(); print(len(R))
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    plan=pickle.load(open(W+'/sheets_plan2.pkl','rb'))
    fg,ax=plt.subplots(figsize=(18,18))
    for g in plan['G'].geoms: ax.fill(*g.exterior.xy,color='#ccc')
    for d in plan['decks']:
        ax.fill(*d['poly'].exterior.xy,color='#fbb')
    col={'motorway':'b','trunk':'c','primary':'r','secondary':'g','tertiary':'y'}
    for r in R: ax.plot(r['xy'][:,0],r['xy'][:,1],col[r['cls']],lw=0.7)
    ax.set_xlim(279480,288600); ax.set_ylim(204480,214080); ax.set_aspect('equal'); ax.grid(True)
    ax.set_xticks(np.arange(279500,288600,500)); ax.set_yticks(np.arange(204500,214080,500)); ax.tick_params(labelsize=6,labelrotation=90)
    fg.savefig(W+'/img/ext.png',dpi=60,bbox_inches='tight')

def graph(R):
    import networkx as nx
    U=shapely.unary_union([sg.LineString(r['xy']) for r in R])
    G=nx.Graph()
    for g in U.geoms:
        a=tuple(np.round(g.coords[0],1)); b=tuple(np.round(g.coords[-1],1))
        if a==b: continue
        if G.has_edge(a,b) and G[a][b]['w']<=g.length: continue
        G.add_edge(a,b,w=g.length,geom=g)
    return G
def spath(G,a,b):
    import networkx as nx
    N=np.array(G.nodes); na=tuple(N[np.argmin(np.hypot(*(N-a).T))]); nb=tuple(N[np.argmin(np.hypot(*(N-b).T))])
    p=nx.shortest_path(G,na,nb,weight='w'); out=[]; edges=[]
    for u,v in zip(p[:-1],p[1:]):
        c=np.array(G[u][v]['geom'].coords)
        if np.hypot(*(c[0]-u))>np.hypot(*(c[-1]-u)): c=c[::-1]
        out.append(c if not out else c[1:]); edges.append((u,v))
    return np.vstack(out),edges
def pair(G,a,b):
    """two carriageways between a and b: returns (path a->b, path b->a) using drive-on-right."""
    a=np.array(a,float); b=np.array(b,float)
    p1,e1=spath(G,a,b); H=G.copy(); H.remove_edges_from(e1)
    try:
        p2,e2=spath(H,a,b)
        if sg.LineString(p2).length>1.25*sg.LineString(p1).length or sg.LineString(p2).hausdorff_distance(sg.LineString(p1))>80: raise Exception('far')
    except Exception as ex:
        return p1,p1[::-1].copy(),False
    # side of p2 relative to p1 at midpoint
    l1=sg.LineString(p1); m=l1.interpolate(0.5,normalized=True); t=np.array(l1.interpolate(0.5,normalized=True).coords[0])-np.array(l1.interpolate(0.49,normalized=True).coords[0])
    q=np.array(sg.LineString(p2).interpolate(sg.LineString(p2).project(m)).coords[0])-np.array(m.coords[0])
    left=t[0]*q[1]-t[1]*q[0]>0
    return (p1,p2[::-1].copy(),True) if left else (p2,p1[::-1].copy(),True)
EXT={ # name: (from, to, classes) ; "from" is the end joining the plan area
 'south':((284871,208141),(282950,206600),('motorway',)),
 'i240e':((286904,210504),(288450,211000),('motorway',)),
 'pate':((286881,210227),(287900,210450),('primary','motorway')),
 'patw':((284565,209694),(283050,209330),('primary','motorway')),
 'north':((285960,212678),(285480,214000),('motorway','trunk')),
}
def build_ext():
    R=load(); out={}
    for k,(a,b,cl) in EXT.items():
        G=graph([r for r in R if r['cls'] in cl])
        fwd,back,dual=pair(G,a,b)     # fwd: a->b (leaving plan area), back: b->a (entering)
        out[k]={'out':fwd,'in':back,'dual':dual}
        print(k,'dual',dual,'len %.0f %.0f'%(sg.LineString(fwd).length,sg.LineString(back).length),'start',fwd[0].round(),'end',fwd[-1].round(),'| in start',back[0].round(),'in end',back[-1].round())
    pickle.dump(out,open(W+'/ext.pkl','wb')); return out

def esnap(G,p):
    best=None
    for u,v,d in G.edges(data=True):
        dist=d['geom'].distance(sg.Point(*p))
        if best is None or dist<best[0]: best=(dist,u,v,d['geom'])
    return best
def sub(geom,s0,s1):
    r=so.substring(geom,min(s0,s1),max(s0,s1)); c=np.array(r.coords)
    return c if s0<=s1 else c[::-1]
def epath(G,a,b):
    import networkx as nx
    _,u1,v1,g1=esnap(G,a); _,u2,v2,g2=esnap(G,b)
    s1=g1.project(sg.Point(*a)); s2=g2.project(sg.Point(*b))
    if g1.equals(g2): return sub(g1,s1,s2)
    def endpt(g,n): return 0.0 if np.hypot(*(np.array(g.coords[0])-n))<0.2 else g.length
    best=None
    for n1 in (u1,v1):
        for n2 in (u2,v2):
            try: L=nx.shortest_path_length(G,n1,n2,weight='w')
            except Exception: continue
            tot=abs(endpt(g1,n1)-s1)+L+abs(endpt(g2,n2)-s2)
            if best is None or tot<best[0]: best=(tot,n1,n2)
    _,n1,n2=best; p=nx.shortest_path(G,n1,n2,weight='w'); out=[sub(g1,s1,endpt(g1,n1))]
    for u,v in zip(p[:-1],p[1:]):
        c=np.array(G[u][v]['geom'].coords)
        if np.hypot(*(c[0]-u))>np.hypot(*(c[-1]-u)): c=c[::-1]
        out.append(c[1:])
    out.append(sub(g2,endpt(g2,n2),s2)[1:])
    return np.vstack(out)
def build_ext2():
    R=load(); out={}
    for k,(a,b,cl) in EXT.items():
        G=graph([r for r in R if r['cls'] in cl]); a=np.array(a,float); b=np.array(b,float)
        t=(b-a)/np.hypot(*(b-a)); rn=np.array([t[1],-t[0]])
        fwd=epath(G,a+rn*9,b+rn*9); back=epath(G,b-rn*9,a-rn*9)
        lf,lb=sg.LineString(fwd).length,sg.LineString(back).length
        if lf>1.12*lb: fwd=back[::-1].copy()
        elif lb>1.12*lf: back=fwd[::-1].copy()
        dual=sg.LineString(fwd).hausdorff_distance(sg.LineString(back))>6
        out[k]={'out':fwd,'in':back,'dual':bool(dual)}
        print(k,'dual',dual,'len %.0f %.0f'%(sg.LineString(fwd).length,sg.LineString(back).length),'out',fwd[0].round(),fwd[-1].round(),'| in',back[0].round(),back[-1].round())
    pickle.dump(out,open(W+'/ext.pkl','wb')); return out
