"""Triangulate sheets, solve road elevations with a membrane model."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely, shapely.geometry as sg, triangle as tr
import scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.spatial import cKDTree
from scipy.ndimage import map_coordinates, gaussian_filter
DEM=np.load(W+'/dem_sp5.npy'); DX0,DY1,DRES=277500.0,215200.0,5.0
DEMs=gaussian_filter(DEM,2.0)
def terr(x,y,smooth=True):
    c=(np.asarray(x)-DX0)/DRES-0.5; r=(DY1-np.asarray(y))/DRES-0.5
    return map_coordinates(DEMs if smooth else DEM,[r,c],order=1,mode='nearest')
def tri_poly(poly,maxarea):
    """triangulate a shapely Polygon (with holes) -> verts, tris"""
    poly=poly.simplify(0.25)
    if poly.is_empty or poly.area<5: return None
    if poly.geom_type!='Polygon': return None
    V=[];S=[];H=[]
    def ring(c):
        c=np.asarray(c)[:-1]
        # drop near-duplicate consecutive points
        keep=[0]
        for i in range(1,len(c)):
            if np.hypot(*(c[i]-c[keep[-1]]))>0.05: keep.append(i)
        c=c[keep]
        if len(c)<3: return False
        n0=len(V)
        for p in c: V.append(p)
        for i in range(len(c)): S.append((n0+i,n0+(i+1)%len(c)))
        return True
    if not ring(poly.exterior.coords): return None
    for h in poly.interiors:
        hp=sg.Polygon(h)
        if hp.area<1: continue
        if ring(h.coords): H.append(hp.representative_point().coords[0])
    d={'vertices':np.array(V),'segments':np.array(S)}
    if H: d['holes']=np.array(H)
    t=tr.triangulate(d,'pq28a%.1f'%maxarea)
    return t['vertices'],t['triangles']
def build(plan,maxarea_g=45,maxarea_d=30):
    sheets=[]   # dict(kind, id, poly, verts, tris)
    G=plan['G']
    for g in (G.geoms if hasattr(G,'geoms') else [G]):
        if g.area<30: continue
        r=tri_poly(g,maxarea_g)
        if r: sheets.append({'kind':'G','poly':g,'v':r[0],'t':r[1]})
    for k,d in enumerate(plan['decks']):
        P=d['poly']
        for g in (P.geoms if hasattr(P,'geoms') else [P]):
            if g.area<30: continue
            r=tri_poly(g,maxarea_d)
            if r: sheets.append({'kind':'D','deck':k,'poly':g,'v':r[0],'t':r[1],'iso':d.get('iso',False),'abut':d['abut'],'meta':d})
    return sheets
def solve(sheets,L=160.0,Lex=45.0,clear=6.8,verbose=True,Gexist=None,bend=25.0):
    off=np.cumsum([0]+[len(s['v']) for s in sheets]); N=off[-1]
    XY=np.vstack([s['v'] for s in sheets])
    sid=np.concatenate([np.full(len(s['v']),i) for i,s in enumerate(sheets)])
    isD=np.array([sheets[i]['kind']=='D' for i in sid])
    zt=terr(XY[:,0],XY[:,1])
    # FEM assembly
    I=[];J=[];Vv=[]; mass=np.zeros(N)
    for i,s in enumerate(sheets):
        v=s['v']; t=s['t']+off[i]
        p=XY[t]   # (nt,3,2)
        e0=p[:,2]-p[:,1]; e1=p[:,0]-p[:,2]; e2=p[:,1]-p[:,0]
        A=0.5*np.abs(e1[:,0]*e2[:,1]-e1[:,1]*e2[:,0])+1e-9
        E=[e0,e1,e2]
        for a in range(3):
            for b in range(3):
                I.append(t[:,a]);J.append(t[:,b]);Vv.append((E[a]*E[b]).sum(1)/(4*A))
            np.add.at(mass,t[:,a],A/3)
    K=sp.coo_matrix((np.concatenate(Vv),(np.concatenate(I),np.concatenate(J))),shape=(N,N)).tocsr()
    Kb=(K.T@sp.diags(1.0/np.maximum(mass,1e-3))@K)*(bend**2)
    K=K+Kb
    kk=np.full(N,1.0/L**2)
    if Gexist is not None and not Gexist.is_empty:
        ex=shapely.contains(Gexist.buffer(0.5),shapely.points(XY)); kk[ex]=1.0/Lex**2
    spring=np.where(isD,0.0,kk)*mass      # at-grade pulled to terrain
    target=zt.copy()
    # ---- links
    gidx=np.where(~isD)[0]; gtree=cKDTree(XY[gidx])
    pen_i=[];pen_j=[];pen_w=[]
    didx_all=np.where(isD)[0]
    dist_abut=np.full(N,1e9)
    for i,s in enumerate(sheets):
        if s['kind']!='D': continue
        idx=np.arange(off[i],off[i+1])
        ab=s['abut']
        if ab.is_empty: continue
        pts=shapely.points(XY[idx]); inab=shapely.contains(ab,pts)
        dist_abut[idx]=shapely.distance(ab,pts)
        ia=idx[inab]
        if len(ia)==0: continue
        dd,jj=gtree.query(XY[ia])
        ok=dd<8.0
        pen_i+=list(ia[ok]); pen_j+=list(gidx[jj[ok]]); pen_w+=[30.0]*ok.sum()
        s['nlink']=int(ok.sum())
    # deck-deck proximity links (same iso class, different sheets, not overlapping)
    dsheets=[i for i,s in enumerate(sheets) if s['kind']=='D']
    for a_ in range(len(dsheets)):
        for b_ in range(a_+1,len(dsheets)):
            i,j=dsheets[a_],dsheets[b_]; A_,B_=sheets[i],sheets[j]
            if A_['iso']!=B_['iso']: continue
            if A_['poly'].distance(B_['poly'])>6: continue
            if A_['poly'].intersection(B_['poly']).area>3 and not A_['iso']: pass
            ta=cKDTree(XY[off[j]:off[j+1]]); dd,jj=ta.query(XY[off[i]:off[i+1]])
            ok=dd<7.0
            if ok.sum():
                pen_i+=list(np.arange(off[i],off[i+1])[ok]); pen_j+=list(off[j]+jj[ok]); pen_w+=[10.0]*ok.sum()
    pen_i=np.array(pen_i,int);pen_j=np.array(pen_j,int);pen_w=np.array(pen_w)
    def assemble(extra_t_idx,extra_t_val,extra_t_w,rel_i,rel_j,rel_off,rel_w):
        d=spring.copy(); rhs=spring*target
        np.add.at(d,extra_t_idx,extra_t_w); np.add.at(rhs,extra_t_idx,extra_t_w*extra_t_val)
        # pair penalties: w*(zi-zj-o)^2
        pi=np.concatenate([pen_i,rel_i]).astype(int); pj=np.concatenate([pen_j,rel_j]).astype(int); pw=np.concatenate([pen_w,rel_w]); po=np.concatenate([np.zeros(len(pen_i)),rel_off])
        Pm=sp.coo_matrix((np.concatenate([pw,pw,-pw,-pw]),(np.concatenate([pi,pj,pi,pj]),np.concatenate([pi,pj,pj,pi]))),shape=(N,N)).tocsr()
        np.add.at(rhs,pi,pw*po); np.add.at(rhs,pj,-pw*po)
        return (K+Pm+sp.diags(d)).tocsc(),rhs
    # stacked relations: deck node above G / above other deck
    upper=[];lower=[]
    gpolys=[(i,s) for i,s in enumerate(sheets) if s['kind']=='G']
    for i,s in enumerate(sheets):
        if s['kind']!='D': continue
        idx=np.arange(off[i],off[i+1]); pts=shapely.points(XY[idx])
        far=dist_abut[idx]>20.0
        for gi,g in gpolys:
            if not g['poly'].intersects(s['poly']): continue
            inside=shapely.contains(g['poly'],pts)&far
            ii=idx[inside]
            if len(ii)==0: continue
            t_=cKDTree(XY[off[gi]:off[gi+1]]); dd,jj=t_.query(XY[ii])
            upper+=list(ii); lower+=list(off[gi]+jj)
        if s['iso']:
            for j,o in enumerate(sheets):
                if o['kind']!='D' or o['iso'] or not o['poly'].intersects(s['poly']): continue
                inside=shapely.contains(o['poly'].buffer(1.0),pts)
                ii=idx[inside]
                if len(ii)==0: continue
                t_=cKDTree(XY[off[j]:off[j+1]]); dd,jj=t_.query(XY[ii])
                upper+=list(ii); lower+=list(off[j]+jj)
    upper=np.array(upper,int); lower=np.array(lower,int)
    # terrain clearance for deck nodes
    need=np.clip(0.12*(dist_abut-8.0),0,7.0); need[~isD]=0; need[dist_abut>1e8]=4.0
    act_rel=np.zeros(len(upper),bool); act_t=np.zeros(N,bool)
    z=None; e=np.zeros(0,int); ef=np.zeros(0)
    relw=6.0*np.maximum(mass[upper],1.0) if len(upper) else np.zeros(0)
    for it in range(25):
        ti=np.where(act_t)[0]
        ri=upper[act_rel]; rt=np.minimum(z[lower[act_rel]]+clear,zt[ri]+28.0) if z is not None and act_rel.any() else np.zeros(0)
        idx=np.concatenate([ti,ri]).astype(int); val=np.concatenate([zt[ti]+need[ti],rt]); wt=np.concatenate([4.0*np.maximum(mass[ti],1.0),relw[act_rel]])
        A_,rhs=assemble(idx,val,wt,e,e,ef,ef)
        zn=spl.spsolve(A_,rhs)
        dz=np.abs(zn-z).max() if z is not None else 99.0
        z=zn
        v1=(z[upper]-z[lower]<clear-0.05)&~act_rel if len(upper) else np.zeros(0,bool)
        v2=isD&(z<zt+need-0.05)&~act_t
        if verbose: print('iter',it,'new rel viol',int(v1.sum()),'new terr viol',int(v2.sum()),'dz %.3f'%dz,flush=True)
        if v1.sum()==0 and v2.sum()==0 and dz<0.03: break
        if len(upper): act_rel|=v1
        act_t|=v2
    return XY,z,sid,off,zt
if __name__=='__main__':
    plan=pickle.load(open(W+'/sheets_plan.pkl','rb'))
    import time; t0=time.time()
    sheets=build(plan); print('sheets',len(sheets),'verts',sum(len(s['v']) for s in sheets),'tris',sum(len(s['t']) for s in sheets),round(time.time()-t0,1),'s')
    XY,z,sid,off,zt=solve(sheets,Gexist=plan.get('Gexist'))
    print('z range',z.min(),z.max(),'max |z-zt| at grade',np.abs(z-zt)[np.array([sheets[i]['kind']=='G' for i in sid])].max())
    for i,s in enumerate(sheets):
        if s['kind']=='D':
            zz=z[off[i]:off[i+1]]; c=s['poly'].representative_point()
            print('deck',s['deck'],'iso' if s['iso'] else '','links',s.get('nlink',0),'z %.1f..%.1f'%(zz.min(),zz.max()),'above terr %.1f..%.1f'%((zz-zt[off[i]:off[i+1]]).min(),(zz-zt[off[i]:off[i+1]]).max()),(round(c.x),round(c.y)))
    pickle.dump({'sheets':[{k:v for k,v in s.items() if k not in('meta',)} for s in sheets],'XY':XY,'z':z,'sid':sid,'off':off,'zt':zt},open(W+'/mesh.pkl','wb'))
