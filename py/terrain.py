"""Conform terrain to at-grade roads; build two-level terrain tiles."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
from surf import Surf
import shapely, shapely.geometry as sg
from scipy.ndimage import distance_transform_edt, map_coordinates, zoom, gaussian_filter
from rasterio.features import rasterize
from rasterio.transform import from_origin
SX0,SY0,SX1,SY1=279480.0,204480.0,288600.0,214080.0     # scene bbox (multiples of 240 from origin)
RES=2.5
TILE=240.0; FINE=6.0; COARSE=24.0
def conform(surf):
    DEM=np.load(W+'/dem_sp5.npy'); DX0,DY1=277500.0,215200.0
    w=int((SX1-SX0)/RES); h=int((SY1-SY0)/RES)
    # upsample DEM onto raster
    cols=(SX0+(np.arange(w)+0.5)*RES-DX0)/5.0-0.5; rows=(DY1-(SY1-(np.arange(h)+0.5)*RES))/5.0-0.5
    rr,cc=np.meshgrid(rows,cols,indexing='ij')
    base=map_coordinates(DEM,[rr,cc],order=1,mode='nearest').astype(np.float32)
    zr=np.full((h,w),np.nan,np.float32)
    import matplotlib.tri as mtri
    for i,s in enumerate(surf.sheets):
        if s['kind']!='G': continue
        tri=np.vstack([s['t'],s.get('tw',np.zeros((0,3),int))]).astype(int)
        if len(tri)==0: continue
        x0,y0,x1,y1=s['poly'].bounds
        c0=max(0,int((x0-SX0)/RES)-1); c1=min(w,int((x1-SX0)/RES)+2); r0=max(0,int((SY1-y1)/RES)-1); r1=min(h,int((SY1-y0)/RES)+2)
        if c1<=c0 or r1<=r0: continue
        xs=SX0+(np.arange(c0,c1)+0.5)*RES; ys=SY1-(np.arange(r0,r1)+0.5)*RES
        X,Y=np.meshgrid(xs,ys)
        itp=mtri.LinearTriInterpolator(mtri.Triangulation(s['v'][:,0],s['v'][:,1],tri),surf.z[surf.off[i]:surf.off[i+1]])
        z=np.ma.filled(itp(X,Y),np.nan)
        sub=zr[r0:r1,c0:c1]; m=np.isfinite(z); sub[m]=z[m]
    mask=np.isfinite(zr)
    d,(ir,ic)=distance_transform_edt(~mask,return_indices=True); d=d*RES
    znear=zr[ir,ic]
    # grade the ground to the road: level beside the pavement, then 2:1 cut and fill slopes out to natural ground
    lim=np.maximum(d-6.0,0)*0.5
    out=np.clip(base,znear-0.45-lim,znear-0.45+lim).astype(np.float32)
    import caps as CP, pickle as _pk
    CAPS=_pk.load(open(W+'/mesh.pkl','rb')).get('caps',[])
    for c in CAPS:       # ground in the corridor under an overpass is cut down with the road
        r=c['wd']/2+300; x0,x1=c['cen'][0]-r,c['cen'][0]+r; y0,y1=c['cen'][1]-r,c['cen'][1]+r
        c0=max(0,int((x0-SX0)/RES)); c1=min(w,int((x1-SX0)/RES)+1); r0=max(0,int((SY1-y1)/RES)); r1=min(h,int((SY1-y0)/RES)+1)
        if c1<=c0 or r1<=r0: continue
        X,Y=np.meshgrid(SX0+(np.arange(c0,c1)+0.5)*RES,SY1-(np.arange(r0,r1)+0.5)*RES)
        out[r0:r1,c0:c1]=np.minimum(out[r0:r1,c0:c1],CP.field([c],X,Y)-0.45)
    out=np.where(mask,out,gaussian_filter(out,1.2)).astype(np.float32)
    for i,s in enumerate(surf.sheets):      # keep ground below bridge decks
        if s['kind']!='D': continue
        x0,y0,x1,y1=s['poly'].bounds
        c0=max(0,int((x0-SX0)/RES)-1); c1=min(w,int((x1-SX0)/RES)+2); r0=max(0,int((SY1-y1)/RES)-1); r1=min(h,int((SY1-y0)/RES)+2)
        if c1<=c0 or r1<=r0: continue
        xs=SX0+(np.arange(c0,c1)+0.5)*RES; ys=SY1-(np.arange(r0,r1)+0.5)*RES
        X,Y=np.meshgrid(xs,ys); z=surf.z_on(i,np.c_[X.ravel(),Y.ravel()]).reshape(X.shape)
        sub=out[r0:r1,c0:c1]; m=np.isfinite(z)&(sub>z-2.0); sub[m]=z[m]-2.0
    # deck mask for later use (e.g. shadows) not needed
    np.save(W+'/terrain_conf.npy',out); np.save(W+'/road_dist.npy',d.astype(np.float32))
    return out,d
def clamp_under_roads(surf,H,margin_g=0.3,margin_d=1.5):
    """Keep the rendered terrain below every road surface.
    The viewer draws terrain as triangles between grid points FINE m apart, so a grid point is limited by the lowest
    pavement anywhere within one grid step of it (minus a margin): then no terrain triangle can rise through a road."""
    import matplotlib.tri as mtri
    from scipy.ndimage import minimum_filter
    h,w=H.shape; RM=np.full((h,w),np.inf,np.float32)
    for i,s in enumerate(surf.sheets):
        t=s['t']
        if len(t)==0: continue
        x0,y0,x1,y1=s['poly'].bounds
        c0=max(0,int((x0-SX0)/RES)-1); c1=min(w,int((x1-SX0)/RES)+2); r0=max(0,int((SY1-y1)/RES)-1); r1=min(h,int((SY1-y0)/RES)+2)
        if c1<=c0 or r1<=r0: continue
        X,Y=np.meshgrid(SX0+(np.arange(c0,c1)+0.5)*RES,SY1-(np.arange(r0,r1)+0.5)*RES)
        zz=surf.z[surf.off[i]:surf.off[i+1]]
        z=np.ma.filled(mtri.LinearTriInterpolator(mtri.Triangulation(s['v'][:,0],s['v'][:,1],t),zz)(X,Y),np.nan)
        # vertices too: the raster cells alone can miss narrow slivers of pavement
        vc=((s['v'][:,0]-SX0)/RES).astype(int); vr=((SY1-s['v'][:,1])/RES).astype(int); ok=(vc>=0)&(vc<w)&(vr>=0)&(vr<h)
        m=margin_d if s['kind']=='D' else margin_g
        sub=RM[r0:r1,c0:c1]; f=np.isfinite(z); sub[f]=np.minimum(sub[f],z[f]-m)
        np.minimum.at(RM,(vr[ok],vc[ok]),(zz[ok]-m).astype(np.float32))
    k=2*int(np.ceil((FINE+RES)/RES))+1
    RM=minimum_filter(RM,size=k,mode='nearest')
    n=int((H>RM).sum()); H=np.minimum(H,RM).astype(np.float32)
    print('clamped terrain cells under roads',n)
    return H
def sample(H,x,y):
    c=(np.asarray(x)-SX0)/RES-0.5; r=(SY1-np.asarray(y))/RES-0.5
    return map_coordinates(H,[r,c],order=1,mode='nearest')
if __name__=='__main__':
    surf=Surf()
    H,d=conform(surf)
    print('conformed',H.shape,H.min(),H.max())
    H=clamp_under_roads(surf,H); np.save(W+'/terrain_conf.npy',H)
    nx=int((SX1-SX0)/TILE); ny=int((SY1-SY0)/TILE)
    # fine tiles: near roads
    plan=pickle.load(open(W+'/sheets_plan2.pkl','rb'))
    allroad=shapely.unary_union([plan['G']]+[dd['poly'] for dd in plan['decks']]).buffer(130).simplify(5)
    fine=np.zeros((ny,nx),bool)
    for j in range(ny):
        for i in range(nx):
            if allroad.intersects(sg.box(SX0+i*TILE,SY0+j*TILE,SX0+(i+1)*TILE,SY0+(j+1)*TILE)): fine[j,i]=True
    print('tiles',nx,ny,'fine',fine.sum())
    # coarse grid heights (whole scene) and fine tile heights
    gx=SX0+np.arange(int((SX1-SX0)/COARSE)+1)*COARSE; gy=SY0+np.arange(int((SY1-SY0)/COARSE)+1)*COARSE
    X,Y=np.meshgrid(gx,gy); Hc=sample(H,X,Y).astype(np.float32)   # [j,i] y ascending
    nf=int(TILE/FINE)
    tiles=[]
    fh=[]
    for j in range(ny):
        for i in range(nx):
            if not fine[j,i]: continue
            xs=SX0+i*TILE+np.arange(nf+1)*FINE; ys=SY0+j*TILE+np.arange(nf+1)*FINE
            X,Y=np.meshgrid(xs,ys); h=sample(H,X,Y).astype(np.float32)
            # edges adjoining coarse tiles: snap to coarse linear interpolation
            k=int(COARSE/FINE)
            def snap(line):
                idx=np.arange(nf+1); a=(idx//k)*k; b=np.minimum(a+k,nf); t=(idx-a)/k
                return line[a]*(1-t)+line[b]*t
            if j==0 or not fine[j-1,i]: h[0,:]=snap(h[0,:])
            if j==ny-1 or not fine[j+1,i]: h[-1,:]=snap(h[-1,:])
            if i==0 or not fine[j,i-1]: h[:,0]=snap(h[:,0])
            if i==nx-1 or not fine[j,i+1]: h[:,-1]=snap(h[:,-1])
            tiles.append((i,j)); fh.append(h)
    np.savez(W+'/terrain_tiles.npz',Hc=Hc,fine=fine,tiles=np.array(tiles),fh=np.array(fh),bbox=np.array([SX0,SY0,SX1,SY1]))
    print('coarse',Hc.shape,'fine tiles',len(tiles))
