import sys, glob, os, math, numpy as np, rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin, from_bounds
from PIL import Image
W=os.environ.get('I26_WORK') or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'work')
R=6378137.0
def build(z,x0,y0,x1,y1,res,out):
    tiles=[tuple(map(int,os.path.basename(f)[:-4].split('_'))) for f in glob.glob(f'{W}/dem/{z}_*.png')]
    xs=sorted(set(t[1] for t in tiles)); ys=sorted(set(t[2] for t in tiles))
    mos=np.zeros((len(ys)*256,len(xs)*256),np.float32)
    for _,tx,ty in tiles:
        a=np.asarray(Image.open(f'{W}/dem/{z}_{tx}_{ty}.png').convert('RGB')).astype(np.float32)
        h=a[...,0]*256+a[...,1]+a[...,2]/256-32768
        mos[(ty-ys[0])*256:(ty-ys[0]+1)*256,(tx-xs[0])*256:(tx-xs[0]+1)*256]=h
    n=2**z; ext=2*math.pi*R
    west=xs[0]/n*ext-ext/2; north=ext/2-ys[0]/n*ext; px=ext/n/256
    src_t=from_origin(west,north,px,px)
    w=int(round((x1-x0)/res)); h=int(round((y1-y0)/res))
    dst=np.zeros((h,w),np.float32)
    reproject(mos,dst,src_transform=src_t,src_crs='EPSG:3857',dst_transform=from_origin(x0,y1,res,res),dst_crs='EPSG:32119',resampling=Resampling.cubic)
    np.save(out,dst); print(out,dst.shape,dst.min(),dst.max())
build(14,277500,203200,289500,215200,5,W+'/dem_sp5.npy')     # 2400x2400
build(11,245000,170000,325000,250000,100,W+'/dem_far.npy')
