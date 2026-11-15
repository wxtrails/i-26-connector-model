import sys, glob, numpy as np, rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from PIL import Image
W=os.environ.get('I26_WORK') or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'work')
def mosaic(x0,y0,x1,y1,res,files=None):
    w=int(round((x1-x0)/res)); h=int(round((y1-y0)/res))
    dst_t=from_origin(x0,y1,res,res)
    out=np.zeros((3,h,w),np.uint8)
    files=files or (sorted(glob.glob(W+'/maxar/*_2020.tif'))+sorted(f for f in glob.glob(W+'/maxar/*.tif') if '_2020' not in f))
    for f in files:
        with rasterio.open(f) as src:
            tmp=np.zeros((3,h,w),np.uint8)
            reproject(rasterio.band(src,[1,2,3]),tmp,dst_transform=dst_t,dst_crs='EPSG:32119',resampling=Resampling.average if res>0.6 else Resampling.bilinear,src_nodata=0,dst_nodata=0)
            m=tmp.sum(0)>0
            out[:,m]=tmp[:,m]
    return out
if __name__=='__main__':
    x0,y0,x1,y1,res=map(float,sys.argv[1:6])
    o=mosaic(x0,y0,x1,y1,res)
    Image.fromarray(np.moveaxis(o,0,-1)).save(sys.argv[6],quality=88)
    print(o.shape,(o.sum(0)==0).mean())
