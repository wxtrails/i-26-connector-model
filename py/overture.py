import sys, os, json, time
for k in ('AWS_SECRET_ACCESS_KEY','AWS_ACCESS_KEY_ID','AWS_SESSION_TOKEN'): os.environ.pop(k,None)
import pyarrow as pa, pyarrow.dataset as ds, pyarrow.compute as pc, pyarrow.fs as fs
import shapely, shapely.wkb
REL='2026-10-21.0'
bbox=(-82.665,35.53,-82.535,35.635)
proxy=os.environ.get('HTTPS_PROXY') or os.environ.get('https_proxy')
s3=fs.S3FileSystem(anonymous=True,region='us-west-2',proxy_options=proxy,connect_timeout=20,request_timeout=120)
def get(theme,typ,cols):
    t0=time.time()
    d=ds.dataset(f'overturemaps-us-west-2/release/{REL}/theme={theme}/type={typ}/',filesystem=s3,format='parquet')
    f=(pc.field('bbox','xmin')<bbox[2])&(pc.field('bbox','xmax')>bbox[0])&(pc.field('bbox','ymin')<bbox[3])&(pc.field('bbox','ymax')>bbox[1])
    cols=[c for c in cols if c in d.schema.names]
    t=d.to_table(filter=f,columns=cols)
    print(theme,typ,t.num_rows,'rows',round(time.time()-t0),'s',flush=True)
    return t
out=sys.argv[1]
which=sys.argv[2:]
spec={'place':('places','place',['id','geometry','names','categories','basic_category','confidence','brand','addresses']),
 'building':('buildings','building',['id','geometry','names','height','num_floors','class','subtype','roof_shape','roof_color','facade_color']),
 'segment':('transportation','segment',['id','geometry','names','class','subtype','subclass','road_flags','level_rules','routes','speed_limits','road_surface']),
 'water':('base','water',['id','geometry','names','class','subtype']),
 'land_use':('base','land_use',['id','geometry','names','class','subtype']),
 'land_cover':('base','land_cover',['id','geometry','subtype']),
 'infrastructure':('base','infrastructure',['id','geometry','names','class','subtype'])}
for w in which:
    th,ty,cols=spec[w]
    t=get(th,ty,cols)
    rows=t.to_pylist()
    feats=[]
    for r in rows:
        g=shapely.wkb.loads(r.pop('geometry'))
        feats.append({'type':'Feature','geometry':shapely.geometry.mapping(g),'properties':r})
    json.dump({'type':'FeatureCollection','features':feats},open(f'{out}/{w}.geojson','w'),default=str)
