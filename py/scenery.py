"""Buildings, storefront signs, trees, far terrain, synthesized markings."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
import shapely, shapely.geometry as sg
from scipy.spatial import cKDTree
from PIL import Image
OX,OY=285000.0,210000.0
SX0,SY0,SX1,SY1=279480.0,204480.0,288600.0,214080.0
rng=np.random.default_rng(7)
def polys_of(f):
    g=f['geometry']; cs=[g['coordinates']] if g['type']=='Polygon' else g['coordinates']
    out=[]
    for c in cs:
        a=np.array(c[0]); x,y=ll2sp(a[:,0],a[:,1]); out.append(np.c_[x,y])
    return out
def run():
    R=pickle.load(open(W+'/routes.pkl','rb')); plan=pickle.load(open(W+'/sheets_plan2.pkl','rb'))
    RL=shapely.unary_union([sg.LineString(r['P']) for r in R]); shapely.prepare(RL)
    PAT=shapely.unary_union([sg.LineString(r['P']) for r in R if r['id'].startswith('pat')])
    U=shapely.unary_union([plan['G']]+[d['poly'] for d in plan['decks']]); Ub=U.buffer(1.5); shapely.prepare(Ub)
    res={}
    # ---------- buildings
    B=json.load(open(W+'/ov/building.geojson'))['features']
    xy=[];cnt=[];hh=[];typ=[];polys=[]
    for f in B:
        p=f['properties']
        for a in polys_of(f):
            if len(a)<4: continue
            c=a.mean(0)
            if not (SX0+50<c[0]<SX1-50 and SY0+50<c[1]<SY1-50): continue
            pg=sg.Polygon(a)
            if not pg.is_valid: pg=pg.buffer(0)
            if pg.is_empty or pg.geom_type!='Polygon' or pg.area<25: continue
            d=RL.distance(pg)
            if d>520: continue
            if Ub.intersects(pg): continue      # removed for the project / under new pavement
            pg=pg.simplify(0.4); pg=sg.polygon.orient(pg,1.0)
            a=np.array(pg.exterior.coords)[:-1]
            if len(a)<3 or len(a)>60: continue
            h=p.get('height'); nf=p.get('num_floors'); area=pg.area; cl=p.get('class')
            if h is None: h=3.2*nf+1 if nf else (4.5 if area<260 else 6.5 if area<1500 else 8.5)
            h=float(np.clip(h,3.0,80.0))
            if cl in('house','detached','static_caravan','residential') or (cl is None and area<230 and h<9): t=0
            elif h>16: t=2
            elif area>1800 and h<11: t=3
            else: t=1
            xy.append(a); cnt.append(len(a)); hh.append(h); typ.append(t); polys.append(pg)
    print('buildings',len(cnt),'types',np.bincount(typ))
    res['bld_xy']=np.vstack([np.c_[a[:,0]-OX,-(a[:,1]-OY)] for a in xy]).astype(np.float32); res['bld_n']=np.array(cnt,np.uint16); res['bld_h']=np.array(hh,np.float32); res['bld_t']=np.array(typ,np.uint8)
    # ---------- storefront signs from places
    P=json.load(open(W+'/ov/place.geojson'))['features']
    tree=shapely.STRtree(polys); signs=[]; used={}
    cand=[]
    for f in P:
        p=f['properties']; n=p.get('names'); n=json.loads(n) if isinstance(n,str) else n
        name=(n or {}).get('primary')
        if not name or len(name)>34 or (p.get('confidence') or 0)<0.6: continue
        x,y=ll2sp(*f['geometry']['coordinates']); pt=sg.Point(x,y); d=RL.distance(pt)
        if d>110: continue
        cand.append((d-25*(PAT.distance(pt)<60),name,p.get('basic_category') or '',x,y))
    cand.sort(key=lambda t:t[0]); seen=set()
    for d,name,cat,x,y in cand:
        if name.lower() in seen or len(signs)>=256: continue
        pt=sg.Point(x,y); k=tree.nearest(pt)
        if polys[k].distance(pt)>18: continue
        if typ[k]==0 and polys[k].area<140: continue
        a=np.array(polys[k].exterior.coords); best=None
        for i in range(len(a)-1):
            p0,p1=a[i],a[i+1]; Lg=np.hypot(*(p1-p0))
            if Lg<5.5: continue
            mid=(p0+p1)/2; e=(p1-p0)/Lg; nrm=np.array([e[1],-e[0]])     # outward for CCW ring
            q=np.array(RL.interpolate(RL.project(sg.Point(*mid))).coords[0]); to=q-mid; dist=np.hypot(*to)
            face=(to@nrm)/(dist+1e-9)
            if face<0.3: continue
            sc=dist-12*face-0.02*Lg
            if best is None or sc<best[0]: best=(sc,i,p0,e,nrm,Lg)
        if best is None: continue
        _,i,p0,e,nrm,Lg=best; slots=max(1,int(Lg//8.0)); u=used.get((k,i),0)
        if u>=slots: continue
        used[(k,i)]=u+1; seen.add(name.lower())
        w=min(7.0,Lg/slots-1.0); c=p0+e*(Lg*(u+0.5)/slots)
        signs.append({'x':round(float(c[0]-OX),2),'z':round(float(-(c[1]-OY)),2),'ex':round(float(e[0]),3),'ez':round(float(-e[1]),3),'nx':round(float(nrm[0]),3),'nz':round(float(-nrm[1]),3),'w':round(float(w),2),'t':name,'c':cat,'b':int(k)})
    print('signs',len(signs)); res['signs']=signs
    # ---------- trees from imagery (dark, textured vegetation) near routes
    im=np.asarray(Image.open(W+'/scene/base.jpg')).astype(np.float32)/255; H_,W_=im.shape[:2]; px=(SX1-SX0)/W_
    lum=im.mean(2); r,g,b=im[...,0],im[...,1],im[...,2]
    veg=(lum<0.40)&(lum>0.10)&(g>=b*0.98)&(r<0.5)
    from rasterio.features import rasterize
    from rasterio.transform import from_origin
    T=from_origin(SX0,SY1,px,px)
    ov=json.load(open(W+'/ov/segment.geojson'))['features']; rl=[]
    for f in ov:
        if f['properties'].get('subtype')!='road' or f['properties'].get('class') in('footway','path','steps','cycleway','bridleway','track'): continue
        gm=f['geometry']
        for c in ([gm['coordinates']] if gm['type']=='LineString' else gm['coordinates']):
            c=np.array(c); x,y=ll2sp(c[:,0],c[:,1]); rl.append(sg.LineString(np.c_[x,y]).buffer(7.0))
    wat=[]
    for f in json.load(open(W+'/ov/water.geojson'))['features']:
        if f['geometry']['type'] in('Polygon','MultiPolygon'):
            for a in polys_of(f): 
                if len(a)>3: wat.append(sg.Polygon(a).buffer(4))
    excl=rasterize([(g,1) for g in [U.buffer(9)]+[p.buffer(4) for p in polys]+rl+wat],out_shape=(H_,W_),transform=T,fill=0,dtype='uint8').astype(bool)
    near=rasterize([(RL.buffer(380),1)],out_shape=(H_,W_),transform=T,fill=0,dtype='uint8').astype(bool)
    ok=veg&~excl&near; ii=np.argwhere(ok); print('veg px',len(ii))
    sel=ii[rng.choice(len(ii),size=min(len(ii),9000),replace=False)]
    tx=SX0+(sel[:,1]+rng.random(len(sel)))*px; ty=SY1-(sel[:,0]+rng.random(len(sel)))*px
    d=shapely.distance(RL,shapely.points(np.c_[tx,ty])); keep=rng.random(len(sel))<np.clip(1.25-d/380,0.25,1)
    tx,ty=tx[keep],ty[keep]; th=rng.uniform(8,19,len(tx)); ev=(rng.random(len(tx))<0.22).astype(np.float32)
    res['trees']=np.c_[tx-OX,-(ty-OY),th,ev].astype(np.float32); print('trees',len(tx))
    # ---------- far terrain
    far=np.load(W+'/dem_far.npy')[::2,::2]; res['far']=np.clip((far-200)*10,0,65535).astype(np.uint16); res['far_meta']={'x0':245000.0,'y1':250000.0,'res':200.0,'w':int(far.shape[1]),'h':int(far.shape[0])}
    # ---------- synthesized lane markings on extension legs
    from marks import Strips
    wh=Strips(); ye=Strips(); done=set()
    def emit(st,P3,dashed,w=0.15):
        if len(P3)<2: return
        if not dashed: st.add(P3[:,:2],P3[:,2],w); return
        s=np.r_[0,np.cumsum(np.hypot(*np.diff(P3[:,:2],axis=0).T))]
        for s0 in np.arange(0,s[-1]-3,12.0):
            ss=np.array([s0,s0+3.0]); q=np.c_[np.interp(ss,s,P3[:,0]),np.interp(ss,s,P3[:,1]),np.interp(ss,s,P3[:,2])]; st.add(q[:,:2],q[:,2],w)
    for r in R:
        Q=r['P']; t=np.gradient(Q,axis=0); t/=np.hypot(t[:,0],t[:,1])[:,None]+1e-9; nr=np.c_[t[:,1],-t[:,0]]
        legs=[l for rid,grp,name,l in __import__('routes').ROUTES if rid==r['id']][0]
        isext=np.array([legs[k][0].startswith('e_') for k in r['lid']]); south=np.array(['south' in legs[k][0] for k in r['lid']])
        two=np.array([legs[k][3] for k in r['lid']])
        key=np.array([(round(p[0]/3),round(p[1]/3)) in done for p in Q]); 
        new=isext&~key
        for p in Q[isext]: done.add((round(p[0]/3),round(p[1]/3)))
        idx=np.where(new)[0]
        if len(idx)<2: continue
        for run in np.split(idx,np.where(np.diff(idx)>1)[0]+1):
            if len(run)<3: continue
            for side,bit in ((-1,1),(1,2)):
                has=(r['mask'][run]&bit)>0
                for val in (True,False):
                    sub=run[has==val]
                    if len(sub)<2: continue
                    for rr in np.split(sub,np.where(np.diff(sub)>1)[0]+1):
                        if len(rr)<2: continue
                        P3=np.c_[Q[rr]+nr[rr]*side*1.8,r['z'][rr]]
                        if val: emit(wh,P3,True)
                        elif not south[rr[0]]:
                            if side<0 and two[rr[0]]:
                                emit(ye,np.c_[Q[rr]+nr[rr]*side*1.95,r['z'][rr]],False); emit(ye,np.c_[Q[rr]+nr[rr]*side*2.25,r['z'][rr]],False)
                            else: emit(ye if side<0 else wh,P3,False)
                        if val and not south[rr[0]]:
                            emit(ye if (side<0 and not two[rr[0]]) else wh,np.c_[Q[rr]+nr[rr]*side*5.4,r['z'][rr]],False)
    res['syn_w']=wh.get(); res['syn_y']=ye.get(); print('synth marks',res['syn_w'][0].shape,res['syn_y'][0].shape)
    pickle.dump(res,open(W+'/scenery.pkl','wb'))
if __name__=='__main__': run()
