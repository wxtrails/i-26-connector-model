"""Export scene binaries + manifest for the web app."""
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from common import *
from surf import Surf
import shapely, shapely.geometry as sg, os
OX,OY=285000.0,210000.0
OUT=W+'/scene'
man={'origin':[OX,OY],'bufs':{},'imgs':{}}
def put(name,arr):
    arr=np.ascontiguousarray(arr); fn=name+'.bin'; arr.tofile(f'{OUT}/{fn}')
    man['bufs'][name]={'dtype':str(arr.dtype),'n':int(arr.size)}
def toscene(x,y,z): return np.c_[np.asarray(x)-OX,np.asarray(z),-(np.asarray(y)-OY)].astype(np.float32)
def stage_terrain():
    T=np.load(W+'/terrain_tiles.npz')
    Hc=T['Hc']; man['terrain']={'bbox':T['bbox'].tolist(),'coarse':24.0,'fine':6.0,'tile':240.0,'cw':int(Hc.shape[1]),'ch':int(Hc.shape[0]),'ntiles':int(len(T['tiles'])),'fn':int(T['fh'].shape[1])}
    q=lambda h: np.clip((h-500.0)*20,0,65535).astype(np.uint16)
    put('terr_coarse',q(Hc)); put('terr_tiles',T['tiles'].astype(np.uint16)); put('terr_fine',q(T['fh']))
def stage_roads(surf):
    m=surf.m; XY=m['XY']; z=m['z']; off=m['off']
    for kind,name in (('G','road_g'),('D','road_d')):
        P=[];I=[];n=0
        for i,s in enumerate(surf.sheets):
            if s['kind']!=kind: continue
            v=s['v']; zz=z[off[i]:off[i+1]]
            P.append(toscene(v[:,0],v[:,1],zz)); 
            t=s['t'].copy()
            # ensure upward-facing winding (counter-clockwise in x,y => after mapping z=-y normal up needs check in app; fix here)
            a=v[t[:,0]];b=v[t[:,1]];c=v[t[:,2]]; cr=(b[:,0]-a[:,0])*(c[:,1]-a[:,1])-(b[:,1]-a[:,1])*(c[:,0]-a[:,0])
            t[cr<0]=t[cr<0][:,::-1]
            I.append(t+n); n+=len(v)
        put(name+'_p',np.vstack(P)); put(name+'_i',np.vstack(I).astype(np.uint32))
    P=[];I=[];n=0      # level breaks between adjacent roadways -> retaining walls
    for i,s in enumerate(surf.sheets):
        tw=s.get('tw')
        if tw is None or len(tw)==0: continue
        v=s['v']; zz=z[off[i]:off[i+1]]
        pa=np.c_[v[tw[:,0]],zz[tw[:,0]]]; pb=np.c_[v[tw[:,1]],zz[tw[:,1]]]; pc=np.c_[v[tw[:,2]],zz[tw[:,2]]]; nr=np.cross(pb-pa,pc-pa)
        tw=tw[np.hypot(nr[:,0],nr[:,1])/(np.abs(nr[:,2])+1e-9)>1.2]      # gentler breaks are left to the graded ground
        if len(tw)==0: continue
        used=np.unique(tw); remap=np.full(len(v),-1); remap[used]=np.arange(len(used))
        P.append(toscene(v[used,0],v[used,1],zz[used])); I.append(remap[tw]+n); n+=len(used)
    put('wall_p',np.vstack(P) if P else np.zeros((0,3),np.float32)); put('wall_i',(np.vstack(I) if I else np.zeros((0,3))).astype(np.uint32))
    # deck sides + parapets + piers
    plan=pickle.load(open(W+'/sheets_plan2.pkl','rb'))
    Gp=plan['G']
    SP=[];SI=[];n=0; BP=[];BI=[];bn=0
    T=np.load(W+'/terrain_conf.npy')
    from terrain import sample
    PIER=[]
    dsheets=[(i,s) for i,s in enumerate(surf.sheets) if s['kind']=='D']
    others=shapely.STRtree([s['poly'] for i,s in dsheets])
    for di,(i,s) in enumerate(dsheets):
        poly=s['poly']
        # neighbours (linked decks): non-stacked ones nearby
        nb=[dsheets[j][1]['poly'] for j in others.query(poly.buffer(7)) if j!=di and dsheets[j][1]['iso']==s['iso'] and dsheets[j][1]['poly'].intersection(poly).area<0.3*min(poly.area,dsheets[j][1]['poly'].area)]
        nbu=shapely.unary_union(nb).buffer(3.5) if nb else sg.Polygon()
        ab=s['abut']
        for ring in [poly.exterior]+list(poly.interiors):
            ls=sg.LineString(ring.coords); L=ls.length; npt=max(4,int(L/4))
            pts=np.array([ls.interpolate(t).coords[0] for t in np.linspace(0,L,npt+1)])
            zz=surf.z_on(i,pts)
            # fill nan by nearest interior sample
            if np.isnan(zz).any():
                c=np.array(poly.representative_point().coords[0]); q=pts+(c-pts)*0.02/np.maximum(1e-6,np.hypot(*(c-pts).T))[:,None]
                z2=surf.z_on(i,pts*0.985+np.array(poly.centroid.coords[0])*0.015)
                zz=np.where(np.isnan(zz),z2,zz)
                ok=np.isfinite(zz)
                if ok.sum()<2: continue
                zz=np.interp(np.arange(len(zz)),np.where(ok)[0],zz[ok])
            mid=(pts[:-1]+pts[1:])/2
            mp=shapely.points(mid)
            open_=shapely.contains(ab,mp) if not ab.is_empty else np.zeros(len(mid),bool)
            if not nbu.is_empty: open_|=shapely.contains(nbu,mp)
            for k in range(len(pts)-1):
                a,b=pts[k],pts[k+1]; za,zb=zz[k],zz[k+1]
                # side skirt (always)
                q=toscene([a[0],b[0],b[0],a[0]],[a[1],b[1],b[1],a[1]],[za+0.02,zb+0.02,zb-1.7,za-1.7])
                SP.append(q); SI.append(np.array([[0,1,2],[0,2,3]])+n); n+=4
                if not open_[k]:
                    # parapet: inset 0.25 m, 0.25 wide, 0.95 high
                    d=b-a; l=np.hypot(*d)+1e-9; nrm=np.array([-d[1],d[0]])/l
                    # inward normal: test
                    if not poly.contains(sg.Point(*(mid[k]+nrm*0.6))): nrm=-nrm
                    a1=a+nrm*0.05;b1=b+nrm*0.05;a2=a+nrm*0.45;b2=b+nrm*0.45
                    xs=[a1[0],b1[0],b1[0],a1[0], a2[0],b2[0],b2[0],a2[0]]; ys=[a1[1],b1[1],b1[1],a1[1],a2[1],b2[1],b2[1],a2[1]]
                    zs=[za,zb,zb+0.95,za+0.95, za,zb,zb+0.95,za+0.95]
                    q=toscene(xs,ys,zs); BP.append(q)
                    BI.append(np.array([[0,1,2],[0,2,3],[5,4,7],[5,7,6],[3,2,6],[3,6,7]])+bn); bn+=8
        # piers on a grid inside the deck
        inner=poly.buffer(-2.5)
        if inner.is_empty: continue
        x0,y0,x1,y1=poly.bounds
        gx=np.arange(x0+5,x1,11.0); gy=np.arange(y0+5,y1,11.0)
        if len(gx)==0 or len(gy)==0: continue
        X,Y=np.meshgrid(gx,gy); pp=np.c_[X.ravel(),Y.ravel()]
        ins=shapely.contains(inner,shapely.points(pp)); pp=pp[ins]
        if len(pp)==0: continue
        zt=sample(T,pp[:,0],pp[:,1]); zd=surf.z_on(i,pp)
        gm=~shapely.contains(Gp.buffer(2.0),shapely.points(pp))
        placed=[]
        for p_,a_,b_,g_ in zip(pp,zt,zd,gm):
            if g_ and np.isfinite(b_) and b_-a_>3.0 and all(np.hypot(*(p_-q))>26.0 for q in placed):
                placed.append(p_); PIER.append([p_[0]-OX,a_-1.0,-(p_[1]-OY),b_-1.6-(a_-1.0)])
    put('deck_side_p',np.vstack(SP)); put('deck_side_i',np.vstack(SI).astype(np.uint32))
    put('deck_par_p',np.vstack(BP)); put('deck_par_i',np.vstack(BI).astype(np.uint32))
    put('piers',np.array(PIER,np.float32)); print('piers',len(PIER))
def stage_extra():
    mk=pickle.load(open(W+'/marks.pkl','rb')); sc=pickle.load(open(W+'/scenery.pkl','rb'))
    put('mk_w_p',mk['white'][0]); put('mk_w_i',mk['white'][1]); put('mk_y_p',mk['yellow'][0]); put('mk_y_i',mk['yellow'][1]); put('conc_p',mk['conc'][0]); put('conc_i',mk['conc'][1])
    put('syn_w_p',sc['syn_w'][0]); put('syn_w_i',sc['syn_w'][1]); put('syn_y_p',sc['syn_y'][0]); put('syn_y_i',sc['syn_y'][1])
    for k in ('bld_xy','bld_n','bld_h','bld_t','trees','far'): put(k,sc[k])
    man['far']=sc['far_meta']; man['signs']=sc['signs']
    man['atlas']=json.load(open(OUT+'/atlas.json'))
    man['imgs']['base']={'file':'base.jpg'}
    for n in range(man['atlas']['n']): man['imgs'][f'atlas{n}']={'file':f'atlas{n}.jpg'}
    import routes as RT
    NAMES={'b0':'I-26','e_south':'I-26 / I-240','e_north':'I-26 / US 19-23','b6':'I-240','e_i240e':'I-240 toward downtown','b10':'Patton Avenue','b23':'Patton Avenue','b8':'Patton Avenue, Bowen Bridge','b9':'Patton Avenue, Bowen Bridge','e_patw':'Patton Avenue, West Asheville','e_pate':'Patton Avenue into downtown',
           'b5':'Ramp: I-240 West to I-26 north','b2':'Ramp: I-26 to I-240 East','b7':'Flyover: I-240 West to I-26 south','b3':'Flyover: I-26 to I-240 East','b22':'Ramp: Patton Avenue to I-26 north','b24':'Ramp: Patton Avenue to I-26 south'}
    R=pickle.load(open(W+'/routes.pkl','rb')); man['routes']=[]
    for r in R:
        legs=[l for rid,grp,name,l in RT.ROUTES if rid==r['id']][0]
        P=r['P']; arr=np.c_[P[:,0]-OX,r['z'],-(P[:,1]-OY),r['v'],r['mask'].astype(np.float32)].astype(np.float32); put('rt_'+r['id'],arr)
        segs=[]; 
        for i,k in enumerate(r['lid']):
            c=legs[k][0]; key=c if c in NAMES else c.rsplit('_',1)[0]; lab=NAMES.get(key,'')
            if not segs or segs[-1][1]!=lab: segs.append([i,lab])
        man['routes'].append({'id':r['id'],'group':r['group'],'name':r['name'],'segs':segs})
    man['labels']=[{'t':'Westgate','x':285650,'y':210150,'k':'pl','h':45,'r':6000},{'t':"Sam's Club",'x':285150,'y':210150,'k':'pl','h':40,'r':6000},{'t':'New Belgium Brewing','x':285958,'y':209719,'k':'pl','h':45,'r':6000},{'t':'I-26','x':285130,'y':209350,'k':'hw','h':45},{'t':'I-26','x':285640,'y':211800,'k':'hw','h':45},{'t':'I-240','x':286520,'y':210260,'k':'hw','h':45},{'t':'Patton Avenue','x':284800,'y':209750,'k':'hw','h':35},
      {'t':'Patton Ave bridge over the river (Bowen Bridge)','x':285970,'y':210010,'h':40,'r':3500},{'t':'French Broad River','x':285830,'y':210330,'k':'rv','h':25},{'t':'New I-26 river bridges','x':285660,'y':210760,'h':60,'r':4000},
      {'t':'I-26 over Patton Avenue','x':285400,'y':209905,'h':55,'r':4000},{'t':'Haywood Road','x':284945,'y':208655,'h':40,'r':5000},{'t':'Downtown Asheville','x':287750,'y':210420,'k':'pl','h':90,'r':20000},{'t':'West Asheville','x':283700,'y':209300,'k':'pl','h':80,'r':20000},
      {'t':'To I-40','x':283000,'y':206800,'k':'pl','h':70,'r':20000},{'t':'To Weaverville (US 19/23)','x':285520,'y':213900,'k':'pl','h':70,'r':20000},{'t':'Riverside Drive','x':286080,'y':210900,'h':30,'r':3000},{'t':'Brevard Road','x':283860,'y':207130,'h':40,'r':5000},{'t':'I-240 flyovers','x':285900,'y':210780,'h':75,'r':4000}]
if __name__=='__main__':
    surf=Surf()
    stage_terrain(); stage_roads(surf); stage_extra()
    json.dump(man,open(OUT+'/manifest.json','w'))
    print({k:v['n'] for k,v in man['bufs'].items()})
