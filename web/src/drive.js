// ---------- routes
const DS=4.0;
A.routes=(M.routes||[]).map(r=>{const P=buf('rt_'+r.id);return Object.assign({},r,{P,n:P.length/5,len:(P.length/5-1)*DS});});
const _v=new T.Vector3(), UP=new T.Vector3(0,1,0);
function rsample(r,s,out){ // position, speed, mask
  const f=Math.min(Math.max(s/DS,0),r.n-1.001), i=Math.floor(f), u=f-i, P=r.P, a=5*i, b=a+5;
  out.x=P[a]+(P[b]-P[a])*u; out.y=P[a+1]+(P[b+1]-P[a+1])*u; out.z=P[a+2]+(P[b+2]-P[a+2])*u; out.v=P[a+3]+(P[b+3]-P[a+3])*u; out.mask=P[a+4]|0; return out;
}
function rtangent(r,s,out){
  const a=rsample(r,s-3,{}), b=rsample(r,s+3,{}); out.set(b.x-a.x,b.y-a.y,b.z-a.z); if(out.lengthSq()<1e-6) out.set(1,0,0); return out.normalize();
}
// ---------- vehicles
function vbox(G,w,h,l,x,y,z,c){ // centred box, y = bottom
  const g=new T.BoxGeometry(w,h,l); g.translate(x,y+h/2,z); const n=g.attributes.position.count, col=new Float32Array(n*3);
  for(let i=0;i<n;i++){col[3*i]=c[0];col[3*i+1]=c[1];col[3*i+2]=c[2];} g.setAttribute('color',new T.BufferAttribute(col,3)); G.push(g);
}
function mergeGeoms(G){
  let n=0,ni=0; for(const g of G){n+=g.attributes.position.count;ni+=g.index.count;}
  const p=new Float32Array(n*3),nr=new Float32Array(n*3),c=new Float32Array(n*3),ix=new Uint32Array(ni); let o=0,oi=0;
  for(const g of G){p.set(g.attributes.position.array,o*3);nr.set(g.attributes.normal.array,o*3);c.set(g.attributes.color.array,o*3);
    for(let k=0;k<g.index.count;k++) ix[oi+k]=g.index.array[k]+o; o+=g.attributes.position.count; oi+=g.index.count;}
  const m=new T.BufferGeometry(); m.setAttribute('position',new T.BufferAttribute(p,3)); m.setAttribute('normal',new T.BufferAttribute(nr,3)); m.setAttribute('color',new T.BufferAttribute(c,3)); m.setIndex(new T.BufferAttribute(ix,1)); return m;
}
const PAINT=[1,1,1], GLASS=[0.10,0.13,0.17], TIRE=[0.04,0.04,0.045], LITE=[1.6,1.5,1.2], TAIL=[0.9,0.05,0.05], GREY=[0.25,0.25,0.26];
function wheels(G,w,l,xs,r){ for(const z of xs) for(const sx of [-1,1]) vbox(G,0.26,2*r,2*r,sx*(w/2-0.10),0,z,TIRE); }
function carLights(G,w,l,y){ vbox(G,w*0.26,0.14,0.06,-w*0.3,y,l/2,LITE); vbox(G,w*0.26,0.14,0.06,w*0.3,y,l/2,LITE); vbox(G,w*0.26,0.14,0.06,-w*0.3,y,-l/2,TAIL); vbox(G,w*0.26,0.14,0.06,w*0.3,y,-l/2,TAIL); }
const VT=[
 {k:'sedan',len:4.6,wd:1.82,p:0.33,mk:G=>{vbox(G,1.82,0.62,4.6,0,0.28,0,PAINT);vbox(G,1.6,0.5,2.3,0,0.9,-0.25,GLASS);vbox(G,1.5,0.06,1.9,0,1.4,-0.25,PAINT);wheels(G,1.82,4.6,[1.4,-1.4],0.33);carLights(G,1.82,4.6,0.62);}},
 {k:'suv',len:4.8,wd:1.95,p:0.26,mk:G=>{vbox(G,1.95,0.85,4.8,0,0.32,0,PAINT);vbox(G,1.78,0.62,3.0,0,1.17,-0.45,GLASS);vbox(G,1.7,0.07,2.8,0,1.79,-0.45,PAINT);wheels(G,1.95,4.8,[1.45,-1.45],0.38);carLights(G,1.95,4.8,0.85);}},
 {k:'pickup',len:5.6,wd:2.0,p:0.15,mk:G=>{vbox(G,2.0,0.8,5.6,0,0.38,0,PAINT);vbox(G,1.82,0.62,1.9,0,1.18,0.55,GLASS);vbox(G,1.74,0.07,1.75,0,1.8,0.55,PAINT);vbox(G,1.7,0.12,2.1,0,1.18,-1.6,GREY);wheels(G,2.0,5.6,[1.75,-1.7],0.42);carLights(G,2.0,5.6,0.9);}},
 {k:'van',len:5.4,wd:2.05,p:0.08,mk:G=>{vbox(G,2.05,1.75,5.4,0,0.35,0,PAINT);vbox(G,2.07,0.6,1.1,0,1.35,2.0,GLASS);vbox(G,1.9,0.55,0.08,0,1.4,2.7,GLASS);wheels(G,2.05,5.4,[1.75,-1.6],0.38);carLights(G,2.05,5.4,0.75);}},
 {k:'box',len:7.6,wd:2.4,p:0.06,mk:G=>{vbox(G,2.2,1.7,2.0,0,0.5,2.8,PAINT);vbox(G,2.22,0.6,0.9,0,1.5,3.2,GLASS);vbox(G,2.4,2.5,5.5,0,0.85,-1.0,[0.93,0.93,0.9]);vbox(G,2.0,0.4,7.2,0,0.45,0,GREY);wheels(G,2.4,7.6,[2.8,-2.3],0.46);carLights(G,2.2,7.6,0.8);}},
 {k:'semi',len:19.5,wd:2.6,p:0.08,mk:G=>{vbox(G,2.5,2.2,2.8,0,0.6,8.3,PAINT);vbox(G,2.52,0.7,1.0,0,1.9,9.1,GLASS);vbox(G,2.5,0.9,1.4,0,2.8,7.7,PAINT);vbox(G,2.6,2.75,14.6,0,1.25,-2.4,[0.92,0.92,0.93]);vbox(G,2.2,0.35,6.5,0,0.6,6.0,GREY);wheels(G,2.6,19,[8.9,5.6,4.4,-7.4,-8.6],0.52);carLights(G,2.5,19.5,1.0);}},
 {k:'bus',len:11.5,wd:2.55,p:0.04,mk:G=>{vbox(G,2.55,2.75,11.5,0,0.4,0,PAINT);vbox(G,2.58,0.95,10.3,0,1.65,0,GLASS);vbox(G,2.3,0.9,0.08,0,1.7,5.75,GLASS);wheels(G,2.55,11.5,[3.9,-3.4],0.5);carLights(G,2.55,11.5,0.9);}},
];
const CAR_COLS=[0xf2f2f2,0xd9dadb,0x9da2a6,0x55585c,0x1c1d1f,0x1c1d1f,0xf2f2f2,0x8c1c1c,0x1f3f7a,0x2c5d8f,0x6b7f5a,0xb9b2a0,0x7a1f2b,0x3a3f44,0xc8cdd2,0xc4622d];
const MAXV=480;
VT.forEach(t=>{const G=[]; t.mk(G); t.geom=mergeGeoms(G); t.mesh=new T.InstancedMesh(t.geom,new T.MeshLambertMaterial({vertexColors:true}),MAXV); t.mesh.count=0; t.mesh.frustumCulled=false; t.mesh.instanceMatrix.setUsage(T.DynamicDrawUsage); scene.add(t.mesh);});
const shadowGeom=new T.PlaneGeometry(1,1); shadowGeom.rotateX(-Math.PI/2);
const shadowMesh=new T.InstancedMesh(shadowGeom,new T.MeshBasicMaterial({color:0x000000,transparent:true,opacity:0.33,depthWrite:false,polygonOffset:true,polygonOffsetFactor:-3,polygonOffsetUnits:-6}),MAXV+1); shadowMesh.frustumCulled=false; shadowMesh.count=0; scene.add(shadowMesh);
// ego car
const ego={route:null,s:0,v:0,vdes:0,limit:1e9,active:false,paused:true,done:false,pos:new T.Vector3(),fwd:new T.Vector3(1,0,0),lookF:new T.Vector3(1,0,0)};
(function(){const G=[]; VT[1].mk(G); const g=mergeGeoms(G); const col=g.attributes.color; for(let i=0;i<col.count;i++){ if(col.getX(i)===1&&col.getY(i)===1&&col.getZ(i)===1) col.setXYZ(i,0.72,0.08,0.08); }
  ego.mesh=new T.Mesh(g,new T.MeshLambertMaterial({vertexColors:true})); ego.mesh.visible=false; ego.mesh.matrixAutoUpdate=false; scene.add(ego.mesh);
  const hg=new T.BoxGeometry(1.6,0.05,0.75); hg.translate(0,0.60,1.38); ego.hood=new T.Mesh(hg,new T.MeshLambertMaterial({color:0xa01414})); ego.hood.visible=false; ego.hood.matrixAutoUpdate=false; scene.add(ego.hood);})();
const _f=new T.Vector3(), _r=new T.Vector3(), _u=new T.Vector3(), _m=new T.Matrix4(), _p=new T.Vector3(), _col=new T.Color(), _s={};
function basis(m,pos,f,scaleX,scaleY,scaleZ){ // model: +Z forward, +Y up
  _r.crossVectors(UP,f).normalize(); _u.crossVectors(f,_r).normalize();
  m.set(_r.x*scaleX,_u.x*scaleY,f.x*scaleZ,pos.x, _r.y*scaleX,_u.y*scaleY,f.y*scaleZ,pos.y, _r.z*scaleX,_u.z*scaleY,f.z*scaleZ,pos.z, 0,0,0,1); return m;
}
// ---------- traffic
const cars=[]; const trng=rnd(2024);
A.traffic={density:1.0};
function pickType(r){ let x=trng(), street=r.group.indexOf('Patton Avenue')===0&&r.id.indexOf('i26')<0; for(let i=0;i<VT.length;i++){ let p=VT[i].p; if(street&&VT[i].k==='semi') p=0.01; x-=p; if(x<=0) return i; } return 0; }
function spawnTraffic(){
  cars.length=0; const R=A.routes; if(!R.length) return;
  const wts=R.map(r=>r.len*(r.id==='i26n'||r.id==='i26s'?1.0:(r.id.indexOf('pat_')===0?0.8:0.45))), tot=wts.reduce((a,b)=>a+b,0);
  const N=Math.min(MAXV,Math.round(230*A.traffic.density));
  for(let k=0;k<N;k++){
    let x=trng()*tot, ri=0; while(ri<R.length-1&&x>wts[ri]){x-=wts[ri];ri++;}
    const r=R[ri], t=pickType(r), s=trng()*r.len, big=VT[t].k==='semi'||VT[t].k==='bus'||VT[t].k==='box';
    const c={r,t,s,v:0,f:(big?0.80:0.86)+trng()*0.2,lane:trng()<0.34?0:(trng()<0.5?-1:1),lat:0,col:big&&VT[t].k!=='bus'?[0xf2f2f2,0x2c5d8f,0x8c1c1c,0x3a3f44,0xc8cdd2][Math.floor(trng()*5)]:CAR_COLS[Math.floor(trng()*CAR_COLS.length)],x:0,y:0,z:0,fx:1,fy:0,fz:0};
    if(big&&c.lane<0) c.lane=0;
    rsample(r,s,_s); c.v=_s.v*c.f; c.lat=c.lane*3.6*(((c.lane<0&&(_s.mask&1))||(c.lane>0&&(_s.mask&2)))?1:0); cars.push(c);
  }
}
const grid=new Map(); const CELL=45;
function gkey(x,z){return (Math.floor(x/CELL)+4096)*8192+(Math.floor(z/CELL)+4096);}
function stepTraffic(dt){
  // positions
  for(const c of cars){
    rsample(c.r,c.s,_s); rtangent(c.r,c.s,_f);
    const ok=c.lane===0||(c.lane<0&&(_s.mask&1))||(c.lane>0&&(_s.mask&2)); const tl=ok?c.lane*3.6:0; c.lat+=Math.sign(tl-c.lat)*Math.min(Math.abs(tl-c.lat),1.1*dt);
    // right of travel = f x up
    const rx=-_f.z, rz=_f.x, rl=Math.hypot(rx,rz)||1;
    c.x=_s.x+rx/rl*c.lat; c.y=_s.y; c.z=_s.z+rz/rl*c.lat; c.fx=_f.x; c.fy=_f.y; c.fz=_f.z; c.vt=_s.v*c.f;
  }
  grid.clear();
  for(const c of cars){const k=gkey(c.x,c.z); let a=grid.get(k); if(!a){a=[];grid.set(k,a);} a.push(c);}
  const E=ego.active?{x:ego.pos.x,y:ego.pos.y,z:ego.pos.z,fx:ego.fwd.x,fy:ego.fwd.y,fz:ego.fwd.z,v:(ego.paused?0:ego.v),t:1,ego:true}:null;
  if(E){const k=gkey(E.x,E.z); let a=grid.get(k); if(!a){a=[];grid.set(k,a);} a.push(E);}
  ego.limit=1e9;
  if(E){ // the viewer keeps a safe gap like everyone else
    const gap=10+1.5*E.v, cx=Math.floor(E.x/CELL), cz=Math.floor(E.z/CELL);
    for(let i=-1;i<=1;i++)for(let j=-1;j<=1;j++){const a=grid.get((cx+i+4096)*8192+(cz+j+4096)); if(!a) continue;
      for(const o of a){ if(o.ego) continue; const dx=o.x-E.x,dz=o.z-E.z; if(Math.abs(o.y-E.y)>3) continue;
        const ah=dx*E.fx+dz*E.fz; if(ah<=0.5||ah>60) continue; if(Math.abs(-dx*E.fz+dz*E.fx)>1.9) continue; if(o.fx*E.fx+o.fz*E.fz<0.4) continue;
        o.boost=ego.vdes*1.05; if(ah<gap+VT[o.t].len*0.5){ const free=ah-VT[o.t].len*0.5-4.5; const lim=free<=0?Math.min(o.v*0.6,2):o.v+free*0.45; if(lim<ego.limit) ego.limit=lim; } }}
  }
  for(const c of cars){
    let vt=c.vt; if(c.boost){ if(c.boost>vt) vt=c.boost; c.boost=0; }
    const gap=9+1.5*c.v+VT[c.t].len*0.5; const cx=Math.floor(c.x/CELL), cz=Math.floor(c.z/CELL);
    for(let i=-1;i<=1;i++)for(let j=-1;j<=1;j++){const a=grid.get((cx+i+4096)*8192+(cz+j+4096)); if(!a) continue;
      for(const o of a){ if(o===c) continue; const dx=o.x-c.x,dz=o.z-c.z,dy=o.y-c.y; if(Math.abs(dy)>3) continue;
        const ah=dx*c.fx+dz*c.fz; if(ah<=0.5||ah>gap+VT[o.t].len*0.5) continue; const lt=Math.abs(-dx*c.fz+dz*c.fx); if(lt>1.9) continue; if(o.fx*c.fx+o.fz*c.fz<0.4) continue;
        const free=ah-VT[o.t].len*0.5-VT[c.t].len*0.5-2.5; const lim=free<=0?Math.min(o.v*0.6,2):o.v+free*0.45; if(lim<vt) vt=lim; }}
    c.v+=Math.max(-7*dt,Math.min(2.6*dt,vt-c.v)); if(c.v<0) c.v=0;
  }
  for(const c of cars){ c.s+=c.v*dt; if(c.s>c.r.len-6){ c.s=4+trng()*30; c.lat=0; rsample(c.r,c.s,_s); c.v=_s.v*c.f; } }
}
function drawTraffic(){
  const counts=VT.map(()=>0); let sc=0;
  for(const c of cars){
    const t=VT[c.t]; _p.set(c.x,c.y+0.03,c.z); _f.set(c.fx,c.fy,c.fz); basis(_m,_p,_f,1,1,1);
    t.mesh.setMatrixAt(counts[c.t],_m); _col.setHex(c.col); t.mesh.setColorAt(counts[c.t],_col); counts[c.t]++;
    _p.y=c.y+0.06; basis(_m,_p,_f,t.wd*1.12,1,t.len*1.06); shadowMesh.setMatrixAt(sc++,_m);
  }
  if(ego.active&&ego.mesh.visible){ _p.copy(ego.pos); _p.y+=0.06; basis(_m,_p,ego.fwd,2.2,1,5.1); shadowMesh.setMatrixAt(sc++,_m); }
  VT.forEach((t,i)=>{t.mesh.count=counts[i]; t.mesh.instanceMatrix.needsUpdate=true; if(t.mesh.instanceColor) t.mesh.instanceColor.needsUpdate=true;});
  shadowMesh.count=sc; shadowMesh.instanceMatrix.needsUpdate=true;
}
// ---------- camera
const view={tx:500,ty:620,tz:-500,dist:2600,az:0.5,pol:1.0};
A.view=view; A.camMode='orbit'; A.timeScale=1; A.follow=true;
function applyOrbit(){
  const sp=Math.sin(view.pol),cp=Math.cos(view.pol);
  camera.position.set(view.tx+view.dist*sp*Math.sin(view.az),view.ty+view.dist*cp,view.tz+view.dist*sp*Math.cos(view.az));
  const g=terr.height(camera.position.x,camera.position.z)+2; if(camera.position.y<g) camera.position.y=g;
  camera.lookAt(view.tx,view.ty,view.tz);
}
(function(){
  let drag=null; const pts=new Map(); let pinch=0;
  canvas.addEventListener('pointerdown',e=>{pts.set(e.pointerId,{x:e.clientX,y:e.clientY}); drag={x:e.clientX,y:e.clientY,b:e.button}; canvas.setPointerCapture(e.pointerId); if(pts.size===2){const a=[...pts.values()]; pinch=Math.hypot(a[0].x-a[1].x,a[0].y-a[1].y);} });
  const up=e=>{pts.delete(e.pointerId); drag=null;}; canvas.addEventListener('pointerup',up); canvas.addEventListener('pointercancel',up);
  canvas.addEventListener('pointermove',e=>{
    if(!pts.has(e.pointerId)) return; pts.set(e.pointerId,{x:e.clientX,y:e.clientY});
    if(pts.size===2){const a=[...pts.values()], d=Math.hypot(a[0].x-a[1].x,a[0].y-a[1].y); if(pinch>0) view.dist=Math.min(30000,Math.max(12,view.dist*pinch/d)); pinch=d; return;}
    if(!drag) return; const dx=e.clientX-drag.x,dy=e.clientY-drag.y; drag.x=e.clientX; drag.y=e.clientY;
    if(A.camMode!=='orbit'){ if(Math.abs(dx)+Math.abs(dy)>3&&A.setCamMode) { if(ego.active){view.tx=ego.pos.x;view.ty=ego.pos.y;view.tz=ego.pos.z;view.dist=Math.min(view.dist,160);} A.setCamMode('orbit'); } return; }
    if(drag.b===2||e.shiftKey||e.ctrlKey){const s=view.dist*0.0013;const c=Math.cos(view.az),sn=Math.sin(view.az);view.tx-=(dx*c+dy*sn)*s;view.tz-=(-dx*sn+dy*c)*s;view.ty=terr.height(view.tx,view.tz); A.follow=false;}
    else{view.az-=dx*0.005;view.pol=Math.min(1.55,Math.max(0.06,view.pol-dy*0.005));}
  });
  canvas.addEventListener('wheel',e=>{e.preventDefault(); if(A.camMode!=='orbit') return; view.dist=Math.min(30000,Math.max(12,view.dist*Math.exp(e.deltaY*0.0011)));},{passive:false});
  canvas.addEventListener('contextmenu',e=>e.preventDefault());
})();
// keyboard camera, held keys act every frame: W/A/S/D slide the orbit centre forward/left/back/right as seen from the camera (like right-drag);
// Q/E turn the view right/left around the same centre and R/F tilt it down (more overhead) / up (toward the horizon), like left-drag
const camKeys=new Set(), CAM_KEYS=['KeyW','KeyA','KeyS','KeyD','KeyQ','KeyE','KeyR','KeyF'];
window.addEventListener('keydown',e=>{ const t=e.target.tagName; if(t==='SELECT'||t==='INPUT'||t==='TEXTAREA'||e.ctrlKey||e.metaKey||e.altKey) return;
  if(CAM_KEYS.includes(e.code)){ camKeys.add(e.code); e.preventDefault(); } });
window.addEventListener('keyup',e=>camKeys.delete(e.code)); window.addEventListener('blur',()=>camKeys.clear());
function keyCam(dt){
  if(!camKeys.size||A.camMode==='manual') return; const k=c=>camKeys.has(c)?1:0;
  const f=k('KeyW')-k('KeyS'), r=k('KeyD')-k('KeyA'), turn=k('KeyQ')-k('KeyE'), tilt=k('KeyF')-k('KeyR'); if(!f&&!r&&!turn&&!tilt) return;
  if(A.camMode!=='orbit'&&A.setCamMode){ if(ego.active){view.tx=ego.pos.x;view.ty=ego.pos.y;view.tz=ego.pos.z;view.dist=Math.min(view.dist,160);} A.setCamMode('orbit'); }
  if(turn) view.az-=turn*1.05*dt; // about 60 degrees a second; Q swings the heading clockwise seen from above, as drag-right does
  if(tilt) view.pol=Math.min(1.55,Math.max(0.06,view.pol+tilt*0.61*dt)); // about 35 degrees a second, within the left-drag limits
  if(f||r){ const s=view.dist*0.6*dt/Math.hypot(f,r), c=Math.cos(view.az), sn=Math.sin(view.az); // about a third of the screen width per second at any zoom
    view.tx+=(r*c-f*sn)*s; view.tz+=(-r*sn-f*c)*s; view.ty=terr.height(view.tx,view.tz); A.follow=false; }
}
function resize(){const w=canvas.clientWidth,h=canvas.clientHeight;renderer.setSize(w,h,false);camera.aspect=w/h;camera.updateProjectionMatrix();}
window.addEventListener('resize',resize); resize();
const _a={}, _b={}, _t=new T.Vector3(), _look=new T.Vector3();
function updateEgo(dt){
  const r=ego.route; if(!r) return;
  if(!ego.paused&&!ego.done){ rsample(r,ego.s,_a); ego.vdes=_a.v; const tv=Math.min(_a.v,ego.limit==null?1e9:ego.limit); ego.v+=Math.max(-7*dt,Math.min(3*dt,tv-ego.v)); if(ego.v<0) ego.v=0; ego.s+=ego.v*dt; if(ego.s>=r.len-2){ego.s=r.len-2; ego.done=true; ego.paused=true; if(A.onDone) A.onDone();} }
  rsample(r,ego.s,_a); ego.pos.set(_a.x,_a.y,_a.z); rtangent(r,ego.s,ego.fwd);
  basis(ego.mesh.matrix,_p.set(_a.x,_a.y+0.03,_a.z),ego.fwd,1,1,1); ego.hood.matrix.copy(ego.mesh.matrix);
}
function updateCamera(dt){
  keyCam(dt); const mode=A.camMode, r=ego.route;
  ego.mesh.visible=ego.active&&mode!=='driver'; ego.hood.visible=ego.active&&mode==='driver';
  if(mode==='orbit'||!ego.active){ if(ego.active&&A.follow){view.tx=ego.pos.x;view.ty=ego.pos.y+1;view.tz=ego.pos.z;} applyOrbit(); return; }
  const k=1-Math.exp(-dt*4.0);
  if(mode==='driver'){
    rsample(r,ego.s+34,_b); _t.set(_b.x-ego.pos.x,_b.y-ego.pos.y,_b.z-ego.pos.z).normalize(); ego.lookF.lerp(_t,k).normalize();
    _r.crossVectors(UP,ego.fwd).normalize();
    camera.position.copy(ego.pos).addScaledVector(UP,1.32).addScaledVector(_r,0.36).addScaledVector(ego.fwd,0.25);
    _look.copy(camera.position).addScaledVector(ego.lookF,30); _look.y-=1.2; camera.lookAt(_look);
  } else if(mode==='chase'){
    rsample(r,ego.s-13,_b); _t.set(_b.x,_b.y+4.6,_b.z); if(ego.s<13) _t.copy(ego.pos).addScaledVector(ego.fwd,-13).setY(ego.pos.y+4.6);
    camera.position.lerp(_t,ego.camInit?1:Math.min(1,k*2)); rsample(r,ego.s+16,_b); camera.lookAt(_b.x,_b.y+1.2,_b.z);
  } else { // helicopter
    _t.copy(ego.pos).addScaledVector(ego.fwd,-95); _t.y=ego.pos.y+85; camera.position.lerp(_t,ego.camInit?1:k*0.6); rsample(r,ego.s+45,_b); camera.lookAt(_b.x,_b.y,_b.z);
  }
  ego.camInit=false;
}
A.startRoute=function(id){
  const r=A.routes.find(x=>x.id===id); if(!r) return; ego.route=r; ego.s=2; ego.v=rsample(r,2,_a).v; ego.vdes=ego.v; ego.limit=1e9; ego.done=false; ego.paused=false; ego.active=true; ego.camInit=true;
  rsample(r,ego.s,_a); ego.pos.set(_a.x,_a.y,_a.z); rtangent(r,ego.s,ego.fwd); ego.lookF.copy(ego.fwd);
  if(A.camMode==='orbit'&&A.setCamMode) A.setCamMode('driver');
};
A.stopRoute=function(){ ego.active=false; ego.route=null; ego.paused=true; ego.mesh.visible=false; ego.hood.visible=false; };
A.ego=ego; A.spawnTraffic=spawnTraffic; A.rsample=rsample;
