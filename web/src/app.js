(function(){
'use strict';
const M=window.MANIFEST, T=THREE;
const TYPES={float32:Float32Array,uint32:Uint32Array,uint16:Uint16Array,uint8:Uint8Array,int16:Int16Array,int32:Int32Array};
function buf(name){
  const el=document.getElementById('buf_'+name), info=M.bufs[name];
  if(!el||!info) return null;
  const bin=atob(el.textContent); el.textContent='';
  const u8=new Uint8Array(bin.length);
  for(let i=0;i<bin.length;i++) u8[i]=bin.charCodeAt(i);
  return new TYPES[info.dtype](u8.buffer,0,info.n);
}
function imgURL(name){
  const el=document.getElementById('img_'+name); if(!el) return null;
  return 'data:'+(el.dataset.mime||'image/jpeg')+';base64,'+el.textContent;
}
const OX=M.origin[0], OY=M.origin[1];
const canvas=document.getElementById('c');
const renderer=new T.WebGLRenderer({canvas,antialias:true,powerPreference:'high-performance'});
renderer.setPixelRatio(Math.min(window.devicePixelRatio||1,2));
renderer.outputColorSpace=T.SRGBColorSpace;
const scene=new T.Scene();
const SKY_TOP=new T.Color(0x3f7fd6), SKY_HOR=new T.Color(0xcfe2f1);
scene.fog=new T.Fog(SKY_HOR.getHex(),1500,42000);
const camera=new T.PerspectiveCamera(55,1,0.5,90000);
// ---------- sky dome
(function(){
  const g=new T.SphereGeometry(80000,32,16);
  const m=new T.ShaderMaterial({side:T.BackSide,depthWrite:false,fog:false,
    uniforms:{top:{value:SKY_TOP},hor:{value:SKY_HOR},sun:{value:new T.Vector3(0.45,0.62,0.64).normalize()}},
    vertexShader:'varying vec3 d;void main(){d=normalize(position);vec4 p=modelViewMatrix*vec4(position,1.0);gl_Position=projectionMatrix*p;gl_Position.z=gl_Position.w;}',
    fragmentShader:'varying vec3 d;uniform vec3 top;uniform vec3 hor;uniform vec3 sun;void main(){float h=clamp(d.y,0.0,1.0);vec3 c=mix(hor,top,pow(h,0.55));float s=max(dot(normalize(d),sun),0.0);c+=vec3(1.0,0.95,0.8)*(pow(s,600.0)*1.5+pow(s,12.0)*0.12);if(d.y<0.0)c=hor;gl_FragColor=vec4(c,1.0);}'});
  const sky=new T.Mesh(g,m); sky.frustumCulled=false; sky.renderOrder=-10; scene.add(sky);
  scene.userData.sky=sky;
})();
const sunDir=new T.Vector3(0.45,0.62,0.64).normalize();
const sun=new T.DirectionalLight(0xfff4e0,2.2); sun.position.copy(sunDir).multiplyScalar(1000); scene.add(sun);
scene.add(new T.HemisphereLight(0xcfe3ff,0x6b6a5a,1.25));

// ---------- terrain
const TR=M.terrain;
const terr={};
(function(){
  const [X0,Y0,X1,Y1]=TR.bbox; const cw=TR.cw,ch=TR.ch,cs=TR.coarse;
  const hc=buf('terr_coarse'), tiles=buf('terr_tiles'), hf=buf('terr_fine');
  const fn=TR.fn, ts=TR.tile, nx=Math.round((X1-X0)/ts), ny=Math.round((Y1-Y0)/ts);
  const fineIdx=new Int32Array(nx*ny).fill(-1);
  for(let k=0;k<TR.ntiles;k++) fineIdx[tiles[2*k+1]*nx+tiles[2*k]]=k;
  terr.hc=hc; terr.hf=hf; terr.fineIdx=fineIdx; terr.nx=nx; terr.ny=ny;
  const pos=new Float32Array(cw*ch*3), uv=new Float32Array(cw*ch*2);
  for(let j=0;j<ch;j++)for(let i=0;i<cw;i++){const k=j*cw+i;const x=X0+i*cs,y=Y0+j*cs;pos[3*k]=x-OX;pos[3*k+1]=hc[k]/20+500;pos[3*k+2]=-(y-OY);uv[2*k]=i/(cw-1);uv[2*k+1]=j/(ch-1);}
  const per=Math.round(ts/cs); const idx=[];
  for(let j=0;j<ch-1;j++)for(let i=0;i<cw-1;i++){
    if(fineIdx[Math.floor(j/per)*nx+Math.floor(i/per)]>=0) continue;
    const a=j*cw+i,b=a+1,c=a+cw,d=c+1; idx.push(a,b,d,a,d,c);
  }
  const g=new T.BufferGeometry(); g.setAttribute('position',new T.BufferAttribute(pos,3)); g.setAttribute('uv',new T.BufferAttribute(uv,2)); g.setIndex(idx); g.computeVertexNormals();
  terr.matCoarse=new T.MeshLambertMaterial({color:0x7d8f63});
  terr.coarse=new T.Mesh(g,terr.matCoarse); scene.add(terr.coarse);
  // fine tiles
  const n1=fn, fs=TR.fine; terr.tiles=[];
  const tidx=[]; for(let j=0;j<n1-1;j++)for(let i=0;i<n1-1;i++){const a=j*n1+i,b=a+1,c=a+n1,d=c+1;tidx.push(a,b,d,a,d,c);}
  const tuv=new Float32Array(n1*n1*2); for(let j=0;j<n1;j++)for(let i=0;i<n1;i++){tuv[2*(j*n1+i)]=i/(n1-1);tuv[2*(j*n1+i)+1]=j/(n1-1);}
  const uvAttr=new T.BufferAttribute(tuv,2), idxAttr=new T.BufferAttribute(new Uint16Array(tidx),1);
  for(let k=0;k<TR.ntiles;k++){
    const ti=tiles[2*k],tj=tiles[2*k+1]; const p=new Float32Array(n1*n1*3);
    for(let j=0;j<n1;j++)for(let i=0;i<n1;i++){const q=j*n1+i;p[3*q]=X0+ti*ts+i*fs-OX;p[3*q+1]=hf[k*n1*n1+q]/20+500;p[3*q+2]=-(Y0+tj*ts+j*fs-OY);}
    const tg=new T.BufferGeometry(); tg.setAttribute('position',new T.BufferAttribute(p,3)); tg.setAttribute('uv',uvAttr); tg.setIndex(idxAttr); tg.computeVertexNormals();
    const mesh=new T.Mesh(tg,terr.matCoarse); scene.add(mesh); terr.tiles.push(mesh);
  }
  // height lookup (scene coords)
  terr.height=function(sx,sz){
    const x=sx+OX,y=-sz+OY; if(x<X0||x>=X1||y<Y0||y>=Y1) return 600;
    const ti=Math.floor((x-X0)/ts),tj=Math.floor((y-Y0)/ts),k=fineIdx[tj*nx+ti];
    if(k>=0){const fx=(x-X0-ti*ts)/fs,fy=(y-Y0-tj*ts)/fs,i=Math.min(n1-2,Math.floor(fx)),j=Math.min(n1-2,Math.floor(fy)),u=fx-i,v=fy-j,b=k*n1*n1+j*n1+i;
      return ((hf[b]*(1-u)+hf[b+1]*u)*(1-v)+(hf[b+n1]*(1-u)+hf[b+n1+1]*u)*v)/20+500;}
    const fx=(x-X0)/cs,fy=(y-Y0)/cs,i=Math.min(cw-2,Math.floor(fx)),j=Math.min(ch-2,Math.floor(fy)),u=fx-i,v=fy-j,b=j*cw+i;
    return ((hc[b]*(1-u)+hc[b+1]*u)*(1-v)+(hc[b+cw]*(1-u)+hc[b+cw+1]*u)*v)/20+500;
  };
})();

// ---------- roads
function meshFrom(pn,inn,mat){
  const p=buf(pn), i=buf(inn); if(!p) return null;
  const g=new T.BufferGeometry(); g.setAttribute('position',new T.BufferAttribute(p,3)); g.setIndex(new T.BufferAttribute(i,1)); g.computeVertexNormals();
  const m=new T.Mesh(g,mat); scene.add(m); return m;
}
const matAsphalt=new T.MeshLambertMaterial({color:0x4a4c50,polygonOffset:true,polygonOffsetFactor:-1,polygonOffsetUnits:-2});
const matDeck=new T.MeshLambertMaterial({color:0x55575b});
const matConc=new T.MeshLambertMaterial({color:0xb9b5aa,side:T.DoubleSide});
meshFrom('road_g_p','road_g_i',matAsphalt);
meshFrom('road_d_p','road_d_i',matDeck);
meshFrom('deck_side_p','deck_side_i',matConc);
meshFrom('deck_par_p','deck_par_i',matConc);
(function(){
  const P=buf('piers'); if(!P) return; const n=P.length/4;
  const g=new T.CylinderGeometry(1.1,1.1,1,10); g.translate(0,0.5,0);
  const im=new T.InstancedMesh(g,matConc,n); const m=new T.Matrix4();
  for(let k=0;k<n;k++){m.makeScale(1,P[4*k+3],1);m.setPosition(P[4*k],P[4*k+1],P[4*k+2]);im.setMatrixAt(k,m);}
  scene.add(im);
})();

// ---------- camera control (orbit)
const view={tx:500,ty:620,tz:-500,dist:2500,az:0.6,pol:1.0};
function applyView(){
  const sp=Math.sin(view.pol),cp=Math.cos(view.pol);
  camera.position.set(view.tx+view.dist*sp*Math.sin(view.az),view.ty+view.dist*cp,view.tz+view.dist*sp*Math.cos(view.az));
  camera.lookAt(view.tx,view.ty,view.tz);
}
(function(){
  let drag=null;
  canvas.addEventListener('pointerdown',e=>{drag={x:e.clientX,y:e.clientY,b:e.button||(e.shiftKey?2:0)};canvas.setPointerCapture(e.pointerId);});
  canvas.addEventListener('pointerup',e=>{drag=null;});
  canvas.addEventListener('pointermove',e=>{
    if(!drag) return; const dx=e.clientX-drag.x,dy=e.clientY-drag.y; drag.x=e.clientX; drag.y=e.clientY;
    if(drag.b===2||e.shiftKey){const s=view.dist*0.0012;const c=Math.cos(view.az),sn=Math.sin(view.az);view.tx-=(dx*c+dy*sn)*s;view.tz-=(-dx*sn+dy*c)*s;view.ty=terr.height(view.tx,view.tz);}
    else{view.az-=dx*0.005;view.pol=Math.min(1.54,Math.max(0.08,view.pol-dy*0.005));}
  });
  canvas.addEventListener('wheel',e=>{e.preventDefault();view.dist=Math.min(30000,Math.max(15,view.dist*Math.exp(e.deltaY*0.001)));},{passive:false});
  canvas.addEventListener('contextmenu',e=>e.preventDefault());
})();
function resize(){const w=canvas.clientWidth,h=canvas.clientHeight;renderer.setSize(w,h,false);camera.aspect=w/h;camera.updateProjectionMatrix();}
window.addEventListener('resize',resize); resize();
let manualCam=false;
function frame(){
  if(!manualCam) applyView();
  scene.userData.sky.position.copy(camera.position);
  renderer.render(scene,camera);
  requestAnimationFrame(frame);
}
document.getElementById('loading').style.display='none';
window.APP={scene,camera,view,terr,renderer,
  setCam(px,py,pz,tx,ty,tz){manualCam=true;camera.position.set(px,py,pz);camera.lookAt(tx,ty,tz);},
  sp(x,y){return [x-OX,-(y-OY)];}};
requestAnimationFrame(frame);
window.READY=true;
})();
