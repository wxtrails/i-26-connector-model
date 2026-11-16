'use strict';
const M=window.MANIFEST, T=THREE, A={};
const TYPES={float32:Float32Array,uint32:Uint32Array,uint16:Uint16Array,uint8:Uint8Array,int16:Int16Array,int32:Int32Array};
function buf(name){
  const el=document.getElementById('buf_'+name), info=M.bufs[name];
  if(!el||!info) return null;
  const bin=atob(el.textContent); el.textContent='';
  const u8=new Uint8Array(bin.length);
  for(let i=0;i<bin.length;i++) u8[i]=bin.charCodeAt(i);
  return new TYPES[info.dtype](u8.buffer,0,info.n);
}
let pendingImgs=0;
function loadTex(name,opts){
  const el=document.getElementById('img_'+name); const tex=new T.Texture();
  tex.colorSpace=T.SRGBColorSpace; tex.anisotropy=8; Object.assign(tex,opts||{});
  if(!el) return tex;
  const img=new Image(); pendingImgs++;
  img.onload=()=>{tex.image=img;tex.needsUpdate=true;pendingImgs--;};
  img.onerror=()=>{pendingImgs--;};
  img.src='data:'+(el.dataset.mime||'image/jpeg')+';base64,'+el.textContent; el.textContent='';
  return tex;
}
const OX=M.origin[0], OY=M.origin[1];
const SP=(x,y)=>[x-OX,-(y-OY)];
const canvas=document.getElementById('c');
const renderer=new T.WebGLRenderer({canvas,antialias:true,powerPreference:'high-performance'});
renderer.setPixelRatio(Math.min(window.devicePixelRatio||1,2));
renderer.outputColorSpace=T.SRGBColorSpace;
const scene=new T.Scene();
const SKY_TOP=new T.Color(0x3a78cf), SKY_HOR=new T.Color(0xc9dcec);
scene.fog=new T.Fog(0xbfd3e6,1200,34000);
const camera=new T.PerspectiveCamera(58,1,0.6,120000);
const sunDir=new T.Vector3(0.42,0.66,0.62).normalize();
// ---------- sky dome with soft clouds
(function(){
  const g=new T.SphereGeometry(90000,40,20);
  const m=new T.ShaderMaterial({side:T.BackSide,depthWrite:false,fog:false,
    uniforms:{top:{value:SKY_TOP},hor:{value:SKY_HOR},sun:{value:sunDir},time:{value:0}},
    vertexShader:'varying vec3 d;void main(){d=normalize(position);vec4 p=modelViewMatrix*vec4(position,1.0);gl_Position=projectionMatrix*p;gl_Position.z=gl_Position.w;}',
    fragmentShader:`varying vec3 d;uniform vec3 top;uniform vec3 hor;uniform vec3 sun;uniform float time;
float h2(vec2 p){return fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453);}
float n2(vec2 p){vec2 i=floor(p),f=fract(p);f=f*f*(3.0-2.0*f);return mix(mix(h2(i),h2(i+vec2(1,0)),f.x),mix(h2(i+vec2(0,1)),h2(i+vec2(1,1)),f.x),f.y);}
float fbm(vec2 p){float a=0.5,s=0.0;for(int i=0;i<5;i++){s+=a*n2(p);p=p*2.03+vec2(1.7,9.2);a*=0.5;}return s;}
void main(){vec3 n=normalize(d);float h=clamp(n.y,0.0,1.0);vec3 c=mix(hor,top,pow(h,0.5));
float s=max(dot(n,sun),0.0);c+=vec3(1.0,0.95,0.82)*(pow(s,900.0)*2.0+pow(s,14.0)*0.10);
if(n.y>0.0){vec2 uv=n.xz/(n.y+0.18)*1.6+vec2(time*0.004,0.0);float f=fbm(uv);float cl=smoothstep(0.50,0.78,f)*smoothstep(0.02,0.22,n.y);
vec3 cc=mix(vec3(0.80,0.84,0.90),vec3(1.0),smoothstep(0.55,0.9,f+0.25*s));c=mix(c,cc,cl*0.85);}
else c=hor;gl_FragColor=vec4(c,1.0);}`});
  const sky=new T.Mesh(g,m); sky.frustumCulled=false; sky.renderOrder=-10; scene.add(sky); A.sky=sky;
})();
const sun=new T.DirectionalLight(0xfff3dc,2.5); sun.position.copy(sunDir).multiplyScalar(1000); scene.add(sun);
scene.add(new T.HemisphereLight(0xd6e6ff,0x7a776a,1.45));

// ---------- terrain
const TR=M.terrain, terr={};
(function(){
  const [X0,Y0,X1,Y1]=TR.bbox; const cw=TR.cw,ch=TR.ch,cs=TR.coarse;
  const hc=buf('terr_coarse'), tiles=buf('terr_tiles'), hf=buf('terr_fine');
  const fn=TR.fn, ts=TR.tile, nx=Math.round((X1-X0)/ts), ny=Math.round((Y1-Y0)/ts);
  const fineIdx=new Int32Array(nx*ny).fill(-1);
  for(let k=0;k<TR.ntiles;k++) fineIdx[tiles[2*k+1]*nx+tiles[2*k]]=k;
  const pos=new Float32Array(cw*ch*3), uv=new Float32Array(cw*ch*2);
  for(let j=0;j<ch;j++)for(let i=0;i<cw;i++){const k=j*cw+i;const x=X0+i*cs,y=Y0+j*cs;pos[3*k]=x-OX;pos[3*k+1]=hc[k]/20+500;pos[3*k+2]=-(y-OY);uv[2*k]=i/(cw-1);uv[2*k+1]=j/(ch-1);}
  const per=Math.round(ts/cs); const idx=[];
  for(let j=0;j<ch-1;j++)for(let i=0;i<cw-1;i++){
    if(fineIdx[Math.floor(j/per)*nx+Math.floor(i/per)]>=0) continue;
    const a=j*cw+i,b=a+1,c=a+cw,d=c+1; idx.push(a,b,d,a,d,c);
  }
  const g=new T.BufferGeometry(); g.setAttribute('position',new T.BufferAttribute(pos,3)); g.setAttribute('uv',new T.BufferAttribute(uv,2)); g.setIndex(idx); g.computeVertexNormals();
  A.baseTex=loadTex('base');
  terr.matCoarse=new T.MeshLambertMaterial({map:A.baseTex,color:0xffffff});
  scene.add(new T.Mesh(g,terr.matCoarse));
  // fine tiles, textured from atlases
  const AT=M.atlas, n1=fn, fs=TR.fine;
  const mats=[]; for(let a=0;a<AT.n;a++) mats.push(new T.MeshLambertMaterial({map:loadTex('atlas'+a)}));
  const tidx=[]; for(let j=0;j<n1-1;j++)for(let i=0;i<n1-1;i++){const a=j*n1+i,b=a+1,c=a+n1,d=c+1;tidx.push(a,b,d,a,d,c);}
  const idxAttr=new T.BufferAttribute(new Uint16Array(tidx),1); const S=AT.cell*AT.per;
  for(let k=0;k<TR.ntiles;k++){
    const ti=tiles[2*k],tj=tiles[2*k+1]; const p=new Float32Array(n1*n1*3), u=new Float32Array(n1*n1*2);
    const q0=k%(AT.per*AT.per), r=Math.floor(q0/AT.per), c=q0%AT.per;
    for(let j=0;j<n1;j++)for(let i=0;i<n1;i++){const q=j*n1+i;p[3*q]=X0+ti*ts+i*fs-OX;p[3*q+1]=hf[k*n1*n1+q]/20+500;p[3*q+2]=-(Y0+tj*ts+j*fs-OY);
      const fx=i/(n1-1),fy=j/(n1-1);u[2*q]=(c*AT.cell+AT.gut+AT.tile*fx)/S;u[2*q+1]=1-(r*AT.cell+AT.gut+AT.tile*(1-fy))/S;}
    const tg=new T.BufferGeometry(); tg.setAttribute('position',new T.BufferAttribute(p,3)); tg.setAttribute('uv',new T.BufferAttribute(u,2)); tg.setIndex(idxAttr); tg.computeVertexNormals();
    scene.add(new T.Mesh(tg,mats[Math.floor(k/(AT.per*AT.per))]));
  }
  terr.height=function(sx,sz){
    const x=sx+OX,y=-sz+OY; if(x<X0||x>=X1||y<Y0||y>=Y1) return 620;
    const ti=Math.floor((x-X0)/ts),tj=Math.floor((y-Y0)/ts),k=fineIdx[tj*nx+ti];
    if(k>=0){const fx=(x-X0-ti*ts)/fs,fy=(y-Y0-tj*ts)/fs,i=Math.min(n1-2,Math.floor(fx)),j=Math.min(n1-2,Math.floor(fy)),u=fx-i,v=fy-j,b=k*n1*n1+j*n1+i;
      return ((hf[b]*(1-u)+hf[b+1]*u)*(1-v)+(hf[b+n1]*(1-u)+hf[b+n1+1]*u)*v)/20+500;}
    const fx=(x-X0)/cs,fy=(y-Y0)/cs,i=Math.min(cw-2,Math.floor(fx)),j=Math.min(ch-2,Math.floor(fy)),u=fx-i,v=fy-j,b=j*cw+i;
    return ((hc[b]*(1-u)+hc[b+1]*u)*(1-v)+(hc[b+cw]*(1-u)+hc[b+cw+1]*u)*v)/20+500;
  };
  terr.baseUV=(sx,sz)=>[(sx+OX-X0)/(X1-X0),(-sz+OY-Y0)/(Y1-Y0)];
  // distant mountains
  const F=buf('far'), fm=M.far;
  if(F){
    const w=fm.w,h=fm.h,res=fm.res; const fp=new Float32Array(w*h*3), fc=new Float32Array(w*h*3);
    for(let j=0;j<h;j++)for(let i=0;i<w;i++){const k=j*w+i,x=fm.x0+(i+0.5)*res,y=fm.y1-(j+0.5)*res,e=F[k]/10+200;
      const inside=(x>X0+150&&x<X1-150&&y>Y0+150&&y<Y1-150);
      fp[3*k]=x-OX;fp[3*k+1]=inside?e-60:e-4;fp[3*k+2]=-(y-OY);
      const t=Math.min(1,Math.max(0,(e-600)/900)), n=((i*7+j*13)%5)/5*0.03;
      fc[3*k]=0.34-0.06*t+n;fc[3*k+1]=0.37-0.02*t+n;fc[3*k+2]=0.27+0.06*t+n;}
    const fi=[]; for(let j=0;j<h-1;j++)for(let i=0;i<w-1;i++){
      const x=fm.x0+(i+0.5)*res,y=fm.y1-(j+0.5)*res; if(x>X0+res&&x+res<X1-res&&y<Y1-res&&y-res>Y0+res) continue;
      const a=j*w+i,b=a+1,c=a+w,d=c+1; fi.push(a,c,d,a,d,b);}
    const fg=new T.BufferGeometry(); fg.setAttribute('position',new T.BufferAttribute(fp,3)); fg.setAttribute('color',new T.BufferAttribute(fc,3)); fg.setIndex(fi); fg.computeVertexNormals();
    scene.add(new T.Mesh(fg,new T.MeshLambertMaterial({vertexColors:true})));
  }
})();

// ---------- roads, bridges, markings
function meshFrom(pn,inn,mat){
  const p=buf(pn), i=buf(inn); if(!p||!p.length) return null;
  const g=new T.BufferGeometry(); g.setAttribute('position',new T.BufferAttribute(p,3)); g.setIndex(new T.BufferAttribute(i,1)); g.computeVertexNormals();
  const m=new T.Mesh(g,mat); scene.add(m); return m;
}
(function(){
  const matAsphalt=new T.MeshLambertMaterial({color:0x4b4d51,polygonOffset:true,polygonOffsetFactor:-1,polygonOffsetUnits:-2});
  const matDeck=new T.MeshLambertMaterial({color:0xb4afa4});   // light concrete, as in the NCDOT renderings, so structures read apart from ground-level roads
  const matConc=new T.MeshLambertMaterial({color:0xcfcbc1,side:T.DoubleSide});
  const matWhite=new T.MeshBasicMaterial({color:0xecece4,polygonOffset:true,polygonOffsetFactor:-4,polygonOffsetUnits:-8});
  const matYellow=new T.MeshBasicMaterial({color:0xe2b21c,polygonOffset:true,polygonOffsetFactor:-4,polygonOffsetUnits:-8});
  meshFrom('road_g_p','road_g_i',matAsphalt);
  meshFrom('road_d_p','road_d_i',matDeck);
  meshFrom('deck_side_p','deck_side_i',matConc);
  meshFrom('deck_par_p','deck_par_i',matConc);
  meshFrom('conc_p','conc_i',matConc); meshFrom('wall_p','wall_i',new T.MeshLambertMaterial({color:0xa9a59b,side:T.DoubleSide}));
  meshFrom('fill_p','fill_i',new T.MeshLambertMaterial({color:0x76895a}));   // graded fill between adjacent roads at different levels
  meshFrom('mk_w_p','mk_w_i',matWhite); meshFrom('syn_w_p','syn_w_i',matWhite); meshFrom('syn_y_p','syn_y_i',matYellow); meshFrom('mk_y_p','mk_y_i',matYellow);
  const P=buf('piers'); if(P){ const n=P.length/4;
    const g=new T.CylinderGeometry(1.1,1.1,1,10); g.translate(0,0.5,0);
    const im=new T.InstancedMesh(g,matConc,n); const m=new T.Matrix4();
    for(let k=0;k<n;k++){m.makeScale(1,P[4*k+3],1);m.setPosition(P[4*k],P[4*k+1],P[4*k+2]);im.setMatrixAt(k,m);}
    scene.add(im);}
})();
