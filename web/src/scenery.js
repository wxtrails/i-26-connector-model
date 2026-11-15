// ---------- procedural textures
function canvasTex(w,h,draw,rep){
  const c=document.createElement('canvas'); c.width=w; c.height=h; draw(c.getContext('2d'),w,h);
  const t=new T.CanvasTexture(c); t.colorSpace=T.SRGBColorSpace; t.anisotropy=8;
  if(rep){t.wrapS=t.wrapT=T.RepeatWrapping;} return t;
}
function rnd(seed){let s=seed>>>0;return()=>{s=(s*1664525+1013904223)>>>0;return s/4294967296;};}
const FACADES=[
 // 0 house siding: 8 m x 6.4 m tile
 (g,w,h)=>{const r=rnd(3);g.fillStyle='#d9d4c8';g.fillRect(0,0,w,h);g.strokeStyle='rgba(0,0,0,0.10)';for(let y=0;y<h;y+=8){g.beginPath();g.moveTo(0,y);g.lineTo(w,y);g.stroke();}
   for(let f=0;f<2;f++)for(let b=0;b<2;b++){const x=b*w/2+w*0.14,y=f*h/2+h*0.11;g.fillStyle='#f4f2ec';g.fillRect(x-4,y-4,w*0.22+8,h*0.24+8);g.fillStyle='#2f3d4c';g.fillRect(x,y,w*0.22,h*0.24);g.fillStyle='rgba(255,255,255,0.18)';g.fillRect(x,y,w*0.22,h*0.09);g.fillStyle='#f4f2ec';g.fillRect(x+w*0.105,y,3,h*0.24);}},
 // 1 brick commercial
 (g,w,h)=>{g.fillStyle='#a0583f';g.fillRect(0,0,w,h);for(let y=0;y<h;y+=6){for(let x=(y/6%2)*8;x<w;x+=16){g.fillStyle='rgba(0,0,0,'+(0.05+0.08*((x*7+y*3)%5)/5)+')';g.fillRect(x,y,15,5);}}
   for(let f=0;f<2;f++)for(let b=0;b<2;b++){const x=b*w/2+w*0.10,y=f*h/2+h*0.10;g.fillStyle='#e6dfd1';g.fillRect(x-4,y-5,w*0.30+8,h*0.28+10);g.fillStyle='#243240';g.fillRect(x,y,w*0.30,h*0.28);g.fillStyle='rgba(180,210,235,0.25)';g.fillRect(x,y,w*0.30,h*0.11);g.fillStyle='#e6dfd1';g.fillRect(x+w*0.148,y,3,h*0.28);}},
 // 2 office glass grid
 (g,w,h)=>{g.fillStyle='#8f979d';g.fillRect(0,0,w,h);for(let f=0;f<2;f++)for(let b=0;b<4;b++){const x=b*w/4+5,y=f*h/2+12;const gr=g.createLinearGradient(0,y,0,y+h/2-30);gr.addColorStop(0,'#7fa6c4');gr.addColorStop(1,'#34506a');g.fillStyle=gr;g.fillRect(x,y,w/4-10,h/2-30);}},
 // 3 warehouse / big box
 (g,w,h)=>{g.fillStyle='#cfcabd';g.fillRect(0,0,w,h);g.strokeStyle='rgba(0,0,0,0.07)';for(let x=0;x<w;x+=10){g.beginPath();g.moveTo(x,0);g.lineTo(x,h);g.stroke();}g.fillStyle='#8a2f2a';g.fillRect(0,0,w,h*0.06);g.fillStyle='#8d8f90';g.fillRect(w*0.62,h*0.62,w*0.2,h*0.38);}
];
A.facadeTex=FACADES.map(f=>canvasTex(256,256,f,true));
// ---------- buildings
(function(){
  const xy=buf('bld_xy'), cnt=buf('bld_n'), hh=buf('bld_h'), ty=buf('bld_t'); if(!xy) return;
  const W=[0,1,2,3].map(()=>({p:[],u:[],c:[],i:[],n:0})), R={p:[],u:[],i:[],n:0};
  const tints=[[1,1,1],[0.93,0.90,0.84],[0.84,0.88,0.93],[0.98,0.94,0.86],[0.80,0.80,0.78],[0.93,0.84,0.78],[0.86,0.92,0.86]];
  A.bldBase=new Float32Array(cnt.length); A.bldTop=new Float32Array(cnt.length);
  let o=0;
  for(let b=0;b<cnt.length;b++){
    const n=cnt[b], t=ty[b], h=hh[b]; let lo=1e9, hi=-1e9;
    for(let k=0;k<n;k++){const e=terr.height(xy[2*(o+k)],xy[2*(o+k)+1]); if(e<lo)lo=e; if(e>hi)hi=e;}
    const top=hi+h, base=lo-1.5; A.bldBase[b]=lo; A.bldTop[b]=top;
    const w=W[t], tint=t===0?tints[(b*7)%tints.length]:(t===1?tints[(b*5)%4]:[1,1,1]); let per=0;
    for(let k=0;k<n;k++){
      const x0=xy[2*(o+k)],z0=xy[2*(o+k)+1],x1=xy[2*(o+(k+1)%n)],z1=xy[2*(o+(k+1)%n)+1]; const L=Math.hypot(x1-x0,z1-z0);
      const u0=per/8,u1=(per+L)/8,v1=(top-lo)/6.4; per+=L;
      w.p.push(x0,base,z0,x1,base,z1,x1,top,z1,x0,top,z0); w.u.push(u0,-1.5/6.4,u1,-1.5/6.4,u1,v1,u0,v1);
      for(let q=0;q<4;q++) w.c.push(tint[0],tint[1],tint[2]);
      w.i.push(w.n,w.n+1,w.n+2,w.n,w.n+2,w.n+3); w.n+=4;
    }
    const contour=[]; for(let k=0;k<n;k++) contour.push(new T.Vector2(xy[2*(o+k)],xy[2*(o+k)+1]));
    const tris=T.ShapeUtils.triangulateShape(contour,[]);
    for(let k=0;k<n;k++){const x=xy[2*(o+k)],z=xy[2*(o+k)+1]; R.p.push(x,top,z); const uv=terr.baseUV(x,z); R.u.push(uv[0],uv[1]);}
    for(const tr of tris) R.i.push(R.n+tr[0],R.n+tr[2],R.n+tr[1]);
    R.n+=n; o+=n;
  }
  A.bld=new T.Group(); scene.add(A.bld);
  W.forEach((w,t)=>{ if(!w.n) return; const g=new T.BufferGeometry();
    g.setAttribute('position',new T.Float32BufferAttribute(w.p,3)); g.setAttribute('uv',new T.Float32BufferAttribute(w.u,2)); g.setAttribute('color',new T.Float32BufferAttribute(w.c,3)); g.setIndex(w.i); g.computeVertexNormals();
    A.bld.add(new T.Mesh(g,new T.MeshLambertMaterial({map:A.facadeTex[t],vertexColors:true,side:T.DoubleSide})));});
  const g=new T.BufferGeometry(); g.setAttribute('position',new T.Float32BufferAttribute(R.p,3)); g.setAttribute('uv',new T.Float32BufferAttribute(R.u,2)); g.setIndex(R.i); g.computeVertexNormals();
  A.bld.add(new T.Mesh(g,new T.MeshLambertMaterial({map:A.baseTex,color:0xdcdcdc,side:T.DoubleSide})));
})();
// ---------- storefronts and signs
(function(){
  const S=M.signs||[]; if(!S.length) return;
  const CW=512,CH=64,COLS=4,ROWS=Math.ceil(S.length/COLS); const c=document.createElement('canvas'); c.width=CW*COLS; c.height=Math.max(64,CH*ROWS); const g=c.getContext('2d');
  const pal=[['#1d3557','#ffffff'],['#9b2226','#fff4d6'],['#f4f1de','#22333b'],['#264653','#e9c46a'],['#fefae0','#bc4749'],['#2b2d42','#edf2f4'],['#386641','#f2e8cf'],['#ffb703','#1b1b1b'],['#ffffff','#0b4f9c'],['#3a0ca3','#ffffff']];
  const sp=[],su=[],si=[]; let sn=0; const fp=[],fu=[],fi=[]; let fnn=0;
  S.forEach((s,k)=>{
    const col=k%COLS,row=Math.floor(k/COLS),x=col*CW,y=row*CH,p=pal[(k*7+s.t.length)%pal.length];
    g.fillStyle=p[0]; g.fillRect(x,y,CW,CH); g.strokeStyle=p[1]; g.lineWidth=3; g.strokeRect(x+4,y+4,CW-8,CH-8);
    g.fillStyle=p[1]; let fs=40; g.textAlign='center'; g.textBaseline='middle';
    do{g.font='700 '+fs+'px "Segoe UI",Arial,sans-serif'; fs-=2;}while(g.measureText(s.t).width>CW-36&&fs>12);
    g.fillText(s.t,x+CW/2,y+CH/2+2);
    const base=terr.height(s.x,s.z), w=s.w, hw=w/2, y0=base+3.25, y1=base+3.25+w*CH/CW*1.0+0.25;
    const ox=s.nx*0.28, oz=s.nz*0.28;
    const ax=s.x-s.ex*hw+ox, az=s.z-s.ez*hw+oz, bx=s.x+s.ex*hw+ox, bz=s.z+s.ez*hw+oz;
    sp.push(ax,y0,az,bx,y0,bz,bx,y1,bz,ax,y1,az);
    const u0=x/c.width,u1=(x+CW)/c.width,v0=1-(y+CH)/c.height,v1=1-y/c.height;
    // text must read left-to-right when seen from outside: viewer looks along -n, their right is n x up
    const flip=(s.ex*s.nz-s.ez*s.nx)<0;
    if(flip) su.push(u1,v0,u0,v0,u0,v1,u1,v1); else su.push(u0,v0,u1,v0,u1,v1,u0,v1);
    si.push(sn,sn+1,sn+2,sn,sn+2,sn+3); sn+=4;
    const fx=s.nx*0.14, fz=s.nz*0.14, e=hw+0.4, kk=k%4;
    fp.push(s.x-s.ex*e+fx,base-0.6,s.z-s.ez*e+fz, s.x+s.ex*e+fx,base-0.6,s.z+s.ez*e+fz, s.x+s.ex*e+fx,base+3.15,s.z+s.ez*e+fz, s.x-s.ex*e+fx,base+3.15,s.z-s.ez*e+fz);
    fu.push(kk/4,0,(kk+1)/4,0,(kk+1)/4,1,kk/4,1); fi.push(fnn,fnn+1,fnn+2,fnn,fnn+2,fnn+3); fnn+=4;
  });
  const tex=new T.CanvasTexture(c); tex.colorSpace=T.SRGBColorSpace; tex.anisotropy=8;
  const sg=new T.BufferGeometry(); sg.setAttribute('position',new T.Float32BufferAttribute(sp,3)); sg.setAttribute('uv',new T.Float32BufferAttribute(su,2)); sg.setIndex(si);
  // storefront strip: 4 variants (awning colour), glass, door
  const ft=canvasTex(1024,256,(g,w,h)=>{const aw=['#b23a48','#2a6f97','#3a5a40','#e09f3e'];
    for(let k=0;k<4;k++){const x=k*256;g.fillStyle='#d8d2c4';g.fillRect(x,0,256,h);g.fillStyle='#1d2a35';g.fillRect(x+10,70,236,150);
      const gr=g.createLinearGradient(0,70,0,220);gr.addColorStop(0,'rgba(170,205,230,0.55)');gr.addColorStop(1,'rgba(60,90,110,0.25)');g.fillStyle=gr;g.fillRect(x+10,70,236,150);
      g.fillStyle='#d8d2c4';for(let m=1;m<4;m++)g.fillRect(x+10+m*59-2,70,4,150);
      g.fillStyle='#3b2f2f';g.fillRect(x+108,96,44,160);g.fillStyle='rgba(170,205,230,0.5)';g.fillRect(x+113,102,34,90);
      g.fillStyle=aw[k];g.beginPath();g.moveTo(x+2,22);g.lineTo(x+254,22);g.lineTo(x+250,66);g.lineTo(x+6,66);g.closePath();g.fill();
      g.fillStyle='rgba(255,255,255,0.35)';for(let m=0;m<8;m++)g.fillRect(x+6+m*31,22,14,44);
      g.fillStyle='#6d6a62';g.fillRect(x,236,256,20);}
  },false);
  const fg=new T.BufferGeometry(); fg.setAttribute('position',new T.Float32BufferAttribute(fp,3)); fg.setAttribute('uv',new T.Float32BufferAttribute(fu,2)); fg.setIndex(fi); fg.computeVertexNormals();
  A.store=new T.Group(); scene.add(A.store);
  A.store.add(new T.Mesh(fg,new T.MeshLambertMaterial({map:ft,side:T.DoubleSide})));
  A.store.add(new T.Mesh(sg,new T.MeshBasicMaterial({map:tex,side:T.DoubleSide})));
})();
// ---------- trees
(function(){
  const TRE=buf('trees'); if(!TRE) return; const n=TRE.length/4;
  const crown=new T.IcosahedronGeometry(1,1), cone=new T.ConeGeometry(1,1,7), trunk=new T.CylinderGeometry(0.5,0.7,1,5); trunk.translate(0,0.5,0); cone.translate(0,0.5,0);
  let nd=0,ne=0; for(let k=0;k<n;k++){ if(TRE[4*k+3]>0.5) ne++; else nd++; }
  const md=new T.InstancedMesh(crown,new T.MeshLambertMaterial({flatShading:true}),nd), me=new T.InstancedMesh(cone,new T.MeshLambertMaterial({flatShading:true}),ne), mt=new T.InstancedMesh(trunk,new T.MeshLambertMaterial({color:0x5b4a3a}),n);
  const m=new T.Matrix4(), col=new T.Color(); let id=0,ie=0; const r=rnd(11);
  for(let k=0;k<n;k++){
    const x=TRE[4*k],z=TRE[4*k+1],h=TRE[4*k+2],y=terr.height(x,z)-0.3;
    m.makeScale(0.22,h*0.45,0.22); m.setPosition(x,y,z); mt.setMatrixAt(k,m);
    if(TRE[4*k+3]>0.5){ m.makeScale(h*0.2,h*0.85,h*0.2); m.setPosition(x,y+h*0.15,z); me.setMatrixAt(ie,m); col.setRGB(0.10+0.05*r(),0.26+0.08*r(),0.13+0.04*r()); me.setColorAt(ie++,col); }
    else{ const s=h*(0.30+0.1*r()); m.makeScale(s,h*0.36,s); m.setPosition(x,y+h*0.66,z); md.setMatrixAt(id,m); col.setRGB(0.20+0.12*r(),0.36+0.14*r(),0.12+0.06*r()); md.setColorAt(id++,col); }
  }
  A.trees=new T.Group(); A.trees.add(md,me,mt); scene.add(A.trees);
})();
