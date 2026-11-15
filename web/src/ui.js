// ---------- UI
const ui=document.getElementById('ui');
const groups=[]; A.routes.forEach(r=>{let g=groups.find(x=>x.name===r.group); if(!g){g={name:r.group,items:[]};groups.push(g);} g.items.push(r);});
ui.innerHTML=`
<div id="panel">
  <div id="phead"><div><b>I-26 Connector</b><span>Asheville, NC &middot; 3D preview of the completed project</span></div><button id="collapse" title="Hide or show panel">&#9776;</button></div>
  <div id="pbody">
    <label class="lbl" for="route">Drive-through route</label>
    <select id="route">${groups.map(g=>`<optgroup label="${g.name}">${g.items.map(r=>`<option value="${r.id}">${r.name} (${(r.len/1609.34).toFixed(1)} mi)</option>`).join('')}</optgroup>`).join('')}</select>
    <div class="row"><button id="start" class="primary">&#9654; Start drive</button><button id="pause" disabled>Pause</button><button id="stop" disabled>Exit</button></div>
    <label class="lbl">Camera</label>
    <div class="row seg" id="cams"><button data-m="driver">Driver</button><button data-m="chase">Chase</button><button data-m="heli">Helicopter</button><button data-m="orbit" class="on">Free orbit</button></div>
    <label class="lbl">Playback speed</label>
    <div class="row seg" id="speeds"><button data-s="0.5">0.5&times;</button><button data-s="1" class="on">1&times;</button><button data-s="2">2&times;</button><button data-s="4">4&times;</button></div>
    <label class="lbl">Traffic <span id="trv">normal</span></label>
    <input type="range" id="traffic" min="0" max="2" step="0.25" value="1">
    <div class="row checks"><label><input type="checkbox" id="ckLabels" checked> Place labels</label><label><input type="checkbox" id="ckBld" checked> Buildings</label><label><input type="checkbox" id="ckTrees" checked> Trees</label></div>
    <label class="lbl">Jump to an overview</label>
    <div class="row wrap" id="views"></div>
    <div class="hint">Free orbit: drag to rotate, right-drag or Shift-drag to pan, scroll to zoom. Keys: Space pause, C camera, 1&ndash;4 speed.</div>
    <button id="about" class="link">About this model and its data sources</button>
  </div>
</div>
<div id="hud" class="hide"><div id="hudroute"></div><div id="hudrow"><span id="hudroad"></span><span id="hudspeed"></span></div><div id="hudbtns"><button id="hPause">Pause</button><button id="hCam">Camera: Driver</button><button id="hSpeed">Speed: 1&times;</button><button id="hRestart">Restart</button><button id="hExit">Exit drive</button></div><input type="range" id="scrub" min="0" max="1000" value="0"><div id="hudtime"></div></div>
<div id="labels"></div>
<div id="done" class="hide"><div class="box"><h2>Drive complete</h2><p id="donemsg"></p><div class="row"><button id="again" class="primary">Replay</button><button id="next">Next route</button><button id="free">Explore freely</button></div></div></div>
<div id="aboutbox" class="hide"><div class="box"><h2>About this model</h2>
<p>This is an unofficial visualization of NCDOT project I-2513 (the I-26 Connector) as it could look when complete. It is meant to help people picture the layout. It is not an engineering drawing and it is not an NCDOT product.</p>
<ul>
<li><b>Road layout:</b> traced from the NCDOT public meeting maps for the section between Haywood Road and Broadway (pavement, bridges, lane lines, barriers and alignments), and from the 2018 NCDOT design file for the I-26/I-240 section toward I-40.</li>
<li><b>Road heights and grades:</b> estimated for this model from ground elevations and bridge clearances. The plans used here do not include NCDOT's vertical profiles, so bridge heights and slopes are approximate.</li>
<li><b>Terrain:</b> Mapzen/AWS Terrain Tiles (USGS 3DEP source data).</li>
<li><b>Aerial imagery:</b> Maxar Open Data Program, image of 9 Nov 2019, &copy; Maxar, CC BY-NC 4.0. It shows the area before construction, so today's roads are visible beside and beneath the new ones.</li>
<li><b>Buildings, streets and business names:</b> Overture Maps Foundation, including &copy; OpenStreetMap contributors (ODbL). Building shapes are simple blocks; facades, storefronts and signs are generated, not photographed.</li>
<li><b>Beyond the plan limits:</b> I-240 east of the river, Patton Avenue toward downtown and West Asheville, and US 19/23 north are simplified from today's road centerlines.</li>
<li><b>Traffic, vehicles, trees, sky:</b> simulated for illustration.</li>
</ul><div class="row"><button id="aboutclose" class="primary">Close</button></div></div></div>`;
const $=id=>document.getElementById(id);
const VIEWS=[['Whole project',285300,210200,5200,0.35,0.95],['Patton Ave interchange',285380,209900,900,0.5,1.0],['River crossings',285780,210700,1100,-0.9,1.05],['Haywood Road',284930,208650,800,0.3,1.0],['Downtown side',286700,210420,1300,-2.2,1.05],['North end',285750,212000,1500,2.9,1.05]];
VIEWS.forEach(v=>{const b=document.createElement('button'); b.textContent=v[0]; b.onclick=()=>{const p=SP(v[1],v[2]); A.follow=false; view.tx=p[0]; view.tz=p[1]; view.ty=terr.height(p[0],p[1]); view.dist=v[3]; view.az=v[4]; view.pol=v[5]; A.setCamMode('orbit');}; $('views').appendChild(b);});
A.setCamMode=function(m){ if(m!=='orbit'&&!ego.active) m='orbit'; if(m==='orbit'&&A.camMode!=='orbit'&&ego.active){ A.follow=true; view.dist=Math.min(Math.max(view.dist,40),220); view.pol=1.15; view.az=Math.atan2(-ego.fwd.x,-ego.fwd.z); }
  A.camMode=m; ego.camInit=true; $('hCam').textContent='Camera: '+({driver:'Driver',chase:'Chase',heli:'Helicopter',orbit:'Free orbit'})[m]; document.querySelectorAll('#cams button').forEach(b=>b.classList.toggle('on',b.dataset.m===m)); };
document.querySelectorAll('#cams button').forEach(b=>b.onclick=()=>{ if(b.dataset.m!=='orbit'&&!ego.active) begin(); A.setCamMode(b.dataset.m); });
function setSpeed(s){A.timeScale=s; $('hSpeed').innerHTML='Speed: '+s+'&times;'; document.querySelectorAll('#speeds button').forEach(b=>b.classList.toggle('on',+b.dataset.s===s));}
document.querySelectorAll('#speeds button').forEach(b=>b.onclick=()=>setSpeed(+b.dataset.s));
function begin(){ const id=$('route').value; const keep=A.camMode; A.startRoute(id); if(keep!=='orbit') A.setCamMode(keep); $('hud').classList.remove('hide'); $('done').classList.add('hide'); $('pause').disabled=false; $('stop').disabled=false; $('pause').textContent=$('hPause').textContent='Pause'; $('start').innerHTML='&#8635; Restart'; if(!A.keepPanel) $('panel').classList.add('min'); $('hudroute').textContent=ego.route.name; }
$('start').onclick=begin;
$('route').onchange=()=>{ if(ego.active) begin(); else { const r=A.routes.find(x=>x.id===$('route').value), i=Math.floor(r.n/2)*5; A.follow=false; view.tx=r.P[i]; view.ty=r.P[i+1]; view.tz=r.P[i+2]; view.dist=Math.max(1500,r.len*0.55); view.pol=0.8; } };
function togglePause(){ if(!ego.active) return; if(ego.done){begin();return;} ego.paused=!ego.paused; $('pause').textContent=$('hPause').textContent=ego.paused?'Resume':'Pause'; }
$('pause').onclick=togglePause;
function exitDrive(){ A.stopRoute(); A.setCamMode('orbit'); A.follow=false; $('hud').classList.add('hide'); $('done').classList.add('hide'); $('pause').disabled=true; $('stop').disabled=true; $('start').innerHTML='&#9654; Start drive'; }
$('stop').onclick=exitDrive; $('hExit').onclick=()=>{exitDrive(); $('panel').classList.remove('min');}; $('hPause').onclick=togglePause; $('hRestart').onclick=begin;
$('hCam').onclick=()=>{const o=['driver','chase','heli','orbit']; A.setCamMode(o[(o.indexOf(A.camMode)+1)%4]);};
$('hSpeed').onclick=()=>{const o=[0.5,1,2,4]; setSpeed(o[(o.indexOf(A.timeScale)+1)%4]);};
$('traffic').oninput=e=>{A.traffic.density=+e.target.value; $('trv').textContent=A.traffic.density===0?'none':A.traffic.density<0.8?'light':A.traffic.density<1.3?'normal':'heavy'; spawnTraffic();};
$('ckBld').onchange=e=>{ if(A.bld) A.bld.visible=e.target.checked; if(A.store) A.store.visible=e.target.checked; };
$('ckTrees').onchange=e=>{ if(A.trees) A.trees.visible=e.target.checked; };
$('ckLabels').onchange=e=>{ $('labels').style.display=e.target.checked?'':'none'; };
$('collapse').onclick=()=>$('panel').classList.toggle('min');
$('about').onclick=()=>$('aboutbox').classList.remove('hide'); $('aboutclose').onclick=()=>$('aboutbox').classList.add('hide');
$('scrub').oninput=e=>{ if(ego.route){ ego.s=Math.max(2,e.target.value/1000*(ego.route.len-4)); ego.done=false; ego.camInit=true; $('done').classList.add('hide'); } };
A.onDone=()=>{ $('donemsg').textContent='You reached the end of: '+ego.route.name+'.'; $('done').classList.remove('hide'); $('pause').textContent='Replay'; };
$('again').onclick=begin;
$('next').onclick=()=>{ const s=$('route'); s.selectedIndex=(s.selectedIndex+1)%s.options.length; begin(); };
$('free').onclick=()=>{ $('done').classList.add('hide'); A.setCamMode('orbit'); };
window.addEventListener('keydown',e=>{ if(e.target.tagName==='SELECT'||e.target.tagName==='INPUT') return;
  if(e.code==='Space'){e.preventDefault(); if(ego.active) togglePause(); else begin();}
  else if(e.key==='c'||e.key==='C'){const o=['driver','chase','heli','orbit']; A.setCamMode(o[(o.indexOf(A.camMode)+1)%4]);}
  else if(e.key>='1'&&e.key<='4') setSpeed([0.5,1,2,4][+e.key-1]); });
// ---------- place labels
const LABELS=(M.labels||[]).map(l=>{const p=SP(l.x,l.y), el=document.createElement('div'); el.className='plabel '+(l.k||''); el.textContent=l.t; $('labels').appendChild(el); return {el,x:p[0],z:p[1],y:terr.height(p[0],p[1])+(l.h||40),t:l.t,r:l.r||9000};});
const _lp=new T.Vector3();
function updateLabels(){
  const w=canvas.clientWidth,h=canvas.clientHeight, driving=ego.active&&A.camMode!=='orbit'&&A.camMode!=='manual', placed=[];
  for(const l of LABELS){ _lp.set(l.x,l.y,l.z); const d=_lp.distanceTo(camera.position); _lp.project(camera);
    const rr=driving?Math.min(l.r,l.r>15000?l.r:1500):l.r;
    if(_lp.z>1||_lp.z<-1||d>rr||Math.abs(_lp.x)>1.05||Math.abs(_lp.y)>1.05){l.el.style.display='none';continue;}
    const sx=(_lp.x+1)/2*w, sy=(1-_lp.y)/2*h, hw=l.t.length*3.9+12; let clash=false;
    for(const q of placed){ if(Math.abs(q[0]-sx)<q[2]+hw&&Math.abs(q[1]-sy)<26){clash=true;break;} }
    if(clash){l.el.style.display='none';continue;} placed.push([sx,sy,hw]);
    l.el.style.display=''; l.el.style.transform=`translate(-50%,-100%) translate(${sx.toFixed(1)}px,${(sy-14).toFixed(1)}px)`; l.el.style.opacity=d>rr*0.8?0.5:1; }
}
function segLabel(r,s){ const i=s/DS; let t=''; for(const sg of r.segs){ if(sg[0]<=i) t=sg[1]; else break; } return t; }
function fmt(t){t=Math.max(0,Math.round(t));return Math.floor(t/60)+':'+String(t%60).padStart(2,'0');}
let hudT=0;
function updateHud(dt){ hudT+=dt; if(hudT<0.15||!ego.active) return; hudT=0; const r=ego.route;
  $('hudroad').textContent=segLabel(r,ego.s); $('hudspeed').textContent=Math.round(ego.v*2.23694)+' mph';
  if(document.activeElement!==$('scrub')) $('scrub').value=Math.round(ego.s/(r.len-4)*1000);
  $('hudtime').textContent=(ego.s/1609.34).toFixed(1)+' of '+(r.len/1609.34).toFixed(1)+' mi'+(A.timeScale!==1?'  ·  '+A.timeScale+'× speed':''); }
// ---------- main loop
let last=performance.now(), tsky=0;
function frame(now){
  let dt=Math.min(0.1,(now-last)/1000); last=now; const sdt=dt*A.timeScale;
  tsky+=dt; A.sky.material.uniforms.time.value=tsky;
  if(ego.active) updateEgo(sdt);
  const steps=sdt>0.06?Math.ceil(sdt/0.05):1; for(let k=0;k<steps;k++) stepTraffic(sdt/steps);
  drawTraffic(); updateCamera(dt); A.sky.position.copy(camera.position);
  renderer.render(scene,camera); updateLabels(); updateHud(dt);
  requestAnimationFrame(frame);
}
spawnTraffic();
(function ready(){ if(pendingImgs>0){ setTimeout(ready,60); return; }
  for(let k=0;k<200;k++) stepTraffic(0.1);
  { const p=SP(285350,210150); view.tx=p[0]; view.tz=p[1]; view.ty=terr.height(p[0],p[1]); view.dist=3400; view.az=0.45; view.pol=0.98; }
  document.getElementById('loading').style.display='none'; requestAnimationFrame(frame); window.READY=true; })();
window.APP=Object.assign(A,{scene,camera,terr,renderer,sp:SP,
  setCam(px,py,pz,tx,ty,tz){A.camMode='manual';camera.position.set(px,py,pz);camera.lookAt(tx,ty,tz);}});
const _uc=updateCamera; updateCamera=function(dt){ if(A.camMode==='manual'){ ego.mesh.visible=ego.active; ego.hood.visible=false; return; } _uc(dt); };
