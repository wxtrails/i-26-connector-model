// usage: node shot.js out.png "px,py,pz,tx,ty,tz" [w h] ; coordinates in state plane x,y and elevation z: "x,y,z,tx,ty,tz"
const { chromium } = require('playwright');
(async()=>{
  const args=process.argv.slice(2);
  const html=args[0];
  const shots=JSON.parse(args[1]);
  const b=await chromium.launch({args:['--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader','--ignore-gpu-blocklist']});
  const p=await b.newPage({viewport:{width:+(args[2]||1280),height:+(args[3]||720)},deviceScaleFactor:1});
  const errs=[]; p.on('console',m=>{if(m.type()==='error'||m.type()==='warning')errs.push(m.type()+': '+m.text());}); p.on('pageerror',e=>errs.push('PAGEERR '+e.message));
  await p.goto('file://'+html);
  try{ await p.waitForFunction('window.READY===true',{timeout:120000}); }catch(e){ console.log('not ready',errs.join('\n')); await b.close(); process.exit(1); }
  for(const s of shots){
    await p.evaluate((s)=>{ if(s.js){ (0,eval)(s.js); } if(s.cam){const c=s.cam; const a=APP.sp(c[0],c[1]), t=APP.sp(c[3],c[4]); APP.setCam(a[0],c[2],a[1],t[0],c[5],t[1]);} },s);
    await p.waitForTimeout(s.wait||1200);
    await p.screenshot({path:s.out,timeout:240000});
  }
  if(errs.length) console.log(errs.slice(0,15).join('\n'));
  await b.close();
})();
