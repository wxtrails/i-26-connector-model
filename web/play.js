const { chromium } = require('playwright');
(async()=>{
  const b=await chromium.launch({args:['--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader','--ignore-gpu-blocklist']});
  const p=await b.newPage({viewport:{width:800,height:450},deviceScaleFactor:1});
  const errs=[]; p.on('console',m=>{if(m.type()==='error')errs.push(m.text());}); p.on('pageerror',e=>errs.push('PAGEERR '+e.message));
  await p.goto('file://'+process.argv[2]); await p.waitForFunction('window.READY===true',{timeout:180000});
  const ids=await p.evaluate(()=>APP.routes.map(r=>[r.id,r.len,r.name]));
  console.log('routes',ids.length);
  await p.click('#start'); await p.waitForTimeout(4000);
  let st=await p.evaluate(()=>({s:APP.ego.s,v:APP.ego.v,mode:APP.camMode,hud:document.getElementById('hudroad').textContent,sp:document.getElementById('hudspeed').textContent,cars:APP.routes.length}));
  console.log('after start',JSON.stringify(st));
  await p.click('#hCam'); await p.click('#hSpeed'); await p.click('#hPause'); await p.waitForTimeout(1500);
  st=await p.evaluate(()=>({paused:APP.ego.paused,mode:APP.camMode,ts:APP.timeScale})); console.log('controls',JSON.stringify(st));
  await p.click('#hPause');
  for(const [id,len] of ids){
    await p.evaluate(id=>{document.getElementById('route').value=id; document.getElementById('route').dispatchEvent(new Event('change')); APP.timeScale=4;},id);
    await p.waitForTimeout(1200);
    const s1=await p.evaluate(()=>APP.ego.s);
    await p.evaluate(()=>{APP.ego.s=APP.ego.route.len-60;});
    await p.waitForFunction('APP.ego.done===true',{timeout:120000}).catch(()=>console.log('  not done',id));
    const r=await p.evaluate(()=>({done:APP.ego.done,box:!document.getElementById('done').classList.contains('hide')}));
    console.log(id,'moved',s1.toFixed(0),JSON.stringify(r));
  }
  await p.click('#free'); await p.click('#hExit'); await p.waitForTimeout(800);
  console.log('exit',await p.evaluate(()=>({active:APP.ego.active,mode:APP.camMode})));
  console.log('errors',errs.length,errs.slice(0,8).join(' | '));
  await b.close();
})();
