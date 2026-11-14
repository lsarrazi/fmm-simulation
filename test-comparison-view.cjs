const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
(async()=>{
const html=fs.readFileSync('comparison.html','utf8'),els={};let callbacks=[],wall=0;async function next(t){const current=callbacks;callbacks=[];for(const callback of current)await callback(t);}
const ctx=new Proxy({createImageData:(w,h)=>({data:new Uint8ClampedArray(w*h*4)})},{get:(o,k)=>o[k]??(()=>{})});
function el(id){return els[id]||={value:'',textContent:'',handlers:{},setAttribute(){},toggleAttribute(){},addEventListener(t,f){this.handlers[t]=f;},getContext(){return ctx;},fire(t){this.handlers[t]();}};}
for(const m of html.matchAll(/\bid="([^"]+)"/g))el(m[1]);
for(const m of html.matchAll(/<input\b([^>]+)>/g)){const id=/\bid="([^"]+)"/.exec(m[1])?.[1];if(id){el(id).value=/\bvalue="([^"]*)"/.exec(m[1])?.[1]||'';el(id).checked=/\bchecked\b/.test(m[1]);}}
for(const m of html.matchAll(/<select\b[^>]*id="([^"]+)"[^>]*>([\s\S]*?)<\/select>/g)){
 const opts=[...m[2].matchAll(/<option\s+value="([^"]+)"([^>]*)>/g)];el(m[1]).value=(opts.find(o=>/selected/.test(o[2]))||opts[0])[1];
}
const s={StellarObservationPolicy:require('./src/observation/stellar-budget'),PlanetarySystems:require('./planetary-systems'),StellarPopulation:require('./stellar-population'),LocalStellarDynamics:require('./local-stellar-dynamics'),...require('./galaxy-comparison.js'),...require('./observation-tools.js'),...require('./stellar-light.js'),PersistentMatter:require('./persistent-matter.js'),console,performance:{now:()=>wall},URLSearchParams,location:{search:''},requestAnimationFrame:f=>{callbacks.push(f);},document:{hidden:false,getElementById:id=>{assert(els[id],id);return els[id];},createElement:()=>({getContext:()=>ctx})}};s.window=s;vm.createContext(s);
el('speedUnit').value='My';el('stepUnit').value='My';el('engine').value='cpu';el('resolution').value='28,25';
el('scenario').value='galaxy';el('softening').value='0.042';
vm.runInContext(fs.readFileSync('comparison-view.js','utf8'),s,{filename:'comparison-view.js'});
for(let i=0;i<20;i++)await Promise.resolve();
assert.equal(el('time').value,'0.000');
el('starBudget').value='5003';el('starOnset').value='0.125';el('applyLightCalibration').fire('click');assert.equal(el('starBudget').value,'5003');assert.match(el('lightCalibrationStatus').textContent,/apparition ×0.125/);el('starBudget').value='12.5';el('applyLightCalibration').fire('click');assert.match(el('lightCalibrationStatus').textContent,/entier/);el('resetLightCalibration').fire('click');assert.equal(el('starBudget').value,'8192');assert.equal(el('starOnset').value,'160');
el('pause').fire('click');await next(wall+=20);assert.equal(el('time').value,'0.020');
el('pause').fire('click');await next(wall+=20);assert.equal(el('time').value,'0.020');
el('step').fire('click');for(let i=0;i<10;i++)await next(wall+=20);assert.equal(el('time').value,'0.120');await next(wall+=20);assert.equal(el('time').value,'0.120');
el('zoomValue').value='10000';el('zoomValue').fire('change');await next(wall+=40);assert.match(el('observationStatus').textContent,/Identités persistantes/);
el('speedUnit').value='s';el('speedUnit').fire('change');el('pause').fire('click');const before=Number.parseFloat(el('physicalTime').textContent);await next(wall+=100);assert.equal(el('time').value,'0.120');
el('pause').fire('click');
el('scenario').value='cold';el('scenario').fire('change');for(let i=0;i<20;i++)await Promise.resolve();assert.equal(el('time').value,'0.000');

console.log('Passed: single continuum view startup, fixed panel, play/pause, physical stepping and scenario reset.');
})().catch(error=>{console.error(error);process.exitCode=1;});
