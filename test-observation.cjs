const assert=require('node:assert/strict');
const{SimulationClock,MODEL_NS,physicalNs,localSample}=require('./observation-tools');
const clock=new SimulationClock();clock.setRate('1','s');for(let i=0;i<100;i++)clock.tick(10);assert.equal(clock.visible(),1000000000n);assert.equal(clock.duration(),0);assert.equal(clock.format('s'),'1 s');
clock.integrated=physicalNs('1000000000','y');clock.debt=0n;clock.tick(1000);assert.equal(clock.visible(),physicalNs('1000000000','y')+1000000000n,'seconds remain visible after a billion years');
clock.setRate('.5'.replace(/^\./,'0.'),'y');assert.equal(clock.rate,physicalNs('0.5','y'));assert.throws(()=>clock.setRate('-1','s'));
clock.setRate('1','Gy');clock.tick(1000);assert(clock.limited);assert(clock.duration()<=.025);clock.consume(.01);assert.equal(clock.integrated,physicalNs('1000000000','y')+MODEL_NS/100n);
assert.equal(physicalNs('1e-3','y'),31557600000000n);
const n=64,g={n,nv:3,extent:1.5,vmax:2,dx:3/n,dv:4/3,mass:new Float64Array(n*n).fill(1/(n*n))};
const o={x:.2,y:.1,zoom:1e10,budget:128,population:10n**21n};
const a=localSample(g,o),b=localSample(g,o);assert.deepEqual(a,b,'stable at a fixed snapshot');assert(a.points.length>0&&a.points.length<30,'only a few stars at deep zoom');assert.equal(new Set(a.points.map(p=>p.id)).size,a.points.length);
const wider=localSample(g,{...o,zoom:5e9,budget:512});for(const p of a.points){const match=wider.points.find(q=>q.id===p.id);assert(match,'same address when zooming out');assert.equal(match.x,p.x);assert.equal(match.y,p.y);}
const global=localSample(g,{zoom:1,budget:128,population:10n**24n});assert.equal(global.points.length,128);assert(global.points.some(p=>p.x<0)&&global.points.some(p=>p.x>0),'bounded sample covers both sides');
const phases=new Map();for(let i=0;i<n*n;i++){const v=new Float32Array(9);v[5]=g.mass[i];phases.set(i,v);}
const before=localSample(g,{...o,phases}),dt=1e-8,v=g.dv;
const moved=localSample(g,{...o,x:o.x+v*dt,phases,offset:dt});assert.equal(before.points.length,moved.points.length);for(const p of before.points){const q=moved.points.find(q=>q.id===p.id);assert(q);assert(Math.abs(q.x-p.x-v*dt)<1e-14,'ballistic local interpolation');}
console.log('Passed: precise physical clock, large-time seconds, rate units, high-rate cap, bounded local population, deterministic deep zoom and ballistic interpolation.');
console.log({deepVisible:a.points.length,deepPopulation:a.estimated.toString(),visited:a.visited,globalDisplayed:global.points.length});
