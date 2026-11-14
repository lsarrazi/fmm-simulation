const assert=require('node:assert/strict'),PersistentMatter=require('./persistent-matter'),{SimulationClock}=require('./observation-tools');
const grid={n:4,nv:3,dx:.75,dv:4/3,extent:1.5,initialPhase:[[(5*3+1)*3+2,1]]};
const zero=new Float64Array(16),ax=new Float64Array(16).fill(.2),ay=new Float64Array(16).fill(-.3);
const m=new PersistentMatter(grid),clone=new PersistentMatter(grid),original=m.evaluate(0,123456789n),id=original.id;
assert.deepEqual(m.resolve(id),original);const bytes=m.data.byteLength+m.forces.byteLength;
let view=m.sample({x:original.x,y:original.y,zoom:1e11,budget:128});assert(view.points.some(p=>p.id===id),'address discoverable locally without global enumeration');
for(let i=0;i<10;i++){m.sample({x:original.x,y:original.y,zoom:1e10,budget:64});m.step(.001,zero,zero);clone.step(.001,zero,zero);}
assert.deepEqual(m.data,clone.data,'observation order does not affect trajectories');
const moved=m.resolve(id);assert.equal(moved.id,id);assert(Math.abs(moved.x-original.x-original.vx*.01)<1e-14);
assert(m.sample({x:moved.x,y:moved.y,zoom:1e11,budget:128}).points.some(p=>p.id===id),'same star after movement and leaving old viewport');
m.setField(ax,ay);const prediction=m.resolve(id,.01);m.step(.01,ax,ay);const actual=m.resolve(id);assert(Math.abs(prediction.x-actual.x)<1e-14&&Math.abs(prediction.vx-actual.vx)<1e-14,'fine-time trajectory continuous at a collective step');
const huge=new PersistentMatter(grid,{population:10n**24n});assert.equal(huge.data.byteLength+huge.forces.byteLength,bytes,'memory independent of virtual population');assert.equal(huge.groups.reduce((s,g)=>s+g.count,0n),10n**24n);
const deform=new PersistentMatter(grid),old=deform.resolve(id),linearX=new Float64Array(16),linearY=new Float64Array(16);for(let i=0;i<16;i++){linearX[i]=-((i%4+.5)*grid.dx-grid.extent);linearY[i]=-(Math.floor(i/4)+.5)*grid.dx+grid.extent;}
deform.step(.01,linearX,linearY);const next=deform.resolve(id);assert(Math.abs(next.x-(old.x+old.vx*.01-.5*old.x*.01**2))<1e-14,'affine deformation follows linear gravitational field');assert(Math.abs(next.vx-(old.vx-old.x*.01))<1e-14);
function motion(unit){const clock=new SimulationClock(),material=new PersistentMatter(grid);clock.setRate('1',unit);const p=material.resolve(id);for(let i=0;i<100;i++){clock.tick(10);const dt=clock.duration();if(dt){material.step(dt,zero,zero);clock.consume(dt);}}return material.resolve(id,clock.fineOffset()).x-p.x;}
const fast=motion('My'),year=motion('y'),second=motion('s');assert(Math.abs(fast/year/1e6-1)<1e-6,'1 y/s slows actual displacement by one million');assert(Math.abs(year/second/31557600-1)<.01,'1 s/s slows actual displacement again by one year in seconds');
console.log('Passed: immutable identity, lazy local discovery, deterministic moving reconstruction, independence from observation, fine-time continuity, population-independent memory and actual motion scaling.');
console.log({identities:m.population.toString(),cohorts:m.groups.length,bytes,displacementMy:fast,displacementYear:year,displacementSecond:second});

const calibrationMatter=new PersistentMatter(grid,{population:1000000n}),largeSample=calibrationMatter.sample({budget:5003});assert.equal(largeSample.points.length,5003,'optical sampling must honor a ceiling beyond 2048');assert.equal(calibrationMatter.sample({budget:0}).points.length,0);assert.equal(calibrationMatter.sample({budget:5003,exact:true}).points.length,5003);console.log('Passed: actual optical population beyond 2048 and a zero display ceiling.');
