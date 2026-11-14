const assert=require('node:assert/strict'),Matter=require('./persistent-matter');
const grid={n:4,nv:3,dx:.75,dv:4/3,extent:1.5,initialPhase:[[(5*3+1)*3+2,1]]},m=new Matter(grid,{population:10n**24n}),zero=new Float64Array(16);m.coordinateCacheLimit=128;
const ids=Array.from({length:64},(_,i)=>m.evaluate(0,BigInt(i)*1234567890000000000n).id);
let calls=0;const original=m.initialPoint.bind(m);m.initialPoint=(...args)=>{calls++;return original(...args);};
for(let k=0;k<10;k++)for(const id of ids)m.resolve(id,k*.001);assert.equal(calls,0,'repeated local positions must not rebuild immutable BigInt addresses');
// Force refinement while keeping the initial material coordinates immutable.
m.data[4]=m.data[9]=2;const expected=ids.map(id=>m.resolve(id));assert(m.refine(0));for(let i=0;i<ids.length;i++){const p=m.resolve(ids[i]);assert(Math.abs(p.x-expected[i].x)<1e-14&&Math.abs(p.y-expected[i].y)<1e-14,'cached identity survives changed cohort indices and subdivisions');}
for(let k=0;k<4;k++){m.step(.001,zero,zero);for(const id of ids){const cached=m.resolve(id,.00001);m.coordinateCache.delete(id);assert.deepEqual(m.resolve(id,.00001),cached,'eviction changes neither identity nor moving reconstruction');}}
for(let i=0;i<512;i++)m.evaluate(0,m.groups[0].start+BigInt(i));assert(m.coordinateCache.size<=128,'memoization remains bounded independently of virtual population');
console.log('Passed: bounded immutable coordinate cache, no repeated BigInt subdivision, exact reconstruction after refinement, movement and eviction.');
