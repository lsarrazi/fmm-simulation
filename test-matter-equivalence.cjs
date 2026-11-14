const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict'),Matter=require('./persistent-matter'),Population=require('./stellar-population');
const sandbox={module:{exports:{}},console};vm.runInNewContext(fs.readFileSync('test-fixtures/persistent-matter-spatial-reference.js','utf8'),sandbox);const Reference=sandbox.module.exports;const populationSandbox={module:{exports:{}},console};vm.runInNewContext(fs.readFileSync('test-fixtures/stellar-population-reference.js','utf8'),populationSandbox);const ReferencePopulation=populationSandbox.module.exports;
const plain=value=>JSON.parse(JSON.stringify(value,(_,v)=>typeof v==='bigint'?v.toString()+'n':v));
const grid={n:8,nv:3,dx:.5,dv:.4,extent:2,initialPhase:[[71,.2],[103,.3],[220,.4],[401,.1]]},zero=new Float64Array(64);let addresses=0,queries=0;
for(const population of [100000000000n,9007199254740991n,9007199254740993n,10n**24n]){
 const a=new Matter(grid,{population}),b=new Reference(grid,{population});
 for(let g=0;g<a.groups.length;g++)for(let k=0;k<300;k++){const count=a.groups[g].count,rank=k<4?BigInt(k):k<8?count-BigInt(k-3):(count-1n)*BigInt(k*7919%100000)/100000n;assert.deepEqual(a.initialPoint(a.groups[g],rank),plain(b.initialPoint(b.groups[g],rank)));addresses++;}
 const anchor=a.evaluate(0,a.groups[0].count/3n);
 for(const zoom of [1,2.75,29,300,1e6])for(const filtered of [false,true]){const options={x:anchor.x+.013/zoom,y:anchor.y-.021/zoom,zoom,budget:37,offset:.0001},x=a.sample({...options,rankFilter:filtered?Population.luminosityFilter(.5):null}),y=b.sample({...options,rankFilter:filtered?ReferencePopulation.luminosityFilter(.5):null});assert.deepEqual(plain(x),plain(y),'same ordered identities, exact coordinates, population estimate and search-limit behaviour');queries++;}
 const zoom=Math.sqrt(Number(population))*100;assert.deepEqual(plain(a.sample({x:anchor.x,y:anchor.y,zoom,budget:128,exact:true})),plain(b.sample({x:anchor.x,y:anchor.y,zoom,budget:128,exact:true})));queries++;
 a.coordinateCacheLimit=b.coordinateCacheLimit=128;a.coordinateCache.clear();b.coordinateCache.clear();for(let k=0;k<700;k++){if(k%17===0){const id=a.coordinateCache.keys().next().value;a.coordinateCache.delete(id);b.coordinateCache.delete(id);}a.evaluate(0,BigInt(k));b.evaluate(0,BigInt(k));assert.deepEqual(Array.from(a.coordinateCache.keys()),Array.from(b.coordinateCache.keys()),'FIFO eviction order stays exact under churn and manual removal');}
 // Moving/refined cohorts must keep the same reconstruction after observation.
 a.data[4]=b.data[4]=2;a.data[9]=b.data[9]=2;a.step(.001,zero,zero);b.step(.001,zero,zero);assert.deepEqual(Array.from(a.data),Array.from(b.data));assert.deepEqual(plain(a.resolve(anchor.id,.0001)),plain(b.resolve(anchor.id,.0001)));
}
console.log('Passed: exact equivalence to the independent scrambled BigInt sampler across '+addresses+' addresses, '+queries+' queries, 53-bit boundaries, 10^24 identities, filtering, refinement and motion.');
