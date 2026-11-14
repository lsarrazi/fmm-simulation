const assert=require('node:assert/strict'),Matter=require('./persistent-matter'),Population=require('./stellar-population');
const grid={n:4,nv:1,dx:1,dv:1,extent:2,initialPhase:[[5,1]]};
const m=new Matter(grid,{population:1048576n,seed:'diagnostic'}),g=m.groups[0];
const filter=Population.luminosityFilter(50.11872336272722);assert.equal(filter.key,8);const f=filter.forGroup(m.universe+'/'+g.key.toString(36)),histX=Array(16).fill(0),histY=Array(16).fill(0),cells=new Set();
for(let j=0n;j<f.count(0n,g.count);j++){const rank=f.at(0n,j),p=m.initialPoint(g,rank),x=Math.floor((p.u+.5)*1024),y=Math.floor((p.v+.5)*1024);histX[x%16]++;histY[y%16]++;cells.add(x+','+y);assert.deepEqual(m.initialPoint(g,rank),p);}
assert.equal(cells.size,4096,'one deterministic address per leaf, no duplicates');
for(const h of [histX,histY]){assert(h.every(n=>n>150&&n<360),'bright stars must occupy every spatial phase, not a repeated lattice');}
// A complete filtered local query must discover exactly the same addresses as a global query.
const all=m.sample({zoom:1,budget:8192,rankFilter:filter}),local=m.sample({x:-.5,y:-.5,zoom:8,budget:8192,rankFilter:filter,exact:true}),expected=all.points.filter(p=>Math.abs(p.x+.5)<=.25&&Math.abs(p.y+.5)<=.25);
assert.deepEqual(local.points.map(p=>p.id).sort(),expected.map(p=>p.id).sort());
console.log('Passed: spatial/luminosity decorrelation, deterministic unique addresses and exact lazy filtered discovery.',{histX,histY,local:local.points.length});
