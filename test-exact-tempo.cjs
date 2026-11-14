const assert=require('node:assert/strict'),{SimulationClock,physicalNs,unitNs}=require('./observation-tools');
for(const u of ['ns','us','ms','s','min','h','d','w','mo','y','ky','My','Gy','Ty'])assert(unitNs(u)>0n);
const c=new SimulationClock();c.setFactor('1/1000000');assert.equal(c.factorRatio(),'1/1000000');assert.equal(c.rate,physicalNs('1','y'));assert.equal(c.quantity('y'),'1');
c.setFactor('1/31557600000000');for(let i=0;i<100;i++)c.tick(10);assert.equal(c.visible(),physicalNs('1','s'));
c.debt=0n;c.setRate('1/3','ns');for(let i=0;i<300;i++)c.tick(10);assert.equal(c.visible(),1n,'rational sub-nanosecond rate keeps its remainder');
c.setRate('0,125','h');assert.equal(c.quantity('min'),'15/2');c.setFactor('1e-6');assert.equal(c.factorRatio(),'1/1000000');assert.throws(()=>c.setFactor('1/0'));assert.throws(()=>c.setFactor('-1'));
console.log('Passed: 14 physical units, exact decimal/scientific/fraction factors, year and second slowdown, rational accumulation and validation.');
const fine=new SimulationClock();fine.setRate('1','y');fine.tick(100);assert.equal(fine.duration(),0);fine.minimumStep=1e-8;assert(fine.duration()>0);assert.equal(fine.fineOffset(),0,'unintegrated local backlog must not be extrapolated as fine display time');const dt=fine.duration();fine.consume(dt);assert(fine.visible()>0n);assert(fine.debt<fine.quantum());
