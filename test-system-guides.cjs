const assert=require('node:assert/strict'),{systemBodyGuide}=require('./observation-tools');
const AU=149597870.7,width=30,large=systemBodyGuide(6371,width),small=systemBodyGuide(100,width);
assert(large.opacity>.5&&small.opacity===0,'planet-size objects precede smaller bodies');
const nearer=systemBodyGuide(100,width/100);assert(nearer.opacity>0,'dark minor bodies acquire a navigation marker on approach');
assert.deepEqual(systemBodyGuide(100,1),systemBodyGuide(1000,10),'equal angular size gives the same guide regardless of mass/category');
let previous=0;for(let k=0;k<=100;k++){const g=systemBodyGuide(1,100*10**(-k/20));assert(g.opacity>=previous&&g.opacity<=1);previous=g.opacity;}
assert.equal(systemBodyGuide(1e-9,1e6,{selected:true}).opacity,1,'an explicitly selected identity remains identifiable');
assert(systemBodyGuide(AU,.5).ringRadius>1000,'at close range the ring follows the physical contour');
console.log('Passed: size-based guide discovery, smooth opacity, scale invariance, selected identity and physical outline.');
