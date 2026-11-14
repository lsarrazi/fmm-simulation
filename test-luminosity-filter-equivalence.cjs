const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict'),Population=require('./stellar-population');
const sandbox={module:{exports:{}}};vm.runInNewContext(fs.readFileSync('test-fixtures/stellar-population-reference.js','utf8'),sandbox);const Reference=sandbox.module.exports;let cases=0;
for(const limit of [0,1e-12,1e-4,.01,.5,1,100,1e4,1e6,1e100])for(const root of ['small','galaxy/root','milkyway:1000000000000000000000000/3lp']){
 const a=Population.luminosityFilter(limit),b=Reference.luminosityFilter(limit);assert.equal(a.key,b.key);assert.equal(a.fraction,b.fraction);const x=a.forGroup(root),y=b.forGroup(root);
 for(const start of [0n,1n,9007199254740991n,10n**21n+111n])for(const count of [0n,1n,31n,2n**30n,10n**12n]){assert.equal(x.count(start,count),y.count(start,count));for(const rank of [0n,1n,4n,999n])assert.equal(x.at(start,rank),y.at(start,rank));cases++;}
}
console.log('Passed: '+cases+' luminosity-filter cases exactly match the previous modulo arithmetic, including empty tails and ranks beyond 10^21.');
