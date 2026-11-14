const assert=require('node:assert/strict'),F=require('./formation-model'),P=require('./stellar-population');
let tested=0,improved=0;const ratios=[];
for(let i=0;i<150;i++){
 const id='clearing-'+i,host=P.properties(id),field=F.build(id,host),control=F.build(id,host,{stellarClearing:false});assert.deepEqual(field.compact,control.compact,'conditioning the mature residual disk preserves the condensed mass field');
 const sum=field.radialCells.reduce((s,c)=>s+c.residual,0);assert(Math.abs(sum-field.residualMass)<Math.max(1e-18,field.residualMass*1e-9));let previous=0;
 for(const c of field.radialCells){assert(c.stellarSurvival>=previous&&c.stellarSurvival<=1);previous=c.stellarSurvival;assert(c.residual<=c.capacity*(1+1e-10),'redistribution respects the density cap');}
 if(field.compact.length&&field.residualMass>0){const r=field.compact[0].a*.5,f=g=>g.radialCells.filter(c=>c.r<r).reduce((s,c)=>s+c.residual,0)/g.residualMass,before=f(control),after=f(field);if(before>1e-8){tested++;if(after<before)improved++;ratios.push(after/before);}}
 assert.deepEqual(field.radialCDF,F.build(id,host).radialCDF,'birth distribution remains deterministic');
}
assert(tested>50&&improved===tested,'every measured inner region is depleted');ratios.sort((a,b)=>a-b);assert(ratios[Math.floor(ratios.length/2)]<.2,'typical central density is substantially reduced');
console.log('Passed: deterministic central birth depletion, smooth stellar survival, conserved mass, unchanged condensates and capped radial redistribution. '+JSON.stringify({systems:150,measured:tested,improved,medianInnerMassRatio:ratios[Math.floor(ratios.length/2)]}));
