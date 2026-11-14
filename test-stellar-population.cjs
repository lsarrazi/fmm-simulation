const assert=require('node:assert/strict'),P=require('./stellar-population'),Persistent=require('./persistent-matter');
const classes={};let low=0,mass=0;for(let i=0;i<65536;i++){
 const p=P.properties('survey/root/'+i);classes[p.type]=(classes[p.type]||0)+1;low+=p.mass<.5;mass+=p.mass;
 assert(p.mass>0&&p.luminosity>0&&p.radius>0&&p.temperature>0&&p.visibleLuminosity<=p.luminosity);
 assert(Math.abs(p.luminosity/(p.radius*p.radius*(p.temperature/P.solarTemperatureK)**4)-1)<1e-12,'Stefan–Boltzmann consistency for every population');
}
assert(low/65536>.65,'volume population dominated by faint low-mass stars');assert(Math.abs(mass/65536-P.meanMass)<.01);
assert(Math.abs(classes['Séquence principale']/65536-4946/5230)<.002);assert(Math.abs(classes['Naine blanche']/65536-264/5230)<.002);assert(Math.abs(((classes['Géante rouge']||0)+(classes['Géante très lumineuse']||0))/65536-20/5230)<.001);
assert(P.meanLuminosity>0&&P.meanVisibleLuminosity>0);assert(Math.abs(P.msLuminosity(1)-10**(-.007))<1e-12);
let previousTail={probability:1,luminosity:P.meanVisibleLuminosity};for(const limit of [0,1e-5,.001,.1,1,10,100,10000,1e10,Infinity]){const tail=P.luminosityTailMoments(limit);assert(tail.probability>=0&&tail.probability<=previousTail.probability+1e-12);assert(tail.luminosity>=0&&tail.luminosity<=previousTail.luminosity+1e-12);previousTail=tail;}assert.equal(previousTail.luminosity,0);
const fake={n:4,nv:3,dx:.75,dv:4/3,extent:1.5,initialPhase:[[49,1]]},m=new Persistent(fake,{population:100000n}),initialBytes=m.data.byteLength,filter=P.luminosityFilter(1000),query=m.sample({zoom:1,budget:2048,rankFilter:filter}),found=new Set(query.points.map(p=>p.id));
assert.equal(query.estimated,m.population,'bright query still reports full physical population');assert(query.visited<100,'rare sources found without walking the stars');assert(query.points.length<500);
let bright=0;for(let i=0;i<100000;i++){const id=m.evaluate(0,BigInt(i)).id;if(P.properties(id).luminosity>=1000){bright++;assert(found.has(id),'no physically bright star excluded by lazy filter');}}
assert(bright>0);const selected=query.points[0],before=P.properties(selected.id);m.data[4]*=8;m.refine(0);assert.deepEqual(P.properties(m.resolve(selected.id).id),before,'properties survive refinement without registry');
const huge=new Persistent(fake,{population:10n**24n}),q=huge.sample({budget:64,rankFilter:P.luminosityFilter(10000)});assert(q.points.length<=64&&q.visited<100);assert.equal(huge.data.byteLength,initialBytes);
const rare=huge.sample({budget:64,rankFilter:P.luminosityFilter(1e6)});assert(rare.points.some(p=>P.properties(p.id).mass>80),'extreme massive-star tail remains accessible beyond a 32-bit probability grid');
console.log('Passed: CNS5 proportions, low-mass dominance, empirical main sequence, physical L/R/T consistency, deterministic marks and lazy rare-bright discovery at 10^24.');console.log({classes,meanMass:P.meanMass,meanLuminosity:P.meanLuminosity,rareStars:bright,visited:query.visited});
