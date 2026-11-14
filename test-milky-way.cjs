const assert=require('node:assert/strict'),{MilkyWay:MW,createMilkyWay,stateMW}=require('./milky-way'),{TVDMatter,DirectBodies}=require('./galaxy-comparison'),Persistent=require('./persistent-matter'),{projectDensityLight,stellarName}=require('./stellar-light');
const b=createMilkyWay(2048),c=createMilkyWay(2048);assert.deepEqual(b,c,'repeatable physical initial state');
assert.equal(MW.massMsun,Number(MW.stars)*require('./stellar-population').meanMass,'collective gravity and inferred stellar population share the same expected mass budget');
assert.equal(MW.stars,100000000000n);assert(Math.max(...b.map(p=>Math.hypot(p.x,p.y)))*MW.lengthKpc<=MW.radiusKpc);assert(Math.abs(b.calibration.solarCircularKmS-230)<1e-8);assert(b.calibration.meanSpeedKmS>180&&b.calibration.meanSpeedKmS<240);assert(MW.numberDensitySolar>.1&&MW.numberDensitySolar<.2);
assert(Math.abs(MW.velocityKmS-9777.9222)<.01);assert(b.physics.gravity<.001,'physical gravity replaces arbitrary .48');assert(b.every(p=>Math.abs(p.vx)<b.physics.vmax&&Math.abs(p.vy)<b.physics.vmax));
const regular={...b.physics,structure:0},ring=Array.from({length:128},(_,i)=>stateMW(.65*Math.cos(i*Math.PI/64),.65*Math.sin(i*Math.PI/64),b.physics));
assert(Math.max(...ring.map(p=>p.disc))/Math.min(...ring.map(p=>p.disc))>3,'initial arms are resolved density structure');
assert(ring.some(p=>Math.abs(p.vr)*MW.velocityKmS>10),'coherent noncircular flow drives physical evolution');
assert.equal(stateMW(.65,0,regular).disc,stateMW(0,.65,regular).disc,'zero perturbation recovers circular disk');
assert.equal(stateMW(.65,0,regular).vr,0);assert.notEqual(stateMW(.03,.03,b.physics).bulge,stateMW(-.03,.03,b.physics).bulge,'central bar breaks circular symmetry');
// Independent Monte Carlo versus analytic angular probability over a broad annulus.
const mc=createMilkyWay(32768),bins=new Float64Array(12),expected=new Float64Array(12);let total=0;
for(const p of mc)if(p.disc&&Math.hypot(p.x,p.y)>.3&&Math.hypot(p.x,p.y)<1.2)bins[Math.floor((Math.atan2(p.y,p.x)+Math.PI)/(2*Math.PI)*12)%12]++;
for(let r=.302;r<1.2;r+=.004)for(let k=0;k<720;k++){const a=(k+.5)/720*2*Math.PI-Math.PI,m=stateMW(r*Math.cos(a),r*Math.sin(a),b.physics).disc*r;expected[Math.floor(k/60)]+=m;total+=m;}
const observed=bins.reduce((a,b)=>a+b,0);assert(bins.reduce((s,v,i)=>s+Math.abs(v/observed-expected[i]/total),0)<.06,'N-body samples the same asymmetric density as the continuum');
const g=new TVDMatter({n:32,nv:33,extent:2,vmax:.06});g.load(b);g.field(b.physics);const withHalo=g.ax.slice();g.field({...b.physics,haloV2:0});for(let i=0;i<g.mass.length;i++){const x=(i%g.n+.5)*g.dx-g.extent,y=(Math.floor(i/g.n)+.5)*g.dx-g.extent;assert(Math.abs(withHalo[i]-g.ax[i]+b.physics.haloV2*x/(x*x+y*y+.04))<1e-12);}
const projection=projectDensityLight(g,{resolution:192});assert(projection.image.every(v=>Number.isFinite(v)&&v>=0));assert(Math.abs(projection.mass-1)<.002,'smooth density reconstruction keeps stellar mass');
const fake={n:4,nv:3,dx:.75,dv:4/3,extent:1.5,initialPhase:[[49,1]]},m=new Persistent(fake),p=m.evaluate(0,123456789n);m.data[4]*=32;m.data[9]*=16;const stretched=m.resolve(p.id);
for(let k=0;k<8;k++)if(!m.refine(0))break;
assert(m.groups.length>1);assert.equal(m.groups.reduce((s,g)=>s+g.count,0n),m.population);const refined=m.resolve(p.id);assert(Math.abs(refined.x-stretched.x)<1e-13&&Math.abs(refined.y-stretched.y)<1e-13,'subdivision preserves the same star position');
for(let i=0;i<m.groups.length;i++)assert(Math.max(Math.hypot(m.data[i*12+4],m.data[i*12+5]),Math.hypot(m.data[i*12+8],m.data[i*12+9]))<=1.75*fake.dx,'no unbounded affine ribbon');
assert.equal(stellarName(p.id),stellarName(refined.id),'name survives adaptive subdivision');
const name=stellarName(p.id);m.step(.01,new Float64Array(16).fill(.01),new Float64Array(16));assert.equal(stellarName(m.resolve(p.id).id),name,'name survives motion');assert.notEqual(stellarName(p.id+'0'),name);assert.equal(new Set(Array.from({length:1000},(_,i)=>stellarName('universe/root/'+i))).size,1000,'useful name diversity');
console.log('Passed: MW physical scale, repeatability, stellar count, radial density, solar rotation, halo consistency, smooth positive projection and identity-preserving subdivision.');console.log({meanSpeedKmS:b.calibration.meanSpeedKmS,solarNumberDensity:MW.numberDensitySolar,adaptiveGroups:m.groups.length});
