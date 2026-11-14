const assert=require('node:assert/strict'),{lightPlan,addPSF,addPhotosphere,stellarProperties,cameraGeometry,sourceOptics}=require('./stellar-light'),P=require('./stellar-population'),Persistent=require('./persistent-matter');
function field(width,expected=100){return{image:new Float64Array([expected/1e11]),resolution:1,width,pixel:width,left:-width/2,bottom:-width/2,expected,mass:expected/1e11};}
const solar={id:'sun',x:0,y:0,photometry:{mass:1,luminosity:1,radius:1,temperature:5772,visibleLuminosity:P.visibleFraction(5772),color:[1,1,1]}},a=cameraGeometry(field(1)),b=cameraGeometry(field(2));
assert(Math.abs(sourceOptics(solar,a).fluxWm2/sourceOptics(solar,b).fluxWm2-4)<1e-12,'double distance quarters received flux');
const sameHeight=field(2*Math.tan(20*Math.PI/180)*a.heightPc/10000),normal=lightPlan(field(1),[solar],{population:10n**11n}),changedLens=lightPlan(sameHeight,[solar],{population:10n**11n,fovDegrees:40});assert(Math.abs(normal.candidates[0].linearFlux/changedLens.candidates[0].linearFlux-1)<1e-12,'lens field of view does not invent received light at constant camera distance');
const lateral=sourceOptics({...solar,x:a.heightPc/10000},a);assert(Math.abs(lateral.fluxWm2/sourceOptics(solar,a).fluxWm2-.5)<1e-12,'off-axis distance uses all three coordinates');
assert.equal(lateral.radiusPixels,sourceOptics(solar,a).radiusPixels,'projected size uses camera Z rather than off-axis Euclidean distance');
assert(Math.abs(sourceOptics(solar,a).radiusPixels/sourceOptics(solar,b).radiusPixels-2)<1e-12,'doubling camera Z halves the physical disk');
assert.equal(sourceOptics({...solar,photometry:{...solar.photometry,mass:100,luminosity:1e5}},a).radiusPixels,sourceOptics(solar,a).radiusPixels,'mass and luminosity never act as radius multipliers');
assert(sourceOptics(solar,a).radiusPixels<1e-6,'photosphere has no minimum fake point size');
const bright={...solar,id:'bright',photometry:{...solar.photometry,luminosity:1e5,visibleLuminosity:5e4}},dim={...solar,id:'dim',photometry:{...solar.photometry,luminosity:.001,visibleLuminosity:.0001}};
const far=lightPlan(field(1),[bright,dim],{population:10n**11n});assert(far.drawn.some(p=>p.id==='bright')&&!far.drawn.some(p=>p.id==='dim'),'luminous star becomes detectable before faint star');
const near=lightPlan(field(1e-5,2),[bright,dim],{population:10n**11n,complete:true});assert.equal(near.drawn.length,2);assert.equal(near.residual,0,'no diffuse duplication when local sources are enumerated');
const less=lightPlan(field(1),[bright],{population:10n**11n,exposureEV:-8});assert(less.drawn.length<far.drawn.length,'exposure controls detectability');
function image(R,sigma=1,px=64.3){const buffer=new Float32Array(128*128*3);addPhotosphere(buffer,128,px,63.7,100,[1,1,1],R,sigma);let total=0,peak=0;for(let i=0;i<128*128;i++){total+=buffer[i*3];peak=Math.max(peak,buffer[i*3]);}return{total,peak};}
for(const R of [0,.2,3,15,30])assert(Math.abs(image(R).total-100)<.5,'PSF and resolved photosphere conserve source flux');
assert(image(30).peak<image(3).peak,'physical disk spreads flux over its real area');assert(image(30,1,1).total<60,'offscreen photosphere flux is not renormalized into the viewport');
const buffer=new Float32Array(64*64*3);addPSF(buffer,64,32.3,31.7,100,[1,1,1],2);assert(Math.abs(buffer.filter((_,i)=>i%3===0).reduce((a,b)=>a+b,0)-100)<1e-4);
assert.deepEqual(stellarProperties('identity'),stellarProperties('identity'));
const grid={n:4,nv:3,dx:.75,dv:4/3,extent:1.5,initialPhase:[[49,1]]},m=new Persistent(grid),star=m.evaluate(0,1234n),before=m.groups[0].count;m.step(1000,new Float64Array(16),new Float64Array(16));assert.equal(m.resolve(star.id).id,star.id);assert.equal(m.groups[0].count,before,'no expiration after leaving visible domain');
console.log('Passed: inverse-square flux, 3D camera distance, luminous-first detection, exposure, physical photosphere radius, limb-darkened disk/PSF flux and persistent properties.');
