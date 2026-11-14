const assert=require('node:assert/strict'),{lightPlan,diffuseCrowdingWeight,addPhotosphere,addStellarLight,projectDensityLight,photosphereResolution}=require('./stellar-light');
const N=10n**11n,size=600;
function field(count,resolution=32){const mass=count/Number(N);return{image:new Float64Array(resolution*resolution).fill(mass/(resolution*resolution)),mass,expected:count,resolution,width:1,pixel:1/resolution,left:-.5,bottom:-.5};}
// Only the expected population changes. No stars need to be instantiated.
let previous=0;for(const count of [0,1,10,100,1000,10000,100000,1e6,1e11]){const p=lightPlan(field(count),[],{population:N});assert(p.residual>=previous-1e-12&&p.residual<=1);previous=p.residual;if(count<=1000)assert(p.residual<.001);if(count>=1e6)assert(p.residual>.999);}
const p=lightPlan(field(20000),[],{population:N}),hi=lightPlan(field(20000,64),[],{population:N,size:1200,sigma:2});assert(Math.abs(p.residual-hi.residual)<1e-12,'same camera/optics at different raster resolution');
const before=lightPlan(field(20000),[],{population:N}),after=lightPlan(field(20000*.9999),[],{population:N});assert(Math.abs(before.residual-after.residual)<.0001,'small zoom changes produce a smooth statistical background');
const lit=lightPlan(field(1e6),[{id:'bright',x:0,y:0,photometry:{radius:1,luminosity:10000,visibleLuminosity:1000,color:[1,1,1]}}],{population:N});assert(lit.sources.length===1);assert(lit.diffuseImage.every(v=>v===lit.diffuseImage[0]),'materializing a bright source cannot carve a local dark halo in a uniform ensemble background');
for(const lambda of [1e-12,1e-6,.000999999,.001000001,1,10,50]){const a=diffuseCrowdingWeight(lambda),b=diffuseCrowdingWeight(lambda*(1+1e-6));assert(Number.isFinite(a)&&a>=0&&b>=a-1e-14);assert(Math.abs(b-a)<1e-5);}
const grid={n:4,dx:1,extent:2,mass:new Float64Array(16).fill(1/16)},wide=projectDensityLight(grid,{zoom:1,resolution:32}),close=projectDensityLight(grid,{zoom:10,resolution:32});assert(close.mass<wide.mass/50,'local population integrates the density footprint');
function profile(radius){const b=new Float32Array(64*64*3);addPhotosphere(b,64,31.3,32.1,100,[1,1,1],radius,1);return b;}
const distant=profile(0),unresolved=profile(.1);assert.deepEqual(unresolved,distant,'unresolved photosphere footprint is optical, not proportional to physical angular size');
assert.equal(photosphereResolution(.1).resolved,false);assert.equal(photosphereResolution(10).resolved,true);
const glow=radius=>{const b=new Float32Array(128*128*3);addStellarLight(b,128,64.2,63.8,100,[1,1,1],radius,1,{glare:.8});return b;};
for(const radius of [.3,.8,1.1])assert.deepEqual(glow(radius),glow(0),'all unresolved stellar radii share the same core and halo profile');
for(const radius of [1.1774100225,2,3]){const a=profile(radius*(1-1e-5)),b=profile(radius*(1+1e-5));let delta=0,flux=0;for(let i=0;i<a.length;i+=3){delta+=Math.abs(a[i]-b[i]);flux+=a[i];}assert(delta/flux<.001,'no jump between point PSF and resolved photosphere');}
console.log('Passed: estimated-density background fade, no hidden-star enumeration, zoom/raster continuity, fixed unresolved optical footprint and continuous photosphere resolution.');
