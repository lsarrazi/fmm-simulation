const assert=require('node:assert/strict'),{addStellarLight}=require('./stellar-light');
function render(flux=100,profile='starburst',angle=0){const b=new Float32Array(256*256*3);addStellarLight(b,256,128.25,127.75,flux,[1,.7,.35],0,1,{profile,angle,glare:.8});return b;}
const a=render(),b=render(200);let total=0,outer=0,error=0;for(let y=0;y<256;y++)for(let x=0;x<256;x++){const i=(y*256+x)*3;total+=a[i];if(Math.hypot(x-128,y-128)>8)outer+=a[i];error+=Math.abs(b[i]-2*a[i]);}
assert(Math.abs(total-100)<1e-4,'glare redistributes source flux');assert(outer>.1,'visible halo/spikes extend beyond the point core');assert(error<1e-6,'brighter sources scale photon energy, not a drawn point radius');
const six=render(100,'six'),rotated=render(100,'six',Math.PI/6),at=(image,x,y)=>image[(y*256+x)*3];assert(at(six,148,128)>at(rotated,148,128),'orientation belongs to the instrument');
const soft=render(100,'soft');assert(at(six,148,128)>at(soft,148,128),'soft optic removes directional diffraction spikes');
console.log('Passed: normalized halo/diffraction optics, linear photon flux, long wings, instrument orientation and optical profiles.');
