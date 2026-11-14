/* WebGPU kinetic solver: phase-space buffers remain on GPU. The GPU FMM uses
 * Orbit's P2M/M2M/M2L/L2L scheme with linear local expansions, on the spatial grid.
 * Only the projected density and diagnostics are read back for the shared renderer.
 */
(function(root){
'use strict';
const shader = `
struct Config { n:u32, nv:u32, size:u32, depth:u32, dx:f32, dv:f32, extent:f32, vmax:f32, gravity:f32, eps2:f32, duration:f32, haloV2:f32 }
struct Stage { axis:u32, level:u32, pad0:u32, pad1:u32 }
struct Cell { mass:f32, mx:f32, my:f32, ax:f32, ay:f32, jxx:f32, jxy:f32, jyy:f32 }
struct State { dt:f32, time:f32, space:f32, velocity:f32 }
struct Seed { index:u32, mass:f32 }
@group(0) @binding(0) var<uniform> cfg:Config;
@group(0) @binding(1) var<uniform> stage:Stage;
@group(0) @binding(2) var<storage,read> src:array<f32>;
@group(0) @binding(3) var<storage,read_write> dst:array<f32>;
@group(0) @binding(4) var<storage,read_write> flux:array<f32>;
@group(0) @binding(5) var<storage,read_write> cells:array<Cell>;
@group(0) @binding(6) var<storage,read_write> state:State;
// Per-transport-line boundary losses, reduced without float atomics.
@group(0) @binding(7) var<storage,read_write> loss:array<f32>;
@group(0) @binding(8) var<storage,read> seeds:array<Seed>;
@group(0) @binding(9) var<storage,read_write> result:array<f32>;
fn offset(level:u32)->u32 { return ((1u << (2u*level))-1u)/3u; }
fn index2(id:vec3<u32>)->u32 {return id.x+id.y*65535u*256u;}
fn velocity(k:u32)->f32 { return (f32(k)-f32(cfg.nv-1u)*0.5)*cfg.dv; }
fn pos(i:u32, side:u32)->vec2<f32> {return (vec2<f32>(f32(i%side),f32(i/side))+0.5)*(2.0*cfg.extent/f32(side))-cfg.extent;}
@compute @workgroup_size(256) fn initialize(@builtin(global_invocation_id) id:vec3<u32>){let i=index2(id);if(i<arrayLength(&seeds)){dst[seeds[i].index]=seeds[i].mass;}}
var<workgroup> partial:array<f32,256>;
var<workgroup> angularPartial:array<f32,256>;
fn velocityPDF(k:u32,ca:f32,sa:f32,meanR:f32,meanT:f32)->vec2<f32>{let vx=velocity(k%cfg.nv);let vy=velocity(k/cfg.nv);let vr=ca*vx+sa*vy-meanR;let vp=-sa*vx+ca*vy-meanT;
 let a=vr*9777.9222/30.0;let b=vp*9777.9222/20.0;let sigma=9777.9222/90.0;return vec2<f32>(exp(-0.5*(a*a+b*b)),exp(-0.5*(vx*vx+vy*vy)*sigma*sigma));
}
@compute @workgroup_size(256) fn initializeMW(@builtin(workgroup_id) group:vec3<u32>,@builtin(local_invocation_index) lane:u32){
 let cell=group.x+group.y*65535u;if(cell>=cfg.n*cfg.n){return;}let p=pos(cell,cfg.n);let r=length(p);let base=98u+4u*cell;let density=vec2<f32>(seeds[base].mass,seeds[base+1u].mass);let meanR=seeds[base+2u].mass;let meanT=seeds[base+3u].mass;let ca=p.x/max(r,1e-20);let sa=p.y/max(r,1e-20);
 var sums=vec2<f32>(0.0);for(var k=lane;k<cfg.nv*cfg.nv;k+=256u){sums+=velocityPDF(k,ca,sa,meanR,meanT);}partial[lane]=sums.x;angularPartial[lane]=sums.y;workgroupBarrier();
 for(var step=128u;step>0u;step/=2u){if(lane<step){partial[lane]+=partial[lane+step];angularPartial[lane]+=angularPartial[lane+step];}workgroupBarrier();}
 for(var k=lane;k<cfg.nv*cfg.nv;k+=256u){let pdf=velocityPDF(k,ca,sa,meanR,meanT);dst[cell*cfg.nv*cfg.nv+k]=cfg.dx*cfg.dx*(density.x*pdf.x/max(partial[0],1e-30)+density.y*pdf.y/max(angularPartial[0],1e-30))/seeds[97u].mass;}
}
@compute @workgroup_size(256) fn project(@builtin(workgroup_id) group:vec3<u32>,@builtin(local_invocation_index) lane:u32){
 let i=group.x+group.y*65535u; if(i>=cfg.n*cfg.n){return;}
 var sum=0.0;var angular=0.0;let position=pos(i,cfg.n);for(var k=lane;k<cfg.nv*cfg.nv;k+=256u){let m=src[i*cfg.nv*cfg.nv+k];sum+=m;angular+=m*(position.x*velocity(k/cfg.nv)-position.y*velocity(k%cfg.nv));}
 partial[lane]=sum;angularPartial[lane]=angular;workgroupBarrier();
 for(var stride=128u;stride>0u;stride/=2u){if(lane<stride){partial[lane]+=partial[lane+stride];angularPartial[lane]+=angularPartial[lane+stride];}workgroupBarrier();}
 if(lane==0u){let p=position;cells[offset(cfg.depth)+i]=Cell(partial[0],partial[0]*p.x,partial[0]*p.y,0.0,0.0,0.0,0.0,0.0);result[8u+cfg.n*cfg.n+i]=angularPartial[0];}
}
@compute @workgroup_size(256) fn upward(@builtin(global_invocation_id) id:vec3<u32>){
 let i=index2(id);let side=1u<<stage.level;if(i>=side*side){return;}
 var c=Cell(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
 for(var y=0u;y<2u;y++){for(var x=0u;x<2u;x++){let a=cells[offset(stage.level+1u)+(i/side*2u+y)*side*2u+i%side*2u+x];c.mass+=a.mass;c.mx+=a.mx;c.my+=a.my;}}
 cells[offset(stage.level)+i]=c;
}
@compute @workgroup_size(256) fn interaction(@builtin(global_invocation_id) id:vec3<u32>){
 let i=index2(id);let side=1u<<stage.level;if(i>=side*side){return;}
 let x=i32(i%side);let y=i32(i/side);let p=pos(i,side);var out=cells[offset(stage.level)+i];
 for(var py=y/2-1;py<=y/2+1;py++){for(var px=x/2-1;px<=x/2+1;px++){
  if(px<0||py<0||px>=i32(side/2u)||py>=i32(side/2u)){continue;}
  for(var cy=0;cy<2;cy++){for(var cx=0;cx<2;cx++){
   let sx=px*2+cx;let sy=py*2+cy;if(abs(sx-x)<=1&&abs(sy-y)<=1){continue;}
   let a=cells[offset(stage.level)+u32(sy)*side+u32(sx)];if(a.mass<=0.0){continue;}
   let r=vec2<f32>(a.mx,a.my)/a.mass-p;let q=dot(r,r)+cfg.eps2;let k=cfg.gravity*a.mass/(q*sqrt(q));
   out.ax+=r.x*k;out.ay+=r.y*k;out.jxx+=k*(-1.0+3.0*r.x*r.x/q);out.jxy+=k*3.0*r.x*r.y/q;out.jyy+=k*(-1.0+3.0*r.y*r.y/q);
  }}
 }}let address=offset(stage.level)+i;cells[address].ax=out.ax;cells[address].ay=out.ay;cells[address].jxx=out.jxx;cells[address].jxy=out.jxy;cells[address].jyy=out.jyy;
}
@compute @workgroup_size(256) fn downward(@builtin(global_invocation_id) id:vec3<u32>){
 let i=index2(id);let side=1u<<stage.level;if(i>=side*side){return;}
 let pi=(i/side/2u)*(side/2u)+i%side/2u;let parent=cells[offset(stage.level-1u)+pi];
 let d=pos(i,side)-pos(pi,side/2u);var c=cells[offset(stage.level)+i];
 c.ax+=parent.ax+parent.jxx*d.x+parent.jxy*d.y;c.ay+=parent.ay+parent.jxy*d.x+parent.jyy*d.y;
 c.jxx+=parent.jxx;c.jxy+=parent.jxy;c.jyy+=parent.jyy;cells[offset(stage.level)+i]=c;
}
@compute @workgroup_size(256) fn nearField(@builtin(global_invocation_id) id:vec3<u32>){
 let i=index2(id);if(i>=cfg.n*cfg.n){return;}let leaf=offset(cfg.depth);var c=cells[leaf+i];
 let x=i32(i%cfg.n);let y=i32(i/cfg.n);let p=pos(i,cfg.n);
 for(var oy=-1;oy<=1;oy++){for(var ox=-1;ox<=1;ox++){
  let sx=x+ox;let sy=y+oy;if(sx<0||sy<0||sx>=i32(cfg.n)||sy>=i32(cfg.n)||(ox==0&&oy==0)){continue;}
  let a=cells[leaf+u32(sy)*cfg.n+u32(sx)];let r=pos(u32(sy)*cfg.n+u32(sx),cfg.n)-p;
  let q=dot(r,r)+cfg.eps2;let k=cfg.gravity*a.mass/(q*sqrt(q));c.ax+=r.x*k;c.ay+=r.y*k;
 }}let halo=cfg.haloV2/(dot(p,p)+0.04);cells[leaf+i].ax=c.ax-p.x*halo;cells[leaf+i].ay=c.ay-p.y*halo;
}
@compute @workgroup_size(1) fn chooseStep(){
 var maximum=0.0;let leaf=offset(cfg.depth);for(var i=0u;i<cfg.n*cfg.n;i++){let c=cells[leaf+i];maximum=max(maximum,max(abs(c.ax),abs(c.ay)));}
 state.dt=min(cfg.duration,0.8*cfg.dx/cfg.vmax);if(maximum>0.0){state.dt=min(state.dt,0.8*cfg.dv/maximum);}state.time+=state.dt;
}
@compute @workgroup_size(256) fn saveField(@builtin(global_invocation_id) id:vec3<u32>){let i=index2(id);if(i>=cfg.n*cfg.n){return;}let c=cells[offset(cfg.depth)+i];result[8u+2u*cfg.n*cfg.n+i]=c.ax;result[8u+3u*cfg.n*cfg.n+i]=c.ay;}
fn strideFor()->u32 {switch stage.axis {case 0u:{return cfg.nv*cfg.nv;}case 1u:{return cfg.n*cfg.nv*cfg.nv;}case 2u:{return 1u;}default:{return cfg.nv;}}}
fn speedFor(i:u32)->f32 {switch stage.axis {case 0u:{return velocity(i%cfg.nv);}case 1u:{return velocity(i/cfg.nv%cfg.nv);}case 2u:{return cells[offset(cfg.depth)+i/(cfg.nv*cfg.nv)].ax;}default:{return cells[offset(cfg.depth)+i/(cfg.nv*cfg.nv)].ay;}}}
@compute @workgroup_size(256) fn makeFlux(@builtin(global_invocation_id) id:vec3<u32>){
 let i=index2(id);if(i>=cfg.size){return;}let stride=strideFor();let count=select(cfg.nv,cfg.n,stage.axis<2u);let coord=i/stride%count;
 let v=src[i];if(v==0.0){flux[i]=0.0;return;}let speed=speedFor(i);let c=abs(speed)*state.dt/select(cfg.dv,cfg.dx,stage.axis<2u);
 var left=v;var right=-v;if(coord>0u){left-=src[i-stride];}if(coord+1u<count){right+=src[i+stride];}
 var slope=0.0;if(left*right>0.0){slope=sign(left)*min(min(2.0*abs(left),2.0*abs(right)),abs(left+right)*0.5);}
 flux[i]=c*(v+0.5*(1.0-c)*select(-1.0,1.0,speed>=0.0)*slope);
}
@compute @workgroup_size(256) fn transport(@builtin(global_invocation_id) id:vec3<u32>){
 let i=index2(id);if(i>=cfg.size){return;}let stride=strideFor();let count=select(cfg.nv,cfg.n,stage.axis<2u);let coord=i/stride%count;let speed=speedFor(i);
 var incoming=0.0;if(speed>=0.0&&coord>0u){incoming=flux[i-stride];}if(speed<0.0&&coord+1u<count){incoming=flux[i+stride];}
 dst[i]=src[i]-flux[i]+incoming;
 if(coord==0u){let line=i/(stride*count)*stride+i%stride;let last=i+(count-1u)*stride;loss[line]=select(flux[i],flux[last],speed>=0.0);}
}
@compute @workgroup_size(256) fn collectLoss(@builtin(local_invocation_index) lane:u32){
 let count=select(cfg.nv,cfg.n,stage.axis<2u);var sum=0.0;for(var i=lane;i<cfg.size/count;i+=256u){sum+=loss[i];}
 partial[lane]=sum;workgroupBarrier();for(var stride=128u;stride>0u;stride/=2u){if(lane<stride){partial[lane]+=partial[lane+stride];}workgroupBarrier();}
 if(lane==0u){if(stage.axis<2u){state.space+=partial[0];}else{state.velocity+=partial[0];}}
}
@compute @workgroup_size(1) fn summary(){
 result[0]=state.dt;result[1]=state.time;result[2]=state.space;result[3]=state.velocity;
 var mass=0.0;var angular=0.0;
 for(var i=0u;i<cfg.n*cfg.n;i++){mass+=cells[offset(cfg.depth)+i].mass;angular+=result[8u+cfg.n*cfg.n+i];}
 result[4]=mass;result[5]=angular;
 for(var i=0u;i<cfg.n*cfg.n;i++){result[8u+i]=cells[offset(cfg.depth)+i].mass;}
}
`;

class GPUKineticMatter {
 static async create({n,nv,bodies,extent=1.5,vmax=2,physics={}}){
  if(!navigator.gpu)throw new Error('WebGPU indisponible : utilisez Chrome ou Edge avec accélération graphique.');
  const adapter=await navigator.gpu.requestAdapter({powerPreference:'high-performance'});if(!adapter)throw new Error('Aucun adaptateur WebGPU disponible.');
  const size=n*n*nv*nv,bytes=size*4;
  if(n<4||(n&(n-1))!==0)throw new Error('La grille GPU doit avoir une taille spatiale puissance de deux.');
  const limit=Math.min(adapter.limits.maxStorageBufferBindingSize,adapter.limits.maxBufferSize);
  if(bytes>limit||bytes*3>1024**3)throw new Error(`Grille trop grande : ${(bytes*3/1048576).toFixed(0)} Mio pour les trois grands buffers, limite par buffer ${(limit/1048576).toFixed(0)} Mio.`);
  const device=await adapter.requestDevice({requiredLimits:{maxStorageBufferBindingSize:Math.max(bytes,134217728),maxBufferSize:Math.max(bytes,268435456)}});
  const engine=new GPUKineticMatter();Object.assign(engine,{device,n,nv,size,extent,vmax,physics,dx:2*extent/n,dv:2*vmax/nv,depth:Math.log2(n),time:0,spatialLoss:0,velocityLoss:0,mass:new Float64Array(n*n),buffers:[],timeMs:0});
  try {await engine.initialize(bodies);return engine;}catch(e){engine.destroy();throw e;}
 }
 buffer(label,size,usage,data){const b=this.device.createBuffer({label,size:Math.max(16,size),usage});this.buffers.push(b);if(data)this.device.queue.writeBuffer(b,0,data);return b;}
 async initialize(bodies){
  const d=this.device,U=GPUBufferUsage;
  d.lost.then(info=>{this.failure=new Error('GPU perdu : '+info.message);});
  d.addEventListener('uncapturederror',event=>{this.failure=event.error;});
  const module=d.createShaderModule({label:'Orbit continuous GPU',code:shader});
  const info=await module.getCompilationInfo();const errors=info.messages.filter(m=>m.type==='error');if(errors.length)throw new Error(errors.map(m=>`${m.lineNum}: ${m.message}`).join('\n'));
  // Explicit shared layout keeps every kernel compatible even with unused resources.
  const entries=Array.from({length:10},(_,binding)=>({binding,visibility:GPUShaderStage.COMPUTE,buffer:{type:binding<2?'uniform':binding===2||binding===8?'read-only-storage':'storage'}}));
  // Ten bindings include eight storage buffers, within WebGPU's default per-stage limit.
  const layout=d.createBindGroupLayout({entries});const pipelineLayout=d.createPipelineLayout({bindGroupLayouts:[layout]});
  this.pipelines={};for(const entryPoint of ['initialize','initializeMW','project','upward','interaction','downward','nearField','chooseStep','saveField','makeFlux','transport','collectLoss','summary'])this.pipelines[entryPoint]=await d.createComputePipelineAsync({layout:pipelineLayout,compute:{module,entryPoint}});
  this.config=this.buffer('kinetic config',48,U.UNIFORM|U.COPY_DST);
  this.states=this.buffer('kinetic time and losses',16,U.STORAGE|U.COPY_DST);
  this.phase=[0,1].map(i=>this.buffer('phase '+i,this.size*4,U.STORAGE|U.COPY_SRC));
  this.flux=this.buffer('TVD flux',this.size*4,U.STORAGE);
  this.tree=this.buffer('spatial multipoles',((4**(this.depth+1)-1)/3)*32,U.STORAGE);
  this.loss=this.buffer('line outflow',this.size/Math.min(this.n,this.nv)*4,U.STORAGE);
  this.output=this.buffer('density, field and diagnostics',(8+4*this.n*this.n)*4,U.STORAGE|U.COPY_SRC);
  this.readback=this.buffer('small readback',(8+4*this.n*this.n)*4,U.MAP_READ|U.COPY_DST);
  this.ax=new Float64Array(this.n*this.n);this.ay=new Float64Array(this.n*this.n);
  const sparse=new Map();
  const w=(value,origin,spacing,count)=>{const u=(value-origin)/spacing-.5,lo=Math.floor(u);if(lo<0||lo+1>=count)throw new Error('Conditions initiales hors domaine.');return [[lo,1-(u-lo)],[lo+1,u-lo]];};
  if(!this.physics.analytic)for(const b of bodies)for(const [x,a]of w(b.x,-this.extent,this.dx,this.n))for(const [y,c]of w(b.y,-this.extent,this.dx,this.n))for(const [u,e]of w(b.vx,-this.vmax,this.dv,this.nv))for(const [v,f]of w(b.vy,-this.vmax,this.dv,this.nv)){
   const i=((y*this.n+x)*this.nv+v)*this.nv+u;sparse.set(i,Math.fround((sparse.get(i)||0)+b.mass*a*c*e*f));
  }
  const seedData=new ArrayBuffer((this.physics.analytic?98+4*this.n*this.n:sparse.size)*8),ints=new Uint32Array(seedData),floats=new Float32Array(seedData);let k=0;for(const[i,m]of sparse){ints[k*2]=i;floats[k*2+1]=m;k++;}
  this.initialPhase=Array.from(sparse.entries());
  if(this.physics.analytic){let norm=0;for(let y=0;y<this.n;y++)for(let x=0;x<this.n;x++){const a=root.stateMW((x+.5)*this.dx-this.extent,(y+.5)*this.dx-this.extent,this.physics),base=98+4*(y*this.n+x);norm+=(a.disc+a.bulge)*this.dx*this.dx;for(const [i,value]of [a.disc,a.bulge,a.vr,a.vt].entries())floats[(base+i)*2+1]=value;}floats[97*2+1]=norm;this.initialCohorts=root.cohortsForMW(this.n,this.extent,this.physics);}
  this.seeds=this.buffer('sparse initial distribution',seedData.byteLength,U.STORAGE|U.COPY_DST,seedData);
  this.groups=new Map();
  this.bind=(axis=0,level=0,input=0)=>{
   const key=`${axis}:${level}:${input}`;if(this.groups.has(key))return this.groups.get(key);
   const stage=this.buffer('axis and level',16,U.UNIFORM|U.COPY_DST,new Uint32Array([axis,level,0,0]));
   const resources=[this.config,stage,this.phase[input],this.phase[1-input],this.flux,this.tree,this.states,this.loss,this.seeds,this.output];
   const group=d.createBindGroup({layout,entries:resources.map((buffer,binding)=>({binding,resource:{buffer}}))});this.groups.set(key,group);return group;
  };
  this.writeConfig(.01,this.physics);
  const encoder=d.createCommandEncoder();if(this.physics.analytic)this.dispatch(encoder,'initializeMW',this.n*this.n,0,0,1,true);else this.dispatch(encoder,'initialize',sparse.size,0,0,1);this.dispatch(encoder,'project',this.n*this.n,0,0,0,true);this.encodeField(encoder);this.dispatch(encoder,'summary',1,0,0,0);encoder.copyBufferToBuffer(this.output,0,this.readback,0,(8+4*this.n*this.n)*4);d.queue.submit([encoder.finish()]);await this.read();
 }
 writeConfig(duration,{gravity=.48,softening=.042,haloV2=0}){const b=new ArrayBuffer(48),u=new Uint32Array(b),f=new Float32Array(b);u.set([this.n,this.nv,this.size,this.depth]);f.set([this.dx,this.dv,this.extent,this.vmax,gravity,softening*softening,duration,haloV2],4);this.device.queue.writeBuffer(this.config,0,b);}
 dispatch(encoder,name,count,axis=0,level=0,input=0,groups=false){const pass=encoder.beginComputePass();pass.setPipeline(this.pipelines[name]);pass.setBindGroup(0,this.bind(axis,level,input));const work=groups?count:Math.ceil(count/256);pass.dispatchWorkgroups(Math.min(work,65535),Math.ceil(work/65535));pass.end();}
 encodeField(e){for(let l=this.depth-1;l>=0;l--)this.dispatch(e,'upward',4**l,0,l);for(let l=2;l<=this.depth;l++)this.dispatch(e,'interaction',4**l,0,l);for(let l=1;l<=this.depth;l++)this.dispatch(e,'downward',4**l,0,l);this.dispatch(e,'nearField',this.n*this.n);this.dispatch(e,'saveField',this.n*this.n);}
 async read(){await this.readback.mapAsync(GPUMapMode.READ);const a=new Float32Array(this.readback.getMappedRange()).slice();this.readback.unmap();if(this.failure)throw this.failure;const count=this.n*this.n;this.spatialLoss=a[2];this.velocityLoss=a[3];this.totalMass=a[4];this.angularMomentum=a[5];this.mass.set(a.subarray(8,8+count));this.ax.set(a.subarray(8+2*count,8+3*count));this.ay.set(a.subarray(8+3*count));return a[0];}
 async advance(duration,parameters={}){
  const started=performance.now(),tolerance=Math.max(1e-12,duration*1e-6);let remaining=duration,steps=0;this.fieldMs=0;this.transportMs=0;
  while(remaining>tolerance&&steps<12){
   if(this.failure)throw this.failure;this.writeConfig(remaining,parameters);const e=this.device.createCommandEncoder();
   this.dispatch(e,'project',this.n*this.n,0,0,0,true);
   this.encodeField(e);this.dispatch(e,'chooseStep',1);
   let input=0;for(const axis of [2,3,0,1]){this.dispatch(e,'makeFlux',this.size,axis,0,input);this.dispatch(e,'transport',this.size,axis,0,input);this.dispatch(e,'collectLoss',1,axis,0,input);input=1-input;}
   this.dispatch(e,'project',this.n*this.n,0,0,0,true);this.dispatch(e,'summary',1);
   e.copyBufferToBuffer(this.output,0,this.readback,0,(8+4*this.n*this.n)*4);this.device.queue.submit([e.finish()]);
   const dt=await this.read();if(!(dt>0))throw new Error('Pas GPU invalide.');
   // Host double precision time avoids float32 accumulation swallowing slow steps.
   this.time+=dt;
   if(this.onTransport)await this.onTransport({dt,ax:this.ax,ay:this.ay});remaining-=dt;steps++;
  }
  // Account for roundoff of the float32 requested timestep without an extra
  // transport whose increments would be below the precision of f itself.
  if(remaining<=tolerance){this.time+=remaining;remaining=0;}
  this.timeMs=performance.now()-started;return{advanced:duration-Math.max(0,remaining),steps};
 }
 diagnostics(){return {mass:this.totalMass,time:this.time,spatialLoss:this.spatialLoss,velocityLoss:this.velocityLoss,angularMomentum:this.angularMomentum};}
 async readLocalPhase(indices){
  if(indices.length>64)throw new Error('Observation locale limitée à 64 cellules spatiales.');
  const block=this.nv*this.nv*4,bytes=indices.length*block;if(!bytes)return new Map();
  if(!this.localReadback||this.localReadCapacity<bytes){this.localReadback?.destroy();this.localReadCapacity=bytes;this.localReadback=this.buffer('local phase observation',bytes,GPUBufferUsage.COPY_DST|GPUBufferUsage.MAP_READ);}
  const e=this.device.createCommandEncoder();indices.forEach((cell,i)=>e.copyBufferToBuffer(this.phase[0],cell*block,this.localReadback,i*block,block));this.device.queue.submit([e.finish()]);
  await this.localReadback.mapAsync(GPUMapMode.READ);const data=new Float32Array(this.localReadback.getMappedRange()).slice();this.localReadback.unmap();
  return new Map(indices.map((cell,i)=>[cell,data.subarray(i*this.nv*this.nv,(i+1)*this.nv*this.nv)]));
 }
 async snapshot(){const b=this.buffer('test phase readback',this.size*4,GPUBufferUsage.COPY_DST|GPUBufferUsage.MAP_READ);const e=this.device.createCommandEncoder();e.copyBufferToBuffer(this.phase[0],0,b,0,this.size*4);this.device.queue.submit([e.finish()]);await b.mapAsync(GPUMapMode.READ);const f=new Float32Array(b.getMappedRange()).slice();b.unmap();b.destroy();return f;}
 destroy(){for(const b of this.buffers||[])b.destroy();this.device?.destroy();}
}
root.GPUKineticMatter=GPUKineticMatter;
})(typeof window!=='undefined'?window:globalThis);
