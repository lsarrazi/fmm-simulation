/* Bounded local refinement of the collective field, using Orbit's shared FMM.
 * Actual stellar masses, complete (unfiltered) guard patch, fine-minus-coarse
 * force correction. Outside the patch only the collective tide is represented.
 * Stored corrections are physical encounter states, never visibility retention.
 */
(function(root){
'use strict';
const createFMM=typeof module!=='undefined'?require('./fmm-core'):root.createOrbitFMMField;
const Population=typeof module!=='undefined'?require('./stellar-population'):root.StellarPopulation;
class LocalStellarDynamics{
 constructor({budget=2048,maxRecords=16384}={}){this.budget=budget;this.maxRecords=maxRecords;this.records=new Map();this.active=[];this.fine=createFMM();this.coarse=createFMM();this.patch=null;this.time=0;this.enabled=true;this.reason='Zoomer pour activer la FMM locale.';this.cost=0;this.steps=0;this.stepLimit=Infinity;}
 async attachGPU(device){this.gpu=await root.createOrbitGPUForceField({device,capacity:8192});this.budget=8192;}
 destroy(){this.gpu?.destroy();this.gpu=null;}
 resolve(materials,id,offset=0){return this.correct(materials.resolve(id,offset),materials.time+offset);}
 correct(p,time){const r=this.records.get(p.id);if(!r)return p;const t=Math.max(0,time-r.time);return{...p,x:p.x+r.dx+r.dvx*t+.5*r.ax*t*t,y:p.y+r.dy+r.dvy*t+.5*r.ay*t*t,vx:p.vx+r.dvx+r.ax*t,vy:p.vy+r.dvy+r.ay*t};}
 visible(materials,camera,extent,offset=0){const half=extent/camera.zoom,points=[];for(const id of this.active){const p=this.resolve(materials,id,offset);if(Math.abs(p.x-camera.x)<=half&&Math.abs(p.y-camera.y)<=half)points.push(p);}return points;}
 settle(r,time){const t=Math.max(0,time-r.time);r.dx+=r.dvx*t+.5*r.ax*t*t;r.dy+=r.dvy*t+.5*r.ay*t*t;r.dvx+=r.ax*t;r.dvy+=r.ay*t;r.time=time;}
 deactivate(time,reason){for(const id of this.active){const r=this.records.get(id);this.settle(r,time);r.ax=r.ay=0;}this.active=[];this.patch=null;this.stepLimit=Infinity;this.reason=reason;}
 prepare(materials,grid,camera,{gravity,totalMass,lengthPc=10000,softening=.015,guard=4,enabled=true}={}){
  this.enabled=enabled;const time=materials.time;
  if(!enabled||!(totalMass>0)){this.deactivate(time,'FMM locale inactive pour cette configuration.');return;}
  // Reuse complete guard membership while the view remains in its full-force
  // interior. Force integration still runs each step; only redundant discovery
  // and the second evaluation of an already-current force are avoided.
  const previous=this.patch;if(previous&&previous.guard===guard&&previous.gravity===gravity&&previous.totalMass===totalMass&&previous.lengthPc===lengthPc&&previous.softening===softening&&previous.view===2*grid.extent/camera.zoom&&time>=previous.checkedTime&&time<previous.refreshAfter&&Math.max(Math.abs(camera.x-previous.x),Math.abs(camera.y-previous.y))+previous.view/2<previous.width*.275){this.updateLimit(materials);this.reason='FMM locale active';return;}
  const ix=Math.max(0,Math.min(grid.n-1,Math.floor((camera.x+grid.extent)/grid.dx))),iy=Math.max(0,Math.min(grid.n-1,Math.floor((camera.y+grid.extent)/grid.dx))),density=grid.mass?Number(materials.population)*grid.mass[iy*grid.n+ix]/grid.dx**2:0;
  const view=2*grid.extent/camera.zoom,minimum=(density>0?Math.min(10/lengthPc,Math.sqrt(384/density)):1/lengthPc)*guard/4,width=Math.max(guard*view,minimum),half=width/2;
  if(density*width*width>this.budget*1.5){this.deactivate(time,'Zone locale trop peuplée : zoomer davantage.');return;}
  // All stars, including faint ones. Reject incomplete discovery instead of
  // assigning real masses to a brightness-selected or weighted sample.
  const query=materials.sample({x:camera.x,y:camera.y,zoom:2*grid.extent/width,budget:this.budget,offset:0,exact:true});
  if(query.truncated||query.estimated>BigInt(this.budget)||BigInt(query.points.length)!==query.estimated){this.deactivate(time,'Zone locale trop peu résolue ou trop peuplée : zoomer davantage.');return;}
  const ids=new Map(query.points.map(p=>[p.id,p]));
  for(const [id]of this.records){const p=this.resolve(materials,id);if(Math.abs(p.x-camera.x)<half&&Math.abs(p.y-camera.y)<half)ids.set(id,p);}
  if(ids.size>this.budget){this.deactivate(time,'Zone de garde trop peuplée pour le budget physique.');return;}
  const additions=Array.from(ids.keys()).filter(id=>!this.records.has(id));if(this.records.size+additions.length>this.maxRecords){this.deactivate(time,'Budget des trajectoires locales atteint ; états conservés.');return;}
  for(const id of this.active)if(!ids.has(id)){const r=this.records.get(id);this.settle(r,time);r.ax=r.ay=0;}
  for(const id of additions)this.records.set(id,{dx:0,dy:0,dvx:0,dvy:0,ax:0,ay:0,time,mass:Population.properties(id).mass});
  this.active=Array.from(ids.keys()).sort();this.patch={x:camera.x,y:camera.y,width,view,guard,gravity,totalMass,lengthPc,softening,checkedTime:time,refreshAfter:time};this.time=time;
  this.reason='FMM locale active';const pending=this.accelerations(materials,grid,0);if(pending?.then)return pending.then(()=>this.updateLimit(materials,true));this.updateLimit(materials,true);
 }
 field(grid,x,y){const u=(x+grid.extent)/grid.dx-.5,v=(y+grid.extent)/grid.dx-.5,ix=Math.max(0,Math.min(grid.n-2,Math.floor(u))),iy=Math.max(0,Math.min(grid.n-2,Math.floor(v))),tx=Math.max(0,Math.min(1,u-ix)),ty=Math.max(0,Math.min(1,v-iy)),i=iy*grid.n+ix;
  const read=a=>(a[i]*(1-tx)+a[i+1]*tx)*(1-ty)+(a[i+grid.n]*(1-tx)+a[i+grid.n+1]*tx)*ty;return[read(grid.ax),read(grid.ay)];
 }
 accelerations(materials,grid,offset){if(!this.patch||!this.active.length)return;const patch=this.patch,mean=Population.meanMass,bases=this.active.map(id=>materials.resolve(id,offset)),states=bases.map(p=>this.correct(p,materials.time+offset));
  let xmin=Infinity,xmax=-Infinity,ymin=Infinity,ymax=-Infinity;for(const p of states){xmin=Math.min(xmin,p.x);xmax=Math.max(xmax,p.x);ymin=Math.min(ymin,p.y);ymax=Math.max(ymax,p.y);}const ox=(xmin+xmax)/2,oy=(ymin+ymax)/2,scale=Math.max((xmax-xmin)/2,(ymax-ymin)/2,4e-4/patch.lengthPc);
  const particles=states.map(p=>{const r=this.records.get(p.id),d=Math.max(Math.abs(p.x-patch.x),Math.abs(p.y-patch.y))/(patch.width/2),u=Math.max(0,Math.min(1,(.98-d)/.43)),window=u*u*(3-2*u);return{id:p.id,x:(p.x-ox)/scale,y:(p.y-oy)/scale,window,mass:r.mass/mean*window};});
  const gravity=patch.gravity*mean/(patch.totalMass*scale**3),eps=1e-4/patch.lengthPc/scale;
  const apply=correction=>{for(let i=0;i<particles.length;i++){const p=particles[i],id=p.id,r=this.records.get(id),base=bases[i],actual=states[i],a=this.field(grid,actual.x,actual.y),b=this.field(grid,base.x,base.y);r.ax=correction[i*2]*scale*p.window+a[0]-b[0];r.ay=correction[i*2+1]*scale*p.window+a[1]-b[1];}};
  if(this.gpu)return this.gpu.compute(particles,{gravity,softening:eps,coarseSoftening:patch.softening/scale}).then(apply);
  const report=this.fine.compute(particles,{gravity,softening:eps,order:1,nearLimit:Infinity,minDepth:0}),fine=particles.map(p=>[p.ax*scale,p.ay*scale]);
  this.coarse.compute(particles,{gravity,softening:patch.softening/scale,order:1,nearLimit:Infinity,minDepth:0});
  for(let i=0;i<particles.length;i++){const p=particles[i],id=p.id,r=this.records.get(id),base=bases[i],actual=states[i],a=this.field(grid,actual.x,actual.y),b=this.field(grid,base.x,base.y);r.ax=(fine[i][0]-p.ax*scale)*p.window+a[0]-b[0];r.ay=(fine[i][1]-p.ay*scale)*p.window+a[1]-b[1];}
  this.operations=report?.operations||0;
 }
 updateLimit(materials,refresh=false){if(!this.patch||!this.active.length){this.stepLimit=Infinity;return;}let vmax=0,absoluteSpeed=0,amax=0;const origin=materials.resolve(this.active[0]),margin=(this.patch.width-this.patch.view)/2;
  for(const id of this.active){const p=this.resolve(materials,id),r=this.records.get(id);absoluteSpeed=Math.max(absoluteSpeed,Math.hypot(p.vx,p.vy));vmax=Math.max(vmax,Math.hypot(p.vx-origin.vx,p.vy-origin.vy));amax=Math.max(amax,Math.hypot(r.ax,r.ay));}
  // Guard crossing and acceleration limits. No wall-clock speedup skips steps.
  const refreshAt=materials.time+Math.min(.01,margin/(32*Math.max(absoluteSpeed,1e-12)),.1*Math.sqrt(margin/Math.max(amax,1e-30)));this.patch.refreshAfter=refresh?refreshAt:Math.min(this.patch.refreshAfter,refreshAt);
  this.stepLimit=Math.min(.005,margin/(8*Math.max(vmax,1e-12)),.03*Math.sqrt((1e-4/this.patch.lengthPc)/Math.max(amax,1e-30)));
 }
 advance(dt,materials,grid){if(this.gpu)return this.advanceGPU(dt,materials,grid);if(!this.patch||!this.active.length)return;const started=performance.now();let t=0;this.steps=0;
  while(t<dt-1e-15){if(this.steps>=128)throw Error('Pas local trop long : réduire la cadence.');const h=Math.min(dt-t,this.stepLimit),now=materials.time+t;
   for(const id of this.active){const r=this.records.get(id);this.settle(r,now);r.dvx+=.5*r.ax*h;r.dvy+=.5*r.ay*h;r.dx+=r.dvx*h;r.dy+=r.dvy*h;r.time=now+h;r.ax=r.ay=0;}
   t+=h;this.accelerations(materials,grid,t);
   for(const id of this.active){const r=this.records.get(id);r.dvx+=.5*r.ax*h;r.dvy+=.5*r.ay*h;}
   this.steps++;
  }this.time=materials.time+dt;this.cost=performance.now()-started;
 }
 async advanceGPU(dt,materials,grid){if(!this.patch||!this.active.length)return;const started=performance.now();let t=0;this.steps=0;while(t<dt-1e-15){if(this.steps>=128)throw Error('Pas local trop long : réduire la cadence.');const h=Math.min(dt-t,this.stepLimit),now=materials.time+t;for(const id of this.active){const r=this.records.get(id);this.settle(r,now);r.dvx+=.5*r.ax*h;r.dvy+=.5*r.ay*h;r.dx+=r.dvx*h;r.dy+=r.dvy*h;r.time=now+h;r.ax=r.ay=0;}t+=h;await this.accelerations(materials,grid,t);for(const id of this.active){const r=this.records.get(id);r.dvx+=.5*r.ax*h;r.dvy+=.5*r.ay*h;}this.steps++;}this.time=materials.time+dt;this.cost=performance.now()-started;}
 status(){return{active:!!this.patch,count:this.active.length,backend:this.gpu?'GPU':'CPU',budget:this.budget,guard:this.patch?.guard||0,widthPc:this.patch?this.patch.width*this.patch.lengthPc:0,records:this.records.size,ms:this.cost,steps:this.steps,reason:this.reason};}
}
root.LocalStellarDynamics=LocalStellarDynamics;if(typeof module!=='undefined')module.exports=LocalStellarDynamics;
})(typeof window!=='undefined'?window:globalThis);
