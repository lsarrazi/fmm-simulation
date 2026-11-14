const assert=require('node:assert/strict');
const {initialBodies,TVDMatter,DirectBodies}=require('./galaxy-comparison.js');
const Base=require('./continuum-solver.js');
function error(a,b){return a.reduce((s,v,i)=>s+Math.abs(v-b[i]),0)/2;}
const initial=initialBodies(1024);
assert.deepEqual(initialBodies(1024),initial,'seeded repeatability');
const angular=b=>b.reduce((s,p)=>s+p.mass*(p.x*p.vy-p.y*p.vx),0);
assert(angular(initial)>.05,'galaxy carries original rotation');
assert.equal(angular(initialBodies(256,'cold')),0);
function run(n,nv,tvd=true){
 const g=new TVDMatter({n,nv});g.load(initial);
 if(!tvd)g.transport=Base.prototype.transport;
 const r=new DirectBodies(initial),initialError=error(g.mass,r.project(g));
 assert(initialError<1e-6,'same projected initial state');
 assert(Math.abs(g.diagnostics().angularMomentum-angular(initial))<1e-6,'same initial angular momentum');
 g.onTransport=({dt})=>r.advance(dt,{softening:.042});
 for(let i=0;i<20;i++)g.advance(.01,{softening:.042});
 assert(Math.abs(g.time-r.time)<1e-12,'synchronized time');
 assert(g.f.every(v=>Number.isFinite(v)&&v>=-1e-12),'positive finite distribution');
 const d=g.diagnostics();assert(Math.abs(d.mass+d.spatialLoss+d.velocityLoss-1)<2e-6,'mass balance');
 assert(Math.abs(angular(r.bodies)-angular(initial))<1e-12,'reference angular momentum conservation');
 return error(g.mass,r.project(g));
}
const coarse=run(20,25),fine=run(32,41),upwind=run(32,41,false);
console.log({time:.2,coarse,fine,upwind});
assert(fine<upwind,'TVD improves agreement relative to upwind at the same resolution');
// A galaxy changes its resolved structure when refining, so L1 at different grid
// scales cannot establish convergence. Verify convergence against an analytic
// translation with a fixed physical profile instead.
function translation(n) {
 const g=new TVDMatter({n,nv:5});g.f.fill(0);
 let total=0;const profile=(x,y)=>Math.exp(-((x+.4)**2+y*y)/.08);
 for(let y=0;y<n;y++)for(let x=0;x<n;x++) {
  const v=profile((x+.5)*g.dx-g.extent,(y+.5)*g.dx-g.extent);
  g.f[((y*n+x)*5+2)*5+3]=v;total+=v;
 }
 for(let i=0;i<g.size;i++)g.f[i]/=total;
 const duration=.3, speed=g.velocities[3];
 for(let i=0;i<30;i++)g.transport(0,duration/30);
 g.project();const exact=new Float64Array(n*n);let norm=0;
 for(let y=0;y<n;y++)for(let x=0;x<n;x++)norm+=exact[y*n+x]=profile((x+.5)*g.dx-g.extent-speed*duration,(y+.5)*g.dx-g.extent);
 for(let i=0;i<exact.length;i++)exact[i]/=norm;
 return error(g.mass,exact);
}
const e24=translation(24),e48=translation(48);assert(e48<e24*.6,'TVD converges on analytic smooth advection');
console.log({analyticAdvection24:e24,analyticAdvection48:e48});
console.log('Passed: common initialization and angular momentum, time synchronization, positivity, mass balance, direct reference torque, refinement and reduced diffusion.');
