const fs=require('node:fs'),assert=require('node:assert/strict');
const createFMM=require('./fmm-core'),{initialBodies,TVDMatter,FMMBodies}=require('./galaxy-comparison');
const source=fs.readFileSync('test-fixtures/fmm-original-reference.js','utf8');
const helpers=source.slice(source.indexOf('  function createLevel('),source.indexOf('  function computeFMM()'));
const kernel=source.slice(source.indexOf('  function gravityKernel('),source.indexOf('  function resetBodies('));
const computation=source.slice(source.indexOf('  function computeFMM()'),source.indexOf('    if (impulse) {',source.indexOf('  function computeFMM()')));
const original=new Function('bodies','options',`const NEAR_CELL_LIMIT=28,BOX_MIN=-1,BOX_SIZE=2;let treeCache=null,boundaryMode='open',torusActive=false;const gravityStrength=options.gravity,gravityExponent=options.exponent,gravitySoftening2=options.softening**2,els={boundary:{value:options.boundary},accuracy:{value:options.order}};${kernel}${helpers}${computation}}computeFMM();`);
for(const boundary of ['open','reflect','torus'])for(const order of [0,1]){
 const a=initialBodies(1024),b=initialBodies(1024),o={gravity:.48,exponent:2,softening:.042,boundary,order};
 original(a,o);createFMM().compute(b,o);
 for(let i=0;i<a.length;i++){assert.equal(b[i].ax,a[i].ax);assert.equal(b[i].ay,a[i].ay);}
}
const oldModule={exports:{}};
new Function('require','module','exports',fs.readFileSync('test-fixtures/galaxy-original-reference.js','utf8'))(require,oldModule,oldModule.exports);
const initial=initialBodies(256),old=new oldModule.exports.TVDMatter({n:16,nv:17}),current=new TVDMatter({n:16,nv:17});old.load(initial);current.load(initial);
for(let i=0;i<old.ax.length;i++){old.ax[i]=current.ax[i]=Math.sin(i)*.2;old.ay[i]=current.ay[i]=Math.cos(i)*.2;}
for(let pass=0;pass<3;pass++)for(const axis of [2,3,0,1]){old.transport(axis,.01);current.transport(axis,.01);assert.deepEqual(current.f,old.f,'identical optimized transport');}
const r=new FMMBodies(initial);r.advance(.01);assert.equal(r.time,.01);assert(r.bodies.every(p=>Number.isFinite(p.x)&&Number.isFinite(p.vx)));
const untouched=new FMMBodies(initial);untouched.advance(.01,{gravity:0});for(let i=0;i<initial.length;i++)assert.equal(untouched.bodies[i].vx,initial[i].vx);
// Warmed timing, same scenario, same resolution. No timing threshold in tests.
function measure(C){const g=new C({n:40,nv:41});g.load(initialBodies(2048));let elapsed=0;for(let i=0;i<6;i++){const start=performance.now();g.advance(.01,{softening:.042});if(i>=3)elapsed+=performance.now()-start;}return elapsed/3;}
const before=measure(oldModule.exports.TVDMatter),after=measure(TVDMatter);
console.log(JSON.stringify({gridBeforeMs:before,gridAfterMs:after,speedup:before/after}));
console.log('Passed: bit-exact shared FMM vs original in all boundaries and orders, bit-exact optimized TVD transport, reference advancement and zero gravity.');
