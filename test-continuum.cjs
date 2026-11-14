const assert = require('node:assert/strict');
const KineticMatter = require('./continuum-solver.js');
function check(s, label) {
  const d = s.diagnostics();
  assert(s.f.every(x => Number.isFinite(x) && x >= 0), label + ': positive finite distribution');
  assert(Math.abs(d.mass + d.spatialLoss + d.velocityLoss - 1) < 3e-6, label + ': mass accounting');
  return d;
}
const cold = new KineticMatter({ n: 16, nv: 15, sigma: 0, fluctuation: 0 });
let d = check(cold, 'initial');
assert(Math.abs(d.mass - 1) < 1e-6);
assert.equal(d.kinetic, 0);
assert.equal(d.angularMomentum, 0);
const initial = cold.f.slice();
cold.advance(0.1, { gravity: 0 });
assert.deepEqual(cold.f, initial, 'cold cloud stays still without force');
function radius2(s) {
  let r = 0;
  for (let y = 0; y < s.n; y++) for (let x = 0; x < s.n; x++) {
    const px = (x + .5) * s.dx - s.extent, py = (y + .5) * s.dx - s.extent;
    r += s.mass[y * s.n + x] * (px * px + py * py);
  }
  return r;
}
const beforeRadius = radius2(cold);
for (let i = 0; i < 30; i++) cold.advance(.01, { gravity: .48, softening: .15 });
d = check(cold, 'gravitational collapse');
assert(radius2(cold) < beforeRadius, 'cold cloud contracts');
assert(Math.abs(d.momentumX) < 1e-6 && Math.abs(d.momentumY) < 1e-6, 'symmetric cloud has no net drift');
assert(Math.abs(d.angularMomentum) < 1e-6, 'no generated global rotation');
const warm = new KineticMatter({ n: 16, nv: 15, sigma: .45, fluctuation: .3 });
const beforeWarm = radius2(warm);
for (let i = 0; i < 30; i++) warm.advance(.02, { gravity: 0 });
check(warm, 'free transport');
assert(radius2(warm) > beforeWarm, 'warm cloud expands without gravity');
// Deliberately narrow velocity domain: a strong field must report truncation losses.
const limited = new KineticMatter({ n: 12, nv: 11, vmax: .2, sigma: 0 });
for (let i = 0; i < 20; i++) limited.advance(.01, { gravity: 1, softening: .15 });
assert(check(limited, 'velocity outflow').velocityLoss > .001);
// Position outflow also contributes to the mass balance, without rescaling f.
const outflow = new KineticMatter({ n: 12, nv: 11, sigma: .6, radius: 1.2 });
for (let i = 0; i < 25; i++) outflow.advance(.04, { gravity: 0 });
assert(check(outflow, 'spatial outflow').spatialLoss > .001);
console.log('Passed: normalization, positivity, no-force rest, contraction, symmetry, free expansion, spatial and velocity outflow.');
const bench = new KineticMatter({ n: 28, nv: 25 });
const start = performance.now(); bench.advance(.01);
console.log(`Default grid: ${bench.size} cells; one 0.01 step: ${(performance.now() - start).toFixed(1)} ms.`);
