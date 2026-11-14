const assert = require('node:assert/strict');
const KineticMatter = require('./continuum-solver.js');
const Lazy = require('./test-fixtures/legacy-identities.js');
const solver = new KineticMatter({ n: 12, nv: 11, sigma: .25 });
const lazy = new Lazy(solver, { population: 10n ** 60n, cacheLimit: 128 });
solver.onTransport = r => lazy.record(r);
const largeId = '987654321012345678909876543210123456789098765432101234567890';
const comparable = p => { const { work, ...state } = p; return state; };
const first = lazy.resolve(largeId);
assert.equal(first.id, largeId, 'exact IDs beyond double precision');
assert.notEqual(lazy.resolve('9007199254740992').id, lazy.resolve('9007199254740993').id);
assert.throws(() => lazy.resolve('-1'));
assert.throws(() => lazy.resolve((10n ** 60n).toString()));
assert.deepEqual(lazy.candidates(20), lazy.candidates(100).slice(0, 20), 'stable candidate prefix');
assert.equal(new Set(lazy.candidates(100)).size, 100, 'unique candidate IDs');
let online;
for (let i = 0; i < 35; i++) { solver.advance(.015); online = lazy.resolve(largeId); }
lazy.forget();
const replay = lazy.resolve(largeId);
assert.deepEqual(comparable(replay), comparable(online), 'replay matches incrementally observed trajectory');
lazy.forget();
let partial = lazy.resolve(largeId, 3);
assert(!partial.ready && partial.work === 3, 'bounded replay work');
while (!partial.ready) partial = lazy.resolve(largeId, 3);
assert.deepEqual(comparable(partial), comparable(online), 'chunking does not change identity');
for (const id of lazy.candidates(200)) lazy.resolve(id);
assert(lazy.cache.size <= 128, 'bounded cache');
assert.deepEqual(comparable(lazy.resolve(largeId)), comparable(online), 'eviction does not change identity');
// Sampler observes copies: replay, cache and population do not perturb physical state.
const physics = solver.f.slice();
lazy.forget(); for (const id of lazy.candidates(50)) lazy.resolve(id);
assert.deepEqual(solver.f, physics);
const cappedSolver = new KineticMatter({ n: 8, nv: 9 });
const cap = new Lazy(cappedSolver, { maxHistoryBytes: 0 });
cappedSolver.onTransport = r => cap.record(r);
const before = cappedSolver.f.slice();
assert.equal(cappedSolver.advance(.1).advanced, 0);
assert.deepEqual(cappedSolver.f, before, 'history budget stops before altering physics');
assert(cap.full);
// Ensemble must follow the same upwind transition probabilities as f, not ballistic guesses.
const free = new KineticMatter({ n: 10, nv: 9, sigma: .4 });
const ensemble = new Lazy(free, { population: 1000000000n, cacheLimit: 32 });
free.onTransport = r => ensemble.record(r);
for (let i = 0; i < 12; i++) free.advance(.025, { gravity: 0 });
let sx = 0, sy = 0, radius2 = 0, alive = 0;
const count = 12000;
for (const id of ensemble.candidates(count)) {
  const p = ensemble.resolve(id);
  if (p.alive) { sx += p.x; sy += p.y; radius2 += p.x * p.x + p.y * p.y; alive++; }
}
let expectedR2 = 0, expectedMass = 0;
for (let y = 0; y < free.n; y++) for (let x = 0; x < free.n; x++) {
  const m = free.mass[y * free.n + x], px = (x + .5) * free.dx - free.extent, py = (y + .5) * free.dx - free.extent;
  expectedMass += m; expectedR2 += m * (px * px + py * py + 2 * (.75 * free.dx) ** 2 / 12);
}
assert(Math.abs(alive / count - expectedMass) < .025, 'survival frequency matches f');
assert(Math.abs(radius2 / count - expectedR2) < .025, 'spatial second moment matches f');
assert(Math.abs(sx / count) < .06 && Math.abs(sy / count) < .06, 'no spurious collective drift');
console.log('Passed: huge exact IDs, stable candidates, replay, chunking, eviction, bounded memory, no added mass, history cap and ensemble consistency.');
