const assert = require("node:assert/strict");
const AnchorEngine = require("./src/poc/anchor-engine");
const {
  AnchorPopulation,
  normalCDF,
  normalQuantile,
  covariance,
} = require("./src/poc/anchor-population");

(async () => {
  for (const p of [0.00001, 0.01, 0.25, 0.5, 0.75, 0.99, 0.99999])
    assert(Math.abs(normalCDF(normalQuantile(p)) - p) < 1e-6);
  const cloud = {
    index: 0,
    x: 0.2,
    y: -0.1,
    vx: 0.01,
    vy: -0.02,
    A: new Float64Array([0.1, 0, 0, 0, 0, 0.1, 0, 0]),
    B: new Float64Array([0, 0, 0.002, 0, 0, 0, 0, 0.002]),
  };
  const model = { seed: 417, count: 1, time: 0, anchors: [cloud] },
    population = new AnchorPopulation(model, { population: 10n ** 24n });
  const id = population.universe + "/0/999999999999999999999999",
    first = population.resolve(id);
  population.cache.clear();
  assert.deepEqual(
    population.resolve(id),
    first,
    "cache eviction must not change an inferred object",
  );
  assert.throws(() =>
    population.resolve(population.universe + "/0/1000000000000000000000000"),
  );
  const before = population.resolve(id);
  cloud.x += 0.3;
  cloud.y -= 0.2;
  const after = population.resolve(id);
  assert(Math.abs(after.x - before.x - 0.3) < 1e-14);
  assert(Math.abs(after.y - before.y + 0.2) < 1e-14);
  const left = population.densityAt(0.499999, -0.3),
    right = population.densityAt(0.500001, -0.3);
  assert(
    Math.abs(left / right - 1) < 1e-6,
    "Gaussian clouds have no spatial step at their centres",
  );
  const projected = population.project({
    x: 0.5,
    y: -0.3,
    width: 2,
    resolution: 128,
  });
  assert(
    Math.abs(projected.mass - 1) < 1e-6,
    "Gaussian projection conserves mass when the cloud is enclosed",
  );
  const sample = population.sample({
    x: 0.51,
    y: -0.3,
    width: 0.04,
    budget: 128,
    maxVisits: 20000,
  });
  assert(
    sample.points.length > 80,
    "lazy local discovery should find the requested neighbourhood",
  );
  for (const p of sample.points) {
    assert(Math.abs(p.x - 0.51) <= 0.02 && Math.abs(p.y + 0.3) <= 0.02);
    assert.deepEqual(population.resolve(p.id), p);
  }
  assert.equal(
    new Set(sample.points.map((p) => p.id)).size,
    sample.points.length,
  );
  assert(
    sample.visited < 25000,
    "query work is bounded independently of theoretical population",
  );
  const moments = Array.from({ length: 4 }, () => ({ sum: 0, squared: 0 }));
  for (let j = 0; j < 4096; j++) {
    const rank = (BigInt(2 * j + 1) * population.population) / 8192n,
      point = population.evaluate(0, rank),
      z = [
        (point.x - cloud.x) / 0.1,
        (point.y - cloud.y) / 0.1,
        (point.vx - cloud.vx) / 0.002,
        (point.vy - cloud.vy) / 0.002,
      ];
    for (let k = 0; k < 4; k++) {
      moments[k].sum += z[k];
      moments[k].squared += z[k] * z[k];
    }
  }
  for (const m of moments) {
    assert(Math.abs(m.sum / 4096) < 0.04);
    assert(
      Math.abs(m.squared / 4096 - 1) < 0.06,
      "material addresses must reproduce the Gaussian phase-space distribution",
    );
  }
  // A source/target direct solver is injected solely as a small independent oracle.
  const force = {
    async compute(p, { gravity, softening }) {
      const result = new Float64Array(p.length * 2);
      for (let i = 0; i < p.length; i++)
        for (let j = 0; j < p.length; j++)
          if (i !== j && p[j].mass) {
            const x = p[j].x - p[i].x,
              y = p[j].y - p[i].y,
              q = x * x + y * y + softening ** 2,
              k = (gravity * p[j].mass) / q ** 1.5;
            result[2 * i] += x * k;
            result[2 * i + 1] += y * k;
          }
      return result;
    },
    destroy() {},
  };
  const engine = new AnchorEngine({
    count: 32,
    fieldSide: 8,
    force,
    maxStep: 0.05,
  });
  await engine.initialize();
  const inferred = new AnchorPopulation(engine),
    address = inferred.universe + "/0/12345";
  const p0 = inferred.resolve(address);
  await engine.advance(0.2);
  const p1 = inferred.resolve(address);
  assert.equal(p1.id, p0.id);
  assert.notEqual(p1.x, p0.x);
  assert(Math.abs(engine.time - 0.2) < 1e-12);
  assert.equal(engine.anchors.length, 32);
  assert(Math.abs(engine.anchors.reduce((s, p) => s + p.mass, 0) - 1) < 1e-12);
  assert(
    engine.anchors.every((p) => Number.isFinite(p.x) && covariance(p).det > 0),
  );
  assert(
    Math.abs(
      (engine.angularMomentum() - engine.initialAngular) /
        engine.initialAngular,
    ) < 1e-8,
    "central halo and pairwise gravity conserve angular momentum with KDK",
  );
  engine.destroy();
  assert.equal(engine.destroyed, true);
  console.log(
    "Passed: persistent Gaussian 4D identities to 10^24, lazy bounded local queries, smooth density, source/target units and KDK conservation.",
  );
})().catch((e) => {
  console.error(e);
  process.exitCode = 1;
});
