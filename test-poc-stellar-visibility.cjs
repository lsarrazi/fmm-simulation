const assert = require("node:assert/strict");
const { renderContext, viewportGeometry } = require("./observation-tools");
const L = require("./stellar-light");
const V = require("./src/poc/stellar-visibility");
const context = renderContext(viewportGeometry(1180, 856), {
  worldShortSpan: 0.004,
  sigmaCss: 1,
  lengthPc: 10000,
});
const optics = V.options({
  population: 100000000000n,
  context,
  exposureEV: -4.5,
  onset: 160,
});
const star = (id, luminosity, x = 0) => ({
  id,
  x,
  y: 0,
  photometry: {
    luminosity,
    visibleLuminosity: luminosity,
    radius: 1,
    color: [1, 1, 1],
  },
});
const target = star("target", 0.03);
const projection = (expected) => ({
  width: context.worldSpan,
  cameraWidth: context.worldShortSpan,
  left: -context.worldSpan / 2,
  bottom: -context.worldSpan / 2,
  centerX: 0,
  centerY: 0,
  resolution: 1,
  image: new Float64Array(1),
  mass: expected / 1e11,
  viewMass: expected / 1e11,
  expected,
});
let peak;
for (const expected of [831642, 1211047, 2650285, 1e9]) {
  const p = projection(expected),
    visible = V.select(p, [target], { ...optics, budget: 8192 });
  assert.equal(
    visible.length,
    1,
    "density increase must not hide a physically detectable star",
  );
  peak ??= visible[0].peak;
  assert.equal(
    visible[0].peak,
    peak,
    "same camera and star, same sensor signal",
  );
  assert.equal(
    L.minimumLuminosity(p, optics),
    L.minimumLuminosity(projection(1), optics),
  );
}
assert.equal(
  V.select(projection(1e6), [target], {
    ...V.options({ ...optics, onset: 80 }),
    budget: 8192,
  }).length,
  0,
  "sensitivity affects detection, not intrinsic brightness",
);
const a = star("retained", 1),
  b = star("new-bright", 2);
assert.equal(
  V.select(projection(1e6), [a, b], {
    ...optics,
    budget: 1,
    priorIds: ["retained"],
  })[0].id,
  "retained",
);
assert.equal(
  V.select(projection(1e6), [a, b], { ...optics, budget: 1, priorIds: [] })[0]
    .id,
  "new-bright",
  "zoom-out can select the brightest without retention",
);
const near = star("near", 1),
  far = star("far", 1, 0.002);
assert.equal(
  V.select(projection(1e6), [far, near], { ...optics, budget: 1 })[0].id,
  "near",
  "rank by received flux, including camera distance",
);
const out = V.history([a], { enabled: true, zoom: 900, previousZoom: 1000 });
assert.deepEqual(out.knownIds, ["retained"]);
assert.deepEqual(
  out.priorIds,
  [],
  "zoom-out removes priority, not knowledge of real candidates",
);
const known = new Map([[a.id, a]]);
const candidates = [b, ...out.knownIds.map((id) => known.get(id))];
assert.deepEqual(
  new Set(
    V.select(projection(1e6), candidates, {
      ...optics,
      budget: 2,
      priorIds: out.priorIds,
    }).map((p) => p.id),
  ),
  new Set(["retained", "new-bright"]),
  "under budget, a still-detectable known source survives an incomplete fresh discovery",
);
assert.equal(
  V.select(projection(1e6), candidates, {
    ...optics,
    budget: 1,
    priorIds: out.priorIds,
  })[0].id,
  "new-bright",
  "under budget pressure, zoom-out does not favour history",
);
assert.equal(
  V.select(projection(1e6), [b, star("retained", 1, 1)], {
    ...optics,
    budget: 2,
    priorIds: out.priorIds,
  }).length,
  1,
  "out-of-frame sources are discarded",
);
assert.equal(
  V.select(projection(1e6), [b, star("retained", 0.000001)], {
    ...optics,
    budget: 2,
    priorIds: out.priorIds,
  }).length,
  1,
  "sources below the detector threshold are discarded",
);
assert.deepEqual(
  V.history([a], { enabled: false, zoom: 1, previousZoom: 2 }).knownIds,
  [],
);
console.log(
  "Passed: density-independent stellar flux and detection, sensitivity, apparent-flux priority and optional source retention.",
);
