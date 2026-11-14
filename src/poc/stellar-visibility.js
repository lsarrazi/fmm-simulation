/* Detector sensitivity is independent of the changing population in the frame.
 * Brightness follows photometry and distance; density only affects discovery
 * workload and the unresolved background, never an individual star's flux.
 */
(function (root) {
  "use strict";
  const light =
    typeof module !== "undefined" ? require("../../stellar-light") : root;
  const observation =
    typeof module !== "undefined" ? require("../../observation-tools") : root;
  function options({ population, context, exposureEV, onset = 160 }) {
    if (!Number.isFinite(onset) || onset <= 0)
      throw Error("Sensibilité d’apparition invalide.");
    return {
      population,
      context,
      exposureEV,
      gain: 1,
      lengthPc: 10000,
      minimumPeak: (0.000125 * 160) / onset,
      detailWeight: 1,
    };
  }
  function select(
    projection,
    points,
    { budget, priorIds = [], selectedId = null, ...optics },
  ) {
    const plan = light.lightPlan(projection, points, {
      ...optics,
      detailWeight: 1,
    });
    const selected = observation.selectObservedSources(
      plan.sources.map((point) => ({ point, flux: point.visibleFluxWm2 })),
      priorIds,
      budget,
    );
    const inspected = selectedId
      ? plan.sources.find((p) => p.id === selectedId)
      : null;
    if (inspected && budget > 0 && !selected.some((p) => p.id === selectedId)) {
      if (selected.length >= budget) selected.pop();
      selected.push(inspected);
    }
    return selected;
  }
  root.POCStellarVisibility = { options, select };
  if (typeof module !== "undefined") module.exports = { options, select };
})(typeof window !== "undefined" ? window : globalThis);
