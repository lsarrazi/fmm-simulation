(function (root) {
  "use strict";
  const { adaptiveStellarBudget } =
    typeof module !== "undefined" ? require("../../observation-tools") : root;

  /** Pure observation policy: no DOM reads, physics mutation or rendering. */
  function stellarSourceBudget(
    field,
    {
      lengthPc,
      population,
      calibration,
      extent,
      adaptive,
      manualBudget,
      sourceLimit,
    },
  ) {
    const widthPc = field.width * lengthPc;
    const density =
      field.expected /
      Math.max(
        widthPc * widthPc * (field.aspectX ?? 1) * (field.aspectY ?? 1),
        1e-300,
      );
    const reference = calibration?.surfaceSolar
      ? (calibration.surfaceSolar * Number(population)) /
        Number(calibration.stars)
      : Number(population) / (2 * extent * lengthPc) ** 2;
    return {
      densityRatio: reference > 0 ? density / reference : 0,
      budget: adaptive
        ? adaptiveStellarBudget(density, reference, sourceLimit)
        : manualBudget,
    };
  }
  const api = { stellarSourceBudget };
  root.StellarObservationPolicy = api;
  if (typeof module !== "undefined") module.exports = api;
})(typeof window !== "undefined" ? window : globalThis);
