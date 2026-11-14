/* Local statistical observations at a fixed kinetic snapshot, not tracked orbits.
 * Counts and addresses use BigInt. Sub-cell positions are procedural assumptions.
 * Numerical time is integrated separately from the fine physical display clock.
 */
(function (root) {
  "use strict";
  const {
    viewportGeometry,
    renderContext,
    coveredCameraView,
    viewportProjection,
    cameraPresentation,
    committedCameraView,
    FramePublication,
  } =
    typeof module !== "undefined"
      ? require("./src/observation/camera-frame")
      : root;

  const NS = 1000000000n,
    YEAR = 31557600n,
    MODEL_NS = 1000000n * YEAR * NS,
    MIN_MODEL_STEP = 1e-6;
  const UNITS = { s: 1n, y: YEAR, My: 1000000n * YEAR, Gy: 1000000000n * YEAR };
  const UNIT_NS = {
    ns: 1n,
    us: 1000n,
    ms: 1000000n,
    s: NS,
    min: 60n * NS,
    h: 3600n * NS,
    d: 86400n * NS,
    w: 604800n * NS,
    mo: (YEAR * NS) / 12n,
    y: YEAR * NS,
    ky: 1000n * YEAR * NS,
    My: 1000000n * YEAR * NS,
    Gy: 1000000000n * YEAR * NS,
    Ty: 1000000000000n * YEAR * NS,
  };
  function unitNs(unit) {
    if (!UNIT_NS[unit]) throw Error("Unité inconnue.");
    return UNIT_NS[unit];
  }
  function gcd(a, b) {
    while (b) {
      const r = a % b;
      a = b;
      b = r;
    }
    return a;
  }
  function decimal(text) {
    text = String(text)
      .replace(/\s/g, "")
      .replace(/,/g, ".")
      .replace(/^\./, "0.");
    if (text.length > 160) throw Error("Valeur trop longue.");
    if (text.includes("/")) {
      const parts = text.split("/");
      if (parts.length !== 2) throw Error("Fraction a/b attendue.");
      const [a, b] = decimal(parts[0]),
        [c, d] = decimal(parts[1]);
      if (c === 0n) throw Error("Division par zéro.");
      return [a * d, b * c];
    }
    const m = text.match(/^(\d+)(?:\.(\d*))?(?:e([+-]?\d+))?$/i);
    if (!m)
      throw Error(
        "Décimal, notation scientifique ou fraction positive attendus.",
      );
    const exponent = Number(m[3] || 0) - (m[2] || "").length;
    if (Math.abs(exponent) > 24) throw Error("Exposant hors limites.");
    const v = BigInt(m[1] + (m[2] || ""));
    return exponent >= 0
      ? [v * 10n ** BigInt(exponent), 1n]
      : [v, 10n ** BigInt(-exponent)];
  }
  class SimulationClock {
    constructor() {
      this.integrated = 0n;
      this.debt = 0n;
      this.remainder = 0n;
      this.rate = MODEL_NS;
      this.rateNumerator = MODEL_NS;
      this.rateDenominator = 1n;
      this.actualRate = 0;
      this.limited = false;
      this.minimumStep = MIN_MODEL_STEP;
    }
    setRate(value, unit) {
      const [a, b] = decimal(value);
      this.configureRate(a * unitNs(unit), b);
    }
    setFactor(value) {
      const [a, b] = decimal(value);
      this.configureRate(a * MODEL_NS, b);
    }
    configureRate(a, b) {
      if (a <= 0n || b <= 0n)
        throw Error("Vitesse strictement positive attendue.");
      const divisor = gcd(a, b);
      this.rateNumerator = a / divisor;
      this.rateDenominator = b / divisor;
      this.rate = this.rateNumerator / this.rateDenominator;
      if (this.debt >= this.quantum()) this.debt = 0n;
      this.remainder = 0n;
    }
    factorRatio() {
      const a = this.rateNumerator,
        b = this.rateDenominator * MODEL_NS,
        g = gcd(a, b);
      return b / g === 1n ? (a / g).toString() : a / g + "/" + b / g;
    }
    rateNumber() {
      return Number(this.rateNumerator) / Number(this.rateDenominator);
    }
    quantity(unit) {
      const a = this.rateNumerator,
        b = this.rateDenominator * unitNs(unit),
        g = gcd(a, b);
      return b / g === 1n ? (a / g).toString() : a / g + "/" + b / g;
    }
    tick(ms) {
      const real = BigInt(Math.max(0, Math.round(ms * 1e6))),
        divisor = this.rateDenominator * NS;
      const numerator = this.rateNumerator * real + this.remainder;
      this.debt += numerator / divisor;
      this.remainder = numerator % divisor;
      const cap = MODEL_NS / 10n;
      this.limited = this.debt > cap;
      if (this.limited) this.debt = cap;
    }
    quantum() {
      return BigInt(
        Math.max(1, Math.round(this.minimumStep * Number(MODEL_NS))),
      );
    }
    duration(maxDuration = 0.025) {
      return this.debt >= this.quantum()
        ? Math.min(maxDuration, Number(this.debt) / Number(MODEL_NS))
        : 0;
    }
    consume(dt) {
      const ns = BigInt(Math.round(dt * Number(MODEL_NS)));
      this.integrated += ns;
      this.debt = this.debt > ns ? this.debt - ns : 0n;
    }
    fineOffset() {
      return this.debt < this.quantum()
        ? Number(this.debt) / Number(MODEL_NS)
        : 0;
    }
    visible() {
      return this.integrated + (this.debt < this.quantum() ? this.debt : 0n);
    }
    format(unit = "My") {
      const divisor = unitNs(unit),
        v = this.visible();
      const whole = v / divisor,
        rem = v % divisor;
      const digits = ((rem * 1000000000n) / divisor)
        .toString()
        .padStart(9, "0")
        .replace(/0+$/, "");
      return whole.toString() + (digits ? "." + digits : "") + " " + unit;
    }
  }
  function hash(a) {
    let x = a >>> 0;
    x = Math.imul(x ^ (x >>> 16), 0x7feb352d);
    x = Math.imul(x ^ (x >>> 15), 0x846ca68b);
    return (x ^ (x >>> 16)) >>> 0;
  }
  function randomID(id, salt) {
    const text = id.toString();
    let h = 2166136261 ^ salt;
    for (let i = 0; i < text.length; i++)
      h = Math.imul(h ^ text.charCodeAt(i), 16777619);
    return (hash(h) + 0.5) / 4294967296;
  }
  function fraction(count, x) {
    return (
      (count *
        BigInt(Math.floor(Math.max(0, Math.min(1, x)) * 9007199254740992))) /
      9007199254740992n
    );
  }
  function children(node) {
    const q = node.count / 4n,
      r = Number(node.count % 4n),
      result = [];
    let start = node.start;
    for (let k = 0; k < 4; k++) {
      const count = q + BigInt(k < r);
      result.push({
        ...node,
        start,
        count,
        x: node.x + ((k % 2) * node.w) / 2,
        y: node.y + (Math.floor(k / 2) * node.w) / 2,
        w: node.w / 2,
      });
      start += count;
    }
    return result;
  }
  function point(node, id) {
    let p = node;
    for (let depth = 0; p.count > 1n && depth < 60; depth++) {
      const a = children(p);
      p = a.find((c) => id >= c.start && id < c.start + c.count);
    }
    return {
      id: id.toString(),
      x: p.x + p.w * randomID(id, 11),
      y: p.y + p.w * randomID(id, 17),
      vx: node.vx || 0,
      vy: node.vy || 0,
    };
  }
  function localSample(
    grid,
    {
      x = 0,
      y = 0,
      zoom = 1,
      budget = 512,
      population = 10n ** 21n,
      phases = null,
      offset = 0,
      boundaryTolerance = null,
    } = {},
  ) {
    budget = Math.max(1, Math.min(2048, Math.floor(budget)));
    population = BigInt(population);
    const width = (2 * grid.extent) / zoom,
      view = {
        l: x - width / 2,
        r: x + width / 2,
        b: y - width / 2,
        t: y + width / 2,
      };
    const groups = [];
    let visited = 0,
      prefix = 0n,
      cumulative = 0;
    const intersects = (p) =>
      p.x < view.r && p.x + p.w > view.l && p.y < view.t && p.y + p.w > view.b;
    const contains = (p) =>
      p.x >= view.l &&
      p.x + p.w <= view.r &&
      p.y >= view.b &&
      p.y + p.w <= view.t;
    const visible = (p) =>
      p.x >= view.l && p.x <= view.r && p.y >= view.b && p.y <= view.t;
    const visit = (node, depth = 0) => {
      if (node.count === 0n || !intersects(node)) return;
      if (++visited > 24000) return;
      if (node.count === 1n) {
        const p = point(node, node.start);
        if (visible(p)) groups.push({ ...node, weight: 1n });
        return;
      }
      if (contains(node)) {
        groups.push({ ...node, weight: node.count });
        return;
      }
      if (node.w < (boundaryTolerance ?? width / 600) || depth >= 55) {
        const area =
          Math.max(
            0,
            Math.min(node.x + node.w, view.r) - Math.max(node.x, view.l),
          ) *
          Math.max(
            0,
            Math.min(node.y + node.w, view.t) - Math.max(node.y, view.b),
          );
        const weight = fraction(node.count, area / (node.w * node.w));
        if (weight) groups.push({ ...node, weight });
        return;
      }
      for (const c of children(node)) visit(c, depth + 1);
    };
    for (let cell = 0; cell < grid.mass.length; cell++) {
      cumulative += Math.max(0, grid.mass[cell]);
      const end = fraction(population, cumulative),
        count = end - prefix;
      const base = {
        start: prefix,
        count,
        x: (cell % grid.n) * grid.dx - grid.extent,
        y: Math.floor(cell / grid.n) * grid.dx - grid.extent,
        w: grid.dx,
      };
      prefix = end;
      // For fine-time interpolation, consider the velocity-expanded cell footprint.
      if (
        !intersects({
          ...base,
          x: base.x - grid.vmax * offset,
          y: base.y - grid.vmax * offset,
          w: base.w + 2 * grid.vmax * offset,
        })
      )
        continue;
      const phase = phases?.get(cell);
      if (phase) {
        let sum = 0;
        for (const m of phase) sum += Math.max(0, m);
        if (sum <= 0) continue;
        let start = base.start,
          acc = 0;
        for (let k = 0; k < phase.length; k++) {
          acc += Math.max(0, phase[k]);
          const stop = base.start + fraction(base.count, acc / sum);
          if (stop > start) {
            const vx = ((k % grid.nv) - (grid.nv - 1) / 2) * grid.dv,
              vy = (Math.floor(k / grid.nv) - (grid.nv - 1) / 2) * grid.dv;
            visit({
              ...base,
              start,
              count: stop - start,
              x: base.x + vx * offset,
              y: base.y + vy * offset,
              vx,
              vy,
            });
          }
          start = stop;
        }
      } else visit(base);
    }
    let total = 0n;
    for (const g of groups) {
      g.prefix = total;
      total += g.weight;
    }
    const points = [],
      seen = new Set(),
      wanted = Number(total < BigInt(budget) ? total : BigInt(budget));
    const attempts = total <= BigInt(budget) ? wanted : budget * 3;
    for (let j = 0; j < attempts && points.length < wanted; j++) {
      const rank =
        total <= BigInt(budget)
          ? BigInt(j)
          : fraction(total, (j * 0.6180339887498949 + 0.5) % 1);
      let lo = 0,
        hi = groups.length - 1;
      while (lo < hi) {
        const m = (lo + hi) >>> 1;
        if (groups[m].prefix + groups[m].weight <= rank) lo = m + 1;
        else hi = m;
      }
      const g = groups[lo];
      if (!g) break;
      const local = rank - g.prefix;
      const id =
        g.start +
        (g.weight === g.count
          ? local
          : fraction(g.count, randomID(local + g.start, 23)));
      if (seen.has(id.toString())) continue;
      const p = point(g, id);
      if (visible(p)) {
        points.push(p);
        seen.add(p.id);
      }
    }
    return { points, estimated: total, visited, truncated: visited > 24000 };
  }
  function physicalNs(value, unit) {
    const [a, b] = decimal(value);
    return (a * unitNs(unit)) / b;
  }
  function selectObservedSources(candidates, priorIds, budget) {
    const prior = new Set(priorIds),
      retained = [],
      fresh = [],
      rank = (a, b) => b.flux - a.flux || a.point.id.localeCompare(b.point.id);
    for (const entry of candidates)
      (prior.has(entry.point.id) ? retained : fresh).push(entry);
    retained.sort(rank);
    if (retained.length >= budget)
      return retained.slice(0, budget).map((entry) => entry.point);
    fresh.sort(rank);
    return retained
      .concat(fresh.slice(0, budget - retained.length))
      .map((entry) => entry.point);
  }
  // Resolve a field only when its expected population approaches the display
  // capacity. The density integral supplies N without enumerating virtual stars.
  // A C2 transition in log population avoids a hard zoom boundary; luminosity
  // and the camera's inverse-square flux still decide which sources appear first.
  function stellarDetailLevel(expected, budget, onset = 32) {
    if (!(budget > 0)) return 0;
    if (!(expected > 0)) return 1;
    const start = Math.log2(onset),
      end = Math.min(-1, start - 1),
      population = Math.log2(expected) - Math.log2(budget);
    const u = Math.max(0, Math.min(1, (start - population) / (start - end)));
    return u * u * u * (10 + u * (-15 + 6 * u));
  }
  // Density comes from the collective PDF, not the number of drawn sources.
  function adaptiveStellarBudget(
    density,
    reference,
    limit = Number.MAX_SAFE_INTEGER,
  ) {
    const ratio = reference > 0 ? Math.max(0, density / reference) : 0,
      target =
        ratio <= 1
          ? 4096 + 4096 * ratio
          : 8192 + 8192 * Math.min(1, (ratio - 1) / 3);
    return Math.max(0, Math.min(limit, Math.round(target)));
  }
  function readLightCalibration(
    budget,
    onset,
    limit = Number.MAX_SAFE_INTEGER,
  ) {
    const numeric = (value) =>
      String(value).trim() === ""
        ? NaN
        : Number(String(value).trim().replace(",", "."));
    const b = numeric(budget),
      o = numeric(onset),
      fail = (message, field) => {
        const error = Error(message);
        error.field = field;
        throw error;
      };
    if (!Number.isSafeInteger(b) || b < 0)
      fail("Le plafond doit être un entier positif ou nul.", "budget");
    if (!Number.isFinite(o) || !(o > 0))
      fail(
        "Le facteur d’apparition doit être un nombre strictement positif.",
        "onset",
      );
    if (b > limit)
      fail(
        "Le rendu GPU courant accepte au plus " +
          limit.toLocaleString("fr-FR") +
          " sources par image.",
        "budget",
      );
    return { budget: b, onset: o };
  }
  class SmoothZoom {
    constructor(value = 1, duration = 220) {
      this.value = this.target = value;
      this.duration = duration;
      this.active = false;
      this.anchor = null;
    }
    to(value, now, anchor = null) {
      this.at(now);
      this.from = Math.log(this.value);
      this.target = value;
      this.started = now;
      this.anchor = anchor;
      this.active = value !== this.value;
    }
    at(now) {
      if (!this.active) return this.value;
      const t = Math.max(0, Math.min(1, (now - this.started) / this.duration));
      this.value =
        t === 1
          ? this.target
          : Math.exp(
              this.from +
                (Math.log(this.target) - this.from) * (1 - (1 - t) ** 3),
            );
      if (t === 1) this.active = false;
      return this.value;
    }
    cancel(value = this.value) {
      this.value = this.target = value;
      this.active = false;
      this.anchor = null;
    }
  }
  // Navigation aid, independent of mass labels, reflected light and exposure.
  // Reference-pixel angular radius keeps discovery consistent across displays.
  function systemBodyGuide(
    radiusKm,
    widthAU,
    { selected = false, context = null } = {},
  ) {
    const radius =
        (radiusKm / 149597870.7 / widthAU) * (context?.cssSpan ?? 600),
      t = Math.max(
        0,
        Math.min(1, Math.log10(Math.max(radius, 1e-30) / 0.0001)),
      ),
      opacity = selected ? 1 : t * t * (3 - 2 * t);
    return {
      radius,
      opacity,
      ringRadius: Math.max(selected ? 6 : 3, radius + 2),
    };
  }
  class FrameRateMeter {
    constructor() {
      this.start = null;
      this.last = null;
      this.frames = 0;
      this.fps = 0;
    }
    tick(now) {
      if (this.start === null) {
        this.start = now;
        this.last = now;
        return;
      }
      this.frames++;
      this.last = now;
      if (now - this.start >= 500) {
        this.fps = (this.frames * 1000) / (now - this.start);
        this.frames = 0;
        this.start = now;
      }
    }
    rate(now) {
      return this.last === null || now - this.last > 1000 ? 0 : this.fps;
    }
    reset() {
      this.start = this.last = null;
      this.frames = 0;
      this.fps = 0;
    }
  }
  const api = {
    SimulationClock,
    localSample,
    MODEL_NS,
    UNITS,
    UNIT_NS,
    unitNs,
    MIN_MODEL_STEP,
    physicalNs,
    selectObservedSources,
    stellarDetailLevel,
    readLightCalibration,
    adaptiveStellarBudget,
    SmoothZoom,
    FrameRateMeter,
    FramePublication,
    viewportGeometry,
    renderContext,
    coveredCameraView,
    viewportProjection,
    cameraPresentation,
    committedCameraView,
    systemBodyGuide,
  };
  Object.assign(root, api);
  if (typeof module !== "undefined") module.exports = api;
})(typeof window !== "undefined" ? window : globalThis);
