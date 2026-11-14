/* Lazy numerical witnesses for the kinetic grid, not additional gravitating matter.
 * An ID selects an initial phase-space cell from f0. Counter-based random transitions
 * reproduce each conservative upwind operator in expectation. Replay depends on the
 * recorded field, timestep order, seed and ID; never on observation or cache order.
 */
(function (root) {
  "use strict";
  function mix(x) {
    x = Math.imul(x ^ x >>> 16, 0x7feb352d);
    x = Math.imul(x ^ x >>> 15, 0x846ca68b);
    return (x ^ x >>> 16) >>> 0;
  }
  function keys(text) {
    let a = 2166136261, b = 3335557771;
    for (let i = 0; i < text.length; i++) {
      a = Math.imul(a ^ text.charCodeAt(i), 16777619);
      b = Math.imul(b ^ text.charCodeAt(i), 2246822519);
    }
    return [mix(a), mix(b)];
  }
  function random(k, tick, axis) {
    return mix(k[0] ^ mix(k[1] ^ Math.imul(tick + 1, 0x9e3779b1)) ^ Math.imul(axis + 7, 0x85ebca6b)) / 4294967296;
  }
  class LazyMatterIdentities {
    constructor(solver, { population = 10n ** 24n, seed = "orbit-1", maxHistoryBytes = 32 * 1024 * 1024, cacheLimit = 2304 } = {}) {
      this.population = BigInt(population);
      if (this.population <= 0n) throw new Error("La population doit être positive.");
      this.n = solver.n; this.nv = solver.nv; this.dx = solver.dx; this.dv = solver.dv;
      this.extent = solver.extent; this.velocities = solver.velocities.slice();
      this.seed = seed; this.cacheLimit = cacheLimit;
      this.history = []; this.cache = new Map(); this.historyBytes = 0;
      this.recordBytes = solver.mass.length * 16 + 8;
      this.maxRecords = Math.floor(maxHistoryBytes / this.recordBytes);
      this.cdf = new Float64Array(solver.size);
      let mass = 0;
      for (let i = 0; i < solver.size; i++) this.cdf[i] = mass += solver.f[i];
      for (let i = 0; i < this.cdf.length; i++) this.cdf[i] /= mass;
      this.full = false;
    }
    canonical(id) {
      const text = String(id).trim();
      if (!/^\d{1,80}$/.test(text)) throw new Error("Identifiant entier positif attendu (80 chiffres maximum).");
      const value = BigInt(text);
      if (value < 0n || value >= this.population) throw new Error("Identifiant en dehors de la population choisie.");
      return value.toString();
    }
    record({ dt, ax, ay }) {
      if (this.history.length >= this.maxRecords) { this.full = true; return false; }
      this.history.push({ dt, ax: ax.slice(), ay: ay.slice() });
      this.historyBytes += this.recordBytes;
      return true;
    }
    create(id) {
      const key = keys(this.seed + ":" + id);
      const u = (Math.floor(random(key, 0, 9) * 0x200000) * 4294967296 + Math.floor(random(key, 0, 10) * 4294967296)) / 9007199254740992;
      let lo = 0, hi = this.cdf.length - 1;
      while (lo < hi) { const mid = (lo + hi) >>> 1; if (this.cdf[mid] <= u) lo = mid + 1; else hi = mid; }
      let q = lo;
      const vx = q % this.nv; q = Math.floor(q / this.nv);
      const vy = q % this.nv; q = Math.floor(q / this.nv);
      return { id, key, x: q % this.n, y: Math.floor(q / this.n), vx, vy, tick: 0, alive: true, exit: null };
    }
    transition(p, record, tick) {
      // Same order as the solver: kick vx, kick vy, drift x, drift y.
      const cell = p.y * this.n + p.x;
      for (const axis of [2, 3, 0, 1]) {
        const speed = axis === 2 ? record.ax[cell] : axis === 3 ? record.ay[cell] : axis === 0 ? this.velocities[p.vx] : this.velocities[p.vy];
        const probability = Math.abs(speed) * record.dt / (axis < 2 ? this.dx : this.dv);
        if (random(p.key, tick, axis) < probability) {
          const coordinate = ["x", "y", "vx", "vy"][axis];
          p[coordinate] += speed > 0 ? 1 : -1;
          if (p[coordinate] < 0 || p[coordinate] >= (axis < 2 ? this.n : this.nv)) {
            p.alive = false; p.exit = axis < 2 ? "espace" : "vitesses"; return;
          }
        }
      }
    }
    resolve(id, maxSteps = Infinity) {
      id = this.canonical(id);
      let p = this.cache.get(id);
      if (!p) p = this.create(id);
      else this.cache.delete(id);
      let work = 0;
      while (p.tick < this.history.length && p.alive && work < maxSteps) {
        this.transition(p, this.history[p.tick], p.tick); p.tick++; work++;
      }
      if (!p.alive) p.tick = this.history.length;
      this.cache.set(id, p);
      if (this.cache.size > this.cacheLimit) this.cache.delete(this.cache.keys().next().value);
      const ready = p.tick === this.history.length;
      const jitterX = random(p.key, 0, 11) - 0.5, jitterY = random(p.key, 0, 12) - 0.5;
      return { id, ready, alive: p.alive, exit: p.exit, work, tick: p.tick,
        x: (p.x + 0.5 + jitterX * 0.75) * this.dx - this.extent,
        y: (p.y + 0.5 + jitterY * 0.75) * this.dx - this.extent,
        vx: this.velocities[p.vx], vy: this.velocities[p.vy] };
    }
    forget() { this.cache.clear(); }
    forgetIdentity(id) { this.cache.delete(this.canonical(id)); }
    candidates(count) {
      // Stable prefix independent of count; generate enough bits even for populations > 2^53.
      const result = [], seen = new Set(), words = Math.ceil(this.population.toString(2).length / 32);
      for (let j = 0; result.length < count && BigInt(result.length) < this.population; j++) {
        const k = keys(this.seed + ":candidate:" + j);
        let value = 0n;
        for (let w = 0; w < words; w++) value = value << 32n | BigInt(Math.floor(random(k, w, 20) * 4294967296));
        // Rejection removes modulo bias for a uniform candidate address.
        const range = 1n << BigInt(words * 32), limit = range - range % this.population;
        if (value >= limit) continue;
        const id = (value % this.population).toString();
        if (!seen.has(id)) { result.push(id); seen.add(id); }
      }
      return result;
    }
  }
  root.LazyMatterIdentities = LazyMatterIdentities;
  if (typeof module !== "undefined") module.exports = LazyMatterIdentities;
})(typeof window !== "undefined" ? window : globalThis);
