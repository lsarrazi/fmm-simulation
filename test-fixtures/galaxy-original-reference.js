/* Controlled comparison: direct softened N-body vs Eulerian Vlasov quadrature.
 * Both start from the SAME seeded weighted phase-space sample. No damping or expansion.
 * The TVD solver intentionally does not use the old first-order identity replay.
 */
(function(root) {
  'use strict';
  const Base = typeof module !== 'undefined' ? require('./continuum-solver.js') : root.KineticMatter;
  function initialBodies(count = 2048, scenario = 'galaxy', seed = 417) {
    let state = seed >>> 0;
    const random = () => { state = (Math.imul(1664525, state) + 1013904223) >>> 0; return (state + .5) / 4294967296; };
    const normal = () => Math.sqrt(-2 * Math.log(random())) * Math.cos(2 * Math.PI * random());
    const bodies = []; let total = 0;
    for (let i = 0; i < count; i++) {
      let x, y, vx, vy;
      if (scenario === 'galaxy') {
        const core = i < Math.max(8, count * .025);
        const r = core ? random() * .06 : Math.pow(random(), .58) * .83 + .035;
        const a = i % 3 * Math.PI * 2 / 3 + r * 7.8 + normal() * (.13 + r * .13);
        x = Math.cos(a) * r + normal() * .012; y = Math.sin(a) * r * .73 + normal() * .012;
        const v = Math.sqrt(.48 * (.16 + Math.min(1, r * 1.35)) / Math.max(r, .035)) * .49;
        vx = -Math.sin(a) * v + normal() * .018; vy = Math.cos(a) * v * .73 + normal() * .018;
      } else {
        const r = Math.sqrt(random()) * .85, a = random() * Math.PI * 2;
        x = r * Math.cos(a); y = r * Math.sin(a);
        vx = scenario === 'cold' ? 0 : normal() * .15; vy = scenario === 'cold' ? 0 : normal() * .15;
      }
      // Exact inverse of the existing default m^-2.35 distribution on [.2, 5].
      const m = Math.pow(Math.pow(.2, -1.35) + random() * (Math.pow(5, -1.35) - Math.pow(.2, -1.35)), -1 / 1.35);
      total += m; bodies.push({ x, y, vx, vy, mass: m });
    }
    for (const b of bodies) b.mass /= total;
    return bodies;
  }
  function weights(value, origin, spacing, count) {
    const u = (value - origin) / spacing - .5;
    const lo = Math.floor(u), t = u - lo;
    if (lo < 0 || lo + 1 >= count) throw new Error('Domaine trop petit pour les conditions initiales.');
    return [[lo, 1 - t], [lo + 1, t]];
  }
  class TVDMatter extends Base {
    constructor(options) { super(options); this.flux = new Float32Array(this.size); }
    load(bodies) {
      this.f.fill(0); this.time = this.spatialLoss = this.velocityLoss = 0;
      for (const b of bodies) {
        const wx = weights(b.x, -this.extent, this.dx, this.n), wy = weights(b.y, -this.extent, this.dx, this.n);
        const wu = weights(b.vx, -this.vmax, this.dv, this.nv), wv = weights(b.vy, -this.vmax, this.dv, this.nv);
        for (const [x, a] of wx) for (const [y, c] of wy) for (const [u, d] of wu) for (const [v, e] of wv)
          this.f[((y * this.n + x) * this.nv + v) * this.nv + u] += b.mass * a * c * d * e;
      }
      this.project();
    }
    transport(axis, dt) {
      const n = this.n, nv = this.nv, input = this.f, output = this.scratch, flux = this.flux;
      const stride = [nv * nv, n * nv * nv, 1, nv][axis], count = axis < 2 ? n : nv;
      let lost = 0;
      for (let i = 0; i < input.length; i++) {
        const coordinate = Math.floor(i / stride) % count;
        const cell = Math.floor(i / (nv * nv));
        const speed = axis === 0 ? this.velocities[i % nv] : axis === 1 ? this.velocities[Math.floor(i / nv) % nv] : axis === 2 ? this.ax[cell] : this.ay[cell];
        const c = Math.abs(speed) * dt / (axis < 2 ? this.dx : this.dv);
        const left = input[i] - (coordinate > 0 ? input[i - stride] : 0);
        const right = (coordinate < count - 1 ? input[i + stride] : 0) - input[i];
        const slope = left * right > 0 ? Math.sign(left) * Math.min(2 * Math.abs(left), 2 * Math.abs(right), Math.abs(left + right) / 2) : 0;
        const direction = speed >= 0 ? 1 : -1;
        flux[i] = c * (input[i] + .5 * (1 - c) * direction * slope);
        if (direction > 0 ? coordinate === count - 1 : coordinate === 0) lost += flux[i];
      }
      for (let i = 0; i < input.length; i++) {
        const coordinate = Math.floor(i / stride) % count, cell = Math.floor(i / (nv * nv));
        const speed = axis === 0 ? this.velocities[i % nv] : axis === 1 ? this.velocities[Math.floor(i / nv) % nv] : axis === 2 ? this.ax[cell] : this.ay[cell];
        const upstream = speed >= 0 ? coordinate > 0 : coordinate < count - 1;
        output[i] = input[i] - flux[i] + (upstream ? flux[i + (speed >= 0 ? -stride : stride)] : 0);
      }
      if (axis < 2) this.spatialLoss += lost; else this.velocityLoss += lost;
      this.f = output; this.scratch = input;
    }
  }
  class DirectBodies {
    constructor(bodies) { this.bodies = bodies.map(b => ({ ...b })); this.time = 0; }
    advance(dt, { gravity = .48, softening = .042 } = {}) {
      const b = this.bodies, ax = new Float64Array(b.length), ay = new Float64Array(b.length);
      for (let i = 0; i < b.length; i++) for (let j = i + 1; j < b.length; j++) {
        const x = b[j].x - b[i].x, y = b[j].y - b[i].y, q = x*x + y*y + softening*softening;
        const k = gravity / (q * Math.sqrt(q));
        ax[i] += b[j].mass * x * k; ay[i] += b[j].mass * y * k;
        ax[j] -= b[i].mass * x * k; ay[j] -= b[i].mass * y * k;
      }
      for (let i = 0; i < b.length; i++) { const p = b[i]; p.vx += ax[i] * dt; p.vy += ay[i] * dt; p.x += p.vx * dt; p.y += p.vy * dt; }
      this.time += dt;
    }
    project(grid) {
      const mass = new Float64Array(grid.n * grid.n);
      for (const b of this.bodies) {
        const ux = (b.x + grid.extent) / grid.dx - .5, uy = (b.y + grid.extent) / grid.dx - .5;
        const ix = Math.floor(ux), iy = Math.floor(uy), tx = ux - ix, ty = uy - iy;
        for (let y = 0; y < 2; y++) for (let x = 0; x < 2; x++) {
          const px = ix + x, py = iy + y;
          if (px >= 0 && py >= 0 && px < grid.n && py < grid.n) mass[py * grid.n + px] += b.mass * (x ? tx : 1 - tx) * (y ? ty : 1 - ty);
        }
      }
      return mass;
    }
  }
  const api = { initialBodies, TVDMatter, DirectBodies };
  Object.assign(root, api); if (typeof module !== 'undefined') module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);
