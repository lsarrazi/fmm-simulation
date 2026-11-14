/* Collisionless kinetic prototype, planar positions with a softened 3D gravity kernel.
 * f stores cell-integrated mass in (x, y, vx, vy), not individual bodies.
 * Spatial and velocity boundaries are absorbing. Lost mass is measured, never renormalized.
 * First-order conservative upwind transport, dimension splitting and adaptive CFL step.
 */
(function (root) {
  "use strict";
  class KineticMatter {
    constructor(options = {}) {
      this.n = options.n ?? 24;
      this.nv = options.nv ?? 21; // Odd: a cold distribution occupies a bin centred at v = 0.
      this.extent = options.extent ?? 1.5;
      this.vmax = options.vmax ?? 2;
      this.dx = 2 * this.extent / this.n;
      this.dv = 2 * this.vmax / this.nv;
      this.size = this.n * this.n * this.nv * this.nv;
      this.f = new Float32Array(this.size);
      this.scratch = new Float32Array(this.size);
      this.mass = new Float64Array(this.n * this.n);
      this.ax = new Float64Array(this.mass.length);
      this.ay = new Float64Array(this.mass.length);
      this.velocities = Float64Array.from({ length: this.nv }, (_, k) => (k - (this.nv - 1) / 2) * this.dv);
      // Species metadata anticipates a later multi-species kinetic model; no EM force here.
      this.species = Object.freeze({ name: "matière neutre", particleMass: 1, charge: 0 });
      this.time = 0;
      this.spatialLoss = 0;
      this.velocityLoss = 0;
      this.reset(options);
    }
    reset({ sigma = 0.12, fluctuation = 0.15, radius = 0.85 } = {}) {
      this.time = this.spatialLoss = this.velocityLoss = 0;
      const weights = this.velocities.map(v => sigma === 0 ? (v === 0 ? 1 : 0) : Math.exp(-0.5 * (v / sigma) ** 2));
      const normalization = weights.reduce((a, b) => a + b, 0);
      for (let k = 0; k < weights.length; k++) weights[k] /= normalization;
      let total = 0;
      for (let y = 0; y < this.n; y++) {
        const py = (y + 0.5) * this.dx - this.extent;
        for (let x = 0; x < this.n; x++) {
          const px = (x + 0.5) * this.dx - this.extent;
          const r = Math.hypot(px, py);
          const edge = Math.max(0, Math.min(1, (radius - r) / (1.5 * this.dx)));
          // Deterministic bounded fluctuations, no galaxy geometry or prescribed orbital motion.
          const noise = (Math.cos(11 * px + 7 * py) + Math.sin(17 * px - 9 * py) + Math.cos(5 * px - 13 * py)) / 3;
          const density = edge * (1 + fluctuation * noise);
          for (let vy = 0; vy < this.nv; vy++) for (let vx = 0; vx < this.nv; vx++) {
            const i = ((y * this.n + x) * this.nv + vy) * this.nv + vx;
            this.f[i] = density * weights[vx] * weights[vy];
            total += this.f[i];
          }
        }
      }
      for (let i = 0; i < this.size; i++) this.f[i] /= total;
      this.project();
    }
    project() {
      const stride = this.nv * this.nv;
      for (let cell = 0; cell < this.mass.length; cell++) {
        let sum = 0;
        for (let i = cell * stride; i < (cell + 1) * stride; i++) sum += this.f[i];
        this.mass[cell] = sum;
      }
    }
    field({ gravity = 0.48, softening = 0.12, exponent = 2 } = {}) {
      this.ax.fill(0); this.ay.fill(0);
      if (gravity === 0) return 0;
      // Equal and opposite pair forces on spatial cells. This is a quadrature of the
      // same softened kernel as Orbit, not a two-dimensional Poisson approximation.
      let maxAcceleration = 0;
      for (let a = 0; a < this.mass.length; a++) {
        const x = a % this.n, y = Math.floor(a / this.n);
        for (let b = a + 1; b < this.mass.length; b++) {
          const rx = (b % this.n - x) * this.dx;
          const ry = (Math.floor(b / this.n) - y) * this.dx;
          const q = rx * rx + ry * ry + softening * softening;
          const k = gravity * Math.pow(q, -(exponent + 1) / 2);
          this.ax[a] += this.mass[b] * k * rx;
          this.ay[a] += this.mass[b] * k * ry;
          this.ax[b] -= this.mass[a] * k * rx;
          this.ay[b] -= this.mass[a] * k * ry;
        }
      }
      for (let i = 0; i < this.mass.length; i++) maxAcceleration = Math.max(maxAcceleration, Math.abs(this.ax[i]), Math.abs(this.ay[i]));
      return maxAcceleration;
    }
    transport(axis, dt) {
      const n = this.n, nv = this.nv, input = this.f, output = this.scratch;
      const stride = [nv * nv, n * nv * nv, 1, nv][axis];
      const spacing = axis < 2 ? this.dx : this.dv;
      let lost = 0;
      for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
        const cell = y * n + x;
        for (let vy = 0; vy < nv; vy++) for (let vx = 0; vx < nv; vx++) {
          const i = (cell * nv + vy) * nv + vx;
          const coordinate = axis === 0 ? x : axis === 1 ? y : axis === 2 ? vx : vy;
          const count = axis < 2 ? n : nv;
          const speed = axis === 0 ? this.velocities[vx] : axis === 1 ? this.velocities[vy] : axis === 2 ? this.ax[cell] : this.ay[cell];
          const c = Math.abs(speed) * dt / spacing;
          const upstream = speed >= 0 ? coordinate > 0 : coordinate < count - 1;
          output[i] = (1 - c) * input[i] + c * (upstream ? input[i + (speed >= 0 ? -stride : stride)] : 0);
          if ((speed > 0 && coordinate === count - 1) || (speed < 0 && coordinate === 0)) lost += c * input[i];
        }
      }
      if (axis < 2) this.spatialLoss += lost; else this.velocityLoss += lost;
      this.f = output; this.scratch = input;
    }
    advance(duration, parameters = {}) {
      let remaining = duration, steps = 0;
      while (remaining > 1e-10 && steps < 12) {
        this.project();
        const acceleration = this.field(parameters);
        const largestVelocity = Math.abs(this.velocities[0]);
        const dt = Math.min(remaining, 0.8 * this.dx / largestVelocity, acceleration > 0 ? 0.8 * this.dv / acceleration : Infinity);
        // Observers record the exact transport inputs; false stops BEFORE changing f.
        if (this.onTransport && this.onTransport({ dt, ax: this.ax, ay: this.ay }) === false) break;
        this.transport(2, dt); this.transport(3, dt);
        this.transport(0, dt); this.transport(1, dt);
        this.time += dt; remaining -= dt; steps++;
      }
      this.project();
      return { advanced: duration - remaining, steps };
    }
    diagnostics() {
      let mass = 0, momentumX = 0, momentumY = 0, angularMomentum = 0, kinetic = 0;
      for (let y = 0; y < this.n; y++) for (let x = 0; x < this.n; x++) {
        const px = (x + 0.5) * this.dx - this.extent, py = (y + 0.5) * this.dx - this.extent;
        for (let vy = 0; vy < this.nv; vy++) for (let vx = 0; vx < this.nv; vx++) {
          const m = this.f[((y * this.n + x) * this.nv + vy) * this.nv + vx];
          const u = this.velocities[vx], v = this.velocities[vy];
          mass += m; momentumX += m * u; momentumY += m * v;
          angularMomentum += m * (px * v - py * u); kinetic += 0.5 * m * (u * u + v * v);
        }
      }
      return { mass, momentumX, momentumY, angularMomentum, kinetic, spatialLoss: this.spatialLoss, velocityLoss: this.velocityLoss, time: this.time };
    }
  }
  root.KineticMatter = KineticMatter;
  if (typeof module !== "undefined") module.exports = KineticMatter;
})(typeof window !== "undefined" ? window : globalThis);
