/* Experimental macro-particle gravity and Gaussian phase-space transport.
 * No kinetic grid is allocated. The small grid below samples force derivatives;
 * it is not a density/velocity transport solver.
 */
(function (root) {
  "use strict";
  const MW = typeof module !== "undefined" ? require("../../milky-way") : root;

  class AnchorEngine {
    static async create(options = {}) {
      const engine = new AnchorEngine(options);
      const adapter = await navigator.gpu.requestAdapter({
        powerPreference: "high-performance",
      });
      if (!adapter) throw Error("WebGPU indisponible pour le POC.");
      const device = await adapter.requestDevice();
      engine.device = device;
      engine.ownsDevice = true;
      try {
        engine.force = await root.createOrbitGPUForceField({
          device,
          capacity: engine.anchors.length + engine.fieldSide ** 2,
        });
        await engine.initialize();
        return engine;
      } catch (error) {
        engine.destroy();
        throw error;
      }
    }

    constructor({
      count = 10000,
      seed = 417,
      structure = 0.85,
      softening = 0.03,
      fieldSide = 32,
      maxStep = 0.1,
      force = null,
    } = {}) {
      if (!Number.isInteger(count) || count < 16 || count > 20000)
        throw Error("Entre 16 et 20 000 ancres.");
      if (!Number.isInteger(fieldSide) || fieldSide < 4 || fieldSide > 128)
        throw Error("Échantillonnage du champ invalide.");
      if (
        !(maxStep > 0) ||
        !Number.isFinite(maxStep) ||
        !(softening > 0) ||
        !Number.isFinite(softening)
      )
        throw Error("Pas de temps ou adoucissement invalide.");
      this.seed = seed;
      this.count = count;
      this.extent = 2;
      this.time = 0;
      this.fieldSide = fieldSide;
      this.maxStep = maxStep;
      this.softening = softening;
      this.force = force;
      const initial = MW.createMilkyWay(count, seed, { structure, softening });
      this.physics = initial.physics;
      this.calibration = initial.calibration;
      this.anchors = initial.map((p, index) => {
        const r = Math.hypot(p.x, p.y),
          ca = p.x / Math.max(r, 1e-12),
          sa = p.y / Math.max(r, 1e-12);
        const flow = MW.stateMW(p.x, p.y, this.physics),
          rho = flow.disc + flow.bulge;
        // Roughly 32 anchor centres within 3 sigma, with bounded resolution.
        const h = Math.max(
          0.008,
          Math.min(
            0.12,
            Math.sqrt(32 / (9 * Math.PI * count * Math.max(rho, 1e-5))),
          ),
        );
        const sv = 5 / MW.MilkyWay.velocityKmS;
        // Reserve part of the initial dispersion for the inferred population.
        const sr = (p.disc ? 30 : 90) / MW.MilkyWay.velocityKmS,
          st = (p.disc ? 20 : 90) / MW.MilkyWay.velocityKmS;
        const mr = p.disc ? flow.vr : 0,
          mt = p.disc ? flow.vt : 0;
        const vr =
          mr + (ca * p.vx + sa * p.vy - mr) * Math.sqrt(1 - (sv / sr) ** 2);
        const vt =
          mt + (-sa * p.vx + ca * p.vy - mt) * Math.sqrt(1 - (sv / st) ** 2);
        return {
          index,
          x: p.x,
          y: p.y,
          vx: ca * vr - sa * vt,
          vy: sa * vr + ca * vt,
          mass: 1 / count,
          ax: 0,
          ay: 0,
          A: new Float64Array([h, 0, 0, 0, 0, h, 0, 0]),
          B: new Float64Array([0, 0, sv, 0, 0, 0, 0, sv]),
          J: new Float64Array(4),
          h,
        };
      });
      this.initialized = false;
      this.lastMs = 0;
      this.steps = 0;
      this.destroyed = false;
    }

    async initialize() {
      await this.updateForces();
      this.initialized = true;
      this.initialAngular = this.angularMomentum();
    }

    async updateForces() {
      const side = this.fieldSide,
        span =
          Math.max(
            this.extent,
            ...this.anchors.map((p) => Math.max(Math.abs(p.x), Math.abs(p.y))),
          ) * 1.1;
      const dx = (2 * span) / (side - 1),
        sources = this.anchors.map((p) => ({
          x: p.x / span,
          y: p.y / span,
          mass: p.mass,
        }));
      for (let y = 0; y < side; y++)
        for (let x = 0; x < side; x++)
          sources.push({
            x: -1 + (2 * x) / (side - 1),
            y: -1 + (2 * y) / (side - 1),
            mass: 0,
          });
      const values = await this.force.compute(sources, {
        gravity: this.physics.gravity / span ** 3,
        softening: this.softening / span,
      });
      const acceleration = new Float64Array(values.length);
      for (let i = 0; i < sources.length; i++) {
        const p = sources[i],
          x = p.x * span,
          y = p.y * span,
          halo = this.physics.haloV2 / (x * x + y * y + 0.04);
        acceleration[2 * i] = values[2 * i] * span - x * halo;
        acceleration[2 * i + 1] = values[2 * i + 1] * span - y * halo;
      }
      const at = (x, y, k) => acceleration[2 * (this.count + y * side + x) + k];
      const gradient = new Float64Array(side * side * 4);
      for (let y = 0; y < side; y++)
        for (let x = 0; x < side; x++) {
          const l = Math.max(0, x - 1),
            r = Math.min(side - 1, x + 1),
            b = Math.max(0, y - 1),
            t = Math.min(side - 1, y + 1),
            o = 4 * (y * side + x);
          gradient[o] = (at(r, y, 0) - at(l, y, 0)) / ((r - l) * dx);
          gradient[o + 3] = (at(x, t, 1) - at(x, b, 1)) / ((t - b) * dx);
          // Symmetric Hessian of the conservative field, at this coarse resolution.
          gradient[o + 1] = gradient[o + 2] =
            0.5 *
            ((at(x, t, 0) - at(x, b, 0)) / ((t - b) * dx) +
              (at(r, y, 1) - at(l, y, 1)) / ((r - l) * dx));
        }
      for (let i = 0; i < this.count; i++) {
        const p = this.anchors[i],
          u = Math.max(0, Math.min(side - 1.000001, (p.x + span) / dx)),
          v = Math.max(0, Math.min(side - 1.000001, (p.y + span) / dx)),
          x = Math.floor(u),
          y = Math.floor(v),
          tx = u - x,
          ty = v - y;
        p.ax = acceleration[2 * i];
        p.ay = acceleration[2 * i + 1];
        for (let k = 0; k < 4; k++)
          p.J[k] =
            gradient[4 * (y * side + x) + k] * (1 - tx) * (1 - ty) +
            gradient[4 * (y * side + x + 1) + k] * tx * (1 - ty) +
            gradient[4 * ((y + 1) * side + x) + k] * (1 - tx) * ty +
            gradient[4 * ((y + 1) * side + x + 1) + k] * tx * ty;
      }
    }

    async advance(duration) {
      if (this.destroyed) throw Error("POC détruit.");
      if (!Number.isFinite(duration) || duration < 0)
        throw Error("Durée invalide.");
      if (!this.initialized) await this.initialize();
      const started = performance.now();
      let remaining = duration,
        steps = 0;
      // A bounded batch lets the UI remain responsive at extreme requested rates.
      while (remaining > 1e-10 && steps < 8) {
        const dt = Math.min(remaining, this.maxStep);
        for (const p of this.anchors) {
          p.vx += 0.5 * p.ax * dt;
          p.vy += 0.5 * p.ay * dt;
          p.x += p.vx * dt;
          p.y += p.vy * dt;
          for (let k = 0; k < 4; k++) {
            const x = p.A[k],
              y = p.A[k + 4];
            p.B[k] += 0.5 * dt * (p.J[0] * x + p.J[1] * y);
            p.B[k + 4] += 0.5 * dt * (p.J[2] * x + p.J[3] * y);
            p.A[k] += dt * p.B[k];
            p.A[k + 4] += dt * p.B[k + 4];
          }
        }
        await this.updateForces();
        for (const p of this.anchors) {
          p.vx += 0.5 * p.ax * dt;
          p.vy += 0.5 * p.ay * dt;
          for (let k = 0; k < 4; k++) {
            const x = p.A[k],
              y = p.A[k + 4];
            p.B[k] += 0.5 * dt * (p.J[0] * x + p.J[1] * y);
            p.B[k + 4] += 0.5 * dt * (p.J[2] * x + p.J[3] * y);
          }
        }
        this.time += dt;
        remaining -= dt;
        steps++;
      }
      this.steps += steps;
      this.lastMs = performance.now() - started;
      return { advanced: duration - remaining, steps };
    }

    angularMomentum() {
      return this.anchors.reduce(
        (s, p) => s + p.mass * (p.x * p.vy - p.y * p.vx),
        0,
      );
    }
    diagnostics() {
      return {
        time: this.time,
        mass: 1,
        angularMomentum: this.angularMomentum(),
        count: this.count,
        lastMs: this.lastMs,
        steps: this.steps,
      };
    }
    destroy() {
      if (this.destroyed) return;
      this.destroyed = true;
      this.force?.destroy();
      if (this.ownsDevice) this.device?.destroy();
    }
  }
  root.AnchorEngine = AnchorEngine;
  if (typeof module !== "undefined") module.exports = AnchorEngine;
})(typeof window !== "undefined" ? window : globalThis);
