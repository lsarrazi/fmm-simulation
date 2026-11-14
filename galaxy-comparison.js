/* Controlled comparison: shared FMM N-body vs Eulerian Vlasov quadrature.
 * Both start from the SAME seeded weighted phase-space sample. No damping or expansion.
 * The TVD solver intentionally does not use the old first-order identity replay.
 */
(function(root) {
  'use strict';
  const Base = typeof module !== 'undefined' ? require('./continuum-solver.js') : root.KineticMatter;
  const createFMM = typeof module !== 'undefined' ? require('./fmm-core.js') : root.createOrbitFMMField;
  const createMW=typeof module!=='undefined'?require('./milky-way.js').createMilkyWay:root.createMilkyWay;
  function initialBodies(count = 2048, scenario = 'galaxy', seed = 417,options={}) {
    if(scenario==='milkyway')return createMW(count,seed,options);
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
    constructor(options) {
      super(options); this.flux = new Float32Array(this.size); this.fmm = createFMM();
      this.fieldPoints = Array.from({length:this.mass.length},(_,i)=>({x:(i%this.n+.5)*this.dx-this.extent,y:(Math.floor(i/this.n)+.5)*this.dx-this.extent,mass:0}));
    }
    field(parameters={}) {
      const start=performance.now();
      for(let i=0;i<this.mass.length;i++)this.fieldPoints[i].mass=this.mass[i];
      this.fmm.compute(this.fieldPoints,parameters);
      let maximum=0;
      for(let i=0;i<this.mass.length;i++){
        const p=this.fieldPoints[i],halo=(parameters.haloV2||0)/(p.x*p.x+p.y*p.y+.04);
        this.ax[i]=p.ax-p.x*halo;this.ay[i]=p.ay-p.y*halo;
        maximum=Math.max(maximum,Math.abs(this.ax[i]),Math.abs(this.ay[i]));
      }
      this.fieldMs=(this.fieldMs||0)+performance.now()-start;
      return maximum;
    }
    advance(duration,parameters={}){this.fieldMs=0;this.transportMs=0;return super.advance(duration,parameters);}
    load(bodies) {
      if(bodies.physics?.analytic){const populate=typeof module!=='undefined'?require('./milky-way.js').populateMilkyWay:root.populateMilkyWay;populate(this,bodies.physics);return;}
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
      const startTime=performance.now();
      const n = this.n, nv = this.nv, input = this.f, output = this.scratch, flux = this.flux;
      const stride = [nv * nv, n * nv * nv, 1, nv][axis], count = axis < 2 ? n : nv;
      let lost = 0;
      // One constant speed per transport line: hoist all index divisions and CFL
      // work out of the inner loop (millions of cells at interactive resolution).
      for(let block=0;block<input.length;block+=stride*count)for(let lane=0;lane<stride;lane++){
        const start=block+lane,cell=Math.floor(start/(nv*nv));
        const speed=axis===0?this.velocities[start%nv]:axis===1?this.velocities[Math.floor(start/nv)%nv]:axis===2?this.ax[cell]:this.ay[cell];
        const c = Math.abs(speed) * dt / (axis < 2 ? this.dx : this.dv);
        const direction=speed>=0?1:-1,correction=.5*(1-c)*direction;
        if(c===0){for(let k=0,i=start;k<count;k++,i+=stride)output[i]=input[i];continue;}
        for(let coordinate=0,i=start;coordinate<count;coordinate++,i+=stride){
        if(input[i]===0){flux[i]=0;continue;}
        const left = input[i] - (coordinate > 0 ? input[i - stride] : 0);
        const right = (coordinate < count - 1 ? input[i + stride] : 0) - input[i];
        const slope = left * right > 0 ? Math.sign(left) * Math.min(2 * Math.abs(left), 2 * Math.abs(right), Math.abs(left + right) / 2) : 0;
        flux[i] = c * (input[i] + correction * slope);
        if (direction > 0 ? coordinate === count - 1 : coordinate === 0) lost += flux[i];
        }
        for(let coordinate=0,i=start;coordinate<count;coordinate++,i+=stride){
        const upstream = speed >= 0 ? coordinate > 0 : coordinate < count - 1;
        output[i] = input[i] - flux[i] + (upstream ? flux[i + (speed >= 0 ? -stride : stride)] : 0);
        }
      }
      if (axis < 2) this.spatialLoss += lost; else this.velocityLoss += lost;
      this.f = output; this.scratch = input;
      this.transportMs=(this.transportMs||0)+performance.now()-startTime;
    }
  }
  class DirectBodies {
    constructor(bodies) { this.bodies = bodies.map(b => ({ ...b })); this.time = 0; }
    advance(dt, { gravity = .48, softening = .042,haloV2=0 } = {}) {
      const b = this.bodies, ax = new Float64Array(b.length), ay = new Float64Array(b.length);
      for (let i = 0; i < b.length; i++) for (let j = i + 1; j < b.length; j++) {
        const x = b[j].x - b[i].x, y = b[j].y - b[i].y, q = x*x + y*y + softening*softening;
        const k = gravity / (q * Math.sqrt(q));
        ax[i] += b[j].mass * x * k; ay[i] += b[j].mass * y * k;
        ax[j] -= b[i].mass * x * k; ay[j] -= b[i].mass * y * k;
      }
      for (let i = 0; i < b.length; i++) { const p = b[i],halo=haloV2/(p.x*p.x+p.y*p.y+.04); p.vx += (ax[i]-p.x*halo) * dt; p.vy += (ay[i]-p.y*halo) * dt; p.x += p.vx * dt; p.y += p.vy * dt; }
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
  class FMMBodies extends DirectBodies {
    constructor(bodies){super(bodies);this.fmm=createFMM();}
    advance(dt,parameters={}){
      this.fmm.compute(this.bodies,parameters);
      for(const p of this.bodies){const halo=(parameters.haloV2||0)/(p.x*p.x+p.y*p.y+.04);p.vx+=(p.ax-p.x*halo)*dt;p.vy+=(p.ay-p.y*halo)*dt;p.x+=p.vx*dt;p.y+=p.vy*dt;}
      this.time+=dt;
    }
  }
  const api = { initialBodies, TVDMatter, DirectBodies, FMMBodies };
  Object.assign(root, api); if (typeof module !== 'undefined') module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);
