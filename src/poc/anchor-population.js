/* Immutable material addresses in overlapping 4D Gaussian clouds.
 * Spatial queries prune rank intervals; no array of the theoretical stars exists.
 * Density is analytic. Discovery has a separately reported work/error budget.
 */
(function (root) {
  "use strict";
  const P =
    typeof module !== "undefined"
      ? require("../../stellar-population")
      : root.StellarPopulation;
  const mix = (x) => {
    x = Math.imul(x ^ (x >>> 16), 0x7feb352d);
    x = Math.imul(x ^ (x >>> 15), 0x846ca68b);
    return (x ^ (x >>> 16)) >>> 0;
  };
  const a = [
      -39.6968302866538, 220.946098424521, -275.928510446969, 138.357751867269,
      -30.6647980661472, 2.50662827745924,
    ],
    b = [
      -54.4760987982241, 161.585836858041, -155.698979859887, 66.8013118877197,
      -13.2806815528857,
    ],
    c = [
      -0.00778489400243029, -0.322396458041136, -2.40075827716184,
      -2.54973253934373, 4.37466414146497, 2.93816398269878,
    ],
    d = [
      0.00778469570904146, 0.32246712907004, 2.445134137143, 3.75440866190742,
    ];
  function normalQuantile(p) {
    p = Math.max(6.661338147750939e-16, Math.min(1 - 6.661338147750939e-16, p));

    if (p < 0.02425 || p > 1 - 0.02425) {
      const q = Math.sqrt(-2 * Math.log(p < 0.5 ? p : 1 - p)),
        v =
          (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q +
            c[5]) /
          ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
      return p < 0.5 ? v : -v;
    }
    const q = p - 0.5,
      r = q * q;
    return (
      ((((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) *
        q) /
      (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)
    );
  }
  function normalCDF(x) {
    const t = 1 / (1 + 0.2316419 * Math.abs(x)),
      p =
        (Math.exp((-x * x) / 2) / Math.sqrt(2 * Math.PI)) *
        t *
        (0.31938153 +
          t *
            (-0.356563782 +
              t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
    return x >= 0 ? 1 - p : p;
  }
  function covariance(p) {
    let xx = 0,
      xy = 0,
      yy = 0;
    for (let k = 0; k < 4; k++) {
      xx += p.A[k] ** 2;
      xy += p.A[k] * p.A[k + 4];
      yy += p.A[k + 4] ** 2;
    }
    return { xx, xy, yy, det: Math.max(xx * yy - xy * xy, 1e-30) };
  }
  const density = (p, q, x, y) =>
    Math.exp(
      (-0.5 *
        (q.yy * (x - p.x) ** 2 -
          2 * q.xy * (x - p.x) * (y - p.y) +
          q.xx * (y - p.y) ** 2)) /
        q.det,
    ) /
    (2 * Math.PI * Math.sqrt(q.det));
  function rectangleProbability(p, q, view) {
    const rx = 8 * Math.sqrt(q.xx),
      ry = 8 * Math.sqrt(q.yy);
    // Tail probability beyond this bound is below floating-point observation
    // precision. Most distant/full-view clouds need no quadrature at all.
    if (
      p.x + rx < view.l ||
      p.x - rx > view.r ||
      p.y + ry < view.b ||
      p.y - ry > view.t
    )
      return 0;
    if (
      p.x - rx >= view.l &&
      p.x + rx <= view.r &&
      p.y - ry >= view.b &&
      p.y + ry <= view.t
    )
      return 1;
    const sx = Math.sqrt(q.xx),
      l = Math.max(-8, (view.l - p.x) / sx),
      r = Math.min(8, (view.r - p.x) / sx);
    if (l >= r) return 0;
    const sy = Math.sqrt(q.det / q.xx),
      slope = q.xy / sx;
    // Simpson integration; the conditional Y distribution handles correlation.
    // Integration is approximate and is used for allocation, never for identity.
    const n = 16,
      dx = (r - l) / n;
    let sum = 0;
    for (let k = 0; k <= n; k++) {
      const z = l + k * dx,
        mean = p.y + slope * z,
        f =
          (Math.exp((-z * z) / 2) / Math.sqrt(2 * Math.PI)) *
          (normalCDF((view.t - mean) / sy) - normalCDF((view.b - mean) / sy));
      sum += (k === 0 || k === n ? 1 : k % 2 ? 4 : 2) * Math.max(0, f);
    }
    return Math.max(0, Math.min(1, (sum * dx) / 3));
  }

  class AnchorPopulation {
    constructor(
      engine,
      { population = 100000000000n, seed = engine.seed } = {},
    ) {
      this.engine = engine;
      this.population = BigInt(population);
      if (this.population <= 0n)
        throw Error("La population doit être positive.");
      this.universe = `anchors:v1:${seed}:${engine.count}:${this.population}`;
      this.cache = new Map();
      this.cacheLimit = 131072;
      this.cacheCursor = null;
      const n = BigInt(engine.count),
        q = this.population / n,
        r = this.population % n;
      this.counts = engine.anchors.map((_, i) => q + BigInt(BigInt(i) < r));
      this.seeds = engine.anchors.map((_, i) =>
        mix(
          Math.floor(
            P.hash(this.universe + "/" + i.toString(36), 907) * 4294967296,
          ),
        ),
      );
    }
    get time() {
      return this.engine.time;
    }
    address(index, rank) {
      const id = this.universe + "/" + index.toString(36) + "/" + rank;
      let z = this.cache.get(id);
      if (z) return { id, z };
      let start = 0n,
        count = this.counts[index],
        seed = this.seeds[index],
        depth = 0;
      const lo = [0, 0, 0, 0],
        hi = [1, 1, 1, 1];
      if (count <= 9007199254740991n) {
        let remaining = Number(count),
          local = Number(rank);
        while (remaining > 1 && depth < 192) {
          const first = Math.ceil(remaining / 2),
            child = local < first ? 0 : 1,
            axis = depth % 4,
            spatial = child ^ (seed & 1),
            mid = (lo[axis] + hi[axis]) / 2;
          if (spatial) lo[axis] = mid;
          else hi[axis] = mid;
          if (child) {
            local -= first;
            remaining -= first;
          } else remaining = first;
          seed = mix(seed ^ Math.imul(child + 1, 0x9e3779b9));
          depth++;
        }
      } else {
        while (count > 1n && depth < 192) {
          const first = (count + 1n) / 2n,
            child = rank < start + first ? 0 : 1,
            axis = depth % 4,
            spatial = child ^ (seed & 1),
            mid = (lo[axis] + hi[axis]) / 2;
          if (spatial) lo[axis] = mid;
          else hi[axis] = mid;
          if (child) {
            start += first;
            count -= first;
          } else count = first;
          seed = mix(seed ^ Math.imul(child + 1, 0x9e3779b9));
          depth++;
        }
      }
      z = Float64Array.from(
        lo.map((v, k) =>
          normalQuantile(v + (hi[k] - v) * P.hash(id, 101 + k * 19)),
        ),
      );
      if (this.cache.size >= this.cacheLimit) {
        this.cacheCursor ??= this.cache.keys();
        let key = this.cacheCursor.next();
        if (key.done) {
          this.cacheCursor = this.cache.keys();
          key = this.cacheCursor.next();
        }
        this.cache.delete(key.value);
      }
      this.cache.set(id, z);
      return { id, z };
    }
    evaluate(index, rank) {
      const { id, z } = this.address(index, rank),
        p = this.engine.anchors[index];
      let x = p.x,
        y = p.y,
        vx = p.vx,
        vy = p.vy;
      for (let k = 0; k < 4; k++) {
        x += p.A[k] * z[k];
        y += p.A[k + 4] * z[k];
        vx += p.B[k] * z[k];
        vy += p.B[k + 4] * z[k];
      }
      return { id, x, y, vx, vy, anchor: index };
    }
    resolve(id) {
      const parts = String(id).split("/");
      if (
        parts.length !== 3 ||
        parts[0] !== this.universe ||
        !/^\d+$/.test(parts[2]) ||
        !/^[0-9a-z]+$/.test(parts[1])
      )
        throw Error("Identité extérieure au POC courant.");
      const index = parseInt(parts[1], 36),
        rank = BigInt(parts[2]);
      if (
        !Number.isSafeInteger(index) ||
        index < 0 ||
        index >= this.counts.length ||
        rank < 0n ||
        rank >= this.counts[index]
      )
        throw Error("Rang hors population.");
      return this.evaluate(index, rank);
    }
    densityAt(x, y) {
      let value = 0;
      for (let i = 0; i < this.counts.length; i++)
        value +=
          (Number(this.counts[i]) / Number(this.population)) *
          density(
            this.engine.anchors[i],
            covariance(this.engine.anchors[i]),
            x,
            y,
          );
      return value;
    }
    gaussians() {
      if (this.gaussianEpoch !== this.engine.time) {
        this.gaussianEpoch = this.engine.time;
        this.gaussianCache = this.engine.anchors.map((anchor, i) => ({
          anchor,
          covariance: covariance(anchor),
          mass: Number(this.counts[i]) / Number(this.population),
        }));
      }
      return this.gaussianCache;
    }
    windowProjection({
      x = 0,
      y = 0,
      width = 4,
      aspectX = 1,
      aspectY = 1,
    } = {}) {
      const view = {
        l: x - (width * aspectX) / 2,
        r: x + (width * aspectX) / 2,
        b: y - (width * aspectY) / 2,
        t: y + (width * aspectY) / 2,
      };
      let mass = 0;
      for (const { anchor, covariance: q, mass: weight } of this.gaussians())
        mass += weight * rectangleProbability(anchor, q, view);
      // Geometry/observation metadata only. This placeholder is never displayed:
      // the direct GPU layer evaluates the actual spatial distribution.
      return {
        image: new Float64Array(1),
        resolution: 1,
        width,
        pixel: width,
        left: x - width / 2,
        bottom: y - width / 2,
        centerX: x,
        centerY: y,
        aspectX,
        aspectY,
        mass,
        viewMass: mass,
        expected: Number(this.population) * mass,
      };
    }
    project({
      x = 0,
      y = 0,
      width = 4,
      resolution = 96,
      aspectX = 1,
      aspectY = 1,
    } = {}) {
      const image = new Float64Array(resolution ** 2),
        pixel = width / resolution,
        left = x - width / 2,
        bottom = y - width / 2;
      for (let i = 0; i < this.counts.length; i++) {
        const p = this.engine.anchors[i],
          q = covariance(p),
          rx = 6 * Math.sqrt(q.xx),
          ry = 6 * Math.sqrt(q.yy);
        const x0 = Math.max(0, Math.floor((p.x - rx - left) / pixel)),
          x1 = Math.min(resolution - 1, Math.ceil((p.x + rx - left) / pixel)),
          y0 = Math.max(0, Math.floor((p.y - ry - bottom) / pixel)),
          y1 = Math.min(resolution - 1, Math.ceil((p.y + ry - bottom) / pixel));
        if (x0 > x1 || y0 > y1) continue;
        const mass = Number(this.counts[i]) / Number(this.population);
        if (Math.max(rx, ry) < pixel) {
          const u = (p.x - left) / pixel - 0.5,
            v = (p.y - bottom) / pixel - 0.5,
            ix = Math.floor(u),
            iy = Math.floor(v),
            tx = u - ix,
            ty = v - iy;
          for (let a = 0; a < 2; a++)
            for (let b = 0; b < 2; b++) {
              const xx = ix + a,
                yy = iy + b;
              if (xx >= 0 && xx < resolution && yy >= 0 && yy < resolution)
                image[yy * resolution + xx] +=
                  mass * (a ? tx : 1 - tx) * (b ? ty : 1 - ty);
            }
          continue;
        }
        for (let yy = y0; yy <= y1; yy++)
          for (let xx = x0; xx <= x1; xx++)
            image[yy * resolution + xx] +=
              mass *
              density(
                p,
                q,
                left + (xx + 0.5) * pixel,
                bottom + (yy + 0.5) * pixel,
              ) *
              pixel ** 2;
      }
      let mass = 0,
        viewMass = 0;
      for (let yy = 0; yy < resolution; yy++)
        for (let xx = 0; xx < resolution; xx++) {
          const value = image[yy * resolution + xx];
          mass += value;
          if (
            Math.abs((xx + 0.5) / resolution - 0.5) <= aspectX / 2 &&
            Math.abs((yy + 0.5) / resolution - 0.5) <= aspectY / 2
          )
            viewMass += value;
        }
      return {
        image,
        resolution,
        width,
        pixel,
        left,
        bottom,
        centerX: x,
        centerY: y,
        aspectX,
        aspectY,
        mass,
        viewMass,
        expected: Number(this.population) * viewMass,
      };
    }

    sample({
      x = 0,
      y = 0,
      width = 4,
      aspectX = 1,
      aspectY = 1,
      budget = 2048,
      rankFilter = null,
      maxVisits = 160000,
    } = {}) {
      if (!Number.isInteger(budget) || budget < 0)
        throw Error("Budget invalide.");
      const view = {
          l: x - (width * aspectX) / 2,
          r: x + (width * aspectX) / 2,
          b: y - (width * aspectY) / 2,
          t: y + (width * aspectY) / 2,
        },
        clouds = [];
      let expected = 0,
        weightSum = 0;
      for (let i = 0; i < this.counts.length; i++) {
        const p = this.engine.anchors[i],
          q = covariance(p),
          prob = rectangleProbability(p, q, view),
          stars = Number(this.counts[i]) * prob;
        expected += stars;
        const weight = stars * (rankFilter?.fraction ?? 1);
        if (weight > 0) {
          weightSum += weight;
          clouds.push({ i, weight, end: weightSum });
        }
      }
      const wanted = Math.min(budget, Math.ceil(weightSum)),
        quotas = new Map();
      let cursor = 0;
      for (let k = 0; k < wanted; k++) {
        const r = ((k + 0.5) / wanted) * weightSum;
        while (cursor < clouds.length - 1 && clouds[cursor].end < r) cursor++;
        const i = clouds[cursor].i;
        quotas.set(i, (quotas.get(i) || 0) + 1);
      }
      const points = [],
        seen = new Set();
      let visited = 0,
        approximate = false;
      for (const [index, quota] of quotas) {
        const p = this.engine.anchors[index],
          filter = rankFilter?.forGroup(
            this.universe + "/" + index.toString(36),
          ),
          cells = [];
        let total = 0,
          work = 0;
        const limit = Math.max(
          32,
          Math.floor((maxVisits * quota) / Math.max(wanted, 1)),
        );
        const visit = (start, count, lo, hi, zlo, zhi, seed, depth) => {
          if (count === 0n) return;
          let l = p.x,
            r = p.x,
            b = p.y,
            t = p.y;
          for (let k = 0; k < 4; k++) {
            const low = zlo[k],
              high = zhi[k],
              a = p.A[k],
              c = p.A[k + 4];
            l += Math.min(a * low, a * high);
            r += Math.max(a * low, a * high);
            b += Math.min(c * low, c * high);
            t += Math.max(c * low, c * high);
          }
          if (r < view.l || l > view.r || t < view.b || b > view.t) return;
          visited++;
          work++;
          const full = l >= view.l && r <= view.r && b >= view.b && t <= view.t;
          if (
            !full &&
            count > 1n &&
            depth < 96 &&
            work < limit &&
            Math.max(r - l, t - b) > width / 256
          ) {
            const axis = depth % 4,
              mid = (lo[axis] + hi[axis]) / 2,
              zmid = normalQuantile(mid),
              first = (count + 1n) / 2n,
              oldLo = lo[axis],
              oldHi = hi[axis],
              oldZlo = zlo[axis],
              oldZhi = zhi[axis];
            for (let child = 0; child < 2; child++) {
              const spatial = child ^ (seed & 1);
              if (spatial) {
                lo[axis] = mid;
                zlo[axis] = zmid;
              } else {
                hi[axis] = mid;
                zhi[axis] = zmid;
              }
              visit(
                start + (child ? first : 0n),
                child ? count - first : first,
                lo,
                hi,
                zlo,
                zhi,
                mix(seed ^ Math.imul(child + 1, 0x9e3779b9)),
                depth + 1,
              );
              lo[axis] = oldLo;
              hi[axis] = oldHi;
              zlo[axis] = oldZlo;
              zhi[axis] = oldZhi;
            }
            return;
          }
          if (!full && count > 1n && work >= limit) approximate = true;
          const eligible = filter ? filter.count(start, count) : count;
          if (!eligible) return;
          const coverage = full
            ? 1
            : (Math.max(0, Math.min(r, view.r) - Math.max(l, view.l)) *
                Math.max(0, Math.min(t, view.t) - Math.max(b, view.b))) /
              Math.max((r - l) * (t - b), 1e-30);
          const weight = Number(eligible) * coverage;
          if (weight > 0) {
            total += weight;
            cells.push({ start, eligible, end: total });
          }
        };
        visit(
          0n,
          this.counts[index],
          [0, 0, 0, 0],
          [1, 1, 1, 1],
          [-8, -8, -8, -8],
          [8, 8, 8, 8],
          this.seeds[index],
          0,
        );
        let made = 0;
        for (let k = 0; k < quota * 24 && made < quota && total > 0; k++) {
          const r = ((k * 0.6180339887498949 + 0.5) % 1) * total;
          let l = 0,
            h = cells.length - 1;
          while (l < h) {
            const m = (l + h) >> 1;
            if (cells[m].end <= r) l = m + 1;
            else h = m;
          }
          const cell = cells[l],
            u = P.hash(index + ":" + k, 23),
            ordinal =
              (cell.eligible * BigInt(Math.floor(u * 9007199254740992))) >> 53n,
            rank = filter
              ? filter.at(cell.start, ordinal)
              : cell.start + ordinal,
            point = this.evaluate(index, rank);
          if (
            point.x >= view.l &&
            point.x <= view.r &&
            point.y >= view.b &&
            point.y <= view.t &&
            !seen.has(point.id)
          ) {
            points.push(point);
            seen.add(point.id);
            made++;
          }
        }
      }
      return {
        points,
        expected,
        visited,
        approximate,
        wanted,
        shortfall: Math.max(0, wanted - points.length),
        clouds: quotas.size,
      };
    }
  }
  const api = {
    AnchorPopulation,
    normalQuantile,
    normalCDF,
    covariance,
    gaussianDensity: density,
    rectangleProbability,
  };
  Object.assign(root, api);
  if (typeof module !== "undefined") module.exports = api;
})(typeof window !== "undefined" ? window : globalThis);
