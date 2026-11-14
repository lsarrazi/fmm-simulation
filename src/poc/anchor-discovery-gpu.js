/* Catalogue discovery, independent of camera history.
 * Traverse immutable rank intervals in parallel. The GPU returns real IDs,
 * never newly sampled positions. CPU evaluation remains authoritative at edges.
 */
(function (root) {
  "use strict";
  const shader = `
struct Cloud { centre:vec4<f32>, ax:vec4<f32>, ay:vec4<f32>, hashes:vec4<u32>, catalogue:vec4<u32> }
struct Config { sizes:vec4<u32>, bounds:vec4<f32> }
struct Node { start:u32, count:u32, seed:u32, depth:u32, lo:vec4<f32>, bounds:vec4<f32> }
struct Output { found:atomic<u32>, visited:atomic<u32>, limited:atomic<u32>, padding:atomic<u32>, items:array<vec2<u32>> }
@group(0) @binding(0) var<storage,read> clouds:array<Cloud>;
@group(0) @binding(1) var<uniform> cfg:Config;
@group(0) @binding(2) var<storage,read_write> out:Output;
@group(0) @binding(3) var<storage,read> quantileBounds:array<vec2<f32>>;
fn mix32(v:u32)->u32 {var x=(v^(v>>16u))*0x7feb352du;x=(x^(x>>15u))*0x846ca68bu;return x^(x>>16u);}
fn quantilePair(p:f32,tail:f32)->f32 {
 if(p<0.02425||tail<0.02425){let q=sqrt(-2.0*log(max(1e-30,min(p,tail))));let v=(((((-0.007784894*q-0.32239646)*q-2.4007583)*q-2.5497325)*q+4.3746643)*q+2.938164)/((((0.0077846957*q+0.32246712)*q+2.4451342)*q+3.7544086)*q+1.0);return select(-v,v,p<0.5);}
 let q=p-0.5;let r=q*q;return (((((-39.69683*r+220.9461)*r-275.9285)*r+138.35776)*r-30.664799)*r+2.5066283)*q/(((((-54.4761*r+161.58583)*r-155.69897)*r+66.801315)*r-13.280682)*r+1.0);
}
// At most 31 binary splits, cycling four axes: all interval endpoints
// lie on the 1/256 lattice. Reuse this tiny table, not logarithms per node.
fn lower(p:f32)->f32 {return quantileBounds[u32(round(p*256.0))].x;}
fn upper(p:f32)->f32 {return quantileBounds[u32(round(p*256.0))].y;}
fn firstIn(start:u32,residue:u32,stride:u32)->u32 {
 if(residue>=start){return residue;}if(stride==0u){return 0xffffffffu;}
 let delta=start-residue;return residue+(delta/stride+select(0u,1u,delta%stride!=0u))*stride;
}
fn append(index:u32,rank:u32){let slot=atomicAdd(&out.found,1u);if(slot<cfg.sizes.y){out.items[slot]=vec2<u32>(index,rank);}}
fn evaluate(index:u32,rank:u32){
 let cloud=clouds[index];var remaining=cloud.catalogue.x;var local=rank;var seed=cloud.catalogue.y;var lo=vec4<f32>(0);var hi=vec4<f32>(1);var depth=0u;
 loop {if(remaining<=1u||depth>=32u){break;}let first=remaining/2u+remaining%2u;let child=select(1u,0u,local<first);let axis=depth%4u;let mid=0.5*(lo[axis]+hi[axis]);if((child^(seed&1u))!=0u){lo[axis]=mid;}else{hi[axis]=mid;}if(child!=0u){local-=first;remaining-=first;}else{remaining=first;}seed=mix32(seed^((child+1u)*0x9e3779b9u));depth++;}
 var digits:array<u32,10>;var n=0u;var number=rank;
 loop {digits[n]=number%10u+48u;n++;number/=10u;if(number==0u){break;}}
 var hashes=cloud.hashes;for(var j=n;j>0u;j--){hashes=(hashes^vec4<u32>(digits[j-1u]))*vec4<u32>(16777619u);}
 var z=vec4<f32>(0);
 for(var k=0u;k<4u;k++){
  let hash=mix32(hashes[k]);let u=(f32(hash)+0.5)/4294967296.0;let complement=(f32(~hash)+0.5)/4294967296.0;let p=lo[k]+(hi[k]-lo[k])*u;let tail=(1.0-hi[k])+(hi[k]-lo[k])*complement;
  z[k]=quantilePair(p,tail);
 }
 let point=cloud.centre.xy+vec2<f32>(dot(cloud.ax,z),dot(cloud.ay,z));let e=cloud.centre.z;
 if((point.x>=cfg.bounds.x-e&&point.x<=cfg.bounds.y+e&&point.y>=cfg.bounds.z-e&&point.y<=cfg.bounds.w+e)){append(index,rank);}
}
@compute @workgroup_size(64) fn discover(@builtin(global_invocation_id) id:vec3<u32>){
 let index=id.x;if(index>=cfg.sizes.x){return;}let cloud=clouds[index];if(cloud.catalogue.x==0u){return;}
 let radius=8.0*vec2<f32>(dot(abs(cloud.ax),vec4<f32>(1)),dot(abs(cloud.ay),vec4<f32>(1)));
 var stack:array<Node,40>;stack[0]=Node(0u,cloud.catalogue.x,cloud.catalogue.y,0u,vec4<f32>(0),vec4<f32>(cloud.centre.xy-radius,cloud.centre.xy+radius));var size=1u;var work=0u;
 loop {
  if(size==0u||atomicLoad(&out.found)>cfg.sizes.y){break;}
  size--;let node=stack[size];if(node.count==0u){continue;}
  let first=firstIn(node.start,cloud.catalogue.z,cloud.catalogue.w);let end=node.start+node.count;if(first>=end){continue;}
  let eligible=select(1u,1u+(end-1u-first)/max(1u,cloud.catalogue.w),cloud.catalogue.w!=0u);
  let l=node.bounds.xy;let r=node.bounds.zw;
  let e=cloud.centre.z;if(r.x<cfg.bounds.x-e||l.x>cfg.bounds.y+e||r.y<cfg.bounds.z-e||l.y>cfg.bounds.w+e){continue;}
  work++;
  // Share the work budget across clouds: a single difficult cloud must not
  // force thinning while the total query is still inexpensive. Claims in
  // blocks avoid a contended global atomic for every visited tree node.
  if((work&63u)==1u&&atomicAdd(&out.visited,64u)>=cfg.sizes.z){atomicAdd(&out.limited,1u);break;}
  if(l.x>=cfg.bounds.x+e&&r.x<=cfg.bounds.y-e&&l.y>=cfg.bounds.z+e&&r.y<=cfg.bounds.w-e){
   var rank=first;loop{if(rank>=end||atomicLoad(&out.found)>cfg.sizes.y){break;}append(index,rank);if(cloud.catalogue.w==0u){break;}rank+=cloud.catalogue.w;}continue;
  }
  if(eligible==1u||node.count==1u){evaluate(index,first);continue;}
  if(size+2u>40u){atomicAdd(&out.limited,1u);break;}
  // Only the split axis changes. Carry the projected bounds instead of
  // recomputing eight inverse-normal values for every node, including rejects.
  // Dyadic interval extent is determined by depth, so no extra stack memory.
  let axis=node.depth%4u;let extent=exp2(-f32(node.depth/4u));let mid=node.lo[axis]+0.5*extent;let half=node.count/2u+node.count%2u;
  let a=vec2<f32>(cloud.ax[axis],cloud.ay[axis]);let low=lower(node.lo[axis]);let high=upper(node.lo[axis]+extent);let ml=lower(mid);let mh=upper(mid);
  let oldMin=min(a*low,a*high);let oldMax=max(a*low,a*high);
  for(var child=0u;child<2u;child++){
   var lo=node.lo;let upperHalf=(child^(node.seed&1u))!=0u;
   if(upperHalf){lo[axis]=mid;}
   let nl=select(low,ml,upperHalf);let nh=select(mh,high,upperHalf);
   let bounds=vec4<f32>(l+(min(a*nl,a*nh)-oldMin),r+(max(a*nl,a*nh)-oldMax));
   stack[size]=Node(node.start+select(0u,half,child!=0u),select(half,node.count-half,child!=0u),mix32(node.seed^((child+1u)*0x9e3779b9u)),node.depth+1u,lo,bounds);size++;
  }
 }
}
`;
  class AnchorDiscoveryGPU {
    static async create({ device, population }) {
      // The current POC has about 10^7 ranks per anchor. Large theoretical
      // catalogues keep their existing CPU path until multiword GPU ranks exist.
      if (population.counts.some((n) => n > 2147483647n)) return null;
      const result = new AnchorDiscoveryGPU(device, population);
      try {
        await result.initialize();
        return result;
      } catch (error) {
        result.destroy();
        throw error;
      }
    }
    constructor(device, population) {
      this.device = device;
      this.population = population;
      this.capacity = 0;
    }
    async initialize() {
      const module = this.device.createShaderModule({
        label: "Immutable 4D catalogue discovery",
        code: shader,
      });
      const errors = (await module.getCompilationInfo()).messages.filter(
        (m) => m.type === "error",
      );
      if (errors.length)
        throw Error(errors.map((m) => m.lineNum + ": " + m.message).join("\n"));
      this.pipeline = await this.device.createComputePipelineAsync({
        layout: "auto",
        compute: { module, entryPoint: "discover" },
      });
      const n = this.population.counts.length;
      this.packed = new ArrayBuffer(n * 80);
      this.floats = new Float32Array(this.packed);
      this.ints = new Uint32Array(this.packed);
      this.clouds = this.device.createBuffer({
        label: "Catalogue clouds",
        size: this.packed.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      });
      this.config = this.device.createBuffer({
        size: 32,
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
      });
      const bounds = new Float32Array(257 * 2);
      for (let i = 0; i <= 256; i++) {
        const p = i / 256;
        bounds[i * 2] = i === 0 ? -8 : root.normalQuantile(p - 1.2e-7) - 1e-4;
        bounds[i * 2 + 1] =
          i === 256 ? 8 : root.normalQuantile(p + 1.2e-7) + 1e-4;
      }
      this.quantiles = this.device.createBuffer({
        size: bounds.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      });
      this.device.queue.writeBuffer(this.quantiles, 0, bounds);
      this.prefixes = this.population.counts.map((_, i) => {
        const prefix = this.population.universe + "/" + i.toString(36) + "/";
        return Array.from({ length: 4 }, (_, k) => {
          let hash = 2166136261 ^ (101 + k * 19);
          for (let j = 0; j < prefix.length; j++)
            hash = Math.imul(hash ^ prefix.charCodeAt(j), 16777619);
          return hash >>> 0;
        });
      });
      this.luminousRanks = this.population.counts.map((_, i) => {
        const key = this.population.universe + "/" + i.toString(36),
          P = root.StellarPopulation;
        const salt =
          (BigInt(Math.floor(P.hash(key, 503) * 4294967296)) << 21n) |
          BigInt(Math.floor(P.hash(key, 509) * 2097152));
        let bits = ((1n << 53n) - 1n) ^ salt,
          reversed = 0n;
        for (let k = 0; k < 53; k++) {
          reversed = (reversed << 1n) | (bits & 1n);
          bits >>= 1n;
        }
        return reversed;
      });
    }
    allocate(capacity) {
      if (this.capacity === capacity) return;
      this.output?.destroy();
      this.readback?.destroy();
      this.capacity = capacity;
      const bytes = 16 + capacity * 8;
      this.output = this.device.createBuffer({
        size: bytes,
        usage:
          GPUBufferUsage.STORAGE |
          GPUBufferUsage.COPY_SRC |
          GPUBufferUsage.COPY_DST,
      });
      this.readback = this.device.createBuffer({
        size: bytes,
        usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
      });
      this.group = this.device.createBindGroup({
        layout: this.pipeline.getBindGroupLayout(0),
        entries: [this.clouds, this.config, this.output, this.quantiles].map(
          (buffer, binding) => ({ binding, resource: { buffer } }),
        ),
      });
    }
    async sample({
      x = 0,
      y = 0,
      width = 4,
      aspectX = 1,
      aspectY = 1,
      budget = 8192,
      rankFilter = null,
    } = {}) {
      const start = performance.now();
      if (!budget || rankFilter?.fraction === 0)
        return {
          points: [],
          complete: true,
          approximate: false,
          visited: 0,
          wanted: 0,
          shortfall: 0,
          backend: "GPU catalogue",
        };
      const pop = this.population,
        n = pop.counts.length,
        capacity = Math.max(32768, Math.min(262144, budget * 4));
      this.allocate(capacity);
      const expected = pop.windowProjection({
        x,
        y,
        width,
        aspectX,
        aspectY,
      }).expected;
      // Start with a uniform luminosity tail when the analytic population
      // estimate already exceeds the return buffer; avoid knowingly overflowing.
      let extra = Math.max(
          0,
          Math.ceil(
            Math.log2(
              (expected * (rankFilter?.fraction ?? 1)) / (capacity * 0.9),
            ),
          ),
        ),
        visits = 0;
      for (let attempt = 0; attempt < 8; attempt++) {
        const strideBase = rankFilter ? 1 / rankFilter.fraction : 1;
        for (let i = 0; i < n; i++) {
          const p = pop.engine.anchors[i],
            o = i * 20,
            filter = rankFilter?.forGroup(pop.universe + "/" + i.toString(36));
          const base = filter ? Number(filter.at(0n, 0n)) : 0,
            stride = strideBase * 2 ** extra;
          // Narrow the same luminosity-tail congruence as the CPU catalogue.
          // It is spatially uniform and favours intrinsically luminous stars.
          const residue = extra
              ? Number(this.luminousRanks[i] & (BigInt(stride) - 1n))
              : base,
            count = Number(pop.counts[i]);
          this.ints[o + 16] = residue >= count ? 0 : count;
          this.ints[o + 17] = pop.seeds[i];
          this.ints[o + 18] = residue < count ? residue : 0;
          this.ints[o + 19] = stride <= 2147483648 ? stride : 0;
          const cx = (p.x - x) / width,
            cy = (p.y - y) / width;
          let span = 0;
          this.floats[o] = cx;
          this.floats[o + 1] = cy;
          for (let k = 0; k < 4; k++) {
            this.floats[o + 4 + k] = p.A[k] / width;
            this.floats[o + 8 + k] = p.A[k + 4] / width;
            span += Math.abs(p.A[k] / width) + Math.abs(p.A[k + 4] / width);
          }
          this.floats[o + 2] =
            0.00001 + 0.001 * span + 0.000002 * (Math.abs(cx) + Math.abs(cy));
          this.ints.set(this.prefixes[i], o + 12);
        }
        this.device.queue.writeBuffer(this.clouds, 0, this.packed);
        const config = new ArrayBuffer(32);
        new Uint32Array(config).set([n, capacity, 8388608, 0]);
        new Float32Array(config).set(
          [-aspectX / 2, aspectX / 2, -aspectY / 2, aspectY / 2],
          4,
        );
        this.device.queue.writeBuffer(this.config, 0, config);
        const e = this.device.createCommandEncoder();
        e.clearBuffer(this.output);
        const pass = e.beginComputePass();
        pass.setPipeline(this.pipeline);
        pass.setBindGroup(0, this.group);
        pass.dispatchWorkgroups(Math.ceil(n / 64));
        pass.end();
        e.copyBufferToBuffer(
          this.output,
          0,
          this.readback,
          0,
          16 + capacity * 8,
        );
        this.device.queue.submit([e.finish()]);
        await this.readback.mapAsync(GPUMapMode.READ);
        const data = new Uint32Array(this.readback.getMappedRange()),
          found = data[0],
          limited = data[2];
        visits += data[1];
        if (found > capacity || limited) {
          this.readback.unmap();
          extra++;
          continue;
        }
        const points = [];
        for (let j = 0; j < found; j++) {
          const p = pop.evaluate(data[4 + j * 2], BigInt(data[5 + j * 2]));
          if (
            Math.abs(p.x - x) <= (width * aspectX) / 2 &&
            Math.abs(p.y - y) <= (width * aspectY) / 2
          )
            points.push(p);
        }
        this.readback.unmap();
        return {
          points,
          expected,
          visited: visits,
          complete: true,
          approximate: extra > 0,
          shortfall: 0,
          wanted: points.length,
          thinning: 2 ** extra,
          backend: "GPU catalogue",
          ms: performance.now() - start,
        };
      }
      // Never publish a spatially biased prefix of an unfinished traversal.
      const fallback = pop.sample({
        x,
        y,
        width,
        aspectX,
        aspectY,
        budget,
        rankFilter,
      });
      return {
        ...fallback,
        complete: false,
        approximate: true,
        backend: "CPU fallback",
        ms: performance.now() - start,
      };
    }
    destroy() {
      this.clouds?.destroy();
      this.config?.destroy();
      this.quantiles?.destroy();
      this.output?.destroy();
      this.readback?.destroy();
    }
  }
  root.AnchorDiscoveryGPU = AnchorDiscoveryGPU;
})(typeof window !== "undefined" ? window : globalThis);
