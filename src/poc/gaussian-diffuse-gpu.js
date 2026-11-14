/* Direct Gaussian projection into the native pixel buffer.
 * The instrument's Gaussian core and an approximation of the pixel box are folded
 * into the covariance, with determinant normalization to preserve total flux.
 * No camera-sized coarse density raster, temporal blending or CPU readback.
 */
(function (root) {
  "use strict";
  const shader = `
struct Gaussian { metric:vec4<f32>, shape:vec4<f32>, bounds:vec4<f32> }
struct Config { camera:vec4<f32>, population:vec4<f32>, color:vec4<f32> }
@group(0) @binding(0) var<storage,read> gaussians:array<Gaussian>;
@group(0) @binding(1) var<uniform> cfg:Config;
struct Vertex { @builtin(position) position:vec4<f32>, @location(0) @interpolate(flat) index:u32 }
@vertex fn vertex(@builtin(vertex_index) vertex:u32,@builtin(instance_index) index:u32)->Vertex {
 let corners=array<vec2<f32>,6>(vec2<f32>(-1,-1),vec2<f32>(1,-1),vec2<f32>(-1,1),vec2<f32>(-1,1),vec2<f32>(1,-1),vec2<f32>(1,1));
 let g=gaussians[index];let pixel=g.metric.xy+corners[vertex]*g.shape.zw;
 var out:Vertex;out.position=vec4<f32>(2.0*pixel.x/cfg.camera.x-1.0,1.0-2.0*pixel.y/cfg.camera.y,0,1);out.index=index;return out;
}
@fragment fn splat(in:Vertex)->@location(0) vec4<f32> {
 let g=gaussians[in.index];let delta=in.position.xy-g.metric.xy;
 let q=g.metric.z*delta.x*delta.x+2.0*g.metric.w*delta.x*delta.y+g.shape.x*delta.y*delta.y;
 if(q>g.bounds.x){discard;}
 return vec4<f32>(g.shape.y*exp(-0.5*max(q,0.0)),0,0,0);
}
`;
  const compositeShader = `
struct Config { camera:vec4<f32>, population:vec4<f32>, color:vec4<f32> }
@group(0) @binding(0) var density:texture_2d<f32>;
@group(0) @binding(1) var<uniform> cfg:Config;
@vertex fn fullscreen(@builtin(vertex_index) i:u32)->@builtin(position) vec4<f32> {
 return vec4<f32>(f32((i<<1u)&2u)*2.0-1.0,f32(i&2u)*2.0-1.0,0,1);
}
fn crowding(lambda:f32)->f32 {
 if(lambda>50.0){return 1.0;}
 if(lambda<0.001){return lambda*lambda*(0.5-lambda/3.0+lambda*lambda/8.0);}
 return clamp(1.0-(1.0+lambda)*exp(-lambda),0.0,1.0);
}
@fragment fn composite(@builtin(position) p:vec4<f32>)->@location(0) vec4<f32> {
 let rho=textureLoad(density,vec2<i32>(p.xy),0).x;
 let weight=max(cfg.population.z,crowding(rho*cfg.population.x*cfg.population.y));
 let delta=(p.xy-0.5*cfg.camera.xy)/cfg.camera.z;
 let light=rho*weight*cfg.population.w*cfg.camera.w/(1.0+dot(delta,delta));
 return vec4<f32>(min(vec3<f32>(10000.0),light*cfg.color.xyz),1);
}
`;

  class GaussianDiffuseGPU {
    static async create(options) {
      const result = new GaussianDiffuseGPU(options);
      try {
        await result.initialize();
        return result;
      } catch (error) {
        result.destroy();
        throw error;
      }
    }
    constructor({ device, population }) {
      this.device = device;
      this.population = population;
      this.capacity = 0;
      this.width = this.height = 0;
      this.work = {};
    }
    async initialize() {
      const d = this.device;
      d.pushErrorScope("validation");
      let scope = true;
      try {
        const module = d.createShaderModule({
          label: "Direct Gaussian density",
          code: shader,
        });
        const composite = d.createShaderModule({
          label: "Unresolved Gaussian light",
          code: compositeShader,
        });
        for (const m of [module, composite]) {
          const errors = (await m.getCompilationInfo()).messages.filter(
            (m) => m.type === "error",
          );
          if (errors.length)
            throw Error(
              errors.map((m) => `${m.lineNum}: ${m.message}`).join("\n"),
            );
        }
        this.pipeline = await d.createRenderPipelineAsync({
          layout: "auto",
          vertex: { module, entryPoint: "vertex" },
          fragment: {
            module,
            entryPoint: "splat",
            targets: [
              {
                format: "r16float",
                blend: {
                  color: {
                    srcFactor: "one",
                    dstFactor: "one",
                    operation: "add",
                  },
                  alpha: {
                    srcFactor: "one",
                    dstFactor: "one",
                    operation: "add",
                  },
                },
              },
            ],
          },
        });
        this.composite = await d.createRenderPipelineAsync({
          layout: "auto",
          vertex: { module: composite, entryPoint: "fullscreen" },
          fragment: {
            module: composite,
            entryPoint: "composite",
            targets: [{ format: "rgba16float" }],
          },
        });
        this.config = d.createBuffer({
          label: "Gaussian projection configuration",
          size: 48,
          usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
        });
        const error = await d.popErrorScope();
        scope = false;
        if (error) throw error;
      } catch (error) {
        if (scope) await d.popErrorScope().catch(() => {});
        throw error;
      }
    }
    resize(width, height) {
      if (this.width === width && this.height === height) return;
      this.density?.destroy();
      this.width = width;
      this.height = height;
      this.density = this.device.createTexture({
        label: "Native Gaussian density",
        size: [width, height],
        format: "r16float",
        usage:
          GPUTextureUsage.RENDER_ATTACHMENT |
          GPUTextureUsage.TEXTURE_BINDING |
          GPUTextureUsage.COPY_SRC,
      });
      this.compositeGroup = this.device.createBindGroup({
        layout: this.composite.getBindGroupLayout(0),
        entries: [
          { binding: 0, resource: this.density.createView() },
          { binding: 1, resource: { buffer: this.config } },
        ],
      });
    }
    encode(
      encoder,
      { target, projection, context, result, gain = 1, exposureEV = 0 },
    ) {
      const start = performance.now(),
        d = this.device;
      const {
        pixelWidth: width,
        pixelHeight: height,
        pixelSpan,
        sigmaPixels,
      } = context;
      this.resize(width, height);
      const clouds = this.population.gaussians();
      if (!this.instances || this.capacity < clouds.length) {
        this.instances?.destroy();
        this.capacity = 2 ** Math.ceil(Math.log2(Math.max(16, clouds.length)));
        this.upload = new Float32Array(this.capacity * 12);
        this.instances = d.createBuffer({
          label: "Direct Gaussian instances",
          size: this.upload.byteLength,
          usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
        });
        this.group = d.createBindGroup({
          layout: this.pipeline.getBindGroupLayout(0),
          entries: [
            { binding: 0, resource: { buffer: this.instances } },
            { binding: 1, resource: { buffer: this.config } },
          ],
        });
      }
      const scale = pixelSpan / projection.width,
        scale2 = scale * scale;
      const pixelFilter = sigmaPixels * sigmaPixels + 1 / 12;
      const colour = root.StellarPopulation.meanColor;
      const lightGain =
        gain *
        context.diffuseScale *
        result.geometry.fieldFactor *
        2 ** exposureEV;
      let count = 0;
      for (const { anchor: p, covariance: q, mass } of clouds) {
        // Subtract camera coordinates in double precision before GPU packing.
        const x = width / 2 + (p.x - projection.centerX) * scale;
        const y = height / 2 - (p.y - projection.centerY) * scale;
        const xx = q.xx * scale2 + pixelFilter,
          xy = -q.xy * scale2,
          yy = q.yy * scale2 + pixelFilter;
        const det =
          q.det * scale2 * scale2 +
          pixelFilter * (q.xx + q.yy) * scale2 +
          pixelFilter * pixelFilter;
        const amplitude = (mass * scale2) / (2 * Math.PI * Math.sqrt(det));
        const cut = Math.max(
          36,
          2 *
            Math.log(
              Math.max(1, (amplitude * lightGain * Math.max(...colour)) / 1e-8),
            ),
        );
        const rx = Math.sqrt(cut * xx),
          ry = Math.sqrt(cut * yy);
        if (x + rx < 0 || x - rx > width || y + ry < 0 || y - ry > height)
          continue;
        const o = count++ * 12;
        this.upload.set(
          [
            x,
            y,
            yy / det,
            -xy / det,
            xx / det,
            amplitude,
            rx,
            ry,
            cut,
            0,
            0,
            0,
          ],
          o,
        );
      }
      if (count)
        d.queue.writeBuffer(
          this.instances,
          0,
          this.upload.buffer,
          0,
          count * 48,
        );
      const policy = result.diffuseParameters;
      d.queue.writeBuffer(
        this.config,
        0,
        new Float32Array([
          width,
          height,
          result.geometry.focalPixels,
          lightGain,
          Number(this.population.population) *
            context.worldPerPixel ** 2 *
            policy.beam,
          policy.unrepresented,
          1 - policy.detailWeight,
          policy.meanRemainder,
          ...colour,
          0,
        ]),
      );
      const splat = encoder.beginRenderPass({
        label: "Gaussian density at native resolution",
        colorAttachments: [
          {
            view: this.density.createView(),
            clearValue: { r: 0, g: 0, b: 0, a: 0 },
            loadOp: "clear",
            storeOp: "store",
          },
        ],
      });
      splat.setPipeline(this.pipeline);
      splat.setBindGroup(0, this.group);
      if (count) splat.draw(6, count);
      splat.end();
      const composite = encoder.beginRenderPass({
        label: "Unresolved light without coarse interpolation",
        colorAttachments: [
          {
            view: target.createView(),
            clearValue: { r: 0, g: 0, b: 0, a: 1 },
            loadOp: "clear",
            storeOp: "store",
          },
        ],
      });
      composite.setPipeline(this.composite);
      composite.setBindGroup(0, this.compositeGroup);
      composite.draw(3);
      composite.end();
      this.work = {
        backend: "direct Gaussian GPU",
        clouds: clouds.length,
        visibleClouds: count,
        width,
        height,
        cpuMs: performance.now() - start,
      };
    }
    destroy() {
      this.instances?.destroy();
      this.config?.destroy();
      this.density?.destroy();
    }
  }
  root.GaussianDiffuseGPU = GaussianDiffuseGPU;
})(typeof window !== "undefined" ? window : globalThis);
