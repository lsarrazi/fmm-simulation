/* GPU optical renderer. Individual kernels, photospheres, their separable PSF,
 * diffuse interpolation, additive HDR composition and the sensor response stay
 * on the shared WebGPU device. No image/particle readback in the render loop.
 * Only the bounded source list and small collective projection are uploaded.
 */
(function(root){
'use strict';
const shader=`
struct Config { screen:vec4<f32>, diffuse:vec4<f32>, color:vec4<f32>, counts:vec4<u32>, viewport:vec4<f32> }
struct Star { center:vec4<f32>, color:vec4<f32>, shape:vec4<f32>, halo:vec4<f32>, tile:vec4<f32> }
struct Prepared { norms:vec4<f32>, blur:vec4<f32> }
struct Tile { start:u32, width:u32, height:u32, star:u32 }
@group(0) @binding(0) var<uniform> cfg:Config;
@group(0) @binding(1) var<storage,read> stars:array<Star>;
@group(0) @binding(2) var<storage,read_write> prepared:array<Prepared>;
@group(0) @binding(3) var<storage,read_write> raw:array<f32>;
@group(0) @binding(4) var<storage,read_write> blurred:array<f32>;
@group(0) @binding(5) var<storage,read> tiles:array<Tile>;
@group(0) @binding(6) var glare:texture_2d<f32>;
@group(0) @binding(7) var density:texture_2d<f32>;
@group(0) @binding(8) var<storage,read> visibleGlare:array<u32>;
var<workgroup> sums:array<vec4<f32>,256>;
fn profile(p:vec2<f32>,s:Star)->f32 {
 var value=0.0;for(var a=0u;a<2u;a++){for(var b=0u;b<2u;b++){
  let d=p+vec2<f32>(0.25+0.5*f32(a),0.25+0.5*f32(b))-s.center.xy;
  let q=dot(d,d)/(s.center.z*s.center.z);if(q<1.0){if(s.tile.w>=1.0){let cosine=clamp(s.tile.w-2.0,-1.0,1.0);let light=vec2<f32>(cos(s.tile.z),sin(s.tile.z))*sqrt(max(0.0,1.0-cosine*cosine));value+=max(0.0,dot(d/s.center.z,light)+sqrt(1.0-q)*cosine);}else{value+=0.4+0.6*sqrt(1.0-q);}}
 }}return value*0.25;
}
fn gaussianRadius(s:Star,width:f32)->f32 {let signal=max(1e-30,s.center.w*max(max(s.color.x,s.color.y),s.color.z));return ceil(width*max(8.0,sqrt(2.0*max(0.0,log(signal)-log(1e-8)))))+1.0;}
fn haloTaper(p:vec2<f32>,s:Star)->f32 {
 if(s.halo.x<0.0){return 1.0;}let u=clamp((s.halo.x-length(p))/max(s.shape.x,0.2*s.halo.x),0.0,1.0);return u*u*(3.0-2.0*u);
}
@compute @workgroup_size(256) fn prepare(@builtin(workgroup_id) group:vec3<u32>,@builtin(local_invocation_index) lane:u32){
 let i=group.x+group.y*65535u;if(i>=cfg.counts.x){return;}let s=stars[i];if(s.tile.w!=0.0){return;}var values=vec4<f32>(0.0);
 let pad=i32(gaussianRadius(s,s.shape.y));for(var k=i32(lane)-pad;k<=pad;k+=256){let d=floor(s.center.xy)+vec2<f32>(f32(k)+0.5)-s.center.xy;values.x+=exp(-0.5*d.x*d.x/(s.shape.y*s.shape.y));values.y+=exp(-0.5*d.y*d.y/(s.shape.y*s.shape.y));}
 let blurPad=i32(gaussianRadius(s,s.shape.x));for(var k=i32(lane)-blurPad;k<=blurPad;k+=256){values.w+=exp(-0.5*f32(k*k)/(s.shape.x*s.shape.x));}
 if(s.shape.z>0.0&&s.center.z<=32.0){let lo=vec2<i32>(floor(s.center.xy-vec2<f32>(s.center.z)));let hi=vec2<i32>(ceil(s.center.xy+vec2<f32>(s.center.z)));let side=hi.x-lo.x+1;let count=side*(hi.y-lo.y+1);
  for(var k=i32(lane);k<count;k+=256){values.z+=profile(vec2<f32>(lo+vec2<i32>(k%side,k/side)),s);}}
 sums[lane]=values;workgroupBarrier();for(var step=128u;step>0u;step/=2u){if(lane<step){sums[lane]+=sums[lane+step];}workgroupBarrier();}
 if(lane==0u){var integral=3.14159265359*s.center.z*s.center.z*0.8;if(s.tile.w>=1.0){let alpha=acos(clamp(s.tile.w-2.0,-1.0,1.0));integral=2.0/3.0*s.center.z*s.center.z*(sin(alpha)+(3.14159265359-alpha)*cos(alpha));}prepared[i].norms=vec4<f32>(sums[0].x*sums[0].y,select(integral,sums[0].z,s.center.z<=32.0),1.0,0.0);prepared[i].blur=vec4<f32>(sums[0].w,f32(blurPad),0.0,0.0);}workgroupBarrier();
 var norm=0.0;if(s.shape.w>0.0&&s.halo.x>=0.0){let side=textureDimensions(glare).x;let origin=f32(side-1u)*0.5;for(var k=lane;k<side*side;k+=256u){let p=vec2<u32>(k%side,k/side);norm+=textureLoad(glare,vec2<i32>(p),0).x*haloTaper(vec2<f32>(p)-vec2<f32>(origin),s);}}
 sums[lane]=vec4<f32>(norm,0.0,0.0,0.0);workgroupBarrier();for(var step=128u;step>0u;step/=2u){if(lane<step){sums[lane]+=sums[lane+step];}workgroupBarrier();}
 if(lane==0u&&s.shape.w>0.0&&s.halo.x>=0.0){prepared[i].norms.z=max(sums[0].x,1e-30);}
}
fn findTile(address:u32)->Tile {var lo=0u;var hi=cfg.counts.y;while(lo+1u<hi){let m=(lo+hi)/2u;if(tiles[m].start<=address){lo=m;}else{hi=m;}}return tiles[lo];}
@compute @workgroup_size(256) fn photospheres(@builtin(global_invocation_id) id:vec3<u32>){
 let address=id.x+id.y*65535u*256u;if(address>=cfg.counts.z){return;}let t=findTile(address);let s=stars[t.star];let k=address-t.start;
 raw[address]=profile(s.tile.xy+vec2<f32>(f32(k%t.width),f32(k/t.width)),s)*s.center.w*(1.0-s.shape.w)*s.shape.z/max(prepared[t.star].norms.y,1e-30);
}
@compute @workgroup_size(256) fn horizontalBlur(@builtin(global_invocation_id) id:vec3<u32>){
 let address=id.x+id.y*65535u*256u;if(address>=cfg.counts.z){return;}let t=findTile(address);let s=stars[t.star];let x=i32((address-t.start)%t.width);let pad=i32(prepared[t.star].blur.y);var value=0.0;
 for(var k=-pad;k<=pad;k++){if(x+k>=0&&x+k<i32(t.width)){value+=raw[u32(i32(address)+k)]*exp(-0.5*f32(k*k)/(s.shape.x*s.shape.x));}}blurred[address]=value/prepared[t.star].blur.x;
}
struct Vertex { @builtin(position) position:vec4<f32>, @location(0) @interpolate(flat) star:u32 }
@vertex fn starVertex(@builtin(vertex_index) vertex:u32,@builtin(instance_index) index:u32)->Vertex {
 let corners=array<vec2<f32>,6>(vec2<f32>(-1.0,-1.0),vec2<f32>(1.0,-1.0),vec2<f32>(-1.0,1.0),vec2<f32>(-1.0,1.0),vec2<f32>(1.0,-1.0),vec2<f32>(1.0,1.0));let s=stars[index];let extent=select(s.color.w,s.tile.z,(cfg.counts.w&2u)!=0u&&s.tile.w==0.0);let p=s.center.xy+corners[vertex]*extent;
 var out:Vertex;out.position=vec4<f32>(p.x/cfg.viewport.x*2.0-1.0,1.0-p.y/cfg.viewport.y*2.0,0.0,1.0);out.star=index;return out;
}
fn haloPixel(p:vec2<i32>,s:Star)->f32 {let side=i32(textureDimensions(glare).x);if(any(p<vec2<i32>(0))||any(p>=vec2<i32>(side))){return 0.0;}return textureLoad(glare,p,0).x*haloTaper(vec2<f32>(p)-vec2<f32>(f32(side-1)*0.5),s);}
fn coreValue(position:vec2<f32>,s:Star,norms:Prepared)->f32 {let d=position-s.center.xy;var value=0.0;
 if(s.shape.z<1.0){let gaussian=exp(-0.5*dot(d,d)/(s.shape.y*s.shape.y))/max(norms.norms.x,1e-30);if(all(abs(d)<=vec2<f32>(gaussianRadius(s,s.shape.y)+1.0))){value+=s.center.w*(1.0-s.shape.w)*(1.0-s.shape.z)*gaussian;}}
 if(s.shape.z>0.0){let p=vec2<i32>(floor(position)-s.tile.xy);let width=i32(s.halo.z);let height=i32(s.halo.w);let pad=i32(norms.blur.y);if(p.x>=0&&p.x<width){for(var k=-pad;k<=pad;k++){if(p.y+k>=0&&p.y+k<height){value+=blurred[u32(i32(s.halo.y)+(p.y+k)*width+p.x)]*exp(-0.5*f32(k*k)/(s.shape.x*s.shape.x))/norms.blur.x;}}}}
 return value;
}
@fragment fn starFragment(in:Vertex)->@location(0) vec4<f32>{
 let s=stars[in.star];let norms=prepared[in.star];let d=in.position.xy-s.center.xy;var value=0.0;
 if(s.tile.w!=0.0){
  let q=dot(d,d)/(s.center.z*s.center.z);let edge=max(1.0/s.center.z,0.00001);let coverage=1.0-smoothstep(1.0-edge,1.0+edge,q);if(coverage<=0.0){discard;}
  let nz=sqrt(max(0.0,1.0-q));var surface=0.4+0.6*nz;var integral=3.14159265359*s.center.z*s.center.z*0.8;
  if(s.tile.w>=1.0){let cosine=clamp(s.tile.w-2.0,-1.0,1.0);let light=vec2<f32>(cos(s.tile.z),sin(s.tile.z))*sqrt(max(0.0,1.0-cosine*cosine));surface=max(0.0,dot(d/s.center.z,light)+nz*cosine);let alpha=acos(cosine);integral=2.0/3.0*s.center.z*s.center.z*(sin(alpha)+(3.14159265359-alpha)*cos(alpha));}
  return vec4<f32>(min(vec3<f32>(10000.0),s.color.xyz*s.center.w*surface*coverage/max(integral,1e-30)),0.0);
 }
 value=coreValue(in.position.xy,s,norms);
 if(s.shape.w>0.0&&(cfg.counts.w&2u)==0u){let origin=f32(textureDimensions(glare).x-1u)*0.5;let p=floor(in.position.xy)-s.center.xy+vec2<f32>(origin);let lo=vec2<i32>(floor(p));let f=fract(p);let w=mix(mix(haloPixel(lo,s),haloPixel(lo+vec2<i32>(1,0),s),f.x),mix(haloPixel(lo+vec2<i32>(0,1),s),haloPixel(lo+vec2<i32>(1,1),s),f.x),f.y)/norms.norms.z;
  if(w*s.center.w*s.shape.w*max(max(s.color.x,s.color.y),s.color.z)>=1e-6){value+=s.center.w*s.shape.w*w;}}
 return vec4<f32>(min(vec3<f32>(10000.0),max(vec3<f32>(0.0),value*s.color.xyz)),0.0);
}
@vertex fn glareVertex(@builtin(vertex_index) vertex:u32,@builtin(instance_index) instance:u32)->Vertex {
 let tileSize=cfg.viewport.z;let side=u32(ceil(f32(textureDimensions(glare).x+2u)/tileSize));let address=visibleGlare[instance];let index=address/(side*side);let tile=address%(side*side);let s=stars[index];
 let corners=array<vec2<f32>,6>(vec2<f32>(0.0,0.0),vec2<f32>(1.0,0.0),vec2<f32>(0.0,1.0),vec2<f32>(0.0,1.0),vec2<f32>(1.0,0.0),vec2<f32>(1.0,1.0));
 let p=floor(s.center.xy)-vec2<f32>(f32(textureDimensions(glare).x-1u)*0.5+1.0)+vec2<f32>(f32(tile%side),f32(tile/side))*tileSize+corners[vertex]*tileSize;
 var out:Vertex;out.position=vec4<f32>(p.x/cfg.viewport.x*2.0-1.0,1.0-p.y/cfg.viewport.y*2.0,0.0,1.0);out.star=index;return out;
}
@fragment fn glareFragment(in:Vertex)->@location(0) vec4<f32>{
 let s=stars[in.star];let origin=f32(textureDimensions(glare).x-1u)*0.5;let p=floor(in.position.xy)-s.center.xy+vec2<f32>(origin);let lo=vec2<i32>(floor(p));let f=fract(p);
 let w=mix(mix(haloPixel(lo,s),haloPixel(lo+vec2<i32>(1,0),s),f.x),mix(haloPixel(lo+vec2<i32>(0,1),s),haloPixel(lo+vec2<i32>(1,1),s),f.x),f.y);
 if(w*s.center.w*s.shape.w*max(max(s.color.x,s.color.y),s.color.z)<1e-6){discard;}
 let halo=w*s.center.w*s.shape.w*s.color.xyz;var core=vec3<f32>(0.0);
 if(all(abs(in.position.xy-s.center.xy)<=vec2<f32>(s.tile.z))){core=coreValue(in.position.xy,s,prepared[in.star])*s.color.xyz;}
 return vec4<f32>(min(vec3<f32>(10000.0),core+halo)-min(vec3<f32>(10000.0),core),0.0);
}
@vertex fn fullscreen(@builtin(vertex_index) i:u32)->@builtin(position) vec4<f32>{return vec4<f32>(f32((i<<1u)&2u)*2.0-1.0,f32(i&2u)*2.0-1.0,0.0,1.0);}
fn densityPixel(p:vec2<i32>)->f32 {let side=i32(textureDimensions(density).x);if(any(p<vec2<i32>(0))||any(p>=vec2<i32>(side))){return 0.0;}return textureLoad(density,p,0).x;}
@fragment fn diffuseFragment(@builtin(position) p:vec4<f32>)->@location(0) vec4<f32>{
 let sensor=p.xy-cfg.viewport.xy*0.5+vec2<f32>(cfg.screen.x*0.5);let uv=vec2<f32>(sensor.x,cfg.screen.x-sensor.y)/cfg.screen.x*cfg.screen.y-vec2<f32>(0.5);let lo=vec2<i32>(floor(uv));let f=fract(uv);let mass=mix(mix(densityPixel(lo),densityPixel(lo+vec2<i32>(1,0)),f.x),mix(densityPixel(lo+vec2<i32>(0,1)),densityPixel(lo+vec2<i32>(1,1)),f.x),f.y);
 let delta=(p.xy-cfg.viewport.xy*0.5)/cfg.screen.z;let value=mass*cfg.diffuse.x*cfg.diffuse.y*cfg.diffuse.z/(1.0+dot(delta,delta));return vec4<f32>(min(vec3<f32>(10000.0),value*cfg.color.xyz),1.0);
}`;
// Generate the instrument kernel on GPU as well. A gather of the rotated
// integer ray samples reproduces the CPU reference's bilinear deposition,
// without building/sorting a JavaScript Map when stars first become visible.
const glareShader=`
struct Config { dimensions:vec4<u32>, style:vec4<f32>, shape:vec4<f32> }
@group(0) @binding(0) var<uniform> cfg:Config;
@group(0) @binding(1) var<storage,read_write> values:array<vec2<f32>>;
@group(0) @binding(2) var<storage,read_write> groups:array<vec2<f32>>;
@group(0) @binding(3) var<storage,read_write> normalization:vec2<f32>;
@group(0) @binding(4) var output:texture_storage_2d<r32float,write>;
@group(0) @binding(5) var<storage,read_write> tileBounds:array<f32>;
var<workgroup> partial:array<vec2<f32>,256>;
@compute @workgroup_size(256) fn generate(@builtin(global_invocation_id) id:vec3<u32>,@builtin(workgroup_id) group:vec3<u32>,@builtin(local_invocation_index) lane:u32){
 let index=id.x+id.y*65535u*256u;let side=cfg.dimensions.x;let sigma=cfg.style.x;var value=vec2<f32>(0.0);
 if(index<side*side){let p=vec2<f32>(f32(index%side),f32(index/side))-vec2<f32>(f32(cfg.dimensions.y));let haloScale=cfg.shape.x/12.0;let r2=dot(p,p)/(sigma*sigma*haloScale*haloScale);value.x=pow(1.0+r2/4.0,-1.7)*exp(-r2/288.0);
  if(cfg.style.z>0.0){let axes=select(6u,3u,cfg.style.w==1.0);let maxT=i32(cfg.dimensions.y);
   for(var a=0u;a<axes;a++){let angle=cfg.style.y+f32(a)*3.141592653589793/f32(axes);let c=cos(angle);let s=sin(angle);let reach=abs(c)+abs(s);let along=dot(p,vec2<f32>(c,s));let across=dot(p,vec2<f32>(-s,c));let short=cfg.style.w==0.0&&a%2u==1u;let falloff=select(1.0,7.0/12.0,short)*cfg.shape.y*sigma;let amplitude=select(1.0,0.5,short);
    let t0=max(-maxT,i32(ceil(along-reach)));let t1=min(maxT,i32(floor(along+reach)));
    for(var t=t0;t<=t1;t++){let width=max(cfg.shape.z*sigma,cfg.shape.w*sigma*(1.0-abs(f32(t))/(cfg.shape.y*56.0/12.0*sigma)));let maxQ=i32(ceil(f32(cfg.dimensions.w)*width));let q0=max(-maxQ,i32(ceil(across-reach)));let q1=min(maxQ,i32(floor(across+reach)));
     for(var q=q0;q<=q1;q++){let position=vec2<f32>(f32(t)*c-f32(q)*s,f32(t)*s+f32(q)*c);let triangle=max(vec2<f32>(0.0),vec2<f32>(1.0)-abs(p-position));let alongScale=f32(t)/falloff;let wing=f32(t)/(3.0*sigma);let acrossScale=f32(q)/width;let w=amplitude*exp(-0.5*alongScale*alongScale)/sqrt(1.0+wing*wing)*exp(-0.5*acrossScale*acrossScale)/width;value.y+=w*triangle.x*triangle.y;}
    }
   }
  }values[index]=value;
 }
 partial[lane]=value;workgroupBarrier();for(var step=128u;step>0u;step/=2u){if(lane<step){partial[lane]+=partial[lane+step];}workgroupBarrier();}if(lane==0u){groups[group.x+group.y*65535u]=partial[0];}
}
@compute @workgroup_size(256) fn reduce(@builtin(local_invocation_index) lane:u32){var sum=vec2<f32>(0.0);for(var i=lane;i<cfg.dimensions.z;i+=256u){sum+=groups[i];}partial[lane]=sum;workgroupBarrier();for(var step=128u;step>0u;step/=2u){if(lane<step){partial[lane]+=partial[lane+step];}workgroupBarrier();}if(lane==0u){normalization=max(partial[0],vec2<f32>(1e-30));}}
@compute @workgroup_size(256) fn normalize(@builtin(global_invocation_id) id:vec3<u32>){let index=id.x+id.y*65535u*256u;let side=cfg.dimensions.x;if(index>=side*side){return;}let v=values[index]/normalization;textureStore(output,vec2<i32>(i32(index%side),i32(index/side)),vec4<f32>(v.x*(1.0-cfg.style.z)+v.y*cfg.style.z,0.0,0.0,0.0));}
@compute @workgroup_size(256) fn bounds(@builtin(workgroup_id) group:vec3<u32>,@builtin(local_invocation_index) lane:u32){
 let side=cfg.dimensions.x;let tiles=(side+17u)/16u;let tile=group.x;let origin=vec2<i32>(i32(tile%tiles)*16-2,i32(tile/tiles)*16-2);var maximum=0.0;
 for(var i=lane;i<324u;i+=256u){let p=origin+vec2<i32>(i32(i%18u),i32(i/18u));if(all(p>=vec2<i32>(0))&&all(p<vec2<i32>(i32(side)))){let v=values[u32(p.y)*side+u32(p.x)]/normalization;maximum=max(maximum,v.x*(1.0-cfg.style.z)+v.y*cfg.style.z);}}
 partial[lane]=vec2<f32>(maximum,0.0);workgroupBarrier();for(var step=128u;step>0u;step/=2u){if(lane<step){partial[lane]=max(partial[lane],partial[lane+step]);}workgroupBarrier();}if(lane==0u){tileBounds[tile]=partial[0].x;}
}
`;
// Visibility uses the actual normalized optical kernel, not a heuristic radius.
// Every retained tile has the identical fragment shader/cutoff as the full quad.
const glareVisibilityShader=`
struct Star { center:vec4<f32>, color:vec4<f32>, shape:vec4<f32>, halo:vec4<f32>, tile:vec4<f32> }
struct Params { count:u32, side:u32, width:u32, origin:u32, height:u32, tileSize:u32, unused1:u32, unused2:u32 }
struct Draw { vertexCount:u32, instanceCount:atomic<u32>, firstVertex:u32, firstInstance:u32 }
@group(0) @binding(0) var<storage,read> stars:array<Star>;
@group(0) @binding(1) var<storage,read> bounds:array<f32>;
@group(0) @binding(2) var<storage,read_write> visible:array<u32>;
@group(0) @binding(3) var<storage,read_write> draw:Draw;
@group(0) @binding(4) var<uniform> cfg:Params;
@compute @workgroup_size(256) fn visibility(@builtin(global_invocation_id) id:vec3<u32>){
 let address=id.x+id.y*65535u*256u;let tiles=cfg.side*cfg.side;
 if(address==0u){draw.vertexCount=6u;draw.firstVertex=0u;draw.firstInstance=0u;}
 if(address>=cfg.count*tiles){return;}let index=address/tiles;let tile=address%tiles;let s=stars[index];
 let signal=s.center.w*s.shape.w*max(max(s.color.x,s.color.y),s.color.z);
 if(s.tile.w!=0.0||signal<0.9999e-6||bounds[tile]*signal<0.9999e-6){return;}
 // Disjoint 16-pixel tiles are aligned with screen pixels even at fractional centres.
 let origin=floor(s.center.xy)-vec2<f32>(f32(cfg.origin));
 let lo=origin+vec2<f32>(f32(tile%cfg.side)*f32(cfg.tileSize),f32(tile/cfg.side)*f32(cfg.tileSize));
 if(any(lo+vec2<f32>(f32(cfg.tileSize))<=vec2<f32>(0.0))||any(lo>=vec2<f32>(f32(cfg.width),f32(cfg.height)))){return;}
 let slot=atomicAdd(&draw.instanceCount,1u);visible[slot]=address;
}
`;

// Expand unresolved optical footprints on the GPU, alongside normalization.
const coarseBoundsShader=`
struct Config { base:u32, side:u32, factor:u32, unused:u32 }
@group(0) @binding(0) var<storage,read> fine:array<f32>;
@group(0) @binding(1) var<storage,read_write> coarse:array<f32>;
@group(0) @binding(2) var<uniform> cfg:Config;
@compute @workgroup_size(256) fn reduce(@builtin(global_invocation_id) id:vec3<u32>){let i=id.x;if(i>=cfg.side*cfg.side){return;}let ox=(i%cfg.side)*cfg.factor;let oy=(i/cfg.side)*cfg.factor;var value=0.0;for(var y=0u;y<cfg.factor;y++){for(var x=0u;x<cfg.factor;x++){if(ox+x<cfg.base&&oy+y<cfg.base){value=max(value,fine[(oy+y)*cfg.base+ox+x]);}}}coarse[i]=value;}
`;
const opticalPackingShader=`
struct Star { center:vec4<f32>, color:vec4<f32>, shape:vec4<f32>, halo:vec4<f32>, tile:vec4<f32> }
struct Params { count:u32, radius:u32, unused:vec2<u32> }
@group(0) @binding(0) var<storage,read_write> stars:array<Star>;
@group(0) @binding(1) var<uniform> cfg:Params;
@compute @workgroup_size(256) fn pack(@builtin(global_invocation_id) id:vec3<u32>){
 let i=id.x+id.y*65535u*256u;if(i>=cfg.count){return;}var s=stars[i];if(s.shape.y>=0.0){return;}
 let sigma=s.shape.x;let fraction=s.shape.z*0.55*s.shape.w;let signal=max(1e-30,s.center.w*max(max(s.color.x,s.color.y),s.color.z));
 let pad=ceil(sigma*max(8.0,sqrt(2.0*max(0.0,log(signal)-log(1e-8)))))+1.0;
 let core=ceil(s.center.z+pad)+2.0;s.color.w=max(core,select(0.0,f32(cfg.radius)+2.0,fraction>0.0));
 s.shape=vec4<f32>(sigma,sigma,0.0,fraction);s.halo=vec4<f32>(-1.0,0.0,0.0,0.0);s.tile=vec4<f32>(0.0,0.0,core,0.0);stars[i]=s;
}
`;
const toneShader=`
struct Config { screen:vec4<f32>, diffuse:vec4<f32>, color:vec4<f32>, counts:vec4<u32>, viewport:vec4<f32> }
struct Presentation { transform:vec4<f32>, viewport:vec4<f32> }
@group(0) @binding(0) var hdr:texture_2d<f32>;
@group(0) @binding(1) var<uniform> cfg:Config;
@group(0) @binding(2) var<uniform> view:Presentation;
@vertex fn fullscreen(@builtin(vertex_index) i:u32)->@builtin(position) vec4<f32>{return vec4<f32>(f32((i<<1u)&2u)*2.0-1.0,f32(i&2u)*2.0-1.0,0.0,1.0);}
fn pixel(p:vec2<i32>)->vec3<f32>{let dimensions=vec2<i32>(textureDimensions(hdr));if(any(p<vec2<i32>(0))||any(p>=dimensions)){return vec3<f32>(0.0);}return min(textureLoad(hdr,p,0).xyz,vec3<f32>(10000.0));}
@fragment fn tone(@builtin(position) p:vec4<f32>)->@location(0) vec4<f32>{
 let source=(p.xy-view.viewport.xy*0.5-view.transform.yz)/view.transform.x+cfg.viewport.xy*0.5;let uv=source-vec2<f32>(0.5);let lo=vec2<i32>(floor(uv));let f=fract(uv);
 var radiance:vec3<f32>;if(view.transform.x==1.0&&all(view.transform.yz==vec2<f32>(0.0))&&all(view.viewport.xy==cfg.viewport.xy)){radiance=pixel(vec2<i32>(p.xy));}else{radiance=mix(mix(pixel(lo),pixel(lo+vec2<i32>(1,0)),f.x),mix(pixel(lo+vec2<i32>(0,1)),pixel(lo+vec2<i32>(1,1)),f.x),f.y);}var rgb=sqrt(max(vec3<f32>(0.0),radiance)/(vec3<f32>(1.0)+radiance));
 if((cfg.counts.w&1u)!=0u){let edge=abs(distance(source,vec2<f32>(cfg.diffuse.w,cfg.color.w))-cfg.screen.w)*view.transform.x;let alpha=clamp(0.5*cfg.viewport.w+0.5-edge,0.0,1.0);rgb=mix(rgb,vec3<f32>(255.0,80.0,80.0)/255.0,alpha);}return vec4<f32>(rgb,1.0);
}
`;
const guideShader=`
struct Guide { shape:vec4<f32> }
@group(0) @binding(0) var<storage,read> guides:array<Guide>;
struct GuideViewport { size:vec4<f32> }
@group(0) @binding(1) var<uniform> viewport:GuideViewport;
struct GuideVertex { @builtin(position) position:vec4<f32>, @location(0) local:vec2<f32>, @location(1) @interpolate(flat) shape:vec2<f32> }
@vertex fn vertex(@builtin(vertex_index) vertex:u32,@builtin(instance_index) index:u32)->GuideVertex {
 let corners=array<vec2<f32>,6>(vec2<f32>(-1.0,-1.0),vec2<f32>(1.0,-1.0),vec2<f32>(-1.0,1.0),vec2<f32>(-1.0,1.0),vec2<f32>(1.0,-1.0),vec2<f32>(1.0,1.0));
 let g=guides[index].shape;let local=corners[vertex]*(g.z+1.0);let p=g.xy+local;var out:GuideVertex;out.position=vec4<f32>(p.x/(0.5*viewport.size.x)-1.0,1.0-p.y/(0.5*viewport.size.y),0.0,1.0);out.local=local;out.shape=g.zw;return out;
}
@fragment fn fragment(in:GuideVertex)->@location(0) vec4<f32> {
 let edge=1.0-smoothstep(0.15,0.85,abs(length(in.local)-in.shape.x));let arc=(atan2(in.local.y,in.local.x)+3.14159265359)*in.shape.x;let phase=fract(arc/4.0)*4.0;let dash=smoothstep(0.0,0.35,phase)*(1.0-smoothstep(1.65,2.0,phase));
 return vec4<f32>(0.612,0.682,0.714,edge*dash*in.shape.y);
}
`;
class GPUStellarLightRenderer {
 static async create(options){const renderer=new GPUStellarLightRenderer(options);try{await renderer.initialize();return renderer;}catch(e){renderer.destroy();throw e;}}
 constructor({device,size=600,capacity=2048}){this.device=device;this.size=size;this.capacity=capacity;this.sourceLimit=Math.min(0xffffffff,Math.floor(Math.min(device.limits.maxStorageBufferBindingSize,device.limits.maxBufferSize)/80));this.backend='GPU';this.buffers=new Set();this.textures=new Set();this.glareCache=new Map();this.atlasCapacity=0;this.sparseGlare=true;this.glareTileCapacity=0;this.gpuPacking=true;this.tileBudgetBytes=Math.min(64*1024*1024,device.limits.maxStorageBufferBindingSize,device.limits.maxBufferSize);this.sourceScratch=new Float32Array(0);}
 buffer(label,size,usage){const b=this.device.createBuffer({label,size:Math.max(16,size),usage});this.buffers.add(b);return b;}
 texture(label,size,format,usage){const t=this.device.createTexture({label,size,format,usage});this.textures.add(t);return t;}
 async initialize(){const d=this.device,U=GPUBufferUsage;d.pushErrorScope('validation');let scope=true;try{
  const module=d.createShaderModule({label:'Stellar optics GPU',code:shader}),tone=d.createShaderModule({label:'Stellar HDR response GPU',code:toneShader}),glareModule=d.createShaderModule({label:'Instrument kernel GPU',code:glareShader});for(const m of [module,tone,glareModule]){const info=await m.getCompilationInfo(),errors=info.messages.filter(x=>x.type==='error');if(errors.length)throw Error(errors.map(x=>x.lineNum+':'+x.linePos+' '+x.message).join('\n'));}
  this.layout=d.createBindGroupLayout({entries:[{binding:0,visibility:GPUShaderStage.COMPUTE|GPUShaderStage.FRAGMENT|GPUShaderStage.VERTEX,buffer:{type:'uniform'}},...Array.from({length:5},(_,i)=>({binding:i+1,visibility:GPUShaderStage.COMPUTE|([1,2,4].includes(i+1)?GPUShaderStage.FRAGMENT:0)|(i===0?GPUShaderStage.VERTEX:0),buffer:{type:[1,5].includes(i+1)?'read-only-storage':'storage'}})),...[6,7].map(binding=>({binding,visibility:GPUShaderStage.COMPUTE|GPUShaderStage.FRAGMENT|GPUShaderStage.VERTEX,texture:{sampleType:'unfilterable-float'}})),{binding:8,visibility:GPUShaderStage.VERTEX,buffer:{type:'read-only-storage'}}]});
  const layout=d.createPipelineLayout({bindGroupLayouts:[this.layout]});this.pipelines={};for(const entryPoint of ['prepare','photospheres','horizontalBlur'])this.pipelines[entryPoint]=await d.createComputePipelineAsync({layout,compute:{module,entryPoint}});
  this.pipelines.diffuse=await d.createRenderPipelineAsync({layout,vertex:{module,entryPoint:'fullscreen'},fragment:{module,entryPoint:'diffuseFragment',targets:[{format:'rgba16float'}]}});
  this.pipelines.stars=await d.createRenderPipelineAsync({layout,vertex:{module,entryPoint:'starVertex'},fragment:{module,entryPoint:'starFragment',targets:[{format:'rgba16float',blend:{color:{srcFactor:'one',dstFactor:'one',operation:'add'},alpha:{srcFactor:'one',dstFactor:'one',operation:'add'}}}]}});
  this.pipelines.glare=await d.createRenderPipelineAsync({layout,vertex:{module,entryPoint:'glareVertex'},fragment:{module,entryPoint:'glareFragment',targets:[{format:'rgba16float',blend:{color:{srcFactor:'one',dstFactor:'one',operation:'add'},alpha:{srcFactor:'one',dstFactor:'one',operation:'add'}}}]}});
  const visibilityModule=d.createShaderModule({label:'Conservative optical tile visibility',code:glareVisibilityShader});this.visibilityPipeline=await d.createComputePipelineAsync({layout:'auto',compute:{module:visibilityModule,entryPoint:'visibility'}});
  const coarseModule=d.createShaderModule({label:'Conservative optical bounds pyramid',code:coarseBoundsShader});this.coarseBoundsPipeline=await d.createComputePipelineAsync({layout:'auto',compute:{module:coarseModule,entryPoint:'reduce'}});
  const packingModule=d.createShaderModule({label:'Unresolved source optical preparation GPU',code:opticalPackingShader});this.packingPipeline=await d.createComputePipelineAsync({layout:'auto',compute:{module:packingModule,entryPoint:'pack'}});this.packingConfig=this.buffer('Optical preparation parameters',16,U.UNIFORM|U.COPY_DST);
  this.glareDraw=this.buffer('Visible optical tiles indirect draw',16,U.STORAGE|U.INDIRECT|U.COPY_DST);this.glareVisibilityConfig=this.buffer('Optical tile visibility parameters',32,U.UNIFORM|U.COPY_DST);this.visibleGlare=this.buffer('Visible optical tiles',4,U.STORAGE);
  const glareLayout=d.createBindGroupLayout({entries:[{binding:0,visibility:GPUShaderStage.COMPUTE,buffer:{type:'uniform'}},...[1,2,3].map(binding=>({binding,visibility:GPUShaderStage.COMPUTE,buffer:{type:'storage'}})),{binding:4,visibility:GPUShaderStage.COMPUTE,storageTexture:{access:'write-only',format:'r32float'}},{binding:5,visibility:GPUShaderStage.COMPUTE,buffer:{type:'storage'}}]});this.glareLayout=glareLayout;const glarePipelineLayout=d.createPipelineLayout({bindGroupLayouts:[glareLayout]});this.glarePipelines={};for(const entryPoint of ['generate','reduce','normalize','bounds'])this.glarePipelines[entryPoint]=await d.createComputePipelineAsync({layout:glarePipelineLayout,compute:{module:glareModule,entryPoint}});
  this.format=navigator.gpu.getPreferredCanvasFormat();this.pipelines.tone=await d.createRenderPipelineAsync({layout:'auto',vertex:{module:tone,entryPoint:'fullscreen'},fragment:{module:tone,entryPoint:'tone',targets:[{format:this.format}]}});
  const guideModule=d.createShaderModule({label:'Size-based minor-body navigation rings',code:guideShader});const guideInfo=await guideModule.getCompilationInfo();if(guideInfo.messages.some(m=>m.type==='error'))throw Error(guideInfo.messages.map(m=>m.message).join('\n'));
  this.pipelines.guides=await d.createRenderPipelineAsync({layout:'auto',vertex:{module:guideModule,entryPoint:'vertex'},fragment:{module:guideModule,entryPoint:'fragment',targets:[{format:this.format,blend:{color:{srcFactor:'src-alpha',dstFactor:'one-minus-src-alpha',operation:'add'},alpha:{srcFactor:'one',dstFactor:'one-minus-src-alpha',operation:'add'}}}]}});
  this.guideViewport=this.buffer('Guide viewport aspect',16,U.UNIFORM|U.COPY_DST);this.config=this.buffer('Stellar optics parameters',80,U.UNIFORM|U.COPY_DST);this.presentConfig=this.buffer('GPU camera interpolation',32,U.UNIFORM|U.COPY_DST);this.ensureSourceCapacity(this.capacity);this.allocateAtlas(1);this.emptyGlare=this.texture('Unused optical kernel',[1,1],'r32float',GPUTextureUsage.TEXTURE_BINDING|GPUTextureUsage.COPY_DST);this.resize(this.size,this.size,true);
  const error=await d.popErrorScope();scope=false;if(error)throw Error(error.message);
 }catch(e){if(scope)await d.popErrorScope().catch(()=>{});throw e;}}
 ensureSourceCapacity(count){
  if(this.stars&&count<=this.capacity)return;
  if(!Number.isSafeInteger(count)||count<0||count>this.sourceLimit)throw Error('Limite matérielle du rendu optique : '+this.sourceLimit+' sources.');
  const capacity=Math.min(this.sourceLimit,2**Math.ceil(Math.log2(Math.max(64,count)))),U=GPUBufferUsage,stars=this.buffer('Bounded optical sources',capacity*80,U.STORAGE|U.COPY_DST),prepared=this.buffer('GPU normalized optical kernels',capacity*32,U.STORAGE),tiles=this.buffer('Resolved photosphere tiles',capacity*16,U.STORAGE|U.COPY_DST);
  for(const buffer of [this.stars,this.prepared,this.tiles])this.retire(buffer);this.stars=stars;this.prepared=prepared;this.tiles=tiles;this.capacity=capacity;
 }
 allocateAtlas(pixels){if(pixels<=this.atlasCapacity)return;this.retire(this.raw);this.retire(this.blur);const bytes=Math.max(1024,2**Math.ceil(Math.log2(pixels))*4);if(bytes>this.device.limits.maxStorageBufferBindingSize)throw Error('Photosphères trop nombreuses pour le tampon optique GPU.');this.raw=this.buffer('GPU photosphere profiles',bytes,GPUBufferUsage.STORAGE);this.blur=this.buffer('GPU convolved photospheres',bytes,GPUBufferUsage.STORAGE);this.atlasCapacity=bytes/4;}
 resize(width,height=width,force=false){width=Math.max(1,Math.round(width));height=Math.max(1,Math.round(height));if(width===this.width&&height===this.height&&!force)return;this.width=width;this.height=height;this.size=Math.max(width,height);this.retire(this.hdr);this.hdr=this.texture('Native rectangular linear starlight',[width,height],'rgba16float',GPUTextureUsage.RENDER_ATTACHMENT|GPUTextureUsage.TEXTURE_BINDING|GPUTextureUsage.COPY_SRC);this.toneGroup=this.device.createBindGroup({layout:this.pipelines.tone.getBindGroupLayout(0),entries:[{binding:0,resource:this.hdr.createView()},{binding:1,resource:{buffer:this.config}},{binding:2,resource:{buffer:this.presentConfig}}]});}
 glareTexture(sigma,profile,angle,signal=1,shape=root.DEFAULT_OPTICAL_SHAPE){
  const support=root.opticalKernelSupport(sigma,signal,shape),key=[sigma.toFixed(3),profile,angle.toFixed(3),support.radius,support.tailSigma,...Object.values(shape)].join(':');if(this.glareCache.has(key))return this.glareCache.get(key);
  const radius=support.radius,side=2*radius+1,count=side*side,groups=Math.ceil(count/256),U=GPUBufferUsage,d=this.device,texture=this.texture('Instrument PSF GPU · '+key,[side,side],'r32float',GPUTextureUsage.TEXTURE_BINDING|GPUTextureUsage.STORAGE_BINDING),config=this.buffer('Instrument kernel settings',48,U.UNIFORM|U.COPY_DST),raw=this.buffer('Instrument kernel samples',count*8,U.STORAGE),sums=this.buffer('Instrument kernel reduction',groups*8,U.STORAGE),norm=this.buffer('Instrument kernel normalizer',16,U.STORAGE),tileSide=Math.ceil((side+2)/16),bounds=this.buffer('Conservative optical tile maxima',tileSide*tileSide*4,U.STORAGE);
  const bytes=new ArrayBuffer(48);new Uint32Array(bytes).set([side,radius,groups,support.tailSigma]);new Float32Array(bytes).set([shape.haloExtent,shape.spikeLength,shape.spikeMin,shape.spikeMax],8);new Float32Array(bytes).set([sigma,angle,profile==='soft'?0:.42,profile==='six'?1:profile==='soft'?2:0],4);d.queue.writeBuffer(config,0,bytes);const group=d.createBindGroup({layout:this.glareLayout,entries:[config,raw,sums,norm].map((buffer,binding)=>({binding,resource:{buffer}})).concat({binding:4,resource:texture.createView()},{binding:5,resource:{buffer:bounds}})}),encoder=d.createCommandEncoder({label:'Instrument PSF entirely on GPU'});
  for(const name of ['generate','reduce','normalize','bounds']){const pass=encoder.beginComputePass();pass.setPipeline(this.glarePipelines[name]);pass.setBindGroup(0,group);pass.dispatchWorkgroups(name==='reduce'?1:name==='bounds'?tileSide*tileSide:groups);pass.end();}d.queue.submit([encoder.finish()]);const release=()=>{for(const b of [config,raw,sums,norm])this.retire(b);};d.queue.onSubmittedWorkDone().then(release,release);
  if(this.glareCache.size>=8){const first=this.glareCache.keys().next().value;this.retire(this.glareCache.get(first).texture);this.retire(this.glareCache.get(first).bounds);for(const level of this.glareCache.get(first).levels.values())this.retire(level.bounds);this.glareCache.delete(first);}const result={texture,radius,bounds,tileSide,tileSize:16,levels:new Map()};this.glareCache.set(key,result);return result;
 }
 adaptiveGlareBounds(kernel,count){if(!kernel.tileSide)return kernel;let tileSize=16,side=kernel.tileSide;while(count*side*side*4>this.tileBudgetBytes&&side>1){tileSize*=2;side=Math.ceil((2*kernel.radius+3)/tileSize);}if(tileSize===16)return kernel;if(kernel.levels.has(tileSize))return kernel.levels.get(tileSize);const d=this.device,bounds=this.buffer('Coarse conservative optical bounds',side*side*4,GPUBufferUsage.STORAGE),config=this.buffer('Optical bounds pyramid parameters',16,GPUBufferUsage.UNIFORM|GPUBufferUsage.COPY_DST);d.queue.writeBuffer(config,0,new Uint32Array([kernel.tileSide,side,tileSize/16,0]));const group=d.createBindGroup({layout:this.coarseBoundsPipeline.getBindGroupLayout(0),entries:[kernel.bounds,bounds,config].map((buffer,binding)=>({binding,resource:{buffer}}))}),encoder=d.createCommandEncoder(),pass=encoder.beginComputePass();pass.setPipeline(this.coarseBoundsPipeline);pass.setBindGroup(0,group);pass.dispatchWorkgroups(Math.ceil(side*side/256));pass.end();d.queue.submit([encoder.finish()]);const release=()=>this.retire(config);d.queue.onSubmittedWorkDone().then(release,release);const result={...kernel,bounds,tileSide:side,tileSize};kernel.levels.set(tileSize,result);return result;}
 render(canvas,projection,points,{population,gain=1,exposureEV=0,sigma=1,complete=false,selected=null,lengthPc=10000,normalizationLengthPc=lengthPc,fovDegrees=50,glare=.65,opticalProfile='starburst',spikeAngle=15,detailWeight=1,opticalShape=root.DEFAULT_OPTICAL_SHAPE,stellarFade=root.DEFAULT_STELLAR_FADE,context=null,diffuseLayer=null,minimumPeak=.000125}={}){
  if(context)sigma=context.sigmaPixels;const selectionScale=context?.pixelRatio??1;
  const size=this.size,width=this.width,height=this.height,d=this.device,exposure=2**exposureEV,result=root.lightPlan(projection,points,{population,size,sigma,complete,gain,exposureEV,lengthPc,normalizationLengthPc,fovDegrees,detailWeight,stellarFade,context,minimumPeak});this.ensureSourceCapacity(result.sources.length);

  let strongest=0,strongestGlare=0;for(const p of result.sources)if(!p.photometry.disc&&p.glareOpacity>0){const signal=p.linearFlux*exposure*Math.max(p.photometry.color[0],p.photometry.color[1],p.photometry.color[2]);strongest=Math.max(strongest,signal);strongestGlare=Math.max(strongestGlare,signal*Math.max(0,Math.min(.8,glare*(p.photometry.glareScale??1)))*.55*p.glareOpacity);}
  const baseKernel=strongestGlare>=.9999e-6&&glare>0?this.glareTexture(sigma,opticalProfile,spikeAngle*Math.PI/180,strongest,opticalShape):{texture:this.emptyGlare,radius:0};const kernel=this.adaptiveGlareBounds(baseKernel,result.sources.length);if(this.sourceScratch.length<Math.max(20,result.sources.length*20))this.sourceScratch=new Float32Array(Math.max(20,this.capacity*20));const packed=this.sourceScratch.subarray(0,Math.max(20,result.sources.length*20)),tileValues=[];let pixels=0,fastPacked=0;
  for(const [i,p]of result.sources.entries()){
   if(p.photometry.disc){const px=width/2+(p.x-result.geometry.centerX)/projection.width*size,py=height/2-(p.y-result.geometry.centerY)/projection.width*size,visibility=p.visibility??1,R=Math.max(p.radiusPixels,1e-6),mode=p.photometry.lightCos===undefined?-1:2+p.photometry.lightCos;packed.set([px,py,R,p.linearFlux*exposure*visibility,...p.photometry.color,R+2,0,0,0,0,0,0,0,0,0,0,p.photometry.lightAngle||0,mode],i*20);continue;}
if(this.gpuPacking&&p.radiusPixels<=Math.max(1,1.1774100225*sigma)){const o=i*20,color=p.photometry.color;packed[o]=width/2+(p.x-result.geometry.centerX)/projection.width*size;packed[o+1]=height/2-(p.y-result.geometry.centerY)/projection.width*size;packed[o+2]=p.radiusPixels;packed[o+3]=p.linearFlux*exposure;packed[o+4]=color[0];packed[o+5]=color[1];packed[o+6]=color[2];packed[o+8]=sigma;packed[o+9]=-1;packed[o+10]=Math.max(0,Math.min(.8,glare*(p.photometry.glareScale??1)));packed[o+11]=p.glareOpacity;fastPacked++;continue;}
const px=width/2+(p.x-result.geometry.centerX)/projection.width*size,py=height/2-(p.y-result.geometry.centerY)/projection.width*size,R=p.radiusPixels,lower=Math.max(1,1.1774100225*sigma),upper=Math.max(3,3*sigma),q=Math.max(0,Math.min(1,(R-lower)/(upper-lower))),mix=q*q*(3-2*q),psf=sigma,u=Math.max(0,Math.min(1,(R-lower)/Math.max(6*sigma-lower,1e-30))),unresolved=1/(1+(Math.max(0,R-lower)/(2*sigma))**4),fraction=Math.max(0,Math.min(.8,glare*(p.photometry.glareScale??1)))*.55*unresolved*(1-u*u*(3-2*u))*p.glareOpacity,signal=p.linearFlux*exposure*Math.max(...p.photometry.color),limit=-1,pad=root.gaussianSupport(sigma,signal),x0=Math.max(-pad,Math.floor(px-R-pad)),x1=Math.min(width+pad,Math.ceil(px+R+pad)),y0=Math.max(-pad,Math.floor(py-R-pad)),y1=Math.min(height+pad,Math.ceil(py+R+pad)),w=Math.max(0,x1-x0),h=Math.max(0,y1-y0),extent=Math.ceil(Math.max(mix<1?root.gaussianSupport(psf,signal):0,R+pad,fraction>0?(limit<0?Math.min(kernel.radius,root.opticalKernelSupport(sigma,signal,opticalShape).radius):Math.min(kernel.radius,limit)):0))+2;
   packed.set([px,py,R,p.linearFlux*exposure,...p.photometry.color,extent,sigma,psf,mix,fraction,limit,pixels,w,h,x0,y0,Math.ceil(Math.max(mix<1?pad:0,R+pad))+2,p.photometry.lightCos===undefined?0:2+p.photometry.lightCos],i*20);if(mix>0&&w*h>0){tileValues.push(pixels,w,h,i);pixels+=w*h;}}
  this.allocateAtlas(Math.max(1,pixels));d.queue.writeBuffer(this.stars,0,packed);if(tileValues.length)d.queue.writeBuffer(this.tiles,0,new Uint32Array(tileValues));
  const rebased=normalizationLengthPc!==lengthPc,densityScale=rebased?gain/projection.pixel**2*(normalizationLengthPc/lengthPc)**2:1;const n=projection.resolution;if(this.densitySize!==n){this.retire(this.density);this.densitySize=n;this.density=this.texture('Small statistical light field',[n,n],'r32float',GPUTextureUsage.TEXTURE_BINDING|GPUTextureUsage.COPY_DST);}d.queue.writeTexture({texture:this.density},Float32Array.from(result.diffuseImage,value=>value*densityScale),{bytesPerRow:n*4},[n,n]);
  const tileCount=kernel.tileSide?kernel.tileSide**2:0,tileSlots=result.sources.length*tileCount,sparse=this.sparseGlare&&tileSlots>0&&tileSlots*4<=Math.min(d.limits.maxStorageBufferBindingSize,d.limits.maxBufferSize);
  if(sparse&&tileSlots>this.glareTileCapacity){this.retire(this.visibleGlare);this.glareTileCapacity=Math.min(Math.floor(Math.min(d.limits.maxStorageBufferBindingSize,d.limits.maxBufferSize)/4),2**Math.ceil(Math.log2(tileSlots)));this.visibleGlare=this.buffer('Compacted visible optical tiles',this.glareTileCapacity*4,GPUBufferUsage.STORAGE);}
  const config=new ArrayBuffer(80),floats=new Float32Array(config),ints=new Uint32Array(config);floats.set([size,n,result.geometry.focalPixels,0,(rebased?1:gain/projection.pixel**2)*(context?.diffuseScale??1),result.geometry.fieldFactor,exposure,0,...root.StellarPopulation.meanColor,Number(sparse)]);floats.set([width,height,kernel.tileSize||16,selectionScale],16);const target=selected?result.candidates.find(p=>p.id===selected):null;if(target){floats[3]=Math.min(size*2,Math.max(7*selectionScale,target.radiusPixels+4*selectionScale));floats[7]=width/2+(target.x-result.geometry.centerX)/projection.width*size;floats[11]=height/2-(target.y-result.geometry.centerY)/projection.width*size;}ints.set([result.sources.length,tileValues.length/4,pixels,Number(!!target)|(Number(sparse)<<1)],12);d.queue.writeBuffer(this.config,0,config);
  const resources=[this.config,this.stars,this.prepared,this.raw,this.blur,this.tiles],group=d.createBindGroup({layout:this.layout,entries:[...resources.map((buffer,binding)=>({binding,resource:{buffer}})),{binding:6,resource:kernel.texture.createView()},{binding:7,resource:this.density.createView()},{binding:8,resource:{buffer:this.visibleGlare}}]}),encoder=d.createCommandEncoder({label:'Stellar optical image on GPU'});
  const compute=(name,count,workgroups=false)=>{if(!count)return;const pass=encoder.beginComputePass();pass.setPipeline(this.pipelines[name]);pass.setBindGroup(0,group);const groups=workgroups?count:Math.ceil(count/256);pass.dispatchWorkgroups(Math.min(groups,65535),Math.ceil(groups/65535));pass.end();};if(fastPacked){d.queue.writeBuffer(this.packingConfig,0,new Uint32Array([result.sources.length,kernel.radius,0,0]));const packing=encoder.beginComputePass();packing.setPipeline(this.packingPipeline);packing.setBindGroup(0,d.createBindGroup({layout:this.packingPipeline.getBindGroupLayout(0),entries:[{binding:0,resource:{buffer:this.stars}},{binding:1,resource:{buffer:this.packingConfig}}]}));const groups=Math.ceil(result.sources.length/256);packing.dispatchWorkgroups(Math.min(groups,65535),Math.ceil(groups/65535));packing.end();}compute('prepare',result.sources.length,true);compute('photospheres',pixels);compute('horizontalBlur',pixels);
  if(sparse){d.queue.writeBuffer(this.glareVisibilityConfig,0,new Uint32Array([result.sources.length,kernel.tileSide,width,kernel.radius+1,height,kernel.tileSize,0,0]));encoder.clearBuffer(this.glareDraw);const visibilityGroup=d.createBindGroup({layout:this.visibilityPipeline.getBindGroupLayout(0),entries:[this.stars,kernel.bounds,this.visibleGlare,this.glareDraw,this.glareVisibilityConfig].map((buffer,binding)=>({binding,resource:{buffer}}))}),visibility=encoder.beginComputePass();visibility.setPipeline(this.visibilityPipeline);visibility.setBindGroup(0,visibilityGroup);const groups=Math.ceil(tileSlots/256);visibility.dispatchWorkgroups(Math.min(groups,65535),Math.ceil(groups/65535));visibility.end();}

  if(diffuseLayer)diffuseLayer.encode(encoder,{target:this.hdr,projection,context,result,gain,exposureEV});
  const pass=encoder.beginRenderPass({colorAttachments:[{view:this.hdr.createView(),clearValue:{r:0,g:0,b:0,a:1},loadOp:diffuseLayer?'load':'clear',storeOp:'store'}]});pass.setBindGroup(0,group);if(!diffuseLayer&&result.residual>0){pass.setPipeline(this.pipelines.diffuse);pass.draw(3);}if(result.sources.length){pass.setPipeline(this.pipelines.stars);pass.draw(6,result.sources.length);}if(sparse){pass.setPipeline(this.pipelines.glare);pass.drawIndirect(this.glareDraw,0);}pass.end();
  d.queue.submit([encoder.finish()]);if(canvas)this.present(canvas);result.work={gpuPreparedSources:fastPacked,sparseGlare:sparse,glareTileSize:kernel.tileSize||0,glareTileSlots:sparse?tileSlots:0,glareTileBytes:this.glareTileCapacity*4};result.backend='GPU';result.atlasPixels=pixels;result.reflectedCount=result.sources.filter(p=>p.photometry.disc&&p.photometry.lightCos!==undefined).length;result.discCount=result.sources.filter(p=>p.photometry.disc).length;if(result.discCount)result.drawn=result.sources.filter(p=>p.photometry.disc||p.peak>=.008);return result;
 }
 present(canvas,{scale=1,x=0,y=0,guides=[],width=this.width,height=this.height,context=null}={}){
  const d=this.device,size=this.size;
  if(this.canvas!==canvas){this.context?.unconfigure();this.canvas=canvas;this.context=canvas.getContext('webgpu');if(!this.context)throw Error('Canvas WebGPU indisponible.');this.context.configure({device:d,format:this.format,alphaMode:'opaque'});}width=Math.max(1,Math.round(width));height=Math.max(1,Math.round(height));if(canvas.width!==width||canvas.height!==height){canvas.width=width;canvas.height=height;}
  if(guides.length){if(!this.guideBuffer||guides.length>this.guideCapacity){this.retire(this.guideBuffer);this.guideCapacity=2**Math.ceil(Math.log2(Math.max(64,guides.length)));this.guideUpload=new Float32Array(this.guideCapacity*4);this.guideBuffer=this.buffer('Minor-body guide instances',this.guideUpload.byteLength,GPUBufferUsage.STORAGE|GPUBufferUsage.COPY_DST);this.guideGroup=d.createBindGroup({layout:this.pipelines.guides.getBindGroupLayout(0),entries:[{binding:0,resource:{buffer:this.guideBuffer}},{binding:1,resource:{buffer:this.guideViewport}}]});}for(let i=0;i<guides.length;i++){const g=guides[i];const o=i*4;this.guideUpload[o]=g.x;this.guideUpload[o+1]=g.y;this.guideUpload[o+2]=g.radius;this.guideUpload[o+3]=g.opacity;}d.queue.writeBuffer(this.guideBuffer,0,this.guideUpload.buffer,0,guides.length*16);}
  d.queue.writeBuffer(this.guideViewport,0,new Float32Array([context?.cssWidth??600*width/size,context?.cssHeight??600*height/size,0,0]));d.queue.writeBuffer(this.presentConfig,0,new Float32Array([scale,x,y,size,width,height,0,0]));const encoder=d.createCommandEncoder({label:'GPU camera and sensor presentation'}),pass=encoder.beginRenderPass({colorAttachments:[{view:this.context.getCurrentTexture().createView(),clearValue:{r:0,g:0,b:0,a:1},loadOp:'clear',storeOp:'store'}]});pass.setPipeline(this.pipelines.tone);pass.setBindGroup(0,this.toneGroup);pass.draw(3);if(guides.length){pass.setPipeline(this.pipelines.guides);pass.setBindGroup(0,this.guideGroup);pass.draw(6,guides.length);}pass.end();d.queue.submit([encoder.finish()]);
 }

 retire(resource){if(!resource)return;resource.destroy();this.buffers.delete(resource);this.textures.delete(resource);}
 destroy(){this.context?.unconfigure();for(const b of this.buffers)b.destroy();for(const t of this.textures)t.destroy();this.buffers.clear();this.textures.clear();this.glareCache.clear();}
}
root.GPUStellarLightRenderer=GPUStellarLightRenderer;
})(typeof window!=='undefined'?window:globalThis);
