/* Perspective camera above a planar stellar distribution. Physical inverse-square
 * flux and photosphere radii; Gaussian instrument PSF; fixed exposure response.
 * No dust, depth distribution, retarded time, or full spectral instrument model.
 */
(function(root){
'use strict';
const Population=typeof module!=='undefined'?require('./stellar-population'):root.StellarPopulation,{hash}=Population;
// Frozen optical birth properties need no re-evaluation between frames. Keep
// only a bounded local LRU; positions, velocities and forces are never cached here.
const stellarPropertyCache=new Map();let stellarPropertyCacheLimit=32768,propertyHead=null,propertyTail=null;
function reserveStellarProperties(count){stellarPropertyCacheLimit=Math.max(stellarPropertyCacheLimit,Math.min(1048576,Math.max(32768,Math.ceil(count)*3)));return stellarPropertyCacheLimit;}
const stellarProperties=id=>{
 let entry=stellarPropertyCache.get(id);
 if(entry){if(entry!==propertyTail){if(entry.prev)entry.prev.next=entry.next;else propertyHead=entry.next;entry.next.prev=entry.prev;entry.prev=propertyTail;entry.next=null;propertyTail.next=entry;propertyTail=entry;}return entry.value;}
 const value=Population.properties(id);Object.freeze(value.color);Object.freeze(value);
 if(stellarPropertyCache.size>=stellarPropertyCacheLimit){const old=propertyHead;propertyHead=old.next;if(propertyHead)propertyHead.prev=null;else propertyTail=null;stellarPropertyCache.delete(old.id);}
 entry={id,value,prev:propertyTail,next:null};if(propertyTail)propertyTail.next=entry;else propertyHead=entry;propertyTail=entry;stellarPropertyCache.set(id,entry);return value;
};
// Fictional designation derived only from the immutable ID, with no name registry.
function stellarName(id){
 const first=['Ael','Bel','Cyr','Dra','El','Fen','Gal','Hel','Iri','Jan','Kel','Lyr','Mor','Ner','Ori','Sol'];
 const middle=['a','ae','an','ar','e','en','eri','i','il','in','o','or','ul','ur','va','y'];
 const last=['lis','na','ria','ron','ris','th','tis','via','n','mir','dor','lia','nus','ra','sel','wen'];
 const pick=(parts,salt)=>parts[Math.floor(hash(id,salt)*parts.length)],code=Math.floor(hash(id,419)*4294967296).toString(36).toUpperCase().padStart(7,'0');
 return pick(first,401)+pick(middle,409)+pick(last,417)+' · '+code;
}
function occupancy(projection,population,size,sigma,context=null){const N=Number(population),beam=4*Math.PI*sigma*sigma,area=(size/projection.resolution)**2;let minimum=Infinity;for(const mass of projection.image)if(mass>0)minimum=Math.min(minimum,N*mass/area*beam);return {minimum:minimum===Infinity?0:minimum,mean:projection.expected/(context?context.pixelWidth*context.pixelHeight:size*size*(projection.aspectX??1)*(projection.aspectY??1))*beam,beam};}
// Probability that a resolution element contains at least two sources. This is
// a statistical representation weight, not absorption of physical stellar light.
function diffuseCrowdingWeight(lambda){if(!(lambda>0))return 0;if(lambda>50)return 1;if(lambda<.001)return lambda*lambda*(.5-lambda/3+lambda*lambda/8);return -Math.expm1(-lambda)-lambda*Math.exp(-lambda);}
function photosphereResolution(radius,sigma=1){const limit=1.1774100225*sigma;return{resolved:radius>=limit,limitRadiusPixels:limit};}
function cameraGeometry(projection,{size=600,lengthPc=10000,fovDegrees=50,context=null}={}){
 const tangent=Math.tan(fovDegrees*Math.PI/360),heightPc=(projection.cameraWidth??projection.width)*lengthPc/(2*tangent);
 return{heightPc,focalPixels:context?.focalPixels??size*(projection.cameraWidth??projection.width)/projection.width/(2*tangent),lengthPc,fieldFactor:(tangent/Math.tan(25*Math.PI/180))**2,centerX:projection.centerX??projection.left+projection.width/2,centerY:projection.centerY??projection.bottom+projection.width/2};
}
// Visibility/ranking only needs the inverse-square visible flux, not a full PSF record.
function sourceVisibleFlux(point,geometry,photometry=point.photometry||stellarProperties(point.id)){
 const x=(point.x-geometry.centerX)*geometry.lengthPc,y=(point.y-geometry.centerY)*geometry.lengthPc,R=photometry.radius*Population.solarRadiusPc*1.000001,distance2=Math.max(x*x+y*y+geometry.heightPc*geometry.heightPc,R*R);return Population.solarFluxAtPc*photometry.visibleLuminosity/distance2;
}
function sourceOptics(point,geometry){const photometry=point.photometry||stellarProperties(point.id),R=photometry.radius*Population.solarRadiusPc;
 const x=(point.x-geometry.centerX)*geometry.lengthPc,y=(point.y-geometry.centerY)*geometry.lengthPc,distance2=x*x+y*y+geometry.heightPc*geometry.heightPc,rawDistance=Number.isFinite(distance2)?Math.sqrt(distance2):Math.hypot(x,y,geometry.heightPc),distancePc=Math.max(rawDistance,R*1.000001);
 const fluxWm2=Population.solarFluxAtPc*photometry.luminosity/(distancePc*distancePc),visibleFluxWm2=Population.solarFluxAtPc*photometry.visibleLuminosity/(distancePc*distancePc);
 return{photometry,distancePc,fluxWm2,visibleFluxWm2,radiusPixels:geometry.focalPixels*R/Math.max(geometry.heightPc,1e-30),insidePhotosphere:rawDistance<=R};
}
const DEFAULT_STELLAR_FADE=Object.freeze({enabled:true,start:1000000,end:1000,curve:10});
function readStellarFade(values={}){const result={...DEFAULT_STELLAR_FADE,...values};for(const key of ['start','end','curve']){result[key]=Number(String(result[key]).trim().replace(',','.'));if(!Number.isFinite(result[key])||!(result[key]>0)){const e=Error('Distances et courbe strictement positives attendues.');e.field='fade'+key[0].toUpperCase()+key.slice(1);throw e;}}if(result.start<=result.end){const e=Error('Le début doit être plus éloigné que la fin de l’atténuation.');e.field='fadeStart';throw e;}result.enabled=!!result.enabled;return result;}
function stellarGlareOpacity(source,fade=DEFAULT_STELLAR_FADE,logRange=Math.log(fade.start/fade.end)){if(!fade.enabled||source.photometry.disc||!(source.photometry.radius>0))return 1;const distance=source.distancePc/(source.photometry.radius*Population.solarRadiusPc);if(distance>=fade.start)return 1;if(distance<=fade.end)return 0;const t=Math.log(distance/fade.end)/logRange;return (t*t*(3-2*t))**fade.curve;}
function minimumLuminosity(projection,{population,size=600,sigma=1,gain=1,exposureEV=0,fovDegrees=50,context=null,minimumPeak=.000125}={}){
 if(context){size=context.exposurePixels;sigma=context.sigmaPixels;}
 const factor=(Math.tan(fovDegrees*Math.PI/360)/Math.tan(25*Math.PI/180))**2;
 return minimumPeak*2*Math.PI*sigma*sigma*Number(population)*Population.meanVisibleLuminosity*(projection.cameraWidth??projection.width)**2/(gain*size*size*2**exposureEV*factor);
}
// Local surface LOD is based on projected physical size and sensor signal,
// never mass labels, zoom alone, or previously-visible source retention.
function surfaceVisibility(radiusPixels,peak,{pixelRatio=1}={}){
 const smooth=(x,lo,hi)=>{const t=Math.max(0,Math.min(1,Math.log2(Math.max(x,1e-30)/lo)/Math.log2(hi/lo)));return t*t*(3-2*t);};
 const angular=smooth(radiusPixels/pixelRatio,.35,1.25),signal=smooth(peak,.00003125,.0005);return angular*signal;
}
function plan(projection,points,{population,size=600,sigma=1,complete=false,gain=1,exposureEV=0,lengthPc=10000,normalizationLengthPc=lengthPc,fovDegrees=50,detailWeight=1,stellarFade=DEFAULT_STELLAR_FADE,context=null,minimumPeak=.000125}={}){
 if(context){size=context.pixelSpan;sigma=context.sigmaPixels;}
 reserveStellarProperties(points.length);
 detailWeight=Math.max(0,Math.min(1,detailWeight));
 const N=Number(population),geometry=cameraGeometry(projection,{size,lengthPc,fovDegrees,context}),fluxGain=gain*(context?.exposurePixels??size)**2*normalizationLengthPc*normalizationLengthPc/(N*Population.meanVisibleLuminosity*Population.solarFluxAtPc*(2*Math.tan(25*Math.PI/180))**2);
 const exposure=2**exposureEV,psfArea=2*Math.PI*sigma*sigma,fadeLogRange=Math.log(stellarFade.start/stellarFade.end),resolvedRadius=1.1774100225*sigma;
 const sources=[],drawn=[],candidates=[];let luminosity=0,leastVisible=Infinity,stellarSourceCount=0;
 for(const p of points){if(Math.abs(p.x-geometry.centerX)>projection.width*(context?.aspectX??projection.aspectX??1)/2||Math.abs(p.y-geometry.centerY)>projection.width*(context?.aspectY??projection.aspectY??1)/2)continue;const ix=Math.floor((p.x-projection.left)/projection.width*projection.resolution),iy=Math.floor((p.y-projection.bottom)/projection.width*projection.resolution);if(ix<0||iy<0||ix>=projection.resolution||iy>=projection.resolution)continue;
  const optics=sourceOptics(p,geometry),linearFlux=optics.visibleFluxWm2*fluxGain*(optics.photometry.isSystemBody?1:detailWeight),peak=linearFlux*exposure/(psfArea+Math.PI*optics.radiusPixels**2*.8),source={...p,photometry:optics.photometry,distancePc:optics.distancePc,fluxWm2:optics.fluxWm2,visibleFluxWm2:optics.visibleFluxWm2,radiusPixels:optics.radiusPixels,insidePhotosphere:optics.insidePhotosphere,linearFlux,peak,photosphereResolved:optics.radiusPixels>=resolvedRadius,glareOpacity:stellarGlareOpacity(optics,stellarFade,fadeLogRange)};source.visibility=p.photometry?.proximityVisibility?surfaceVisibility(optics.radiusPixels,peak,{pixelRatio:context?.pixelRatio??1}):1;candidates.push(source);
  if(source.visibility>0&&(complete||peak>=minimumPeak||p.photometry?.proximityVisibility)){sources.push(source);if(!source.photometry.isSystemBody){stellarSourceCount++;leastVisible=Math.min(leastVisible,source.photometry.visibleLuminosity);luminosity+=source.photometry.visibleLuminosity;}}
  if(peak*source.visibility>=.008)drawn.push(source);
 }
 const diffuseImage=new Float64Array(projection.image.length),beam=4*Math.PI*sigma*sigma,pixelArea=(size/projection.resolution)**2,mass=projection.mass??projection.image.reduce((a,b)=>a+b,0),observedMass=projection.viewMass??mass,unrepresented=Math.max(0,1-stellarSourceCount/Math.max(N*observedMass,1e-300));let remaining=0,coverage=0;
 const tail=Population.luminosityTailMoments(leastVisible),represented=Math.min(1,stellarSourceCount/Math.max(N*observedMass*tail.probability,1e-300)),meanRemainder=Math.max(0,1-detailWeight*represented*tail.luminosity/Population.meanVisibleLuminosity);
 for(let i=0;i<diffuseImage.length;i++){const cellMass=projection.image[i],weight=Math.max(1-detailWeight,diffuseCrowdingWeight(N*cellMass*unrepresented*beam/pixelArea));coverage+=cellMass*weight;diffuseImage[i]=complete?(1-detailWeight)*cellMass:weight*cellMass*meanRemainder;remaining+=diffuseImage[i];}
 return{drawn,sources,candidates,diffuseImage,diffuseParameters:{unrepresented,meanRemainder,detailWeight,beam},residual:remaining/Math.max(mass,1e-300),diffuseCoverage:coverage/Math.max(mass,1e-300),estimatedPopulation:N*observedMass,luminosity,complete,geometry,fluxGain,occupancy:occupancy(projection,population,size,sigma,context)};
}
// Tail support follows the exposed image signal, never the world-coordinate
// magnitude. Bounds are rounded into cacheable levels; omitted light stays below
// the linear sensor floor even for severely overexposed nearby sources.
function gaussianSupport(sigma,signal=1){return Math.ceil(sigma*Math.max(8,Math.sqrt(2*Math.max(0,Math.log(Math.max(1e-30,signal)/1e-8)))))+1;}
const DEFAULT_OPTICAL_SHAPE=Object.freeze({haloExtent:4,spikeLength:9,spikeMin:.35,spikeMax:1.2});
function readOpticalShape(values={}){const result={...DEFAULT_OPTICAL_SHAPE};for(const key of Object.keys(result)){if(values[key]!==undefined){const n=Number(String(values[key]).trim().replace(',','.'));if(!Number.isFinite(n)||!(n>0)){const e=Error('Valeurs optiques strictement positives attendues.');e.field=key;throw e;}result[key]=n;}}if(result.spikeMin>result.spikeMax){const e=Error('L’épaisseur au bout doit être inférieure ou égale à celle au centre.');e.field='spikeMin';throw e;}return result;}
function opticalKernelSupport(sigma,signal=1,shape=DEFAULT_OPTICAL_SHAPE){const tail=Math.max(4,Math.ceil(Math.sqrt(2*Math.max(0,Math.log(Math.max(1e-30,signal)*64/1e-8))))),reach=Math.max(shape.haloExtent,shape.spikeLength,shape.spikeMax),level=Math.max(2,Math.ceil((reach*tail+4)/32));return{radius:Math.ceil(level*32*sigma)+3,tailSigma:Math.ceil(level*32/reach)};}

function addPSF(buffer,size,px,py,flux,rgb,sigma,signal=flux*Math.max(...rgb)){
 const radius=gaussianSupport(sigma,signal),cx=Math.floor(px),cy=Math.floor(py),wx=[],wy=[];let sx=0,sy=0;
 for(let k=-radius;k<=radius;k++){const x=Math.exp(-.5*((cx+k+.5-px)/sigma)**2),y=Math.exp(-.5*((cy+k+.5-py)/sigma)**2);wx.push(x);wy.push(y);sx+=x;sy+=y;}
 for(let y=-radius;y<=radius;y++)for(let x=-radius;x<=radius;x++){const ix=cx+x,iy=cy+y;if(ix<0||iy<0||ix>=size||iy>=(buffer.height??size))continue;const weight=flux*wx[x+radius]*wy[y+radius]/(sx*sy),out=(iy*size+ix)*3;for(let c=0;c<3;c++)buffer[out+c]+=weight*rgb[c];}
}
// Instrument-only cache: no observed identity or visibility is retained here.
const glareKernels=new Map();
function glareKernel(sigma,profile,angle,signal=1,shape=DEFAULT_OPTICAL_SHAPE){const support=opticalKernelSupport(sigma,signal,shape),key=[sigma.toFixed(3),profile,angle.toFixed(3),support.radius,support.tailSigma,...Object.values(shape)].join(':');if(glareKernels.has(key))return glareKernels.get(key);
 const map=new Map(),put=(x,y,w)=>{const ix=Math.floor(x),iy=Math.floor(y),tx=x-ix,ty=y-iy;for(let a=0;a<2;a++)for(let b=0;b<2;b++){const value=w*(a?tx:1-tx)*(b?ty:1-ty);if(!value)continue;const k=(ix+a)+':'+(iy+b);map.set(k,(map.get(k)||0)+value);}},halo=[];let hs=0;
 const radius=support.radius;for(let y=-radius;y<=radius;y++)for(let x=-radius;x<=radius;x++){const r2=(x*x+y*y)/(sigma*sigma*(shape.haloExtent/12)**2),w=(1+r2/4)**(-1.7)*Math.exp(-r2/288);if(w){halo.push([x,y,w]);hs+=w;}}
 const spikeFraction=profile==='soft'?0:.42;for(const [x,y,w]of halo)put(x,y,w/hs*(1-spikeFraction));
 if(spikeFraction){const axes=profile==='six'?3:6,rays=[];let total=0;
  for(let a=0;a<axes;a++){const theta=angle+a*Math.PI/axes,c=Math.cos(theta),s=Math.sin(theta),long=(profile==='starburst'&&a%2?7/12:1)*shape.spikeLength*sigma,amplitude=profile==='starburst'&&a%2?.5:1;
   for(let t=-radius;t<=radius;t++){const width=Math.max(shape.spikeMin*sigma,shape.spikeMax*sigma*(1-Math.abs(t)/(shape.spikeLength*56/12*sigma)));for(let q=-Math.ceil(support.tailSigma*width);q<=Math.ceil(support.tailSigma*width);q++){const w=amplitude*Math.exp(-.5*(t/long)**2)/Math.sqrt(1+(t/(3*sigma))**2)*Math.exp(-.5*(q/width)**2)/width;rays.push([t*c-q*s,t*s+q*c,w]);total+=w;}}
  }
  for(const [x,y,w]of rays)put(x,y,w/total*spikeFraction);
 }
 const result=[];for(const [k,w]of map){const [x,y]=k.split(':').map(Number);result.push([x,y,w]);}result.sort((a,b)=>b[2]-a[2]);if(glareKernels.size>=8)glareKernels.delete(glareKernels.keys().next().value);glareKernels.set(key,result);return result;
}
function addOpticalGlare(buffer,size,px,py,flux,rgb,sigma,{profile='starburst',angle=0,exposure=1,cutoff=0,radiusLimit=Infinity,shape=DEFAULT_OPTICAL_SHAPE}={}){
 const x0=Math.floor(px),y0=Math.floor(py),tx=px-x0,ty=py-y0,kernel=glareKernel(sigma,profile,angle,flux*exposure*Math.max(...rgb),shape),taper=(x,y)=>{if(!Number.isFinite(radiusLimit))return 1;const u=Math.max(0,Math.min(1,(radiusLimit-Math.hypot(x,y))/Math.max(sigma,.2*radiusLimit)));return u*u*(3-2*u);};let norm=1;if(Number.isFinite(radiusLimit)){norm=0;for(const [x,y,w]of kernel)norm+=w*taper(x,y);}if(!(norm>0))return;const peak=Math.max(...rgb)*flux*exposure/norm;
 for(const [dx,dy,w]of kernel){if(w*peak<cutoff)break;const weight=w*taper(dx,dy)/norm;if(!weight)continue;for(let a=0;a<2;a++)for(let b=0;b<2;b++){const x=x0+dx+a,y=y0+dy+b;if(x<0||y<0||x>=size||y>=(buffer.height??size))continue;const value=flux*weight*(a?tx:1-tx)*(b?ty:1-ty),out=(y*size+x)*3;for(let c=0;c<3;c++)buffer[out+c]+=value*rgb[c];}}
}
function addStellarLight(buffer,size,px,py,flux,rgb,R,sigma,{glare=.65,profile='starburst',angle=0,exposure=1,cutoff=0,shape=DEFAULT_OPTICAL_SHAPE,glareOpacity=1}={}){
 const lower=Math.max(1,photosphereResolution(R,sigma).limitRadiusPixels),u=Math.max(0,Math.min(1,(R-lower)/Math.max(6*sigma-lower,1e-30))),unresolved=1/(1+(Math.max(0,R-lower)/(2*sigma))**4),fraction=Math.max(0,Math.min(.8,glare))*.55*unresolved*(1-u*u*(3-2*u))*glareOpacity,radiusLimit=Infinity;
 addPhotosphere(buffer,size,px,py,flux*(1-fraction),rgb,R,sigma,flux*exposure*Math.max(...rgb));
 if(fraction>0)addOpticalGlare(buffer,size,px,py,flux*fraction,rgb,sigma,{profile,angle,exposure,cutoff,radiusLimit,shape});
}
// Physical limb-darkened photosphere, convolved with a normalized instrument PSF.
// Unresolved stars share the instrumental footprint; only resolved disks grow.
function addPhotosphere(buffer,size,px,py,flux,rgb,R,sigma,signal=flux*Math.max(...rgb)){
 // Unresolved stars have the instrument's footprint. Blend its fixed profile
 // to a physical limb-darkened disk continuously as its diameter becomes resolved.
 const lower=Math.max(1,photosphereResolution(R,sigma).limitRadiusPixels),upper=Math.max(3,3*sigma),q=Math.max(0,Math.min(1,(R-lower)/(upper-lower))),mix=q*q*(3-2*q);
 if(mix<1)addPSF(buffer,size,px,py,flux*(1-mix),rgb,sigma,signal);
 if(mix===0)return;addResolvedPhotosphere(buffer,size,px,py,flux*mix,rgb,R,sigma,signal);
}
function addResolvedPhotosphere(buffer,size,px,py,flux,rgb,R,sigma,signal=flux*Math.max(...rgb)){
 const pad=gaussianSupport(sigma,signal),x0=Math.max(-pad,Math.floor(px-R-pad)),x1=Math.min(size+pad,Math.ceil(px+R+pad)),y0=Math.max(-pad,Math.floor(py-R-pad)),y1=Math.min((buffer.height??size)+pad,Math.ceil(py+R+pad)),w=x1-x0,h=y1-y0;if(w<=0||h<=0)return;
 const raw=new Float32Array(w*h),temp=new Float32Array(w*h),kernel=[];let norm=0;for(let k=-pad;k<=pad;k++){const v=Math.exp(-.5*(k/sigma)**2);kernel.push(v);norm+=v;}for(let i=0;i<kernel.length;i++)kernel[i]/=norm;
 const profile=(x,y)=>{let value=0;for(const a of [.25,.75])for(const b of [.25,.75]){const q=((x+a-px)**2+(y+b-py)**2)/(R*R);if(q<1)value+=.4+.6*Math.sqrt(1-q);}return value*.25;};
 let integral=Math.PI*R*R*.8;if(R<=32){integral=0;for(let y=Math.floor(py-R);y<=Math.ceil(py+R);y++)for(let x=Math.floor(px-R);x<=Math.ceil(px+R);x++)integral+=profile(x,y);}const factor=flux/Math.max(integral,1e-30);
 for(let y=0;y<h;y++)for(let x=0;x<w;x++)raw[y*w+x]=profile(x+x0,y+y0)*factor;
 for(let y=0;y<h;y++)for(let x=0;x<w;x++){let value=0;for(let k=-pad;k<=pad;k++)if(x+k>=0&&x+k<w)value+=raw[y*w+x+k]*kernel[k+pad];temp[y*w+x]=value;}
 for(let y=Math.max(0,-y0);y<Math.min(h,(buffer.height??size)-y0);y++)for(let x=Math.max(0,-x0);x<Math.min(w,size-x0);x++){let value=0;for(let k=-pad;k<=pad;k++)if(y+k>=0&&y+k<h)value+=temp[(y+k)*w+x]*kernel[k+pad];const out=((y+y0)*size+x+x0)*3;for(let c=0;c<3;c++)buffer[out+c]+=value*rgb[c];}
}
class StellarLightRenderer{
 constructor(size=600){this.resize(size);}
 resize(width,height=width){if(width!==this.width||height!==this.height){this.width=width;this.height=height;this.size=Math.max(width,height);this.linear=new Float32Array(width*height*3);this.linear.height=height;}}
 render(canvas,projection,points,{population,gain=1,exposureEV=0,sigma=1,complete=false,selected=null,lengthPc=10000,normalizationLengthPc=lengthPc,fovDegrees=50,glare=.65,opticalProfile='starburst',spikeAngle=15,detailWeight=1,opticalShape=DEFAULT_OPTICAL_SHAPE,stellarFade=DEFAULT_STELLAR_FADE,context=null,minimumPeak=.000125}={}){
  const size=this.size,width=this.width,height=this.height,selectionScale=context?.pixelRatio??1,ctx=canvas.getContext('2d'),buffer=this.linear;if(context)sigma=context.sigmaPixels;buffer.fill(0);if(canvas.width!==width||canvas.height!==height){canvas.width=width;canvas.height=height;}
  const result=plan(projection,points,{population,size,sigma,complete,gain,exposureEV,lengthPc,normalizationLengthPc,fovDegrees,detailWeight,stellarFade,context,minimumPeak}),exposure=2**exposureEV,n=projection.resolution,image=ctx.createImageData(width,height);
  const diffuseGain=gain/projection.pixel**2*(normalizationLengthPc/lengthPc)**2*(context?.diffuseScale??1);
  // Bilinear interpolation of the mean photon field, not an opaque density palette.
  if(result.residual>0){const indices=new Int32Array(width),fractions=new Float64Array(width),radii=new Float64Array(width),field=result.diffuseImage,color=Population.meanColor;
   for(let x=0;x<width;x++){const u=(x+.5-width/2+size/2)/size*n-.5;indices[x]=Math.floor(u);fractions[x]=u-indices[x];radii[x]=((x+.5-width/2)/result.geometry.focalPixels)**2;}
   for(let y=0;y<height;y++){const v=(size/2-y-.5+height/2)/size*n-.5,iy=Math.floor(v),ty=v-iy,row=iy*n,next=row+n,dy2=((y+.5-height/2)/result.geometry.focalPixels)**2;
    for(let x=0;x<width;x++){const ix=indices[x],tx=fractions[x];let low=0,high=0;if(iy>=0&&iy<n){if(ix>=0)low+=field[row+ix]*(1-tx);if(ix+1<n)low+=field[row+ix+1]*tx;}if(iy+1>=0&&iy+1<n){if(ix>=0)high+=field[next+ix]*(1-tx);if(ix+1<n)high+=field[next+ix+1]*tx;}
     const radiance=(low*(1-ty)+high*ty)*diffuseGain*result.geometry.fieldFactor/(1+radii[x]+dy2),out=(y*width+x)*3;buffer[out]=radiance*color[0];buffer[out+1]=radiance*color[1];buffer[out+2]=radiance*color[2];
    }
   }
  }
  for(const p of result.sources)addStellarLight(buffer,width,width/2+(p.x-result.geometry.centerX)/projection.width*size,height/2-(p.y-result.geometry.centerY)/projection.width*size,p.linearFlux,p.photometry.color,p.radiusPixels,sigma,{glare:glare*(p.photometry.glareScale??1),profile:opticalProfile,angle:spikeAngle*Math.PI/180,exposure,cutoff:1e-6,shape:opticalShape,glareOpacity:p.glareOpacity});
  // Fixed exposure, shared nonlinear camera response; never normalize to each frame.
  for(let i=0;i<width*height;i++){const o=i*3,R=buffer[o]*exposure,G=buffer[o+1]*exposure,B=buffer[o+2]*exposure;image.data[i*4]=255*Math.sqrt(R/(1+R));image.data[i*4+1]=255*Math.sqrt(G/(1+G));image.data[i*4+2]=255*Math.sqrt(B/(1+B));image.data[i*4+3]=255;}
  ctx.putImageData(image,0,0);
  if(selected){const p=result.candidates.find(p=>p.id===selected);if(p){ctx.strokeStyle='#ff5050';ctx.lineWidth=selectionScale;ctx.beginPath();ctx.arc(width/2+(p.x-result.geometry.centerX)/projection.width*size,height/2-(p.y-result.geometry.centerY)/projection.width*size,Math.min(size*2,Math.max(7*selectionScale,p.radiusPixels+4*selectionScale)),0,Math.PI*2);ctx.stroke();}}
  return result;
 }
}
// Positive cubic B-spline reconstruction of the collective density. Unlike
// stretched affine polygons, this field cannot draw spurious ballistic ribbons.
function projectDensityLight(grid,{x=0,y=0,zoom=1,resolution=192}={}){
 const width=2*grid.extent/zoom,pixel=width/resolution,left=x-width/2,bottom=y-width/2,image=new Float64Array(resolution*resolution),n=grid.n;
 const weights=t=>[(1-t)**3/6,(3*t**3-6*t*t+4)/6,(-3*t**3+3*t*t+3*t+1)/6,t**3/6];
 const xs=Array.from({length:resolution},(_,i)=>{const u=(left+(i+.5)*pixel+grid.extent)/grid.dx-.5;return{index:Math.floor(u)-1,w:weights(u-Math.floor(u))};});
 for(let y=0;y<resolution;y++){const v=(bottom+(y+.5)*pixel+grid.extent)/grid.dx-.5,iy=Math.floor(v)-1,wy=weights(v-Math.floor(v));for(let x=0;x<resolution;x++){const p=xs[x];let mass=0;for(let a=0;a<4;a++)for(let b=0;b<4;b++){const xx=p.index+a,yy=iy+b;if(xx>=0&&yy>=0&&xx<n&&yy<n)mass+=grid.mass[yy*n+xx]*p.w[a]*wy[b];}image[y*resolution+x]=mass*pixel*pixel/(grid.dx*grid.dx);}}
 let mass=0;for(const m of image)mass+=m;return{image,resolution,width,pixel,left,bottom,centerX:x,centerY:y,expected:Number(grid.population||0)*mass,mass};
}
const api={DEFAULT_STELLAR_FADE,readStellarFade,stellarGlareOpacity,DEFAULT_OPTICAL_SHAPE,readOpticalShape,opticalGlareKernel:glareKernel,gaussianSupport,opticalKernelSupport,StellarLightRenderer,stellarProperties,reserveStellarProperties,stellarName,lightPlan:plan,opticalOccupancy:occupancy,addPSF,addPhotosphere,addOpticalGlare,addStellarLight,projectDensityLight,cameraGeometry,sourceOptics,sourceVisibleFlux,surfaceVisibility,minimumLuminosity,diffuseCrowdingWeight,photosphereResolution};Object.assign(root,api);if(typeof module!=='undefined')module.exports=api;
})(typeof window!=='undefined'?window:globalThis);
