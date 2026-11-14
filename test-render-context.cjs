const assert=require('node:assert/strict'),O=require('./observation-tools'),L=require('./stellar-light');
const context=(w,h,dpr=1)=>O.renderContext(O.viewportGeometry(w,h,{pixelRatio:dpr}),{worldShortSpan:1,sigmaCss:1,lengthPc:1});
let firstPeak,firstThreshold;
for(const [w,h,dpr]of [[600,600,1],[1200,600,1],[600,1200,1],[1120,754,1],[1120,754,2]]){
 const c=context(w,h,dpr),p={width:c.worldSpan,cameraWidth:c.worldShortSpan,left:-c.worldSpan/2,bottom:-c.worldSpan/2,centerX:0,centerY:0,resolution:16,image:new Float64Array(256),mass:0,expected:0,aspectX:c.aspectX,aspectY:c.aspectY};
 assert(Object.isFrozen(c));assert.equal(c.sigmaPixels/c.pixelRatio,1);assert(c.worldWidth>=1&&c.worldHeight>=1);assert(Math.abs(Math.min(c.worldWidth,c.worldHeight)-1)<1e-12);
 const point={id:'context-star',x:0,y:0,photometry:{radius:0,luminosity:1,visibleLuminosity:1,color:[1,1,1]}},options={population:100n,gain:.01,lengthPc:1,context:c},r=L.lightPlan(p,[point],options),threshold=L.minimumLuminosity(p,options);firstPeak??=r.sources[0].peak;firstThreshold??=threshold;assert(Math.abs(r.sources[0].peak/firstPeak-1)<1e-12,'same optical peak at fixed camera/exposure, independently of viewport/DPR');assert(Math.abs(threshold/firstThreshold-1)<1e-12);
 assert.equal(O.systemBodyGuide(1,c.worldSpan,{context:c,selected:true}).ringRadius,6);assert.equal(L.surfaceVisibility(.5*dpr,1,{pixelRatio:c.pixelRatio}),L.surfaceVisibility(.5,1));assert(Math.abs(c.observationTolerance*c.cssShort-1)<1e-12);
}
const c=context(1200,600),image={originX:0,originY:0,modelLengthPc:1,unitPc:1,space:'galaxy',x:0,y:0,width:c.worldSpan,cameraWidth:1,zoom:1,follow:null};
for(const requested of [{...image,x:.01},{...image,width:image.width*2},{...image,y:-.4}]){const selected=O.coveredCameraView(image,requested,c,c);assert.equal(selected.pending,true);assert.equal(selected.view,image,'uncomputed pan/dezoom must keep the complete existing frame');}
const zoomed={...image,width:image.width/2,cameraWidth:.5,zoom:2};assert.equal(O.coveredCameraView(image,zoomed,c,c).pending,false,'zoom-in can reuse a covered region');
const portrait=context(357,518),same=O.coveredCameraView(image,image,c,c);assert.equal(same.pending,false);assert.equal(portrait.pixelSpan,518);assert.equal(Math.max(portrait.pixelWidth,portrait.pixelHeight),518);
console.log('Passed: immutable viewport/world/optical context, stable CSS PSF and exposure, full rectangular FOV, fixed CSS guides/LOD, explicit sampling error, and coverage-only frame publication.');
