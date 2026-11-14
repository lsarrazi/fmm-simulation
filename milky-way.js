/* Fiducial MW-like planar model, not a reconstruction of the actual Milky Way.
 * Length = 10 kpc, time = 1 Myr, stellar mass = count × population mean mass.
 * Exponential stellar disc + barred bulge + spiral perturbations + smooth halo.
 * The perturbations define an initial condition, not a prescribed animation.
 */
(function(root){
'use strict';
const population=typeof module!=='undefined'?require('./stellar-population'):root.StellarPopulation;
const MW={stars:100000000000n,lengthKpc:10,timeMyr:1,massMsun:100000000000*population.meanMass,meanStarMass:population.meanMass,radiusKpc:15.33,scaleKpc:2.6,heightPc:300,solarKpc:8.2,solarSpeedKmS:230,velocityKmS:3.0856775814913673e17/3.15576e13};
MW.gravity=4.30091e-6*MW.massMsun/(MW.lengthKpc*MW.velocityKmS**2);
MW.surfaceSolar=Number(MW.stars)*.85/(2*Math.PI*(MW.scaleKpc*1000)**2*(1-Math.exp(-MW.radiusKpc/MW.scaleKpc)*(1+MW.radiusKpc/MW.scaleKpc)))*Math.exp(-MW.solarKpc/MW.scaleKpc);
MW.numberDensitySolar=MW.surfaceSolar/(2*MW.heightPc);
function surfaceComponents(r){const rd=.26,rm=1.533,a=.045,rb=.2;return[r<rm?.85*Math.exp(-r/rd)/(2*Math.PI*rd*rd*(1-Math.exp(-rm/rd)*(1+rm/rd))):0,r<rb?.15/(Math.PI*a*a*(rb*rb/(rb*rb+a*a)))/(1+r*r/(a*a))**2:0];}
function circularSpeed(r,p){const dr=1.7/96,t=Math.max(0,Math.min(95.999,r/dr)),i=Math.floor(t),u=t-i,force=p.radialForces[i]*(1-u)+p.radialForces[i+1]*u;return Math.sqrt(Math.max(0,r*force+p.haloV2*r*r/(r*r+.04)));}
// One shared spatial/velocity prescription for Monte Carlo, CPU, GPU and identities.
// Arms have zero azimuthal mean; the area-preserving bar keeps the bulge mass.
function profileMW(x,y,s){
 const r=Math.hypot(x,y),theta=Math.atan2(y,x),taper=(1-Math.exp(-((r/.18)**2)))*Math.exp(-((r/1.65)**4));
 const phase=2*(theta-.45-Math.log(Math.max(r,.035)/.35)/Math.tan(.32));
 const contrast=1+s*taper*(.8*Math.cos(phase)+.16*Math.cos(theta+.8+.7*Math.log(Math.max(r,.035)/.35)));
 const stretch=1+.85*s,c=Math.cos(.45),a=Math.sin(.45),barRadius=Math.hypot((c*x+a*y)/stretch,(-a*x+c*y)*stretch);
 const disc=surfaceComponents(r)[0]*contrast,bulge=surfaceComponents(barRadius)[1];
 return{disc,bulge,taper,phase,r};
}
function stateMW(x,y,p){const s=p.structure??.85,{disc,bulge,taper,phase,r}=profileMW(x,y,s);
 return{disc,bulge,vr:s*taper*22*Math.sin(phase)/MW.velocityKmS,vt:circularSpeed(r,p)+s*taper*14*Math.cos(phase)/MW.velocityKmS};
}
function cohortsForMW(n,extent,p){const side=Math.min(n,64),dx=2*extent/side,cohorts=[];let mass=0,key=0;
 for(let y=0;y<side;y++)for(let x=0;x<side;x++){const px=(x+.5)*dx-extent,py=(y+.5)*dx-extent,r=Math.hypot(px,py),ca=r?px/r:1,sa=r?py/r:0,{disc,bulge,vr:meanR,vt:meanT}=stateMW(px,py,p);
  for(const[type,density]of [[0,disc],[1,bulge]])if(density>0){const sr=(type?90:30)/MW.velocityKmS,st=(type?90:20)/MW.velocityKmS;
   for(const[vr,vp,w]of [[0,0,.5],[2*sr,0,.125],[-2*sr,0,.125],[0,2*st,.125],[0,-2*st,.125]]){const tangent=vp+(type?0:meanT),radial=vr+(type?0:meanR),m=density*dx*dx*w;cohorts.push({key:key++,x:px,y:py,vx:ca*radial-sa*tangent,vy:sa*radial+ca*tangent,width:dx,mass:m});mass+=m;}
  }
 }for(const c of cohorts)c.mass/=mass;return cohorts;
}
function populateMilkyWay(grid,p){grid.f.fill(0);grid.time=grid.spatialLoss=grid.velocityLoss=0;const nv=grid.nv,discPDF=new Float64Array(nv*nv),bulgePDF=new Float64Array(nv*nv);let total=0;
 for(let y=0;y<grid.n;y++)for(let x=0;x<grid.n;x++){const px=(x+.5)*grid.dx-grid.extent,py=(y+.5)*grid.dx-grid.extent,r=Math.hypot(px,py),ca=r?px/r:1,sa=r?py/r:0,{disc,bulge,vr:meanR,vt:meanT}=stateMW(px,py,p);if(disc+bulge===0)continue;let sd=0,sb=0;
  for(let v=0;v<nv;v++)for(let u=0;u<nv;u++){const vx=(u-(nv-1)/2)*grid.dv,vy=(v-(nv-1)/2)*grid.dv,vr=ca*vx+sa*vy-meanR,vp=-sa*vx+ca*vy-meanT,k=v*nv+u;sd+=discPDF[k]=Math.exp(-.5*((vr*MW.velocityKmS/30)**2+(vp*MW.velocityKmS/20)**2));sb+=bulgePDF[k]=Math.exp(-.5*(vx*vx+vy*vy)*(MW.velocityKmS/90)**2);}
  const cell=y*grid.n+x;for(let k=0;k<nv*nv;k++){const value=grid.dx*grid.dx*(disc*discPDF[k]/sd+bulge*bulgePDF[k]/sb);grid.f[cell*nv*nv+k]=value;total+=value;}
 }for(let i=0;i<grid.f.length;i++)grid.f[i]/=total;grid.project();grid.initialCohorts=cohortsForMW(grid.n,grid.extent,p);
}
function createMilkyWay(count=2048,seed=417,{softening=.015,structure=.85}={}){
 structure=Math.max(0,Math.min(1,structure));
 let state=seed>>>0;const random=()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return(state+.5)/4294967296;},normal=()=>Math.sqrt(-2*Math.log(random()))*Math.cos(2*Math.PI*random());
 const bodies=[];
 for(let i=0;i<count;i++){
  const disc=i>=Math.floor(count*.15);let r;
  if(disc){do{r=-MW.scaleKpc*Math.log(random()*random());}while(r>MW.radiusKpc);}
  else {do{const u=random();r=.45*Math.sqrt(u/(1-u));}while(r>2);}
  let angle=random()*2*Math.PI,x,y;
  if(disc){
   // Rejection from the exponential disk gives precisely the analytic arm PDF.
   while(true){x=Math.cos(angle)*r/10;y=Math.sin(angle)*r/10;const shape=profileMW(x,y,structure);
    if(random()<(shape.disc/surfaceComponents(r/10)[0])/(1+.96*structure))break;
    do{r=-MW.scaleKpc*Math.log(random()*random());}while(r>MW.radiusKpc);angle=random()*2*Math.PI;
   }
  }else{const u=Math.cos(angle)*r/10*(1+.85*structure),v=Math.sin(angle)*r/10/(1+.85*structure);x=Math.cos(.45)*u-Math.sin(.45)*v;y=Math.sin(.45)*u+Math.cos(.45)*v;}
  bodies.push({x,y,vx:0,vy:0,mass:1/count,disc});
 }
 // Azimuthal quadrature of the SAME softened stellar mass used by both solvers.
 const radial=new Float64Array(97),dr=1.7/96;
 const quadrature=[];let qm=0;for(let j=0;j<256;j++){const r=(j+.5)*1.533/256,[disc,bulge]=surfaceComponents(r),mass=2*Math.PI*r*(1.533/256)*(disc+bulge)/128;for(let a=0;a<128;a++){quadrature.push({x:r*Math.cos((a+.5)*2*Math.PI/128),y:r*Math.sin((a+.5)*2*Math.PI/128),mass});qm+=mass;}}
 for(let k=1;k<radial.length;k++){const r=k*dr;let inward=0;for(const b of quadrature){const dx=b.x-r,dy=b.y,q=dx*dx+dy*dy+softening*softening;inward-=MW.gravity*b.mass/qm*dx/(q*Math.sqrt(q));}radial[k]=Math.max(0,inward);}
 const stellarForce=r=>{const t=Math.max(0,Math.min(95.999,r/dr)),i=Math.floor(t),u=t-i;return radial[i]*(1-u)+radial[i+1]*u;};
 const rs=MW.solarKpc/10,haloV2=Math.max(0,(MW.solarSpeedKmS/MW.velocityKmS)**2-rs*stellarForce(rs))*(rs*rs+.04)/(rs*rs);
 const physics={gravity:MW.gravity,softening,haloV2,haloCore:.2,extent:2,vmax:.06,analytic:true,radialForces:radial,structure};
 let speedSum=0,discSum=0,discCount=0;
 for(const b of bodies){const r=Math.hypot(b.x,b.y),ca=b.x/Math.max(r,1e-12),sa=b.y/Math.max(r,1e-12);
  if(b.disc){const flow=stateMW(b.x,b.y,physics),vr=flow.vr+normal()*30/MW.velocityKmS,vphi=flow.vt+normal()*20/MW.velocityKmS;b.vx=ca*vr-sa*vphi;b.vy=sa*vr+ca*vphi;discSum+=Math.hypot(b.vx,b.vy)*MW.velocityKmS;discCount++;}
  else{b.vx=normal()*90/MW.velocityKmS;b.vy=normal()*90/MW.velocityKmS;}
  speedSum+=Math.hypot(b.vx,b.vy)*MW.velocityKmS;
 }
 bodies.physics=physics;
 bodies.calibration={...MW,meanSpeedKmS:speedSum/count,meanDiscSpeedKmS:discSum/discCount,solarCircularKmS:Math.sqrt(rs*stellarForce(rs)+haloV2*rs*rs/(rs*rs+.04))*MW.velocityKmS};
 return bodies;
}
const api={MilkyWay:MW,createMilkyWay,surfaceMW:surfaceComponents,stateMW,cohortsForMW,populateMilkyWay};Object.assign(root,api);if(typeof module!=='undefined')module.exports=api;
})(typeof window!=='undefined'?window:globalThis);
