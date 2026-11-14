/* A deterministic coarse birth model, not a time integration of accretion.
 * One mass field is partitioned into compact condensates and residual matter.
 * The threshold is an effective surface-density stability proxy, not an IAU
 * definition or a resolved gas/streaming-instability calculation.
 */
(function(root){
'use strict';
const P=typeof module!=='undefined'?require('./stellar-population'):root.StellarPopulation;
const G=39.47841760435743,MSUN=1.98847e30,AU_KM=149597870.7;
const rand=(id,s)=>P.hash(id,s),normal=(id,s)=>Math.sqrt(-2*Math.log(rand(id,s)))*Math.cos(2*Math.PI*rand(id,s+1));
function build(starId,host,{densityScale=1,stellarClearing=true}={}){
 if(!(densityScale>0&&Number.isFinite(densityScale)))throw Error('Densité initiale strictement positive attendue.');
 const id=starId+'/primordial/v1',nr=32,nt=16,dtheta=2*Math.PI/nt;
 const inner=Math.max(4*host.radius*P.solarRadiusM/1000/AU_KM,.035*Math.cbrt(host.mass));
 const outer=Math.max(inner*30,10**(.65+1.1*rand(id,1))*Math.sqrt(host.mass)),dl=Math.log(outer/inner)/nr;
 const initialMass=host.mass*10**(-4.3+1.4*rand(id,2))*densityScale;
 const dispersion=2e-5+6e-5*rand(id,3),lumps=Array.from({length:4+Math.floor(rand(id,4)*5)},(_,i)=>({r:rand(id,10+i*7)*nr,t:rand(id,11+i*7)*nt,wr:1.1+2*rand(id,12+i*7),wt:.7+1.4*rand(id,13+i*7),amp:3+12*rand(id,14+i*7)}));
 const cells=[];let weightedArea=0;
 for(let i=0;i<nr;i++)for(let j=0;j<nt;j++){
  const r0=inner*Math.exp(i*dl),r1=r0*Math.exp(dl),r=Math.sqrt(r0*r1),theta=(j+.5)*dtheta,area=.5*(r1*r1-r0*r0)*dtheta;
  let weight=.05*(r/inner)**.15;
  for(const p of lumps){const dt=Math.min(Math.abs(j+.5-p.t),nt-Math.abs(j+.5-p.t));weight+=p.amp*Math.exp(-.5*((i+.5-p.r)/p.wr)**2-.5*(dt/p.wt)**2);}
  // Effective cold-disk support sigma_r and stellar shear Omega. Threshold
  // mass = area * sigma_r*Omega/(pi*G); finite zones define its resolution.
  const omega=Math.sqrt(G*host.mass/r**3),sigma=dispersion*Math.sqrt(G*host.mass/r);
  cells.push({i,j,r0,r1,r,theta,area,surfaceWeight:weight,capacity:area*sigma*omega/(Math.PI*G)});weightedArea+=weight*area;
 }
 // The prior is a surface density, not a probability per logarithmic bin.
 // Integrate it over each sector area (the polar Jacobian) before normalising.
 for(const c of cells){c.initialMass=initialMass*c.surfaceWeight*c.area/weightedArea;c.condensed=c.initialMass>c.capacity?c.initialMass:0;c.residual=c.initialMass-c.condensed;}
 // Connected overdense zones become one compact body. No independent planet
 // count, giant-planet draw or second mass budget exists.
 const visited=new Set(),compact=[];
 for(let index=0;index<cells.length;index++){
  if(visited.has(index)||!(cells[index].condensed>0))continue;
  const queue=[index];visited.add(index);
  for(let q=0;q<queue.length;q++){const k=queue[q],c=cells[k];
   for(const [di,dj]of [[1,0],[-1,0],[0,1],[0,-1]]){const ni=c.i+di,nj=(c.j+dj+nt)%nt,nk=ni*nt+nj;if(ni>=0&&ni<nr&&!visited.has(nk)&&cells[nk].condensed>0){visited.add(nk);queue.push(nk);}}
  }
  compact.push(...fragmentComponent(queue,cells,host.mass));
 }
 compact.sort((a,b)=>a.a-b.a);
 const residualMass=cells.reduce((s,c)=>s+c.residual,0),compactMass=compact.reduce((s,c)=>s+c.mass,0),cdf=new Float64Array(cells.length+1);
 for(let i=0;i<cells.length;i++)cdf[i+1]=cdf[i]+(residualMass>0?cells[i].residual/residualMass:0);cdf[cells.length]=residualMass>0?1:0;
 // Normalise the size spectrum by the actual remaining mass. The inferred
 // number can be enormous; per-identity masses use a stable integral interval.
 const b=2.5,rmin=.01,rmax=400,c0=rmax**(-b),c1=rmin**(-b),q=1-3/b;
 const kgPerKm3=4*Math.PI/3*2.4e12,meanKg=kgPerKm3*(c1**q-c0**q)/(q*(c1-c0));
 const count=residualMass>0?BigInt(Math.max(1,Math.ceil(residualMass*MSUN/meanKg))):0n;
 const field={version:6,nr,nt,inner,outer,dispersion,initialMass,residualMass,compactMass,cells,compact,cdf,count,spectrum:{b,c0,c1,q}};
 Object.assign(field,radialDistribution(field,host,field.compact,{stellarClearing}));return field;
}
// A connected overdensity spanning distant orbits cannot accrete into a single
// coherent object. Partition its radial feeding domain using a Hill-scale proxy,
// rather than drawing an independent planet count or changing the mass budget.
function fragmentComponent(zones,cells,hostMass){
 const rows=new Map();for(const index of zones){const c=cells[index];let row=rows.get(c.i);if(!row){row={r:c.r,mass:0,zones:[]};rows.set(c.i,row);}row.mass+=c.condensed;row.zones.push(index);}
 const result=[];let group=[],mass=0,inner=0;
 const flush=()=>{if(!group.length)return;let angular=0,x=0,y=0;for(const k of group){const c=cells[k];angular+=c.condensed*Math.sqrt(c.r);x+=c.condensed*Math.cos(c.theta);y+=c.condensed*Math.sin(c.theta);}result.push({mass,a:(angular/mass)**2,e:0,omega:0,mean:Math.atan2(y,x),zones:group});group=[];mass=0;};
 for(const row of [...rows.values()].sort((a,b)=>a.r-b.r)){const coherence=10*Math.cbrt((mass+row.mass)/(3*hostMass));if(mass&&Math.log(row.r/inner)>coherence)flush();if(!mass)inner=row.r;mass+=row.mass;group.push(...row.zones.sort((a,b)=>a-b));}flush();return result;
}
// The primordial sectors determine condensation, not the angular appearance of
// the mature residual disk. Phase mixing marginalises them before observation.
// Smooth radial density is then conditioned only on planetary orbital radii and
// masses. No planet phase or sector boundary survives this operation.
// Effective post-formation feeding domain of the central star. This is a birth
// prior conditioned on the neighbouring condensates, not radiative removal of
// kilometre-scale bodies or a resolved integration of stellar accretion.
function stellarClearingRadius(field,host,planets){
 const first=Math.min(...planets.filter(p=>p.mass>0&&p.a>0).map(p=>p.a*(1-(p.e||0)))),photosphere=(host.radius||0)*P.solarRadiusM/1000/AU_KM;
 const thermal=.0344*Math.sqrt(Math.max(0,host.luminosity||0));
 return Math.max(6*photosphere,thermal,Number.isFinite(first)?.65*first:.15*field.outer);
}
function radialDistribution(field,host,planets=field.compact,{stellarClearing=true}={}){
 const n=512,smoothing=Math.max(.12,1.25*Math.log(field.outer/field.inner)/field.nr),samplingOuter=field.outer*Math.exp(8*smoothing),step=Math.log(samplingOuter/field.inner)/n,rings=[];
 const primordial=Array.from({length:field.nr},()=>({mass:0,logR:0}));
 for(const c of field.cells){const p=primordial[c.i];p.mass+=c.residual;p.logR=Math.log(c.r);}
 const centralClearingRadius=stellarClearing?stellarClearingRadius(field,host,planets):0;
 const logInner=Math.log(field.inner),normalisation=step/(Math.sqrt(2*Math.PI)*smoothing);
 for(let i=0;i<n;i++){
  const logR=logInner+(i+.5)*step,r=Math.exp(logR),r0=Math.exp(logInner+i*step),r1=Math.exp(logInner+(i+1)*step);let prior=0;
  for(const p of primordial){const kernel=Math.exp(-.5*((logR-p.logR)/smoothing)**2)*normalisation;prior+=p.mass*kernel;}
  const stellarSurvival=Math.exp(-.5*(centralClearingRadius/r)**4);let survival=stellarSurvival;
  for(const p of planets){if(!(p.mass>0&&p.a>0))continue;const hill=Math.cbrt(p.mass/(3*host.mass)),width=Math.max(.035,3*hill+(p.e||0)),distance=(logR-Math.log(p.a))/width,depletion=.985*p.mass/(p.mass+host.mass*1e-7);survival*=1-depletion*Math.exp(-.5*distance*distance);
   const strength=p.mass/(p.mass+host.mass*.0001);
   for(const ratio of [1/3,2/5,1/2]){const resonance=Math.log(p.a*ratio**(2/3)),d=(logR-resonance)/(.012+.015*strength);survival*=1-.75*strength*Math.exp(-.5*d*d);}
  }
  // Cold-disk cap evaluated at the outer edge is conservative over this ring.
  const area=Math.PI*(r1*r1-r0*r0),omega=Math.sqrt(G*host.mass/r1**3),sigma=field.dispersion*Math.sqrt(G*host.mass/r1);
  rings.push({r0,r1,r,area,capacity:area*sigma*omega/(Math.PI*G),weight:Math.max(prior*survival*(-Math.expm1(-.5*((logR-logInner)/smoothing)**2)),field.residualMass*1e-30),survival,stellarSurvival,residual:0});
 }
 // Redistribution changes probability, never the remaining mass budget. Water
 // filling normalises the profile without pushing any ring above its cap.
 const capacity=rings.reduce((s,c)=>s+c.capacity,0);if(field.residualMass>capacity*(1+1e-10))throw Error('Capacité radiale insuffisante pour le reliquat.');
 const total=lambda=>rings.reduce((s,c)=>s+Math.min(c.capacity,lambda*c.weight),0);let lo=0,hi=1;while(total(hi)<field.residualMass&&hi<1e30)hi*=2;
 for(let i=0;i<80;i++){const mid=(lo+hi)/2;if(total(mid)<field.residualMass)lo=mid;else hi=mid;}
 const radialCDF=new Float64Array(n+1);for(let i=0;i<n;i++){rings[i].residual=Math.min(rings[i].capacity,hi*rings[i].weight);radialCDF[i+1]=radialCDF[i]+(field.residualMass>0?rings[i].residual/field.residualMass:0);}radialCDF[n]=field.residualMass>0?1:0;
 return{radialCells:rings,radialCDF,samplingOuter,centralClearingRadius};
}
function massBetween(field,start,end){const {c0,c1,q}=field.spectrum,N=Number(field.count),t=c0+Number(start)/N*(c1-c0),delta=Number(end-start)/N*(c1-c0);return field.residualMass*t**q*Math.expm1(q*Math.log1p(delta/t))/(c1**q-c0**q);}
function residualBody(field,id,rank){
 const mass=massBetween(field,rank,rank+1n);
 const draw=rand(id,220);let lo=0,hi=field.radialCells.length;while(lo+1<hi){const mid=(lo+hi)>>1;if(field.radialCDF[mid]<=draw)lo=mid;else hi=mid;}
 const c=field.radialCells[lo],r=Math.sqrt(c.r0*c.r0+rand(id,221)*(c.r1*c.r1-c.r0*c.r0)),theta=rand(id,222)*2*Math.PI;
 // Residual bodies have a modest velocity dispersion after condensation.
 // Convert that exact initial position/velocity into bound orbital elements.
 const radial=.025*normal(id,223),tangent=1+.015*normal(id,225),energy=(radial*radial+tangent*tangent)/2-1,a=-r/(2*energy);
 const ex=tangent*tangent-1,ey=-radial*tangent,e=Math.hypot(ex,ey),localOmega=Math.atan2(ey,ex),omega=theta+localOmega;
 const trueAnomaly=-localOmega,E=Math.atan2(Math.sqrt(1-e*e)*Math.sin(trueAnomaly),e+Math.cos(trueAnomaly));
 return{mass,radiusKm:Math.cbrt(mass*MSUN/(4*Math.PI/3*2.4e12)),a,e,omega,mean:E-e*Math.sin(E),birthRing:lo};
}
const api={build,residualBody,massBetween,radialDistribution,stellarClearingRadius};root.FormationModel=api;if(typeof module!=='undefined')module.exports=api;
})(typeof window!=='undefined'?window:globalThis);
