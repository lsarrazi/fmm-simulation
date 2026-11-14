/* Exactly two levels: galactic star -> bounded, star-centred planetary patch.
 * Birth properties are pure ID functions. Checkpoints retain computed encounters;
 * inactive time uses osculating Kepler propagation, not unseen N-body integration.
 */
(function(root){
'use strict';
const P=typeof module!=='undefined'?require('./stellar-population'):root.StellarPopulation;
const Formation=typeof module!=='undefined'?require('./formation-model'):root.FormationModel;
const C={AU_KM:149597870.7,AU_PC:1/206264.80624709636,YEAR_NS:31557600000000000n,G:39.47841760435743,EARTH_MSUN:3.0034896e-6,EARTH_RADIUS_KM:6371,FAMILIAR_NAME_EARTH_MASS:.001,AU_YEAR_KMS:4.7404704635};
const random=(id,salt)=>P.hash(id,salt);
function describe(starId,host=P.properties(starId)){
 const prefix=starId+'/system/v2',formation=Formation.build(starId,host),seeds=formation.compact.map((p,index)=>{const id=prefix+'/body/'+index,earthMass=p.mass/C.EARTH_MSUN,radiusEarth=earthMass<2?earthMass**.28:earthMass<100?1.25*(earthMass/2)**.5:11.2*(earthMass/318)**(-.04);return{...p,id,index,kind:classifyMass(p.mass),name:objectName(id,p.mass),radiusKm:radiusEarth*C.EARTH_RADIUS_KM,albedo:.1+.45*random(id,26),color:earthMass>=95?[.8,.68,.48]:earthMass>2?[.36,.62,.85]:[.6+.3*random(id,27),.43+.3*random(id,28),.3+.35*random(id,29)]};});
 return{starId,prefix,host,seeds,formation,smallCount:formation.count+BigInt(seeds.length),defaultCount:Math.min(Number(formation.count)+seeds.length,Math.round(6500+8500*random(prefix,41))),extent:Math.max(1,formation.outer*1.3,...seeds.map(p=>p.a)),rotation:0};
}
// Display categories are derived only from mass; birth mixtures do not encode a type.
function classifyMass(mass){const earth=mass/C.EARTH_MSUN;return earth<.01?'Astéroïde':earth<2?'Planète rocheuse':earth<95?'Planète intermédiaire':'Géante gazeuse';}
function objectName(id,mass=Infinity){const first=['Ael','Cor','Vel','Ner','Eri','Tal','Ily','Or','Ser','Lun','Mer','Syl','Ar','Cael','Dor','Vyr'],last=['ora','ion','ea','alis','une','ara','en','is','ela','os','ia','elle','an','eth','ys','on'];const parts=id.split('/'),rank=BigInt(parts.at(-1)),host=parts.slice(0,-4),tag=(host.at(-2)||'0').toUpperCase()+'-'+(/^[0-9]+$/.test(host.at(-1)||'')?BigInt(host.at(-1)).toString(36):String(host.at(-1)||'0')).toUpperCase();if(mass<C.EARTH_MSUN*C.FAMILIAR_NAME_EARTH_MASS)return 'OBJ · '+tag+'-'+rank.toString(36).toUpperCase();return first[Math.floor(random(id,701)*16)]+last[Math.floor(random(id,709)*16)]+' · '+tag+'-'+rank.toString(36).toUpperCase();}
function formatMass(mass){const kg=mass*1.98847e30,earth=mass/C.EARTH_MSUN;return earth<.01?kg.toExponential(3)+' kg':earth>=31.78?(earth/317.828).toPrecision(4)+' M♃':earth.toPrecision(4)+' M⊕';}
function body(system,id){
 if(id===system.starId)return{id:system.starId,name:typeof root.stellarName==='function'?root.stellarName(system.starId):system.starId,kind:'Étoile',mass:system.host.mass,radiusKm:system.host.radius*P.solarRadiusM/1000,a:0,e:0,omega:0,mean:0};
 const planet=system.seeds.find(p=>p.id===id);if(planet)return planet;
 const prefix=system.prefix+'/body/';if(!String(id).startsWith(prefix)||!/^\d+$/.test(String(id).slice(prefix.length)))throw Error('Identité inconnue dans ce système.');const rank=BigInt(String(id).slice(prefix.length));if(String(id)!==prefix+rank.toString())throw Error('Rang canonique attendu.');if(rank>=system.smallCount)throw Error('Rang hors de la population de petits corps.');
 const data=Formation.residualBody(system.formation,id,rank-BigInt(system.seeds.length));
 return{id,rank,...data,kind:classifyMass(data.mass),name:objectName(id,data.mass),albedo:.04+.2*random(id,49),color:[.62,.58,.5]};
}
function phase(timeNs,period){const periodNs=BigInt(Math.max(1,Math.round(period*Number(C.YEAR_NS))));return Number((timeNs%periodNs+periodNs)%periodNs)/Number(periodNs)*2*Math.PI;}
function orbit(elements,mass,timeNs=0n){
 if(!elements.a)return{x:0,y:0,vx:0,vy:0};const {a,e,omega}=elements,period=Math.sqrt(a**3/mass),M=elements.mean+phase(timeNs,period);let E=M;for(let i=0;i<12;i++){const step=(E-e*Math.sin(E)-M)/(1-e*Math.cos(E));E-=step;if(Math.abs(step)<1e-13)break;}
 const c=Math.cos(omega),s=Math.sin(omega),root=Math.sqrt(1-e*e),n=2*Math.PI/period,k=n*a/(1-e*Math.cos(E)),x=a*(Math.cos(E)-e),y=a*root*Math.sin(E),vx=-k*Math.sin(E),vy=k*root*Math.cos(E);return{x:x*c-y*s,y:x*s+y*c,vx:vx*c-vy*s,vy:vx*s+vy*c};
}
function propagate(state,mass,elapsedNs){
 if(elapsedNs===0n)return{...state};const r=Math.hypot(state.x,state.y),v2=state.vx**2+state.vy**2,mu=C.G*mass,h=state.x*state.vy-state.y*state.vx,energy=v2/2-mu/r,a=-mu/(2*energy),ex=state.vy*h/mu-state.x/r,ey=-state.vx*h/mu-state.y/r,e=Math.hypot(ex,ey);
 if(!(a>0&&e<.999999&&h>0)){const t=Number(elapsedNs)/Number(C.YEAR_NS);return{...state,x:state.x+state.vx*t,y:state.y+state.vy*t,ballistic:true};}
 const omega=e>1e-10?Math.atan2(ey,ex):Math.atan2(state.y,state.x),c=Math.cos(omega),s=Math.sin(omega),x=state.x*c+state.y*s,y=-state.x*s+state.y*c,E=Math.atan2(y/(a*Math.sqrt(1-e*e)),x/a+e);return{...state,...orbit({a,e,omega,mean:E-e*Math.sin(E)},mass,elapsedNs)};
}
class Manager{
 constructor(){this.device=null;this.active=null;this.histories=new Map();this.maxSystems=8;this.token=0;}
 setDevice(device){this.device=device;}
 deactivate(){this.token++;if(this.active){this.active.engine?.destroy();this.active.engine=null;this.active=null;}}
 destroy(){this.deactivate();this.histories.clear();this.device=null;}
 async activate(starId,timeNs,budget=null,extraId=null){
  if(!this.device)throw Error('Le système détaillé demande WebGPU.');let record=this.histories.get(starId);
  if(!record){if(this.histories.size>=this.maxSystems)throw Error('Huit systèmes mémorisés : réinitialiser pour explorer un nouvel univers.');record={system:describe(starId),records:new Map(),epochNs:timeNs,states:[],engine:null,signature:''};this.histories.set(starId,record);}
  budget=budget===null?record.system.defaultCount:Math.max(record.system.seeds.length,Math.min(30000,Number(record.system.smallCount),Math.floor(budget)));const signature=budget+':'+(extraId||'');if(this.active===record&&record.signature===signature&&record.engine&&!record.rebase)return record;const ids=[starId];for(let i=0;i<budget;i++)ids.push(record.system.prefix+'/body/'+i);if(extraId&&!ids.includes(extraId)){body(record.system,extraId);ids.push(extraId);}
  if(this.active===record&&record.signature===signature&&record.engine&&!record.rebase)return record;
  const rebuilding=this.active===record;this.deactivate();if(rebuilding)this.active=record;const token=this.token;record.engine?.destroy();const states=ids.map(id=>{const meta=body(record.system,id),old=record.records.get(id),p=id===starId?{x:0,y:0,vx:0,vy:0}:old?propagate(old.state,record.system.host.mass+meta.mass,timeNs-old.epochNs):orbit(meta,record.system.host.mass+meta.mass,timeNs);return{...meta,...p};});
  const engine=await root.createOrbitGPUOrbitSystem({device:this.device,bodies:states,lengthScale:Math.max(record.system.extent,...states.map(p=>Math.hypot(p.x,p.y)))*1.4});if(token!==this.token){engine.destroy();return null;}
  record.states=states;record.byId=new Map(states.map(p=>[p.id,p]));record.epochNs=timeNs;record.signature=signature;record.engine=engine;record.stepYear=Math.min(.05,...states.slice(1).map(p=>Math.sqrt(Math.max(.001,p.a)**3/record.system.host.mass)*(1-p.e)**1.5/128));record.lastMs=0;record.rebase=false;this.active=record;this.save();return record;
 }
 save(){const r=this.active;if(!r)return;for(const p of r.states){if(!r.records.has(p.id)&&r.records.size>=65536)throw Error('Budget des rencontres matérialisées de ce système atteint.');let entry=r.records.get(p.id);if(!entry){entry={};r.records.set(p.id,entry);}entry.state=p;entry.epochNs=r.epochNs;}}
 async advance(years){const r=this.active;if(!r||!(years>0))return{advanced:0,steps:0};const limit=Math.max(1,Math.min(64,Math.floor(6/Math.max(.1,(r.lastMs||6)/(r.lastSteps||8))))),steps=Math.max(1,Math.min(limit,Math.ceil(years/r.stepYear))),dt=Math.min(r.stepYear,years/steps),result=await r.engine.advance(dt,steps);if(this.active!==r)return{advanced:0,steps:0};r.states=result.bodies;r.byId=new Map(r.states.map(p=>[p.id,p]));r.epochNs+=BigInt(Math.round(result.advanced*Number(C.YEAR_NS)));r.lastMs=r.engine.lastMs;r.lastSteps=steps;r.rebase=r.states.some(p=>Math.max(Math.abs(p.x),Math.abs(p.y))>r.engine.scale*.8);this.save();return result;}
 resolve(id,offsetYears=0){const r=this.active;if(!r)throw Error('Aucun système actif.');const found=r.byId.get(id);if(found){if(id===r.system.starId)return{...found,x:0,y:0,vx:0,vy:0};const ax=(found.ax||0)-(r.states[0].ax||0),ay=(found.ay||0)-(r.states[0].ay||0);return{...found,x:found.x+found.vx*offsetYears+.5*ax*offsetYears**2,y:found.y+found.vy*offsetYears+.5*ay*offsetYears**2,vx:found.vx+ax*offsetYears,vy:found.vy+ay*offsetYears};}const meta=body(r.system,id),old=r.records.get(id);return{...meta,...(old?propagate(old.state,r.system.host.mass+meta.mass,r.epochNs-old.epochNs+BigInt(Math.round(offsetYears*Number(C.YEAR_NS)))):orbit(meta,r.system.host.mass+meta.mass,r.epochNs+BigInt(Math.round(offsetYears*Number(C.YEAR_NS)))))};}
 points(offsetYears=0){return this.active?this.active.states.map(p=>this.resolve(p.id,offsetYears)):[];}
}
const api={...C,describe,body,orbit,propagate,classifyMass,objectName,formatMass,Manager};root.PlanetarySystems=api;if(typeof module!=='undefined')module.exports=api;
})(typeof window!=='undefined'?window:globalThis);
