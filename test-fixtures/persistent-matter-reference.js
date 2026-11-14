/* Lagrangian affine cohorts. Every virtual identity has immutable initial
 * coordinates. Only cohort maps evolve; no state or trajectory per virtual star.
 * One-way characteristic approximation driven by the collective kinetic field.
 */
(function(root){
'use strict';
function ratio(n,x){return n*BigInt(Math.floor(Math.max(0,Math.min(1,x))*9007199254740992))/9007199254740992n;}
function random(text,salt){let h=2166136261^salt;for(let i=0;i<text.length;i++)h=Math.imul(h^text.charCodeAt(i),16777619);h=Math.imul(h^(h>>>16),0x7feb352d);h=Math.imul(h^(h>>>15),0x846ca68b);return(((h^(h>>>16))>>>0)+.5)/4294967296;}
function split(p){let start=p.start;const q=p.count/4n,r=Number(p.count%4n);return Array.from({length:4},(_,k)=>{const count=q+BigInt(k<r),child={start,count,x:p.x+k%2*p.w/2,y:p.y+Math.floor(k/2)*p.w/2,w:p.w/2};start+=count;return child;});}
class PersistentMatter{
 constructor(grid,{population=10n**21n,seed='417'}={}){
  this.n=grid.n;this.nv=grid.nv;this.dx=grid.dx;this.dv=grid.dv;this.extent=grid.extent;this.seed=String(seed);this.time=0;this.coordinateCache=new Map();this.coordinateCacheLimit=32768;
  let entries=grid.initialCohorts?grid.initialCohorts.map(c=>[c.key,c.mass]):grid.initialPhase;if(!entries){entries=[];for(let i=0;i<grid.f.length;i++)if(grid.f[i]>0)entries.push([i,grid.f[i]]);}
  this.groups=entries.filter(([,m])=>m>0).sort((a,b)=>a[0]-b[0]).map(([key,mass])=>({key,mass,count:0n}));
  this.lookup=new Map(this.groups.map((g,i)=>[g.key,i]));this.data=new Float64Array(this.groups.length*12);this.forces=new Float64Array(this.groups.length*6);
  for(let i=0;i<this.groups.length;i++){const key=this.groups[i].key,cell=Math.floor(key/(this.nv*this.nv)),v=key%(this.nv*this.nv),o=i*12;
   const cohort=grid.initialCohorts?.[i];this.data[o]=cohort?cohort.x:(cell%this.n+.5)*this.dx-this.extent;this.data[o+1]=cohort?cohort.y:(Math.floor(cell/this.n)+.5)*this.dx-this.extent;
   this.data[o+2]=cohort?cohort.vx:(v%this.nv-(this.nv-1)/2)*this.dv;this.data[o+3]=cohort?cohort.vy:(Math.floor(v/this.nv)-(this.nv-1)/2)*this.dv;this.data[o+4]=this.data[o+9]=cohort?cohort.width:this.dx;
  }this.setPopulation(population);
  this.roots=this.groups.map((g,index)=>{const node={index,start:0n,count:g.count,children:null,depth:0};g.start=0n;g.node=node;return node;});this.maxGroups=Math.max(65536,this.groups.length);this.refinementLimited=false;
 }
 setPopulation(population){if(this.roots){if(BigInt(population)===this.population)return;throw Error('Changer la population exige de réinitialiser cette expérience.');}this.population=BigInt(population);this.universe=`${this.seed}:${this.n}x${this.nv}:${this.population}`;const total=this.groups.reduce((s,g)=>s+g.mass,0);let sum=0,previous=0n;
  for(let i=0;i<this.groups.length;i++){const g=this.groups[i];sum+=g.mass;const end=i===this.groups.length-1?this.population:ratio(this.population,sum/total);g.count=end-previous;previous=end;}
 }
 setField(ax,ay){const n=this.n,dx=this.dx,d=this.data;
  for(let o=0;o<d.length;o+=12){const x=d[o],y=d[o+1];let fx=0,fy=0,jxx=0,jxy=0,jyx=0,jyy=0;
   if(Math.abs(x)<=this.extent&&Math.abs(y)<=this.extent){const u=(x+this.extent)/dx-.5,v=(y+this.extent)/dx-.5,ix=Math.max(0,Math.min(n-2,Math.floor(u))),iy=Math.max(0,Math.min(n-2,Math.floor(v))),tx=Math.max(0,Math.min(1,u-ix)),ty=Math.max(0,Math.min(1,v-iy)),i=iy*n+ix;
    const a=ax[i],b=ax[i+1],c=ax[i+n],e=ax[i+n+1],f=ay[i],g=ay[i+1],h=ay[i+n],k=ay[i+n+1];
    fx=(a*(1-tx)+b*tx)*(1-ty)+(c*(1-tx)+e*tx)*ty;fy=(f*(1-tx)+g*tx)*(1-ty)+(h*(1-tx)+k*tx)*ty;
    jxx=((b-a)*(1-ty)+(e-c)*ty)/dx;jxy=((c-a)*(1-tx)+(e-b)*tx)/dx;jyx=((g-f)*(1-ty)+(k-h)*ty)/dx;jyy=((h-f)*(1-tx)+(k-g)*tx)/dx;
   }
   const f=o/2;this.forces[f]=fx;this.forces[f+1]=fy;this.forces[f+2]=jxx;this.forces[f+3]=jxy;this.forces[f+4]=jyx;this.forces[f+5]=jyy;
  }
 }
 refine(dt){const old=this.data,f=this.forces,needMask=new Uint8Array(this.groups.length);let needs=false;
  for(let index=0;index<this.groups.length;index++){const g=this.groups[index],o=index*12,j=index*6;let span=0;
   for(const k of [4,8]){const px=old[o+k]+old[o+k+2]*dt+.5*(f[j+2]*old[o+k]+f[j+3]*old[o+k+1])*dt*dt,py=old[o+k+1]+old[o+k+3]*dt+.5*(f[j+4]*old[o+k]+f[j+5]*old[o+k+1])*dt*dt;span=Math.max(span,Math.hypot(px,py));}
   if(span>1.75*this.dx&&g.count>=4n&&g.node.depth<12){needMask[index]=1;needs=true;}}
  if(!needs)return false;if(this.groups.length+3>this.maxGroups){this.refinementLimited=true;return false;}
  const groups=[],values=[];let changed=false;for(let index=0;index<this.groups.length;index++){
   const g=this.groups[index],o=index*12,need=needMask[index],can=groups.length+this.groups.length-index+3<=this.maxGroups;
   if(need&&can){changed=true;const quad=split({start:g.start,count:g.count,x:-.5,y:-.5,w:1});g.node.children=[];
    for(let k=0;k<4;k++){const u=(k%2-.5)*.5,v=(Math.floor(k/2)-.5)*.5,node={index:groups.length,start:quad[k].start,count:quad[k].count,children:null,depth:g.node.depth+1};g.node.children.push(node);
     groups.push({...g,mass:g.mass*Number(node.count)/Number(g.count),count:node.count,start:node.start,node});const data=Array.from(old.subarray(o,o+12));for(let q=0;q<4;q++)data[q]+=old[o+4+q]*u+old[o+8+q]*v;for(let q=4;q<12;q++)data[q]*=.5;values.push(...data);
    }
   }else{if(need&&!can)this.refinementLimited=true;g.node.index=groups.length;groups.push(g);values.push(...old.subarray(o,o+12));}
  }
  if(changed){this.groups=groups;this.data=Float64Array.from(values);this.forces=new Float64Array(groups.length*6);}return changed;
 }
 step(dt,ax,ay){this.setField(ax,ay);for(let i=0;i<3;i++){if(!this.refine(dt))break;this.setField(ax,ay);}const d=this.data,f=this.forces;
  for(let o=0;o<d.length;o+=12){const j=o/2;d[o]+=d[o+2]*dt+.5*f[j]*dt*dt;d[o+1]+=d[o+3]*dt+.5*f[j+1]*dt*dt;d[o+2]+=f[j]*dt;d[o+3]+=f[j+1]*dt;
   for(const k of [4,8]){const ax=f[j+2]*d[o+k]+f[j+3]*d[o+k+1],ay=f[j+4]*d[o+k]+f[j+5]*d[o+k+1];d[o+k]+=d[o+k+2]*dt+.5*ax*dt*dt;d[o+k+1]+=d[o+k+3]*dt+.5*ay*dt*dt;d[o+k+2]+=ax*dt;d[o+k+3]+=ay*dt;}
  }this.time+=dt;
 }
 initialPoint(group,rank,node=null){let p=node||{start:group.start||0n,count:group.count,x:-.5,y:-.5,w:1};for(let depth=0;p.count>1n&&depth<60;depth++)p=split(p).find(c=>rank>=c.start&&rank<c.start+c.count);
  const key=this.universe+'/'+group.key+'/'+rank;return{u:p.x+p.w*random(key,11),v:p.y+p.w*random(key,17)};
 }
 // Cache only immutable material addresses. Positions/velocities are always
 // evaluated from the current cohort, and a subdivision invalidates the entry.
 evaluate(index,rank,offset=0,node=null){const g=this.groups[index];if(rank<g.start||rank>=g.start+g.count)throw Error('Identité hors population.');const id=this.universe+'/'+g.key.toString(36)+'/'+rank;let address=this.coordinateCache.get(id);if(!address||address.node!==g.node){const point=this.initialPoint(g,rank,node);address={...point,rank,node:g.node};if(this.coordinateCache.size>=this.coordinateCacheLimit&&!this.coordinateCache.has(id))this.coordinateCache.delete(this.coordinateCache.keys().next().value);this.coordinateCache.set(id,address);}const {u,v}=address,d=this.data,o=index*12,f=this.forces,j=index*6,vx=d[o+2]+d[o+6]*u+d[o+10]*v,vy=d[o+3]+d[o+7]*u+d[o+11]*v;
  const dx=d[o+4]*u+d[o+8]*v,dy=d[o+5]*u+d[o+9]*v,ax=f[j]+f[j+2]*dx+f[j+3]*dy,ay=f[j+1]+f[j+4]*dx+f[j+5]*dy;
  return{id,x:d[o]+dx+vx*offset+.5*ax*offset*offset,y:d[o+1]+dy+vy*offset+.5*ay*offset*offset,vx:vx+ax*offset,vy:vy+ay*offset};
 }
 resolve(id,offset=0){const cached=this.coordinateCache.get(id);if(cached){let node=cached.node;while(node.children)node=node.children.find(c=>cached.rank>=c.start&&cached.rank<c.start+c.count);return this.evaluate(node.index,cached.rank,offset);}const parts=String(id).split('/');if(parts.length!==3||parts[0]!==this.universe||!/^\d+$/.test(parts[2]))throw Error('Identité inconnue dans cette expérience.');const key=parseInt(parts[1],36),index=this.lookup.get(key),rank=BigInt(parts[2]);if(index===undefined)throw Error('Groupe inconnu.');let node=this.roots[index];while(node.children){node=node.children.find(c=>rank>=c.start&&rank<c.start+c.count);if(!node)throw Error('Identité hors population.');}return this.evaluate(node.index,rank,offset);}
 projectLight({x=0,y=0,zoom=1,offset=0,resolution=64}={}){
  const image=new Float64Array(resolution*resolution),width=2*this.extent/zoom,pixel=width/resolution,left=x-width/2,bottom=y-width/2,N=Number(this.population),d=this.data,f=this.forces;
  const deposit=(px,py,mass)=>{const u=(px-left)/pixel-.5,v=(py-bottom)/pixel-.5,ix=Math.floor(u),iy=Math.floor(v),tx=u-ix,ty=v-iy;for(let a=0;a<2;a++)for(let b=0;b<2;b++){const xx=ix+a,yy=iy+b;if(xx>=0&&yy>=0&&xx<resolution&&yy<resolution)image[yy*resolution+xx]+=mass*(a?tx:1-tx)*(b?ty:1-ty);}};
  for(let k=0;k<this.groups.length;k++){
   const g=this.groups[k];if(!g.count)continue;const o=k*12,j=k*6,h=.5*offset*offset,cx=d[o]+d[o+2]*offset+f[j]*h,cy=d[o+1]+d[o+3]*offset+f[j+1]*h;
   const au=d[o+4]+d[o+6]*offset+(f[j+2]*d[o+4]+f[j+3]*d[o+5])*h,av=d[o+8]+d[o+10]*offset+(f[j+2]*d[o+8]+f[j+3]*d[o+9])*h,bu=d[o+5]+d[o+7]*offset+(f[j+4]*d[o+4]+f[j+5]*d[o+5])*h,bv=d[o+9]+d[o+11]*offset+(f[j+4]*d[o+8]+f[j+5]*d[o+9])*h;
   const rx=(Math.abs(au)+Math.abs(av))/2,ry=(Math.abs(bu)+Math.abs(bv))/2,det=au*bv-av*bu,mass=Number(g.count)/N;
   if(cx+rx<left||cx-rx>left+width||cy+ry<bottom||cy-ry>bottom+width)continue;
   if(rx*2<pixel&&ry*2<pixel){deposit(cx,cy,mass);continue;}
   if(Math.abs(det)<1e-26){for(let a=0;a<4;a++)for(let b=0;b<4;b++){const u=(a+.5)/4-.5,v=(b+.5)/4-.5;deposit(cx+au*u+av*v,cy+bu*u+bv*v,mass/16);}continue;}
   const x0=Math.max(0,Math.floor((cx-rx-left)/pixel)),x1=Math.min(resolution-1,Math.floor((cx+rx-left)/pixel)),y0=Math.max(0,Math.floor((cy-ry-bottom)/pixel)),y1=Math.min(resolution-1,Math.floor((cy+ry-bottom)/pixel)),value=mass/Math.abs(det)*pixel*pixel;
   for(let iy=y0;iy<=y1;iy++)for(let ix=x0;ix<=x1;ix++){let coverage=0;for(let a=0;a<2;a++)for(let b=0;b<2;b++){const px=left+(ix+(a+.5)/2)*pixel-cx,py=bottom+(iy+(b+.5)/2)*pixel-cy,u=(bv*px-av*py)/det,v=(-bu*px+au*py)/det;if(Math.abs(u)<=.5&&Math.abs(v)<=.5)coverage++;}image[iy*resolution+ix]+=value*coverage/4;}
  }
  let mass=0;for(const m of image)mass+=m;
  return{image,resolution,width,pixel,left,bottom,expected:Number(this.population)*mass,mass};
 }
 sample({x=0,y=0,zoom=1,budget=512,offset=0,rankFilter=null,exact=false}={}){
  budget=Math.floor(Number(budget));if(!Number.isSafeInteger(budget)||budget<0)throw Error('Plafond de sources invalide.');if(exact)budget=Math.min(8192,budget);const width=2*this.extent/zoom,view={l:x-width/2,r:x+width/2,b:y-width/2,t:y+width/2};
  const groups=[];let visited=0,total=0n,fullTotal=0n;const visible=p=>p.x>=view.l&&p.x<=view.r&&p.y>=view.b&&p.y<=view.t;
  for(let index=0;index<this.groups.length;index++){
   const g=this.groups[index];if(!g.count)continue;const filter=rankFilter?.forGroup(this.universe+'/'+g.key.toString(36)),o=index*12,d=this.data,f=this.forces,j=index*6,h=.5*offset*offset,cx=d[o]+d[o+2]*offset+f[j]*h,cy=d[o+1]+d[o+3]*offset+f[j+1]*h,au=d[o+4]+d[o+6]*offset+(f[j+2]*d[o+4]+f[j+3]*d[o+5])*h,av=d[o+8]+d[o+10]*offset+(f[j+2]*d[o+8]+f[j+3]*d[o+9])*h,bu=d[o+5]+d[o+7]*offset+(f[j+4]*d[o+4]+f[j+5]*d[o+5])*h,bv=d[o+9]+d[o+11]*offset+(f[j+4]*d[o+8]+f[j+5]*d[o+9])*h;
   const visit=(node,depth=0)=>{if(!node.count||visited>=250000)return;const mx=node.x+node.w/2,my=node.y+node.w/2,px=cx+au*mx+av*my,py=cy+bu*mx+bv*my,rx=(Math.abs(au)+Math.abs(av))*node.w/2,ry=(Math.abs(bu)+Math.abs(bv))*node.w/2;
    if(px+rx<view.l||px-rx>view.r||py+ry<view.b||py-ry>view.t)return;visited++;
    if(node.count===1n){const p=this.evaluate(index,node.start,offset,node);if(visible(p)){fullTotal++;if(!filter||filter.count(node.start,1n)){groups.push({index,node,eligible:1n,filter,weight:1n,prefix:total});total++;}}return;}
    const eligible=filter?filter.count(node.start,node.count):node.count;let weight=eligible,fullWeight=node.count;
    if(!(px-rx>=view.l&&px+rx<=view.r&&py-ry>=view.b&&py+ry<=view.t)){
     if((exact||Math.max(rx,ry)*2>width/600)&&depth<60){for(const child of split(node))visit(child,depth+1);return;}
     const area=Math.max(0,Math.min(px+rx,view.r)-Math.max(px-rx,view.l))*Math.max(0,Math.min(py+ry,view.t)-Math.max(py-ry,view.b)),coverage=area/Math.max(4*rx*ry,1e-300);weight=ratio(eligible,coverage);fullWeight=ratio(node.count,coverage);
    }
    fullTotal+=fullWeight;if(weight){groups.push({index,node,eligible,filter,weight,prefix:total});total+=weight;}
   };visit({start:g.start,count:g.count,x:-.5,y:-.5,w:1});
  }
  const points=[],seen=new Set(),wanted=Number(total<BigInt(budget)?total:BigInt(budget)),attempts=total<=BigInt(budget)?wanted:budget*3;
  for(let j=0;j<attempts&&points.length<wanted;j++){const rank=total<=BigInt(budget)?BigInt(j):ratio(total,(j*.6180339887498949+.5)%1);let lo=0,hi=groups.length-1;
   while(lo<hi){const m=(lo+hi)>>>1;if(groups[m].prefix+groups[m].weight<=rank)lo=m+1;else hi=m;}const g=groups[lo];if(!g)break;
   const local=rank-g.prefix,ordinal=g.weight===g.eligible?local:ratio(g.eligible,random(g.index+':'+local,23)),id=g.filter?g.filter.at(g.node.start,ordinal):g.node.start+ordinal,p=this.evaluate(g.index,id,offset,g.node);
   if(visible(p)&&!seen.has(p.id)){points.push(p);seen.add(p.id);}
  }
  return{points,estimated:fullTotal,eligible:total,visited,truncated:visited>=250000};
 }
}
root.PersistentMatter=PersistentMatter;if(typeof module!=='undefined')module.exports=PersistentMatter;
})(typeof window!=='undefined'?window:globalThis);
