/* Frozen present-day field-star population. No per-star state or random registry.
 * Number fractions: CNS5 (2023); low-mass slopes: Kroupa (2001), high-mass
 * present-day slope: Bovy (2017). MS calibrations: Eker et al. (2018).
 * Giant luminosity function and WD cooling/mass spread below are approximations,
 * not a fitted Gaia HR diagram or a universal stellar population synthesis model.
 */
(function(root){
'use strict';
const C={solarLuminosityW:3.828e26,solarRadiusM:6.957e8,parsecM:3.0856775814913673e16,solarTemperatureK:5772};
C.solarRadiusPc=C.solarRadiusM/C.parsecM;C.solarFluxAtPc=C.solarLuminosityW/(4*Math.PI*C.parsecM*C.parsecM);
function hash(id,salt=0){let h=2166136261^salt;for(let i=0;i<id.length;i++)h=Math.imul(h^id.charCodeAt(i),16777619);h=Math.imul(h^(h>>>16),0x7feb352d);h=Math.imul(h^(h>>>15),0x846ca68b);return(((h^(h>>>16))>>>0)+.5)/4294967296;}
function reverse32(x){x=((x>>>1)&0x55555555)|((x&0x55555555)<<1);x=((x>>>2)&0x33333333)|((x&0x33333333)<<2);x=((x>>>4)&0x0f0f0f0f)|((x&0x0f0f0f0f)<<4);x=((x>>>8)&0x00ff00ff)|((x&0x00ff00ff)<<8);return((x>>>16)|(x<<16))>>>0;}
const MASK53=(1n<<53n)-1n;
function reverse53(rank){return(BigInt(reverse32(Number(rank&4294967295n)))<<21n)|BigInt(reverse32(Number((rank>>32n)&2097151n))>>>11);}
function rootSalt(id){return(BigInt(Math.floor(hash(id,503)*4294967296))<<21n)|BigInt(Math.floor(hash(id,509)*2097152));}
function variate(id){const cut=id.lastIndexOf('/'),rank=id.slice(cut+1);if(cut<0||!/^\d+$/.test(rank))return hash(id,101);return(Number(reverse53(BigInt(rank)&MASK53)^rootSalt(id.slice(0,cut)))+.5)/9007199254740992;}
function power(a,b,alpha,u){const p=1-alpha;return (a**p+u*(b**p-a**p))**(1/p);}
function integral(a,b,alpha){return(b**(1-alpha)-a**(1-alpha))/(1-alpha);}
function msLuminosity(m){if(m<.179)return msLuminosity(.179)*(m/.179)**3;const [a,b]=m<.45?[2.028,-.976]:m<.72?[4.572,-.102]:m<1.05?[5.743,-.007]:m<2.4?[4.329,.010]:m<7?[3.967,.093]:[2.865,1.105];return 10**b*m**a;}
function msRadius(m,L){if(m<.179)return(.438*.179**2+.479*.179+.075)*(m/.179);if(m<=1.5)return .438*m*m+.479*m+.075;const log=Math.log10(m),T=10**(-.170*log*log+.888*log+3.671);return Math.sqrt(L)/(T/C.solarTemperatureK)**2;}
function wdRadius(m){return .0112*Math.sqrt((1.44/m)**(2/3)-(m/1.44)**(2/3));}
function visibleFraction(T){const a=.01438777/(700e-9*T),b=.01438777/(400e-9*T),n=16,h=(b-a)/n;let sum=0;for(let i=0;i<=n;i++){const x=a+i*h;sum+=(i===0||i===n?1:i%2?4:2)*x*x*x/Math.expm1(x);}return sum*h/3*15/Math.PI**4;}
function color(T){const rgb=[650e-9,550e-9,450e-9].map(w=>Math.expm1(.01438777/(w*6500))/Math.expm1(.01438777/(w*T))),Y=.2126*rgb[0]+.7152*rgb[1]+.0722*rgb[2];return rgb.map(c=>c/Y);}
const bins=[];
// Split at every empirical mass-law boundary so each bin has one power law.
const edges=[.08,.179,.45,.5,.72,1,1.05,1.5,2.4,7,31,100];for(let i=1;i<192;i++)edges.push(.08*(100/.08)**(i/192));edges.sort((a,b)=>a-b);
let norm=0;for(let i=0;i<edges.length-1;i++){const a=edges[i],b=edges[i+1],alpha=a>=1?4.7:a>=.5?2.3:1.3,w=(a>=.5?.5:1)*integral(a,b,alpha);bins.push({kind:'ms',a,b,alpha,weight:w,maxL:Math.max(msLuminosity(a),msLuminosity(b*(1-1e-12)))*10**.15});norm+=w;}
for(const bin of bins)bin.weight*=4946/5230/norm;
// Smooth red-clump + RGB approximation; only its total number fraction is CNS5.
const evolved=[];let gn=0,wn=0;for(let i=0;i<64;i++){
 const a=5*(1e5/5)**(i/64),b=5*(1e5/5)**((i+1)/64),mid=Math.sqrt(a*b),dlog=Math.log(b/a);
 const w=(.75*Math.exp(-.5*(Math.log(mid/50)/.4)**2)/(.4*Math.sqrt(2*Math.PI))+.25*(mid/5)**(-.55)*.55)*dlog;
 evolved.push({kind:'giant',a,b,weight:w,maxL:b});gn+=w;
 const ta=4000*10**(i/64),tb=4000*10**((i+1)/64),ww=integral(ta,tb,3.8);
 evolved.push({kind:'wd',a:ta,b:tb,weight:ww,maxL:wdRadius(.45)**2*(tb/C.solarTemperatureK)**4});wn+=ww;
}
for(const b of evolved){b.weight*=b.kind==='giant'?20/5230/gn:264/5230/wn;bins.push(b);}bins.sort((a,b)=>a.maxL-b.maxL);
let cumulative=0;for(const b of bins){b.lower=cumulative;cumulative+=b.weight;b.upper=cumulative;}bins[bins.length-1].upper=1;
function fromBin(b,u,v,w){let mass,L,R,T,type;
 if(b.kind==='ms'){mass=power(b.a,b.b,b.alpha,u);L=msLuminosity(mass)*10**((v-.5)*.3);R=msRadius(mass,L)*10**((w-.5)*.08);T=C.solarTemperatureK*(L/(R*R))**.25;type='Séquence principale';}
 else if(b.kind==='giant'){L=b.a*(b.b/b.a)**u;mass=Math.max(power(.8,8,3.7,v),Math.min(18,(L/1000)**.4));T=Math.max(3100,4900-420*Math.log10(L/5))*(.97+.06*w);R=Math.sqrt(L)/(T/C.solarTemperatureK)**2;type=L>1e4?'Géante très lumineuse':'Géante rouge';}
 else{mass=.45+.65*(1-Math.sqrt(1-v))**1.7;T=power(b.a,b.b,3.8,u);R=wdRadius(mass);L=R*R*(T/C.solarTemperatureK)**4;type='Naine blanche';}
 const visible=visibleFraction(T);return{mass,luminosity:L,radius:R,temperature:T,type,visibleLuminosity:L*visible,color:color(T)};
}
function properties(id){const u=variate(id);let lo=0,hi=bins.length-1;while(lo<hi){const mid=(lo+hi)>>>1;if(bins[mid].upper<=u)lo=mid+1;else hi=mid;}const b=bins[lo];return fromBin(b,hash(id,601),hash(id,607),hash(id,613));}
// Constant-size moments, including luminous rare tails (not a small Monte Carlo).
const luminositySamples=[],moments={mass:0,luminosity:0,visibleLuminosity:0,color:[0,0,0]};for(const b of bins)for(let i=0;i<16;i++){
 const p=fromBin(b,(i+.5)/16,((i*7)%16+.5)/16,((i*11)%16+.5)/16),w=b.weight/16;luminositySamples.push({L:p.visibleLuminosity,w});moments.mass+=p.mass*w;moments.luminosity+=p.luminosity*w;moments.visibleLuminosity+=p.visibleLuminosity*w;for(let k=0;k<3;k++)moments.color[k]+=p.color[k]*p.visibleLuminosity*w;
}moments.color=moments.color.map(v=>v/moments.visibleLuminosity);
// Smooth number/light CDF of the same bounded population quadrature. It lets
// resolved sources replace their ensemble mean without digging local dark holes.
luminositySamples.sort((a,b)=>a.L-b.L);const lightCDF=[{x:Math.log(luminositySamples[0].L/2),p:0,l:0}];let cp=0,cl=0;
for(const s of luminositySamples){lightCDF.push({x:Math.log(s.L),p:cp+s.w/2,l:cl+s.w*s.L/2});cp+=s.w;cl+=s.w*s.L;}
lightCDF.push({x:Math.log(luminositySamples.at(-1).L*2),p:1,l:moments.visibleLuminosity});
function luminosityTailMoments(limit){if(!(limit>0))return{probability:1,luminosity:moments.visibleLuminosity};const x=Math.log(limit);if(x<=lightCDF[0].x)return{probability:1,luminosity:moments.visibleLuminosity};if(x>=lightCDF.at(-1).x)return{probability:0,luminosity:0};let lo=0,hi=lightCDF.length-1;while(hi-lo>1){const m=(lo+hi)>>>1;if(lightCDF[m].x<x)lo=m;else hi=m;}const a=lightCDF[lo],b=lightCDF[hi],u=(x-a.x)/(b.x-a.x);return{probability:Math.max(0,1-a.p-u*(b.p-a.p)),luminosity:Math.max(0,moments.visibleLuminosity-a.l-u*(b.l-a.l))};}
// Dyadic luminosity tail: one congruence class of ranks, not a scan. The tail
// conservatively includes all eligible bins and at most twice their probability.
function luminosityFilter(minL){const index=bins.findIndex(b=>b.maxL>=minL),fraction=index<0?0:1-bins[index].lower,k=fraction<=0?54:Math.min(53,Math.max(0,Math.floor(-Math.log2(fraction)))),cache=new Map();
 return{key:k,fraction:k>53?0:2**(-k),forGroup(rootKey){if(cache.has(rootKey))return cache.get(rootKey);const bits=Math.min(k,53),shift=BigInt(bits),mod=1n<<shift,lowMask=mod-1n,mask=MASK53^((1n<<BigInt(53-bits))-1n),residue=reverse53(mask^rootSalt(rootKey))&lowMask;
  const first=start=>start+((residue-start)&lowMask);
  const filter={count(start,count){const a=first(start),end=start+count;return k>53||a>=end?0n:1n+((end-1n-a)>>shift);},at(start,index){return first(start)+(index<<shift);}};cache.set(rootKey,filter);return filter;}};
}
const api={...C,hash,properties,variate,bins,moments,luminosityFilter,luminosityTailMoments,msLuminosity,visibleFraction,meanMass:moments.mass,meanLuminosity:moments.luminosity,meanVisibleLuminosity:moments.visibleLuminosity,meanColor:moments.color};root.StellarPopulation=api;if(typeof module!=='undefined')module.exports=api;
})(typeof window!=='undefined'?window:globalThis);
