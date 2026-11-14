const assert=require('node:assert/strict'),{selectObservedSources}=require('./observation-tools');
const a={point:{id:'a',x:2,vx:3},flux:1},b={point:{id:'b',x:4,vx:5},flux:2},bright={point:{id:'new',x:6,vx:7},flux:100};
for(const policy of ['steady','zoom-in']){
 const selected=selectObservedSources([bright,b,a],['a','b'],2);assert.deepEqual(new Set(selected.map(p=>p.id)),new Set(['a','b']),policy+' must not evict an existing source when a brighter candidate is discovered');assert.equal(selected.find(p=>p.id==='a').x,2,'preservation uses updated position, not a frozen snapshot');
}
const zoomOut=selectObservedSources([a,b,bright],[],2);assert.deepEqual(zoomOut.map(p=>p.id),['new','b'],'zoom-out ranks fresh candidates without historical priority');
assert.deepEqual(selectObservedSources([bright,a],['b'],2).map(p=>p.id),['new','a'],'departed sources are not fabricated');assert.deepEqual(selectObservedSources([a,b,bright],['a','b'],1).map(p=>p.id),['b'],'an explicitly reduced budget still bounds source count');
console.log('Passed: existing sources protected at constant zoom and zoom-in, fresh ranking on zoom-out, updated positions and bounded selection.');
