const assert = require('node:assert/strict');
const { stellarSourceBudget } = require('./src/observation/stellar-budget');

const options = Object.freeze({lengthPc:10000,population:100000000000n,extent:2,adaptive:true,manualBudget:8192,sourceLimit:65536,calibration:{surfaceSolar:30,stars:100000000000n}});
const field = Object.freeze({width:.01,aspectX:1,aspectY:.5,expected:150000});
const base = stellarSourceBudget(field,options);
assert.equal(base.densityRatio,1);
assert.equal(stellarSourceBudget({...field,width:.02,expected:600000},options).budget,base.budget,'same physical density, same budget');
assert.equal(stellarSourceBudget({...field,aspectX:.5,aspectY:1},options).budget,base.budget,'portrait and landscape are equivalent');
const poor = stellarSourceBudget({...field,expected:0},options);
const dense = stellarSourceBudget({...field,expected:1e9},options);
assert.equal(poor.budget,4096);
assert.equal(dense.budget,16384);
assert(base.budget>=poor.budget && base.budget<=dense.budget);
assert.equal(stellarSourceBudget(field,{...options,sourceLimit:5000}).budget,5000);
assert.equal(stellarSourceBudget(field,{...options,adaptive:false,manualBudget:12345}).budget,12345);
const fallback = stellarSourceBudget({width:4,expected:1e11}, {...options,calibration:null});
assert.equal(fallback.densityRatio,1,'uncalibrated scenarios use the mean domain density');
console.log('Passed: density-based budgets, rectangular invariance, manual settings and hardware bound.');
