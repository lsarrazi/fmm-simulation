const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const html = fs.readFileSync('test-fixtures/legacy-continuum.html', 'utf8');
const elements = {};
const context2d = new Proxy({ createImageData: (w, h) => ({ data: new Uint8ClampedArray(w * h * 4) }) }, { get: (target, key) => target[key] ?? (() => {}) });
function element(id) {
  return elements[id] ||= { id, value: '', textContent: '', checked: false, hidden: false, listeners: {},
    addEventListener(type, fn) { (this.listeners[type] ||= []).push(fn); },
    setAttribute() {}, getContext() { return context2d; }, getBoundingClientRect() { return { width: 900, height: 700, left: 0, top: 0 }; },
    dispatch(type, event = {}) { for (const fn of this.listeners[type] || []) fn(event); } };
}
for (const match of html.matchAll(/\bid="([^"]+)"/g)) element(match[1]);
for (const match of html.matchAll(/<input\b([^>]+)>/g)) {
  const id = /\bid="([^"]+)"/.exec(match[1])?.[1]; if (!id) continue;
  element(id).value = /\bvalue="([^"]*)"/.exec(match[1])?.[1] ?? '';
  element(id).checked = /\bchecked\b/.test(match[1]);
}
for (const match of html.matchAll(/<select\b[^>]*id="([^"]+)"[^>]*>([\s\S]*?)<\/select>/g)) {
  const options = [...match[2].matchAll(/<option\s+value="([^"]+)"([^>]*)>/g)];
  element(match[1]).value = (options.find(o => /selected/.test(o[2])) || options[0])[1];
}
let queued;
const sandbox = { console, performance: { now: () => 0 }, devicePixelRatio: 1, URLSearchParams, location: { search: '' },
  requestAnimationFrame(fn) { queued = fn; },
  document: { hidden: false, getElementById(id) { assert(elements[id], 'missing HTML element: ' + id); return elements[id]; }, createElement() { return element('imageCanvas'); }, addEventListener() {} } };
sandbox.window = sandbox;
vm.createContext(sandbox);
for (const file of ['continuum-solver.js', 'test-fixtures/legacy-identities.js', 'test-fixtures/legacy-view.js']) vm.runInContext(fs.readFileSync(file, 'utf8'), sandbox, { filename: file });
assert(queued);
for (let i = 1; i <= 4; i++) queued(i * 400);
assert.match(element('sampleStatus').textContent, /1000000000000000000000000/);
element('identityInput').value = '9007199254740993';
element('findIdentity').dispatch('click');
assert.match(element('identityDetail').textContent, /9007199254740993/);
element('pause').dispatch('click');
const detail = element('identityDetail').textContent;
element('rebuildIdentity').dispatch('click');
assert.equal(element('identityDetail').textContent, detail, 'same state after reconstructing while paused');
element('forget').dispatch('click'); queued(2000);
assert.equal(element('identityDetail').textContent, detail);
element('matter').dispatch('pointermove', { clientX: 450, clientY: 350 });
element('matter').dispatch('pointerdown', { clientX: 100, clientY: 100 });
element('identityInput').value = '-1'; element('findIdentity').dispatch('click');
assert.match(element('identityDetail').textContent, /entier/);
element('population').value = '60'; element('reset').dispatch('click'); queued(2400);
assert.match(element('sampleStatus').textContent, /1000000000000000000000000000000000000000000000000000000000000/);
element('initial').value = 'cold'; element('initial').dispatch('change');
assert(element('sigma').disabled);
element('reset').dispatch('click'); queued(2800);
assert.match(element('status').textContent, /Bilan de masse/);
console.log('Passed: view startup, sampling status, exact ID input, paused reconstruction, cache clearing, pointer controls, invalid input, population reset and cold mode.');
