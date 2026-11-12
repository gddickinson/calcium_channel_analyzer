// Shared helper for the engine tests: deep comparison with a relative tolerance of 1e-9 for numbers.
// (Named *.test.mjs so the static build leaves it out; running it alone does nothing.)
import assert from 'node:assert/strict';

export function close(a, b, path = '') {
  if (typeof a === 'number' && typeof b === 'number') {
    assert.ok(Math.abs(a - b) <= 1e-9 * Math.max(1, Math.abs(b)), `${path}: ${a} vs ${b}`); return;
  }
  if (Array.isArray(b)) { assert.equal(a.length, b.length, path + ' length'); b.forEach((v, i) => close(a[i], v, `${path}[${i}]`)); return; }
  if (b && typeof b === 'object') {
    assert.deepEqual(Object.keys(a).sort(), Object.keys(b).sort(), path + ' keys');
    for (const k of Object.keys(b)) close(a[k], b[k], `${path}.${k}`); return;
  }
  assert.equal(a, b, path);
}
