// node seq_engine.test.mjs: the JS sequence analyses reproduce the Python results in tests/fixtures/engine_fixture.json
// (regenerate with python scripts/make_engine_fixture.py).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { close } from './compare.test.mjs';
import { classifyBundle, domainsBundle, fmt1, pySum, hydropathyBundle, parseBundle, parseFasta } from './seq_engine.js';

const fx = JSON.parse(readFileSync(new URL('../../tests/fixtures/engine_fixture.json', import.meta.url)));
const cfg = fx.config;
let n = 0;
for (const [name, c] of Object.entries(fx.cases)) {
  close(parseBundle(c.fasta), c.parse, name + '.parse'); n++;
  close(domainsBundle(c.fasta, cfg), c.domains, name + '.domains'); n++;
  close(classifyBundle(c.fasta, cfg), c.classify, name + '.classify'); n++;
  c.hydropathy.forEach((h, i) => { close(hydropathyBundle(c.fasta, cfg, i), h, `${name}.hydropathy[${i}]`); n++; });
  c.hydropathy_w11_t1.forEach((h, i) => { close(hydropathyBundle(c.fasta, cfg, i, 11, 1.0), h, `${name}.hyd11[${i}]`); n++; });
}
// Python's round-half-even on exact ties, and plain rounding otherwise
assert.equal(fmt1(51.25), '51.2'); assert.equal(fmt1(51.75), '51.8'); assert.equal(fmt1(0.15), '0.1'); assert.equal(fmt1(75), '75.0');
// CPython 3.12 sum() is compensated: sum([0.1] * 10) == 1.0 exactly
assert.equal(pySum(Array(10).fill(0.1)), 1.0);
assert.throws(() => parseFasta('MKV'), /before the first/);
assert.throws(() => parseFasta('  '), /No sequences/);
console.log(`seq_engine: ${n} fixture comparisons ok`);
