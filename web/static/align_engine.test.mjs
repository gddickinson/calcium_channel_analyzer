// node align_engine.test.mjs: alignment scores, identity distances and UPGMA/NJ trees match Biopython
// (tests/fixtures/engine_fixture.json, made by scripts/make_engine_fixture.py).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { alignBundle, globalAlign, scoreAlignment, treeBundle } from './align_engine.js';
import { close } from './compare.test.mjs';

const fx = JSON.parse(readFileSync(new URL('../../tests/fixtures/engine_fixture.json', import.meta.url)));
const cfg = fx.config, sch = { match: 1, mismatch: 0, open: -0.5, extend: -0.1 };
let n = 0;
for (const [name, c] of Object.entries(fx.cases)) {
  const js = alignBundle(c.fasta, cfg);
  assert.deepEqual(js.ids, c.align.ids);
  c.align.pairs.forEach((p, k) => {
    const q = js.pairs[k];
    close(q.score, p.score, `${name} pair ${k} score`);
    // the JS alignment is a valid optimal one: it re-scores to the optimum and spells both inputs
    close(scoreAlignment(q.aligned1, q.aligned2, sch), p.score, `${name} pair ${k} rescore`);
    assert.equal(q.aligned1.replace(/-/g, ''), p.aligned1.replace(/-/g, ''));
    assert.equal(q.aligned2.replace(/-/g, ''), p.aligned2.replace(/-/g, ''));
    close(scoreAlignment(p.aligned1, p.aligned2, sch), p.score, `${name} pair ${k} python rescore`);
    n++;
  });
  for (const [key, m] of [['tree_upgma', 'UPGMA'], ['tree_nj', 'NJ']]) {
    const t = treeBundle(c.fasta, m), py = c[key];
    close(t.distances, py.distances, `${name} ${m} distances`);
    close(t.tree, py.tree, `${name} ${m} tree`);
    n++;
  }
}
assert.equal(globalAlign('ACDEFGHIK', 'ACDGHIKKK', sch).score.toFixed(6), '5.800000');
assert.equal(globalAlign('', 'AAA', sch).aligned1, '---');
console.log(`align_engine: ${n} fixture comparisons ok`);
