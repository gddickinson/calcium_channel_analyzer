// Pairwise alignment and distance trees ported from alignment_tree.py (pure functions, no DOM).
// - globalAlign: Gotoh global alignment with affine gaps, the scores of the app's Bio.Align.PairwiseAligner
//   (match 1, mismatch 0, gap open -0.5, gap extend -0.1; a gap of length L scores open + (L-1)*extend, end gaps too).
//   The optimal score matches Biopython; when several alignments are optimal the one shown may differ from Biopython's.
// - identityDistances: Bio.Phylo DistanceCalculator('identity') on '-'-padded sequences ('-' == '-' counts as a match).
// - upgma / nj: line-by-line ports of Bio.Phylo DistanceTreeConstructor (same pair choice, ties, branch lengths).
// align_engine.test.mjs checks these against tests/fixtures/engine_fixture.json.
import { parseFasta } from './seq_engine.js';

const NEG = -Infinity;

export function globalAlign(a, b, { match = 1, mismatch = 0, open = -0.5, extend = -0.1 } = {}) {
  const n = a.length, m = b.length, W = m + 1, size = (n + 1) * W;
  const M = new Float64Array(size).fill(NEG), X = new Float64Array(size).fill(NEG), Y = new Float64Array(size).fill(NEG);
  // traceback: previous state (0 = M, 1 = X (gap in b), 2 = Y (gap in a)) for each state
  const tM = new Uint8Array(size), tX = new Uint8Array(size), tY = new Uint8Array(size);
  M[0] = 0;
  for (let i = 1; i <= n; i++) { X[i * W] = open + (i - 1) * extend; tX[i * W] = i === 1 ? 0 : 1; }
  for (let j = 1; j <= m; j++) { Y[j] = open + (j - 1) * extend; tY[j] = j === 1 ? 0 : 2; }
  const best3 = (vM, vX, vY) => (vM >= vX && vM >= vY) ? 0 : (vX >= vY ? 1 : 2);
  for (let i = 1; i <= n; i++) {
    for (let j = 1; j <= m; j++) {
      const c = i * W + j, d = c - W - 1, u = c - W, l = c - 1;
      const s = a[i - 1] === b[j - 1] ? match : mismatch;
      let k = best3(M[d], X[d], Y[d]); M[c] = [M[d], X[d], Y[d]][k] + s; tM[c] = k;
      const xo = [M[u] + open, X[u] + extend, Y[u] + open]; k = best3(...xo); X[c] = xo[k]; tX[c] = k;
      const yo = [M[l] + open, X[l] + open, Y[l] + extend]; k = best3(...yo); Y[c] = yo[k]; tY[c] = k;
    }
  }
  const end = n * W + m;
  let state = best3(M[end], X[end], Y[end]);
  const score = [M[end], X[end], Y[end]][state];
  let i = n, j = m, top = [], bot = [];
  while (i > 0 || j > 0) {
    const c = i * W + j;
    if (state === 0) { top.push(a[i - 1]); bot.push(b[j - 1]); state = tM[c]; i--; j--; }
    else if (state === 1) { top.push(a[i - 1]); bot.push('-'); state = tX[c]; i--; }
    else { top.push('-'); bot.push(b[j - 1]); state = tY[c]; j--; }
  }
  return { score, aligned1: top.reverse().join(''), aligned2: bot.reverse().join('') };
}

// Score an existing gapped pair with the same scheme (used by the tests and as a self-check).
export function scoreAlignment(t, u, { match = 1, mismatch = 0, open = -0.5, extend = -0.1 } = {}) {
  let s = 0, prev = 0;
  for (let k = 0; k < t.length; k++) {
    const st = t[k] === '-' ? 2 : u[k] === '-' ? 1 : 0;
    s += st === 0 ? (t[k] === u[k] ? match : mismatch) : (st === prev ? extend : open);
    prev = st;
  }
  return s;
}

export function pairwiseAlignments(ids, seqs, cfg) {
  const p = { match: cfg.match, mismatch: cfg.mismatch, open: cfg.gap_open, extend: cfg.gap_extend }, out = [];
  for (let i = 0; i < seqs.length; i++) for (let j = i + 1; j < seqs.length; j++) {
    const r = globalAlign(seqs[i].slice(0, cfg.prefix), seqs[j].slice(0, cfg.prefix), p);
    let same = 0;
    for (let k = 0; k < r.aligned1.length; k++) if (r.aligned1[k] === r.aligned2[k] && r.aligned1[k] !== '-') same++;
    out.push({ i, j, seq1: ids[i], seq2: ids[j], score: r.score, aligned1: r.aligned1, aligned2: r.aligned2,
      identity: r.aligned1.length ? same / r.aligned1.length : 0 });
  }
  return out;
}

export function identityDistances(seqs) {
  const L = Math.max(...seqs.map(s => s.length)), pad = seqs.map(s => s + '-'.repeat(L - s.length));
  return pad.map(a => pad.map(b => {
    let same = 0;
    for (let k = 0; k < L; k++) if (a[k] === b[k]) same++;
    return L === 0 ? 1 : 1 - same / L;
  }));
}

const leaf = name => ({ name, branch_length: 0, children: [] });
const height = c => c.children.length ? Math.max(...c.children.map(k => height(k) + k.branch_length)) : 0;
const dropIndex = (dm, i) => { dm.splice(i, 1); dm.forEach(row => row.splice(i, 1)); };

export function upgma(ids, dist) {
  const dm = dist.map(r => r.slice()), clades = ids.map(leaf);
  let mi = 0, mj = 0, count = 0, inner = null;
  while (dm.length > 1) {
    let min = dm[1][0];
    for (let i = 1; i < dm.length; i++) for (let j = 0; j < i; j++) if (min >= dm[i][j]) { min = dm[i][j]; mi = i; mj = j; }
    const c1 = clades[mi], c2 = clades[mj];
    inner = { name: 'Inner' + (++count), branch_length: 0, children: [c1, c2] };
    c1.branch_length = min * 1.0 / 2 - height(c1);
    c2.branch_length = min * 1.0 / 2 - height(c2);
    clades[mj] = inner; clades.splice(mi, 1);
    for (let k = 0; k < dm.length; k++) if (k !== mi && k !== mj) { const v = (dm[mi][k] + dm[mj][k]) / 2; dm[mj][k] = v; dm[k][mj] = v; }
    dropIndex(dm, mi);
  }
  inner.branch_length = 0;
  return inner;
}

export function nj(ids, dist) {
  const dm = dist.map(r => r.slice()), clades = ids.map(leaf);
  if (dm.length === 1) return clades[0];
  if (dm.length === 2) {
    const c1 = clades[1], c2 = clades[0];
    c1.branch_length = dm[1][0] / 2.0; c2.branch_length = dm[1][0] - c1.branch_length;
    return { name: 'Inner', branch_length: 0, children: [c1, c2] };
  }
  let count = 0, inner = null;
  while (dm.length > 2) {
    const n = dm.length, nd = dm.map(row => { let s = 0; for (let j = 0; j < n; j++) s += row[j]; return s / (n - 2); });
    let min = dm[1][0] - nd[1] - nd[0], mi = 0, mj = 1;
    for (let i = 1; i < n; i++) for (let j = 0; j < i; j++) {
      const t = dm[i][j] - nd[i] - nd[j];
      if (min > t) { min = t; mi = i; mj = j; }
    }
    const c1 = clades[mi], c2 = clades[mj];
    inner = { name: 'Inner' + (++count), branch_length: 0, children: [c1, c2] };
    c1.branch_length = (dm[mi][mj] + nd[mi] - nd[mj]) / 2.0;
    c2.branch_length = dm[mi][mj] - c1.branch_length;
    clades[mj] = inner; clades.splice(mi, 1);
    for (let k = 0; k < n; k++) if (k !== mi && k !== mj) {
      const v = (dm[mi][k] + dm[mj][k] - dm[mi][mj]) / 2.0; dm[mj][k] = v; dm[k][mj] = v;
    }
    dropIndex(dm, mi);
  }
  if (clades[0] === inner) {
    clades[0].branch_length = 0; clades[1].branch_length = dm[1][0]; clades[0].children.push(clades[1]); return clades[0];
  }
  clades[0].branch_length = dm[1][0]; clades[1].branch_length = 0; clades[1].children.push(clades[0]); return clades[1];
}

// Newick in Bio.Phylo's style (branch lengths %1.5f, inner names kept).
export function newick(c, root = true) {
  const name = /[\s(),:;]/.test(c.name) ? `'${c.name}'` : c.name;
  const body = c.children.length ? `(${c.children.map(k => newick(k, false)).join(',')})${name}` : name;
  return body + ':' + c.branch_length.toFixed(5) + (root ? ';' : '');
}

// --- bundles with the same shape as web/analysis.py -----------------------------------------------------------
export function alignBundle(fasta, cfg) {
  const { records, warnings } = parseFasta(fasta), a = cfg.align;
  if (records.length < 2) throw new Error('Need at least 2 sequences for alignment');
  if (records.length > a.max_sequences) throw new Error(`At most ${a.max_sequences} sequences for all-pairs alignment`);
  const ids = records.map(r => r.id);
  return { ids, prefix: a.prefix, pairs: pairwiseAlignments(ids, records.map(r => r.sequence), a), warnings };
}
export function treeBundle(fasta, method = 'UPGMA') {
  const { records, warnings } = parseFasta(fasta);
  method = String(method).toUpperCase();
  if (!['UPGMA', 'NJ'].includes(method)) throw new Error('method must be UPGMA or NJ');
  if (records.length < 3) throw new Error('Need at least 3 sequences for tree construction');
  const ids = records.map(r => r.id);
  if (new Set(ids).size !== ids.length) throw new Error('Sequence ids must be unique to build a tree');
  const distances = identityDistances(records.map(r => r.sequence));
  const tree = method === 'UPGMA' ? upgma(ids, distances) : nj(ids, distances);
  return { ids, method, distances, tree, newick: newick(tree), warnings };
}
