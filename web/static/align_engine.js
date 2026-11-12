// Pairwise alignment and distance trees ported from alignment_tree.py (pure functions, no DOM).
// - globalAlign: Gotoh global alignment with affine gaps (a gap of length L scores open + (L-1)*extend, end gaps too),
//   with either match/mismatch scores (the app's alignment page: 1 / 0 / -0.5 / -0.1) or a substitution matrix
//   (the tree: BLOSUM62, -10 / -0.5). As in Biopython, scores within epsilon 1e-6 are ties and ties go to M, then X,
//   then Y, which reproduces Biopython's first optimal alignment (alignments[0]); the fixture tests compare the aligned strings exactly.
//   Scores are kept in rolling rows; the traceback is one byte per cell (2 bits per state).
// - alignmentDistances: 1 - identity (identical columns / alignment length) for every pair.
// - upgma / nj: line-by-line ports of Bio.Phylo DistanceTreeConstructor (same pair choice, ties, branch lengths).
// align_engine.test.mjs checks these against tests/fixtures/engine_fixture.json.
import { parseFasta } from './seq_engine.js';

const NEG = -Infinity;
export const MAX_CELLS = 40000000;
const EPS = 1e-6;  // Biopython PairwiseAligner.epsilon: scores within this are ties

// Scoring table: flat K x K Float64Array plus a letter -> index map (letters not in the alphabet map to 'X').
export function matrixScorer(matrix) {
  const K = matrix.alphabet.length, S = new Float64Array(K * K), idx = new Map();
  [...matrix.alphabet].forEach((c, i) => idx.set(c, i));
  matrix.scores.forEach((row, i) => row.forEach((v, j) => { S[i * K + j] = v; }));
  const x = idx.get('X');
  return { K, S, code: seq => Int32Array.from(seq, c => idx.has(c) ? idx.get(c) : x), letters: matrix.alphabet };
}
function matchScorer(seqA, seqB, match, mismatch) {
  const letters = [...new Set(seqA + seqB)].join('') || 'A', K = letters.length, S = new Float64Array(K * K).fill(mismatch);
  for (let i = 0; i < K; i++) S[i * K + i] = match;
  const idx = new Map([...letters].map((c, i) => [c, i]));
  return { K, S, code: seq => Int32Array.from(seq, c => idx.get(c)), letters };
}

export function globalAlign(a, b, { match = 1, mismatch = 0, open = -0.5, extend = -0.1, scorer = null } = {}) {
  const sc = scorer || matchScorer(a, b, match, mismatch);
  const n = a.length, m = b.length, W = m + 1;
  if (n * m > MAX_CELLS) throw new Error(`Sequences too long to align (${n} x ${m} > ${MAX_CELLS} cells)`);
  const A = sc.code(a), B = sc.code(b), K = sc.K, S = sc.S;
  let pM = new Float64Array(W).fill(NEG), pX = new Float64Array(W).fill(NEG), pY = new Float64Array(W).fill(NEG);
  let cM = new Float64Array(W), cX = new Float64Array(W), cY = new Float64Array(W);
  const tb = new Uint8Array((n + 1) * W);  // bits 0-1: M's previous state, 2-3: X's, 4-5: Y's (0 M, 1 X, 2 Y)
  pM[0] = 0;
  for (let j = 1; j <= m; j++) { pY[j] = open + (j - 1) * extend; tb[j] = (j === 1 ? 0 : 2) << 4; }
  for (let i = 1; i <= n; i++) {
    cM[0] = NEG; cY[0] = NEG; cX[0] = open + (i - 1) * extend; tb[i * W] = (i === 1 ? 0 : 1) << 2;
    const row = A[i - 1] * K;
    for (let j = 1; j <= m; j++) {
      let v, k, t;
      // M: diagonal from (i-1, j-1)
      v = pM[j - 1]; k = 0;
      if (pX[j - 1] > v + EPS) { v = pX[j - 1]; k = 1; }
      if (pY[j - 1] > v + EPS) { v = pY[j - 1]; k = 2; }
      cM[j] = v + S[row + B[j - 1]]; t = k;
      // X: gap in b, from (i-1, j)
      v = pM[j] + open; k = 0;
      if (pX[j] + extend > v + EPS) { v = pX[j] + extend; k = 1; }
      if (pY[j] + open > v + EPS) { v = pY[j] + open; k = 2; }
      cX[j] = v; t |= k << 2;
      // Y: gap in a, from (i, j-1)
      v = cM[j - 1] + open; k = 0;
      if (cX[j - 1] + open > v + EPS) { v = cX[j - 1] + open; k = 1; }
      if (cY[j - 1] + extend > v + EPS) { v = cY[j - 1] + extend; k = 2; }
      cY[j] = v; tb[i * W + j] = t | (k << 4);
    }
    [pM, cM] = [cM, pM]; [pX, cX] = [cX, pX]; [pY, cY] = [cY, pY];
  }
  const ends = [pM[m], pX[m], pY[m]];
  let state = 0;
  if (n === 0 && m === 0) state = 0;
  else { if (ends[1] > ends[state] + EPS) state = 1; if (ends[2] > ends[state] + EPS) state = 2; }
  const score = ends[state];
  let i = n, j = m;
  const top = [], bot = [];
  while (i > 0 || j > 0) {
    const t = tb[i * W + j];
    if (state === 0) { top.push(a[i - 1]); bot.push(b[j - 1]); state = t & 3; i--; j--; }
    else if (state === 1) { top.push(a[i - 1]); bot.push('-'); state = (t >> 2) & 3; i--; }
    else { top.push('-'); bot.push(b[j - 1]); state = (t >> 4) & 3; j--; }
  }
  return { score, aligned1: top.reverse().join(''), aligned2: bot.reverse().join('') };
}

const identityOf = (t, u) => {
  let same = 0;
  for (let k = 0; k < t.length; k++) if (t[k] === u[k] && t[k] !== '-') same++;
  return t.length ? same / t.length : 0;
};

// alignment_tree.alignment_distances: residues missing from the matrix are aligned as X (as in Python).
export function alignmentDistances(ids, seqs, tcfg) {
  const sc = matrixScorer(tcfg.matrix), n = seqs.length, d = ids.map(() => Array(n).fill(0)), pairs = [];
  const clean = seqs.map(s => [...s].map(c => tcfg.matrix.alphabet.includes(c) ? c : 'X').join(''));
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
    if (clean[i].length * clean[j].length > MAX_CELLS) throw new Error(`${ids[i]} x ${ids[j]} is too long to align (limit ${MAX_CELLS} cells)`);
    const r = globalAlign(clean[i], clean[j], { open: tcfg.gap_open, extend: tcfg.gap_extend, scorer: sc });
    const identity = identityOf(r.aligned1, r.aligned2);
    d[i][j] = d[j][i] = 1.0 - identity;
    pairs.push({ i, j, score: r.score, identity, length: r.aligned1.length, aligned1: r.aligned1, aligned2: r.aligned2 });
  }
  return { distances: d, pairs };
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
    out.push({ i, j, seq1: ids[i], seq2: ids[j], score: r.score, aligned1: r.aligned1, aligned2: r.aligned2,
      identity: identityOf(r.aligned1, r.aligned2) });
  }
  return out;
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
let lastTree = null;  // the alignments are the slow part: reuse them when only the method changes
export function treeBundle(fasta, cfg, method = 'UPGMA') {
  const { records, warnings } = parseFasta(fasta), t = cfg.tree;
  method = String(method).toUpperCase();
  if (!['UPGMA', 'NJ'].includes(method)) throw new Error('method must be UPGMA or NJ');
  if (records.length > t.max_sequences) throw new Error(`At most ${t.max_sequences} sequences for a tree`);
  if (records.length < 3) throw new Error('Need at least 3 sequences for tree construction');
  const ids = records.map(r => r.id);
  if (new Set(ids).size !== ids.length) throw new Error('Sequence ids must be unique to build a tree');
  const key = JSON.stringify([ids, records.map(r => r.sequence)]);
  if (!lastTree || lastTree.key !== key) lastTree = { key, ...alignmentDistances(ids, records.map(r => r.sequence), t) };
  const { distances, pairs } = lastTree;
  const tree = method === 'UPGMA' ? upgma(ids, distances) : nj(ids, distances);
  return { ids, method, distances, pairs: pairs.map(({ i, j, score, identity, length }) => ({ i, j, score, identity, length })),
    tree, newick: newick(tree), warnings };
}
