// Sequence analyses ported from the repo's Python (pure functions, no DOM): FASTA parsing (sequence_analysis.py),
// domain detection and scoring (domain_detection.py), ML features (ml_features.py), Kyte-Doolittle hydropathy and TM
// regions (sequence_analysis.py), rule-based classification and electrophysiology labels (channel_rules.py).
// Constants (domain patterns, KD scale, amino acids) come from config.py via api('/config'): pass that object as cfg.
// seq_engine.test.mjs checks every function against tests/fixtures/engine_fixture.json (made by the Python code).

export const MAX_SEQUENCES = 200, MAX_RESIDUES = 500000;

// Python's f"{x:.1f}": round half to even on exact binary ties (only x = k.25 or k.75 can tie), else toFixed.
export function fmt1(x) {
  if (Number.isInteger(x * 4) && !Number.isInteger(x * 2)) {
    const s = Math.floor(x * 10);                    // k.25 -> k2, k.75 -> k7 (negative values not used)
    return ((s % 2 === 0 ? s : s + 1) / 10).toFixed(1);
  }
  return x.toFixed(1);
}

// Python 3.12+ sum() of floats: Neumaier-compensated (CPython builtin_sum_impl), so window sums match the Python
// code to the last bit (a plain running sum differs in the last digit and can flip a "score > threshold" test).
export function pySum(values) {
  let f = 0, c = 0;
  for (const x of values) {
    const t = f + x;
    c += Math.abs(f) >= Math.abs(x) ? (f - t) + x : (x - t) + f;
    f = t;
  }
  return (c && Number.isFinite(c)) ? f + c : f;
}

// sequence_analysis.parse_fasta -> { records, warnings }; throws Error with the same messages.
export function parseFasta(text) {
  const records = [], warnings = [];
  let header = null, chunks = [];
  const lines = String(text || '').replace(/\r/g, '').split('\n');
  if (!lines.some(l => l.trim())) throw new Error("No sequences: paste FASTA text starting with a '>' header line.");
  const finish = () => {
    if (header === null) return;
    const raw = chunks.join('').replace(/[\s\d]/g, '').toUpperCase().replace(/\*+$/, '');
    const bad = [...new Set(raw.replace(/[A-Z]/g, ''))].sort();
    const seq = raw.replace(/[^A-Z]/g, '');
    const words = header.split(/\s+/).filter(Boolean);
    const id = words.length ? words[0] : `seq${records.length + 1}`;
    if (bad.length) warnings.push(`${id}: removed characters ${bad.join(' ')}`);
    if (!seq) { warnings.push(`${id}: empty sequence skipped`); return; }
    records.push({ id, description: header, sequence: seq });
  };
  lines.forEach((l, k) => {
    const line = l.trim();
    if (!line) return;
    if (line.startsWith('>')) { finish(); header = line.slice(1).trim(); chunks = []; }
    else if (header === null) throw new Error(`Line ${k + 1}: sequence data before the first '>' header line.`);
    else chunks.push(line);
  });
  finish();
  if (!records.length) throw new Error('No usable sequences found.');
  if (records.length > MAX_SEQUENCES) throw new Error(`Too many sequences (${records.length}); the limit is ${MAX_SEQUENCES}.`);
  const total = records.reduce((a, r) => a + r.sequence.length, 0);
  if (total > MAX_RESIDUES) throw new Error(`Too many residues (${total}); the limit is ${MAX_RESIDUES}.`);
  const seen = new Map();
  records.forEach(r => seen.set(r.id, (seen.get(r.id) || 0) + 1));
  const dups = [...seen].filter(([, v]) => v > 1).map(([k]) => k);
  if (dups.length) warnings.push('Duplicate ids: ' + dups.join(', '));
  return { records, warnings };
}

// domain_detection.calculate_domain_score
export function domainScore(match, name, patterns) {
  const lengthScore = Math.min(match.length / patterns[name].min_length, 1.0) * 50;
  let hydro = 25;
  if (name === 'TM_helix') {
    let c = 0;
    for (const aa of match) if ('AILVMFGW'.includes(aa)) c++;
    hydro = (c / match.length) * 50;
  }
  return lengthScore + hydro;
}

// domain_detection.detect_domains_in_sequence: every non-overlapping regex match, unless minLength or
// usePatternMinLength asks for a minimum; positions 1-based and inclusive.
export function detectDomains(seq, id, patterns, { minLength = null, usePatternMinLength = false } = {}) {
  const out = [];
  for (const [name, info] of Object.entries(patterns)) {
    for (const m of seq.matchAll(new RegExp(info.pattern, 'g'))) {
      const s = m[0];
      if (minLength !== null && s.length < minLength) continue;
      if (usePatternMinLength && s.length < (info.min_length ?? 1)) continue;
      const start = m.index + 1, end = m.index + s.length;
      out.push({ sequence: id, domain: name, position: `${start}-${end}`, start, end, match: s,
        score: fmt1(domainScore(s, name, patterns)) });
    }
  }
  return out;
}

const FAMILY_KEYWORDS = [
  ['Cav1 (L-type)', ['cav1', 'l-type', 'cacna1s', 'cacna1c', 'cacna1d', 'cacna1f', 'dihydropyridine']],
  ['Cav2 (P/Q/N/R-type)', ['cav2', 'p/q-type', 'n-type', 'r-type', 'cacna1a', 'cacna1b', 'cacna1e']],
  ['Cav3 (T-type)', ['cav3', 't-type', 'cacna1g', 'cacna1h', 'cacna1i']],
  ['RyR', ['ryanodine', 'ryr']], ['IP3R', ['ip3', 'inositol', 'itpr']],
  ['CNGC', ['cyclic nucleotide', 'cngc', 'cng channel']], ['GLR', ['glutamate receptor-like', 'glr', 'iglur']],
  ['TPC', ['two-pore', 'tpc', 'tpcn']], ['TRP', ['transient receptor', 'trp', 'trpc', 'trpv', 'trpm']],
  ['Orai', ['orai', 'crac', 'stim']], ['P2X', ['p2x', 'purinergic']], ['MSL', ['mechanosensitive', 'msl', 'mscs']],
  ['MCA', ['mid1-complementing', 'mca']],
];
// domain_detection.predict_family_from_description
export function familyFromDescription(desc) {
  const d = desc.toLowerCase();
  for (const [fam, kws] of FAMILY_KEYWORDS) if (kws.some(k => d.includes(k))) return fam;
  return 'Unknown';
}

// ml_features.extract_ml_features: [length, 20 composition fractions, mean KD hydropathy, domain presence flags]
export function mlFeatures(seq, cfg) {
  const n = seq.length, f = [n];
  for (const aa of cfg.amino_acids) f.push(n > 0 ? (seq.split(aa).length - 1) / n : 0.0);
  const h = pySum([...seq].map(aa => cfg.kd_hydrophobicity[aa] ?? 0));
  f.push(n > 0 ? h / n : 0.0);
  for (const info of Object.values(cfg.domain_patterns)) f.push(new RegExp(info.pattern).test(seq) ? 1.0 : 0.0);
  return f;
}

// sequence_analysis.hydropathy_profile / predict_tm_regions / tm_prediction
export function hydropathyProfile(seq, window, kd) {
  const positions = [], scores = [];
  for (let i = 0; i + window <= seq.length; i++) {
    const w = new Array(window);
    for (let k = 0; k < window; k++) w[k] = kd[seq[i + k]] ?? 0;
    scores.push(pySum(w) / window); positions.push(i + Math.floor(window / 2));
  }
  return { positions, scores };
}
export function tmRegions(positions, scores, threshold) {
  const regions = [];
  let inTm = false, start = 0;
  scores.forEach((s, i) => {
    if (s > threshold && !inTm) { start = positions[i]; inTm = true; }
    else if (s <= threshold && inTm) { regions.push([start, positions[i]]); inTm = false; }
  });
  if (inTm) regions.push([start, positions[positions.length - 1]]);
  return regions;
}
export function tmPrediction(seq, window, threshold, kd) {
  const { positions, scores } = hydropathyProfile(seq, window, kd);
  return { length: seq.length, window, threshold, positions, scores,
    regions: tmRegions(positions, scores, threshold).map(([s, e]) => ({ start: s, end: e, length: e - s, sequence: seq.slice(s, e) })) };
}

// channel_rules.classify_channel / predict_electrophysiology
const has = (name, seq, pats) => new RegExp(pats[name].pattern).test(seq);
const any = (d, terms) => terms.some(t => d.includes(t));
export function markers(seq, pats) {
  return { voltage_sensor: has('Voltage_sensor', seq, pats), ca_selectivity: has('Ca_selectivity', seq, pats),
    cnbd: has('CNBD', seq, pats), spry: has('SPRY', seq, pats), ip3_binding: has('IP3_binding', seq, pats) };
}
export function classifyChannel(desc, seq, pats) {
  const d = desc.toLowerCase(), m = markers(seq, pats);
  const rules = [
    [any(d, ['cav1', 'l-type']), 'Voltage-gated (Cav1/L-type)', 'High'],
    [any(d, ['cav2', 'p/q', 'n-type']), 'Voltage-gated (Cav2/P-Q-N-type)', 'High'],
    [any(d, ['cav3', 't-type']), 'Voltage-gated (Cav3/T-type)', 'High'],
    [any(d, ['ryanodine', 'ryr']), 'Ligand-gated (RyR)', 'High'],
    [d.includes('ip3'), 'Ligand-gated (IP3R)', 'High'],
    [d.includes('cngc'), 'Plant CNGC', 'High'],
    [any(d, ['tpc', 'two-pore']), 'Two-pore channel (TPC)', 'High'],
    [m.voltage_sensor && m.ca_selectivity, 'Voltage-gated Ca2+ channel', 'Medium'],
    [m.spry, 'Ligand-gated (likely RyR)', 'Medium'],
    [m.ip3_binding, 'Ligand-gated (likely IP3R)', 'Medium'],
    [m.cnbd, 'CNGC-like channel', 'Medium'],
  ];
  const hit = rules.find(r => r[0]) || [true, 'Unknown', 'Low'];
  return { classification: hit[1], confidence: hit[2], length: seq.length, markers: m };
}
export function electrophysiology(desc, seq, pats) {
  const d = desc.toLowerCase();
  const [gating, activation] =
    has('Voltage_sensor', seq, pats) ? ['Voltage-gated', 'Depolarization-activated (-40 to +20 mV typical)']
    : has('CNBD', seq, pats) ? ['Ligand-gated (cAMP/cGMP)', 'cAMP/cGMP binding']
    : has('GLR_ligand', seq, pats) ? ['Ligand-gated (Glutamate-like)', 'Glutamate or analog binding']
    : any(d, ['ryr', 'ryanodine']) ? ['Ligand-gated (Ca²⁺, cADPR, ryanodine)', 'Ca²⁺-induced Ca²⁺ release (CICR)']
    : d.includes('ip3') ? ['Ligand-gated (IP3)', 'IP3 binding']
    : ['Unknown/Mixed', 'Multiple mechanisms possible'];
  const [selectivity, permeability] =
    has('Ca_selectivity', seq, pats) ? ['Ca²⁺-selective (PCa/PNa > 100)', 'High Ca²⁺, Low Na⁺/K⁺']
    : any(d, ['tpc', 'two-pore']) ? ['Non-selective cation (Ca²⁺, Na⁺, K⁺)', 'PCa:PNa:PK ≈ 1:1:1']
    : ['Ca²⁺-permeable', 'Variable Ca²⁺/monovalent selectivity'];
  const blockers =
    any(d, ['cav1', 'l-type']) ? 'Dihydropyridines (nifedipine), phenylalkylamines (verapamil), benzothiazepines (diltiazem)'
    : any(d, ['cav2.1', 'p/q']) ? 'ω-agatoxin IVA'
    : any(d, ['cav2.2', 'n-type']) ? 'ω-conotoxin GVIA'
    : any(d, ['cav3', 't-type']) ? 'Mibefradil, TTA-P2'
    : d.includes('ryr') ? 'Ryanodine (high dose), dantrolene'
    : 'Family-specific blockers required';
  return { gating, activation, selectivity, permeability, blockers };
}

// --- bundles with the same shape as web/analysis.py -----------------------------------------------------------
export function parseBundle(fasta) {
  const { records, warnings } = parseFasta(fasta);
  records.forEach(r => { r.length = r.sequence.length; r.family = familyFromDescription(r.description); });
  return { records, warnings };
}
export function domainsBundle(fasta, cfg, useMinLength = false) {
  const { records, warnings } = parseFasta(fasta);
  return { use_min_length: useMinLength, lengths: Object.fromEntries(records.map(r => [r.id, r.sequence.length])), order: records.map(r => r.id),
    hits: records.flatMap(r => detectDomains(r.sequence, r.id, cfg.domain_patterns, { usePatternMinLength: useMinLength })), warnings };
}
export function hydropathyBundle(fasta, cfg, index = 0, window = cfg.tm.window, threshold = cfg.tm.threshold) {
  const { records, warnings } = parseFasta(fasta);
  if (!(index >= 0 && index < records.length)) throw new Error(`index must be between 0 and ${records.length - 1}`);
  if (!(window >= 3 && window <= 51)) throw new Error('window must be between 3 and 51');
  const r = records[index];
  return { id: r.id, description: r.description, ...tmPrediction(r.sequence, window, threshold, cfg.kd_hydrophobicity), warnings };
}
export function classifyBundle(fasta, cfg) {
  const { records, warnings } = parseFasta(fasta), p = cfg.domain_patterns;
  const rows = records.map(r => ({ id: r.id, description: r.description,
    family_from_description: familyFromDescription(r.description), ...classifyChannel(r.description, r.sequence, p),
    electrophysiology: electrophysiology(r.description, r.sequence, p), features: mlFeatures(r.sequence, cfg) }));
  return { rows, feature_names: cfg.feature_names, warnings };
}
