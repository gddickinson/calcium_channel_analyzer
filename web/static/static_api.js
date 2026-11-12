// Static-site data layer: GET paths come from the JSON exported by scripts/build_static.py; the analyses (POST paths)
// run in the browser with seq_engine.js / align_engine.js, which node tests check against the Python code.
// Fetching from NCBI / UniProt exists only in the local app.
import { alignBundle, treeBundle } from './align_engine.js';
import { classifyBundle, domainsBundle, hydropathyBundle, parseBundle } from './seq_engine.js';

const files = new Map();
async function file(rel) {
  if (!files.has(rel)) files.set(rel, fetch('api/' + rel + '.json').then(r => { if (!r.ok) throw new Error(`${r.status} ${rel}`); return r.json(); }));
  return files.get(rel);
}

export const NEEDS_SERVER = 'needs the full app running locally (Python server): see the README, "Run the web app"';

export async function staticApi(path, body) {
  const clean = path.split('?')[0].replace(/^\/+/, '');
  if (!body) return file(clean);
  const cfg = await file('config');
  switch (clean) {
    case 'parse': return parseBundle(body.fasta);
    case 'domains': return domainsBundle(body.fasta, cfg, !!body.use_min_length);
    case 'hydropathy': return hydropathyBundle(body.fasta, cfg, body.index ?? 0, body.window ?? cfg.tm.window, body.threshold ?? cfg.tm.threshold);
    case 'classify': return classifyBundle(body.fasta, cfg);
    case 'align': return alignBundle(body.fasta, cfg);
    case 'tree': return treeBundle(body.fasta, body.method);
    case 'fetch': throw new Error('Fetching from NCBI or UniProt ' + NEEDS_SERVER + '.');
    default: throw new Error('Unknown API path ' + clean);
  }
}
