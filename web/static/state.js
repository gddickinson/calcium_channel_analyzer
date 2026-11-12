// The working sequence set shared by all pages: FASTA text, defaulting to the repo's example file.
// Kept for this browser tab in sessionStorage (optional: everything works without it).
import { api, esc, isStatic } from './core.js';

const KEY = 'cca-fasta';
let fasta = null, source = null;

function load() {
  try { const v = sessionStorage.getItem(KEY); if (v) { const o = JSON.parse(v); fasta = o.fasta; source = o.source; } } catch (e) { /* storage off */ }
}
load();

export async function examples() { return api('/examples'); }

export async function getFasta() {
  if (fasta) return fasta;
  const ex = await examples();
  fasta = ex.fasta; source = `example file ${ex.file}`;
  return fasta;
}
export const getSource = () => source || 'example file';

export function setFasta(text, src) {
  fasta = text; source = src;
  try { sessionStorage.setItem(KEY, JSON.stringify({ fasta, source })); } catch (e) { /* storage off */ }
}
export async function resetFasta() {
  try { sessionStorage.removeItem(KEY); } catch (e) { /* storage off */ }
  fasta = null; source = null;
  return getFasta();
}

// Run an analysis on the current set. Local app: Python on this computer; static site: the JS port in the browser.
export async function analyse(path, extra = {}) {
  return api('/' + path, { body: { fasta: await getFasta(), ...extra } });
}

export async function seqbar(count) {
  await getFasta();
  const where = isStatic() ? 'computed in your browser (JS port, tested against the Python code)' : 'computed by the Python code on this computer';
  return `<div class="seqbar"><span>Sequence set: <b>${esc(getSource())}</b>${count != null ? ` · ${count} sequences` : ''}</span>
    ${getSource().startsWith('example') ? '<span>(placeholder sequences, <a href="#/about">see About</a>)</span>' : ''}
    <a href="#/sequences">change</a><span>· ${where}</span></div>`;
}

export const warningsHtml = w => (w && w.length) ? `<div class="warnings">${w.map(esc).join('<br>')}</div>` : '';
