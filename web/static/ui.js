// Shared page helpers: domain colours (the desktop app's plot_domain_map palette), number formats, small SVG bits.
import { esc } from './core.js';

export const DOMAIN_COLORS = {
  EF_hand: '#FF6B6B', CNBD: '#4ECDC4', Pore_loop: '#45B7D1', SPRY: '#FFA07A', IP3_binding: '#98D8C8',
  Voltage_sensor: '#F7DC6F', TM_helix: '#BB8FCE', Ca_selectivity: '#85C1E2', Calmodulin_binding: '#F8B88B', GLR_ligand: '#FAD7A0',
};
export const color = d => DOMAIN_COLORS[d] || '#cccccc';
export const f = (x, d = 2) => Number.isFinite(x) ? x.toFixed(d) : String(x);
export const short = (id, n = 22) => id.length > n ? id.slice(0, n - 1) + '…' : id;

// Short display name for UniProt-style ids such as sp|Q13936|CACNA1C_HUMAN -> CACNA1C_HUMAN
export const label = id => { const p = id.split('|'); return p.length === 3 && p[2] ? p[2] : id; };

export function svg(w, h, body, aria) {
  return `<svg viewBox="0 0 ${w} ${h}" role="img" aria-label="${esc(aria)}" preserveAspectRatio="xMinYMin meet">${body}</svg>`;
}

// Diverging-free sequential ramp for heat cells: t in [0, 1] -> background colour on the dark theme.
export function heat(t) {
  const a = Math.max(0, Math.min(1, t));
  return `background:rgba(55,232,255,${(0.06 + 0.5 * a).toFixed(3)})`;
}

// Fixed-width blocks of an alignment with a match line.
export function alignmentText(a, b, idA, idB, width = 60) {
  const la = label(idA).padEnd(14).slice(0, 14), lb = label(idB).padEnd(14).slice(0, 14), out = [];
  for (let k = 0; k < a.length; k += width) {
    const x = a.slice(k, k + width), y = b.slice(k, k + width);
    const mid = [...x].map((c, i) => c === y[i] && c !== '-' ? '|' : (c === '-' || y[i] === '-') ? ' ' : '.').join('');
    out.push(`${la} ${x}\n${' '.repeat(15)}${mid}\n${lb} ${y}`);
  }
  return out.join('\n\n');
}
