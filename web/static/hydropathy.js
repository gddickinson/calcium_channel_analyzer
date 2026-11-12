// Hydropathy & TM: the desktop app's "Transmembrane Prediction" (Kyte-Doolittle, window 19, threshold 1.6 by default),
// sequence_analysis.tm_prediction, with an SVG profile and the predicted regions.
import { $, api, app, esc, route } from './core.js';
import { analyse, seqbar, warningsHtml } from './state.js';
import { f, label, svg } from './ui.js';

function plot(r) {
  const W = 1000, H = 300, l = 50, rt = 12, t = 12, b = 36;
  if (!r.scores.length) return `<p class="mute">The sequence (${r.length} aa) is shorter than the window (${r.window}): no profile.</p>`;
  const lo = Math.min(-4.5, ...r.scores), hi = Math.max(4.5, ...r.scores), xmax = Math.max(r.length, 1);
  const sx = x => l + (x / xmax) * (W - l - rt), sy = y => t + (hi - y) / (hi - lo) * (H - t - b);
  let body = '';
  for (const g of r.regions) body += `<rect x="${sx(g.start)}" y="${t}" width="${Math.max(1, sx(g.end) - sx(g.start))}" height="${H - t - b}" fill="rgba(255,177,61,.18)"/>`;
  for (let v = Math.ceil(lo); v <= Math.floor(hi); v++) if (v % 2 === 0) {
    body += `<line x1="${l}" x2="${W - rt}" y1="${sy(v)}" y2="${sy(v)}" stroke="#212c4c"/><text x="${l - 6}" y="${sy(v) + 4}" fill="#94a3c6" font-size="11" text-anchor="end">${v}</text>`;
  }
  for (let k = 0; k <= 5; k++) { const v = Math.round(xmax * k / 5); body += `<text x="${sx(v)}" y="${H - b + 16}" fill="#94a3c6" font-size="11" text-anchor="middle">${v}</text>`; }
  body += `<line x1="${l}" x2="${W - rt}" y1="${sy(r.threshold)}" y2="${sy(r.threshold)}" stroke="#ff5c7a" stroke-dasharray="6 4"/>`;
  body += `<polyline fill="none" stroke="#37e8ff" stroke-width="1.6" points="${r.positions.map((p, i) => `${sx(p).toFixed(1)},${sy(r.scores[i]).toFixed(1)}`).join(' ')}"/>`;
  body += `<text x="${W - rt}" y="${H - 4}" fill="#94a3c6" font-size="11" text-anchor="end">position (window centre, aa)</text>`;
  body += `<text x="12" y="${t + 10}" fill="#94a3c6" font-size="11">KD</text>`;
  return svg(W, H, body, `Kyte-Doolittle hydropathy of ${r.id}; dashed line = threshold, shaded = predicted TM regions`);
}

route('hydropathy', async (args, q) => {
  const [cfg, p] = await Promise.all([api('/config'), analyse('parse')]);
  const idx = Math.min(Math.max(0, +(q.i || 0) || 0), p.records.length - 1);
  app().innerHTML = `
  <h1>Hydropathy &amp; transmembrane prediction</h1>
  ${await seqbar(p.records.length)}
  <p class="lead">Mean Kyte-Doolittle hydropathy over a sliding window; a predicted TM region runs from the first window
    centre above the threshold to the first one at or below it. The desktop app uses window ${cfg.tm.window} and threshold ${cfg.tm.threshold}.</p>
  ${warningsHtml(p.warnings)}
  <div class="panel"><div class="ctlrow">
    <div><label for="sq">Sequence</label><select id="sq">${p.records.map((r, i) => `<option value="${i}" ${i === idx ? 'selected' : ''}>${esc(label(r.id))} (${r.length} aa)</option>`).join('')}</select></div>
    <div><label for="win">Window (odd, 3-51)</label><input id="win" type="number" min="3" max="51" step="2" value="${cfg.tm.window}"></div>
    <div><label for="thr">Threshold</label><input id="thr" type="number" step="0.1" value="${cfg.tm.threshold}"></div>
    <div class="ctl-btn"><button class="alt" id="def">Defaults</button></div></div>
    <div class="plot" id="plot"></div></div>
  <h2>Predicted TM regions <span class="small mute" id="nr"></span></h2><div class="scroll" id="regions"></div>`;

  const run = async () => {
    const window = Math.round(+$('#win').value), threshold = +$('#thr').value;
    try {
      const r = await analyse('hydropathy', { index: +$('#sq').value, window, threshold });
      $('#plot').innerHTML = plot(r);
      $('#nr').textContent = `(${r.regions.length} for ${label(r.id)})`;
      $('#regions').innerHTML = r.regions.length ? `<table><thead><tr><th>#</th><th>Start-end (1-based)</th><th>Length</th><th>Residues</th></tr></thead><tbody>${
        r.regions.map((g, i) => `<tr><td>TM${i + 1}</td><td>${g.start + 1}-${g.end}</td><td>${g.length}</td><td class="seqcell">${esc(g.sequence)}</td></tr>`).join('')}</tbody></table>`
        : '<p class="mute">No window scores above the threshold.</p>';
      const sc = r.scores; if (sc.length) $('#nr').textContent += ` · profile range ${f(Math.min(...sc))} to ${f(Math.max(...sc))}`;
    } catch (e) { $('#plot').innerHTML = `<p class="err">${esc(e.message)}</p>`; $('#regions').innerHTML = ''; }
  };
  $('#sq').onchange = run; $('#win').onchange = run; $('#thr').onchange = run;
  $('#def').onclick = () => { $('#win').value = cfg.tm.window; $('#thr').value = cfg.tm.threshold; run(); };
  run();
});
