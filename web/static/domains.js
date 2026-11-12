// Domain motifs: domain_detection.detect_domains_in_sequence on every sequence, drawn as a domain map, with the
// patterns from config.py and a per-sequence count matrix.
import { $, $$, api, app, esc, route } from './core.js';
import { analyse, seqbar, warningsHtml } from './state.js';
import { color, label, svg } from './ui.js';

function domainMap(r, on) {
  const W = 1000, rowH = 34, left = 170, maxLen = Math.max(...Object.values(r.lengths));
  const sx = x => left + (x / maxLen) * (W - left - 10);
  let body = '';
  r.order.forEach((id, k) => {
    const y = 16 + k * rowH, L = r.lengths[id];
    body += `<text x="${left - 8}" y="${y + 4}" text-anchor="end" fill="#c6d4ee" font-size="12">${esc(label(id))}</text>`;
    body += `<line x1="${sx(0)}" x2="${sx(L)}" y1="${y}" y2="${y}" stroke="#94a3c6" stroke-width="2"/>`;
    for (const h of r.hits) if (h.sequence === id && on.has(h.domain)) {
      body += `<rect x="${sx(h.start)}" y="${y - 8}" width="${Math.max(1.2, sx(h.end) - sx(h.start))}" height="16" fill="${color(h.domain)}"
        fill-opacity=".8"><title>${esc(h.domain)} ${h.start + 1}-${h.end}: ${esc(h.match)} (score ${h.score})</title></rect>`;
    }
  });
  const ay = 16 + r.order.length * rowH;
  for (let t = 0; t <= 4; t++) {
    const v = Math.round(maxLen * t / 4);
    body += `<text x="${sx(v)}" y="${ay}" fill="#94a3c6" font-size="11" text-anchor="middle">${v}</text>`;
  }
  body += `<text x="${W - 10}" y="${ay + 16}" fill="#94a3c6" font-size="11" text-anchor="end">position (aa)</text>`;
  return svg(W, ay + 22, body, 'Domain map: motif matches along each sequence');
}

route('domains', async () => {
  const [cfg, r] = await Promise.all([api('/config'), analyse('domains')]);
  const names = Object.keys(cfg.domain_patterns);
  const on = new Set(names);
  const counts = id => names.map(n => r.hits.filter(h => h.sequence === id && h.domain === n).length);
  app().innerHTML = `
  <h1>Domain motifs</h1>
  ${await seqbar(r.order.length)}
  <p class="lead">Regex signatures from <code>config.py</code>, matched with <code>detect_domains_in_sequence</code>
    (non-overlapping matches, kept when at least the pattern's minimum length). Short, permissive patterns match often in
    any protein: read a hit as a motif match, not a confirmed domain.</p>
  ${warningsHtml(r.warnings)}
  <div class="panel"><div class="legend" id="lg">${names.map(n => `<label><input type="checkbox" data-d="${n}" checked>
    <span class="sw" style="background:${color(n)}"></span>${esc(n)}</label>`).join('')}</div>
    <div class="dmap" id="map"></div><p class="small mute">Hover a block for its position (1-based), residues and score.</p></div>
  <h2>Matches per sequence</h2>
  <div class="scroll"><table class="heat"><thead><tr><th>Sequence</th>${names.map(n => `<th class="mini">${esc(n)}</th>`).join('')}</tr></thead>
    <tbody>${r.order.map(id => `<tr><td>${esc(label(id))}</td>${counts(id).map(c => `<td class="v">${c || '·'}</td>`).join('')}</tr>`).join('')}</tbody></table></div>
  <h2>Patterns</h2>
  <div class="scroll"><table><thead><tr><th>Domain</th><th>Pattern</th><th>Min length</th><th>Description</th><th>Channels</th></tr></thead>
    <tbody>${names.map(n => { const p = cfg.domain_patterns[n]; return `<tr><td><span class="sw" style="background:${color(n)}"></span> ${esc(n)}</td>
      <td class="mini">${esc(p.pattern)}</td><td>${p.min_length}</td><td>${esc(p.description)}</td><td class="small">${esc(p.channels.join(', '))}</td></tr>`; }).join('')}</tbody></table></div>
  <h2>All matches <span class="small mute" id="nh"></span></h2>
  <div class="row"><div><label for="flt">Show domain</label><select id="flt"><option value="">all</option>${names.map(n => `<option>${esc(n)}</option>`).join('')}</select></div></div>
  <div class="scroll" id="hits"></div>`;

  const draw = () => { $('#map').innerHTML = domainMap(r, on); };
  const table = () => {
    const d = $('#flt').value, rows = r.hits.filter(h => !d || h.domain === d), shown = rows.slice(0, 400);
    $('#nh').textContent = `(${rows.length}${rows.length > shown.length ? `, first ${shown.length} shown` : ''})`;
    $('#hits').innerHTML = `<table><thead><tr><th>Sequence</th><th>Domain</th><th>Position (1-based)</th><th>Match</th><th>Score</th></tr></thead><tbody>${
      shown.map(h => `<tr><td>${esc(label(h.sequence))}</td><td>${esc(h.domain)}</td><td>${h.start + 1}-${h.end}</td>
      <td class="seqcell">${esc(h.match)}</td><td>${h.score}</td></tr>`).join('')}</tbody></table>`;
  };
  $$('#lg input').forEach(cb => { cb.onchange = () => { if (cb.checked) on.add(cb.dataset.d); else on.delete(cb.dataset.d); draw(); }; });
  $('#flt').onchange = table;
  draw(); table();
});
