// Classification & features: channel_rules.classify_channel / predict_electrophysiology and
// ml_features.extract_ml_features for every sequence.
import { $, app, esc, route } from './core.js';
import { analyse, seqbar, warningsHtml } from './state.js';
import { f, heat, label } from './ui.js';

const yn = v => v ? '<span class="yes">yes</span>' : '<span class="no">no</span>';

route('classify', async () => {
  const r = await analyse('classify');
  const mk = Object.keys(r.rows[0].markers);
  app().innerHTML = `
  <h1>Classification &amp; features</h1>
  ${await seqbar(r.rows.length)}
  <p class="lead">The desktop app's rule-based classifier: keywords in the FASTA header decide first ("High" confidence);
    otherwise motif markers decide ("Medium"). Electrophysiology labels come from the same kind of rules. These are the
    program's heuristics, not predictions from structure or measurements.</p>
  ${warningsHtml(r.warnings)}
  <h2>Family</h2>
  <div class="scroll"><table><thead><tr><th>Sequence</th><th>Classification</th><th>Confidence</th><th>Family keyword</th><th>Length</th>
    ${mk.map(m => `<th class="mini">${esc(m)}</th>`).join('')}</tr></thead><tbody>${r.rows.map(x => `<tr><td title="${esc(x.description)}">${esc(label(x.id))}</td>
    <td>${esc(x.classification)}</td><td><span class="tag ${x.confidence === 'High' ? 'a' : x.confidence === 'Medium' ? 'b' : ''}">${esc(x.confidence)}</span></td>
    <td>${esc(x.family_from_description)}</td><td>${x.length}</td>${mk.map(m => `<td>${yn(x.markers[m])}</td>`).join('')}</tr>`).join('')}</tbody></table></div>
  <h2>Electrophysiology labels</h2>
  <div class="scroll"><table><thead><tr><th>Sequence</th><th>Gating</th><th>Activation</th><th>Selectivity</th><th>Permeability</th><th>Blockers</th></tr></thead>
    <tbody>${r.rows.map(x => { const e = x.electrophysiology; return `<tr><td>${esc(label(x.id))}</td><td>${esc(e.gating)}</td><td class="small">${esc(e.activation)}</td>
    <td>${esc(e.selectivity)}</td><td class="small">${esc(e.permeability)}</td><td class="small">${esc(e.blockers)}</td></tr>`; }).join('')}</tbody></table></div>
  <h2>ML feature vectors (${r.feature_names.length} dimensions)</h2>
  <p class="small mute">The inputs the desktop app's Random Forest / Gradient Boosting classifiers train on: length, the fraction of each
    amino acid, mean Kyte-Doolittle hydropathy and one presence flag per domain pattern. Training needs labelled data and
    stays in the desktop app. Cell shading is relative within each row.</p>
  <div class="scroll" id="feat"></div>`;

  const rows = r.feature_names.map((name, k) => {
    const vals = r.rows.map(x => x.features[k]), lo = Math.min(...vals), hi = Math.max(...vals);
    const fmt = v => name === 'length' ? String(v) : name.startsWith('has_') ? (v ? '1' : '0') : f(v, 3);
    return `<tr><td class="mini">${esc(name)}</td>${vals.map(v => `<td class="v" style="${heat(hi > lo ? (v - lo) / (hi - lo) : 0)}">${fmt(v)}</td>`).join('')}</tr>`;
  });
  $('#feat').innerHTML = `<table class="heat"><thead><tr><th>Feature</th>${r.rows.map(x => `<th class="mini">${esc(label(x.id))}</th>`).join('')}</tr></thead>
    <tbody>${rows.join('')}</tbody></table>`;
});
