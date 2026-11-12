// Pairwise alignment: the desktop app's "Multiple Alignment" action (global pairwise alignment of every pair of the
// first 500 residues; match 1, mismatch 0, gap open -0.5, extend -0.1), with a score matrix and the alignments.
import { $, $$, api, app, esc, isStatic, route } from './core.js';
import { analyse, seqbar, warningsHtml } from './state.js';
import { alignmentText, f, heat, label } from './ui.js';

route('align', async () => {
  const cfg = (await api('/config')).align;
  app().innerHTML = `<h1>Pairwise alignment</h1>${await seqbar()}<p class="mute" id="wait">Aligning every pair...</p><div id="out"></div>`;
  let r;
  try { r = await analyse('align'); } catch (e) { $('#wait').innerHTML = `<span class="err">${esc(e.message)}</span>`; return; }
  const n = r.ids.length, S = Array.from({ length: n }, () => Array(n).fill(null));
  r.pairs.forEach((p, k) => { S[p.i][p.j] = S[p.j][p.i] = k; });
  const scores = r.pairs.map(p => p.score), lo = Math.min(...scores), hi = Math.max(...scores);
  $('#wait').remove();
  $('#out').innerHTML = `
  <p class="lead">Global alignment (Gotoh, affine gaps) of each pair's first ${r.prefix} residues, scored match ${cfg.match},
    mismatch ${cfg.mismatch}, gap open ${cfg.gap_open}, gap extend ${cfg.gap_extend}. The desktop app notes this is a demonstration;
    for a real multiple alignment use MUSCLE, MAFFT or Clustal Omega.</p>
  ${warningsHtml(r.warnings)}
  <h2>Alignment scores</h2><p class="small mute">Click a cell to show that alignment.</p>
  <div class="scroll"><table class="heat"><thead><tr><th></th>${r.ids.map(id => `<th class="mini">${esc(label(id))}</th>`).join('')}</tr></thead><tbody>${
    r.ids.map((id, i) => `<tr><td class="mini">${esc(label(id))}</td>${r.ids.map((_, j) => S[i][j] === null ? '<td class="v">·</td>'
      : `<td class="v click" tabindex="0" data-k="${S[i][j]}" style="${heat(hi > lo ? (r.pairs[S[i][j]].score - lo) / (hi - lo) : 0.5)}">${f(r.pairs[S[i][j]].score, 1)}</td>`).join('')}</tr>`).join('')}</tbody></table></div>
  <h2 id="ah">Alignment</h2><div id="aln"></div>
  ${isStatic() ? '<p class="small mute">The browser computes the same optimal scores as Biopython; where several alignments share the optimal score, the one shown may differ from the desktop app\'s.</p>' : ''}`;

  const show = k => {
    const p = r.pairs[k];
    $$('td.sel').forEach(td => td.classList.remove('sel'));
    $$(`td[data-k="${k}"]`).forEach(td => td.classList.add('sel'));
    $('#ah').textContent = `${label(p.seq1)} vs ${label(p.seq2)}`;
    $('#aln').innerHTML = `<p class="small">Score <b>${f(p.score, 1)}</b> · aligned length ${p.aligned1.length} · identical positions
      ${(100 * p.identity).toFixed(1)}% of the aligned length</p><pre class="aln">${esc(alignmentText(p.aligned1, p.aligned2, p.seq1, p.seq2))}</pre>`;
  };
  $$('td.click').forEach(td => {
    td.onclick = () => show(+td.dataset.k);
    td.onkeydown = e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); show(+td.dataset.k); } };
  });
  let best = 0; r.pairs.forEach((p, k) => { if (p.score > r.pairs[best].score) best = k; });
  show(best);
});
