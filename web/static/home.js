// Home: hero, real counts from /api/info, cards to the analyses.
import { api, app, esc, route } from './core.js';

const CARDS = [
  ['sequences', 'Sequences', 'The example FASTA from the repo, or paste / open your own.'],
  ['domains', 'Domain motifs', 'Ten regex domain signatures (EF-hand, pore loop, S4 voltage sensor...) mapped along each sequence.'],
  ['hydropathy', 'Hydropathy & TM', 'Kyte-Doolittle sliding-window hydropathy and predicted transmembrane segments.'],
  ['classify', 'Classification & features', 'The rule-based family classifier, electrophysiology labels and the 32 ML features.'],
  ['align', 'Pairwise alignment', 'Global alignment of every pair (affine gaps), with a score matrix and the alignments.'],
  ['tree', 'Distance tree', 'UPGMA or neighbour-joining trees from identity distances, as the desktop app builds them.'],
];

route('home', async () => {
  const i = await api('/info');
  app().innerHTML = `
  <section class="hero">
    <p class="kicker">Bioinformatics · calcium channels</p>
    <h1 class="glow">Calcium Channel Analyzer</h1>
    <p class="lead">Sequence analyses from the Calcium Channel Analyzer desktop app (Dr. George Dickinson, UC Irvine):
      conserved domain motifs, hydropathy and transmembrane prediction, rule-based family classification,
      pairwise alignment and distance trees, on the example sequences or your own FASTA.</p>
    <a class="btn" href="#/domains">Map domain motifs</a><a class="btn alt" href="#/sequences">Use your own sequences</a>
  </section>
  <div class="tiles">
    <div class="tile"><b>${i.n_examples}</b><span>example sequences (${i.example_residues} aa)</span></div>
    <div class="tile"><b>${i.n_domain_patterns}</b><span>domain motif patterns</span></div>
    <div class="tile"><b>${i.n_families}</b><span>channel families in config.py</span></div>
    <div class="tile"><b>${i.n_feature_dims}</b><span>ML feature dimensions</span></div>
  </div>
  <div class="grid">${CARDS.map(([h, t, d]) => `<a class="card" href="#/${h}"><b>${esc(t)}</b><span>${esc(d)}</span></a>`).join('')}</div>
  <p class="small mute" style="margin-top:22px">The analyses here are the app's pattern-based heuristics, as in the desktop program:
    regex motifs and keyword rules, not curated domain annotations. See <a href="#/about">About</a> for what the web
    version covers and what stays in the desktop app.</p>`;
});
