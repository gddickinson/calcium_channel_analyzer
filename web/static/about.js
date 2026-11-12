// About: what the web version covers, how to run the desktop and web apps, credits (from the README).
import { api, app, esc, isStatic, route } from './core.js';

route('about', async () => {
  const ex = await api('/examples');
  const lens = ex.records.map(r => r.length);
  app().innerHTML = `
  <h1>About</h1>
  <p class="lead">The Calcium Channel Analyzer is a Tkinter desktop application for calcium channel sequence analysis
    by Dr. George Dickinson (UC Irvine). This web version shows its sequence analyses; the network, structure and
    machine-learning tools stay in the desktop app.</p>
  <h2>In the web version</h2>
  <ul>
    <li><b>Domain motifs</b>: <code>domain_detection.detect_domains_in_sequence</code> with the ten patterns in <code>config.py</code>.</li>
    <li><b>Hydropathy &amp; TM</b>: the desktop app's Kyte-Doolittle transmembrane prediction (<code>sequence_analysis.py</code>).</li>
    <li><b>Classification &amp; features</b>: the rule-based classifier and electrophysiology labels (<code>channel_rules.py</code>) and
      the ML feature vectors (<code>ml_features.py</code>).</li>
    <li><b>Pairwise alignment</b> and <b>distance tree</b>: the desktop app's Biopython settings (<code>alignment_tree.py</code>).</li>
  </ul>
  <p>${isStatic() ? 'On this static site every analysis runs in your browser, in JavaScript ports of the Python functions that are checked against the Python results by automated tests. Fetching sequences from NCBI or UniProt needs the local app.'
    : 'In the local app every analysis runs in the Python code on this computer, and sequences can be fetched from NCBI protein or UniProt.'}</p>
  <h2>Only in the desktop app</h2>
  <p>BLAST and DIAMOND searches, genome-wide search, orthologs, HMM profiles, AlphaFold structure download and the 3D viewer,
    training and applying the ML classifiers, expression, interaction and pathway queries, and the advanced plot types.</p>
  <h2>The example sequences</h2>
  <p>The example file <code>${esc(ex.file)}</code> holds ${ex.records.length} sequences with UniProt-style headers
    (human Cav1.2, Cav2.1, RyR2, IP3R1; Arabidopsis CNGC2, TPC1, GLR1.1). They are ${Math.min(...lens)}-${Math.max(...lens)} residues long and are
    not the real proteins: they are placeholder sequences built largely from homopolymer and short repeat runs (for example
    long stretches of G, K or GQ), much shorter than the full-length channels named in the headers. They exercise every
    analysis, but the results on them say nothing about those channels. For real analyses fetch the entries (local app) or paste your own FASTA.</p>
  <h2>Run it yourself</h2>
  <pre>git clone https://github.com/gddickinson/calcium_channel_analyzer
cd calcium_channel_analyzer
pip install -r requirements.txt
python calcium_channel_analyzer.py      # desktop app (Tkinter)

pip install -r requirements-web.txt
python run_web.py                       # web app on http://127.0.0.1:8135</pre>
  <h2>Credits</h2>
  <p>Built with Biopython, scikit-learn, Matplotlib/Seaborn, Pandas/NumPy and Tkinter (desktop app); FastAPI for the web app.
    Methods to cite when using the tool (from the README): AlphaFold (Jumper et al., Nature 2021; Abramson et al., Nature 2024),
    DIAMOND (Buchfink et al., Nature Methods 2015), Biopython (Cock et al., Bioinformatics 2009). Built with AI assistance from
    Claude (Anthropic).</p>`;
});
