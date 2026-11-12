// Sequences: the working set (example FASTA or your own), parsed with the repo's FASTA rules; local app can fetch
// from NCBI protein or UniProt.
import { $, api, app, esc, isStatic, needsServer, route } from './core.js';
import { examples, getFasta, getSource, resetFasta, setFasta, warningsHtml } from './state.js';
import { label } from './ui.js';

route('sequences', async () => {
  const ex = await examples();
  app().innerHTML = `
  <h1>Sequences</h1>
  <p class="lead">Every analysis runs on this set. It starts as the example file shipped with the repo
    (<code>${esc(ex.file)}</code>, ${ex.records.length} sequences); paste or open your own protein FASTA to analyse that instead.</p>
  <div class="panel">
    <label for="fa">FASTA (protein, one or more records; current set: <b id="src"></b>)</label>
    <textarea id="fa" rows="10" spellcheck="false"></textarea>
    <div class="row"><div><label for="file">Open a FASTA file (read in your browser, not uploaded)</label>
      <input type="file" id="file" accept=".fasta,.fa,.faa,.txt"></div></div>
    <button id="use">Use this set</button><button class="alt" id="reset">Back to the example</button>
    <div id="msg"></div>
  </div>
  <div class="panel" id="fetchp"><h3>Fetch from NCBI or UniProt</h3>${isStatic() ? needsServer('Fetching sequences from NCBI or UniProt') : `
    <div class="ctlrow"><div><label for="srcsel">Database</label><select id="srcsel"><option value="uniprot">UniProt</option><option value="ncbi">NCBI protein</option></select></div>
    <div><label for="acc">Accession(s), comma separated</label><input id="acc" placeholder="e.g. Q13936"></div>
    <div class="ctl-btn"><button id="fetch">Fetch</button></div></div>
    <p class="small mute">Fetched records are added to the text box; press "Use this set" to analyse them. Needs the network.</p>
    <div id="fmsg"></div>`}</div>
  <h2>Parsed records</h2><div id="tbl"></div>`;
  const fa = $('#fa');
  fa.value = await getFasta();
  $('#src').textContent = getSource();

  const show = async () => {
    try {
      const p = await api('/parse', { body: { fasta: fa.value } });
      $('#tbl').innerHTML = warningsHtml(p.warnings) + `<div class="scroll"><table><thead><tr><th>Id</th><th>Description</th>
        <th>Family (from description)</th><th>Length</th><th>Start of sequence</th></tr></thead><tbody>${p.records.map(r => `
        <tr><td title="${esc(r.id)}">${esc(label(r.id))}</td><td class="small">${esc(r.description)}</td><td>${esc(r.family)}</td>
        <td>${r.length}</td><td class="seqcell">${esc(r.sequence.slice(0, 40))}${r.length > 40 ? '…' : ''}</td></tr>`).join('')}</tbody></table></div>`;
      return p;
    } catch (e) { $('#tbl').innerHTML = `<p class="err">${esc(e.message)}</p>`; return null; }
  };
  await show();

  $('#use').onclick = async () => {
    const p = await show();
    if (!p) { $('#msg').innerHTML = '<p class="err">Not changed: fix the FASTA first.</p>'; return; }
    setFasta(fa.value, 'your sequences');
    $('#src').textContent = getSource();
    $('#msg').innerHTML = `<p class="ok">Using ${p.records.length} sequences. Open an analysis: <a href="#/domains">domains</a>,
      <a href="#/hydropathy">hydropathy</a>, <a href="#/classify">classification</a>, <a href="#/align">alignment</a>, <a href="#/tree">tree</a>.</p>`;
  };
  $('#reset').onclick = async () => { fa.value = await resetFasta(); $('#src').textContent = getSource(); $('#msg').innerHTML = ''; show(); };
  $('#file').onchange = async e => { const f = e.target.files[0]; if (f) { fa.value = await f.text(); show(); } };
  const fb = $('#fetch');
  if (fb) fb.onclick = async () => {
    fb.disabled = true; $('#fmsg').innerHTML = '<p class="mute">Fetching...</p>';
    try {
      const r = await api('/fetch', { body: { source: $('#srcsel').value, accessions: $('#acc').value } });
      if (r.fasta) fa.value = (fa.value.trim() ? fa.value.trim() + '\n' : '') + r.fasta;
      $('#fmsg').innerHTML = (r.errors.length ? `<p class="err">${r.errors.map(esc).join('<br>')}</p>` : '') +
        (r.fasta ? '<p class="ok">Added to the text box.</p>' : '');
      show();
    } catch (e) { $('#fmsg').innerHTML = `<p class="err">${esc(e.message)}</p>`; }
    fb.disabled = false;
  };
});
