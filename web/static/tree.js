// Distance tree: 1 - identity of pairwise global alignments (BLOSUM62), Bio.Phylo UPGMA or NJ (alignment_tree.py),
// drawn as a phylogram, with the distance matrix and Newick text.
import { $, $$, app, download, esc, route } from './core.js';
import { analyse, seqbar, warningsHtml } from './state.js';
import { f, heat, label, svg } from './ui.js';

function layout(root) {
  const leaves = [], nodes = [];
  const walk = (c, x, depth) => {
    const nx = x + (depth ? c.branch_length : 0);
    const kids = c.children.map(k => walk(k, nx, depth + 1));
    const y = kids.length ? (Math.min(...kids.map(k => k.y)) + Math.max(...kids.map(k => k.y))) / 2 : leaves.push(c) - 1;
    const node = { c, x: nx, px: x, y, kids };
    nodes.push(node);
    return node;
  };
  walk(root, 0, 0);
  return { leaves, nodes };
}

function drawTree(root) {
  const { leaves, nodes } = layout(root);
  const W = 1000, rowH = 30, left = 20, labelW = 190, H = leaves.length * rowH + 40;
  const xs = nodes.map(n => n.x).concat(nodes.map(n => n.px)), x0 = Math.min(0, ...xs), x1 = Math.max(...xs) || 1;
  const sx = x => left + (x - x0) / (x1 - x0) * (W - left - labelW), sy = y => 16 + y * rowH;
  let body = '';
  for (const n of nodes) {
    body += `<line x1="${sx(n.px)}" x2="${sx(n.x)}" y1="${sy(n.y)}" y2="${sy(n.y)}" stroke="#37e8ff" stroke-width="2"/>`;
    if (n.kids.length) {
      const ys = n.kids.map(k => sy(k.y));
      body += `<line x1="${sx(n.x)}" x2="${sx(n.x)}" y1="${Math.min(...ys)}" y2="${Math.max(...ys)}" stroke="#37e8ff" stroke-width="2"/>`;
    } else {
      body += `<circle cx="${sx(n.x)}" cy="${sy(n.y)}" r="3.5" fill="#ffb13d"/><text x="${sx(n.x) + 8}" y="${sy(n.y) + 4}" fill="#e8f1ff" font-size="13">${esc(label(n.c.name))}</text>`;
    }
  }
  const sb = Math.pow(10, Math.floor(Math.log10((x1 - x0) / 3))) * 2;
  const by = H - 14;
  body += `<line x1="${left}" x2="${left + sb / (x1 - x0) * (W - left - labelW)}" y1="${by}" y2="${by}" stroke="#94a3c6" stroke-width="2"/>`;
  body += `<text x="${left}" y="${by - 6}" fill="#94a3c6" font-size="11">${sb.toPrecision(1)} distance</text>`;
  return svg(W, H, body, 'Distance tree of the sequence set');
}

route('tree', async (args, q) => {
  let method = (q.m || 'UPGMA').toUpperCase() === 'NJ' ? 'NJ' : 'UPGMA';
  app().innerHTML = `
  <h1>Distance tree</h1>${await seqbar()}
  <p class="lead">Every pair is aligned end to end (Needleman-Wunsch with affine gaps: BLOSUM62, gap open −10, extend −0.5);
    the distance is 1 − identity, identity being identical columns / alignment length. The tree is built from these distances
    with Bio.Phylo's UPGMA or neighbour joining. Pairwise distances are a quick guide, not a phylogeny from a multiple
    alignment and a substitution model. (The desktop app's tree compares unaligned sequences position by position instead.)</p>
  <div class="seg" role="group" aria-label="Tree method"><button data-m="UPGMA">UPGMA</button><button data-m="NJ">Neighbour joining</button></div>
  <div id="out"></div>`;
  const run = async () => {
    $$('.seg button').forEach(b => b.classList.toggle('on', b.dataset.m === method));
    $('#out').innerHTML = '<p class="mute">Aligning every pair (full length; a few seconds for long channels)...</p>';
    await new Promise(res => setTimeout(res, 30));
    let r;
    try { r = await analyse('tree', { method }); } catch (e) { $('#out').innerHTML = `<p class="err">${esc(e.message)}</p>`; return; }
    const flat = r.distances.flat().filter(v => v > 0), lo = Math.min(...flat), hi = Math.max(...flat);
    $('#out').innerHTML = `${warningsHtml(r.warnings)}
      <div class="panel plot">${drawTree(r.tree)}</div>
      <p class="small mute">${method === 'NJ' ? 'Neighbour-joining trees are unrooted; the drawing hangs them from the last join, as Bio.Phylo returns them.'
        : 'Bio.Phylo\'s UPGMA updates distances as the plain average of the two merged rows.'}</p>
      <h2>Distances (1 − identity)</h2>
      <div class="scroll"><table class="heat"><thead><tr><th></th>${r.ids.map(id => `<th class="mini">${esc(label(id))}</th>`).join('')}</tr></thead><tbody>${
        r.distances.map((row, i) => `<tr><td class="mini">${esc(label(r.ids[i]))}</td>${row.map((v, j) => `<td class="v" style="${i === j ? '' : heat(hi > lo ? 1 - (v - lo) / (hi - lo) : 0.5)}">${f(v, 3)}</td>`).join('')}</tr>`).join('')}</tbody></table></div>
      <h2>Newick</h2><pre id="nw">${esc(r.newick)}</pre><button class="alt" id="dl">Download .nwk</button>`;
    $('#dl').onclick = () => download(`tree_${method.toLowerCase()}.nwk`, r.newick + '\n');
  };
  $$('.seg button').forEach(b => { b.onclick = () => { method = b.dataset.m; run(); }; });
  run();
});
