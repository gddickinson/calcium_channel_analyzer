// Entry point: lazy page modules, the mode label, a subtle drifting-ion background, then the router.
import { $, app, isStatic, lazy, route, start } from './core.js';

const PAGES = {
  'home.js': ['home'], 'sequences.js': ['sequences'], 'domains.js': ['domains'], 'hydropathy.js': ['hydropathy'],
  'classify.js': ['classify'], 'align.js': ['align'], 'tree.js': ['tree'], 'about.js': ['about'],
};
lazy(Object.fromEntries(Object.entries(PAGES).flatMap(([file, names]) => names.map(n => [n, () => import('./' + file)]))));
route('notfound', () => { app().innerHTML = '<h1>Not found</h1><p><a href="#/">Back to the start</a></p>'; });

$('#mode').textContent = isStatic()
  ? 'static site: analyses run in your browser (JS ports tested against the Python code)'
  : 'local app: analyses run in Python on this computer';

// Background: a few faint amber "ions" drifting down through two cyan membrane lines. Off under reduced motion.
function background(canvas) {
  if (!canvas || matchMedia('(prefers-reduced-motion: reduce)').matches) return;
  const g = canvas.getContext('2d');
  const fit = () => { canvas.width = innerWidth; canvas.height = innerHeight; };
  fit(); addEventListener('resize', fit);
  const ions = Array.from({ length: 24 }, () => ({ x: Math.random(), y: Math.random(), v: 0.0004 + Math.random() * 0.0008 }));
  let last = 0;
  const tick = t => {
    requestAnimationFrame(tick);
    if (document.hidden || t - last < 50) return;
    last = t;
    const W = canvas.width, H = canvas.height;
    g.clearRect(0, 0, W, H);
    g.strokeStyle = 'rgba(55,232,255,.35)'; g.lineWidth = 1;
    for (const y of [0.42, 0.46]) { g.beginPath(); g.moveTo(0, y * H); g.lineTo(W, y * H); g.stroke(); }
    g.fillStyle = '#ffb13d';
    for (const p of ions) {
      p.y += p.v; if (p.y > 1) { p.y = 0; p.x = Math.random(); }
      g.globalAlpha = 0.5; g.beginPath(); g.arc(p.x * W, p.y * H, 2.2, 0, 2 * Math.PI); g.fill();
    }
    g.globalAlpha = 1;
  };
  requestAnimationFrame(tick);
}
background($('#bg'));
start();
