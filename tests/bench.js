// Frame-time probe: ms per frame for selected elements at a given size.
import { ELEMENTS, ORBITALS } from '../js/elements.js';
import { createRenderer } from '../js/renderer.js';
import { loadTables } from '../js/tables.js';

const params = new URLSearchParams(location.search);
const zs = (params.get('z') || '1,10,18,36').split(',').map(Number)
  .filter(z => z >= 1 && z <= ELEMENTS.length);
const size = +params.get('size') || 800;
const frames = +params.get('frames') || 20;
const color = params.get('color') === 'blend' ? 'blend' : 'nearest';   // app default: nearest
const cut = params.get('cut') === '1';

const tables = await loadTables();
const canvas = document.createElement('canvas');
canvas.width = canvas.height = size;
const r = createRenderer(canvas, tables);
const px = new Uint8Array(4);
const lines = [`size ${size}x${size}, ${frames} frames each, color ${color}${cut ? ', cut' : ''}`];
const results = {};

for (const z of zs) {
  const el = ELEMENTS[z - 1];
  const state = {
    element: el, orientation: [0, 0, 0, 1], zoom: 1, seed: 1, density: 1,
    equalize: false, kMax: 4, gamma: 0.5, mode: 'dark', fixedScale: false,
    visible: ORBITALS.map(() => true), color, cut,
  };
  r.render(state);                                    // warm-up
  r.gl.readPixels(0, 0, 1, 1, r.gl.RGBA, r.gl.UNSIGNED_BYTE, px);
  const t0 = performance.now();
  for (let f = 0; f < frames; f++) {
    state.seed = f + 2;
    r.render(state);
    r.gl.readPixels(0, 0, 1, 1, r.gl.RGBA, r.gl.UNSIGNED_BYTE, px);   // wait for the GPU
  }
  const ms = (performance.now() - t0) / frames;
  const occupied = el.occupancy.filter(o => o > 0).length;
  results[el.symbol] = ms;
  lines.push(`${el.symbol.padEnd(3)} Z=${String(z).padStart(2)}  orbitals ${String(occupied).padStart(2)}  ${ms.toFixed(1)} ms/frame  (${(1000 / ms).toFixed(0)} fps)`);
}
document.getElementById('out').textContent = lines.join('\n');
document.title = 'BENCH ' + Object.entries(results).map(([k, v]) => `${k}=${v.toFixed(1)}ms`).join(' ');
console.log('BENCH ' + lines.join(' | '));
