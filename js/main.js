// App bootstrap and animation loop.
//
// URL parameters (handy for debugging / screenshots):
//   ?z=6          element (1..10)
//   ?debug=1      seed = frame number (reproducible frames)
//   ?mode=white   background mode ('dark' | 'white')
//   ?paused=1     start paused
//   ?density=0.2  exposure (mean detections per pixel; peak/20 if equalized)
//   ?equalize=0   show true relative densities (default: equalized)
//   ?scale=fixed  same view scale for every element (compare sizes)

import { elementByZ, ELEMENTS, ORBITALS } from './elements.js';
import { isValence } from './orbitals.js';
import { qAxisAngle, qMul } from './math.js';
import { createRenderer } from './renderer.js';
import { loadTables } from './tables.js';
import { attachControls } from './controls.js';
import { buildUI, syncControls, updateElement, updateStatus, drawOverlay } from './ui.js';
import { createRadialPlot } from './radial.js';

const params = new URLSearchParams(location.search);

// Default orientation: tilted so x, y and z lobes are all distinguishable.
function defaultOrientation() {
  return qMul(qAxisAngle([1, 0, 0], 0.45), qAxisAngle([0, 1, 0], -0.6));
}

const state = {
  element: elementByZ(clampInt(params.get('z'), 1, ELEMENTS.length, 6)),
  orientation: defaultOrientation(),
  zoom: 1,
  density: +params.get('density') || 0.15,
  equalize: params.get('equalize') !== '0',
  fixedScale: params.get('scale') === 'fixed',
  kMax: +params.get('kmax') || 4,
  gamma: +params.get('gamma') || 0.5,
  mode: params.get('mode') === 'white' ? 'white' : 'dark',
  visible: ORBITALS.map(() => true),
  paused: params.get('paused') === '1',
  debugSeed: params.get('debug') === '1',
  overlay: true,
  frame: 0,
  seed: 0,
};

function clampInt(v, lo, hi, dflt) {
  const n = parseInt(v, 10);
  return Number.isFinite(n) ? Math.min(hi, Math.max(lo, n)) : dflt;
}

const glCanvas = document.getElementById('gl');
const overlay = document.getElementById('overlay');
let renderer, tables;
try {
  tables = await loadTables();
  renderer = createRenderer(glCanvas, tables);
} catch (err) {
  document.getElementById('error').textContent = err.message;
  document.getElementById('error').hidden = false;
  throw err;
}

let dirty = true;          // re-render while paused when something changes
let lastGeo = null;

const app = {
  state,
  setElement(z) {
    state.element = elementByZ(z);
    state.visible = ORBITALS.map(() => true);
    updateElement(app);
    dirty = true;
  },
  setDensity(d) { state.density = d; syncControls(state); dirty = true; },
  setKMax(k) { state.kMax = k; syncControls(state); dirty = true; },
  setGamma(g) { state.gamma = g; syncControls(state); dirty = true; },
  setEqualize(b) { state.equalize = b; syncControls(state); dirty = true; },
  setFixedScale(b) { state.fixedScale = b; syncControls(state); dirty = true; },
  setMode(m) { state.mode = m; syncControls(state); dirty = true; },
  setPaused(p) { state.paused = p; syncControls(state); },
  setDebugSeed(b) { state.debugSeed = b; syncControls(state); dirty = true; },
  setOverlay(b) { state.overlay = b; syncControls(state); dirty = true; },
  setVisible(i, b) { state.visible[i] = b; updateElement(app); dirty = true; },
  setVisibleMany(idx, b) { for (const i of idx) state.visible[i] = b; updateElement(app); dirty = true; },
  showAll() { app.setVisibleMany(ORBITALS.map((o, i) => i), true); },
  showValence() {
    state.visible = ORBITALS.map(o => isValence(o, state.element));
    updateElement(app);
    dirty = true;
  },
  expanded: new Set(),
  toggleExpanded(id) {
    if (app.expanded.has(id)) app.expanded.delete(id); else app.expanded.add(id);
    updateElement(app);
  },
  resetView() { state.orientation = defaultOrientation(); state.zoom = 1; dirty = true; },
  step() { advance(); dirty = true; },
};

attachControls(glCanvas, {
  getOrientation: () => state.orientation,
  setOrientation: q => { state.orientation = q; dirty = true; },
  getZoom: () => state.zoom,
  setZoom: z => { state.zoom = z; dirty = true; },
});

buildUI(app);
const radial = createRadialPlot(
  document.getElementById('radial'),
  document.getElementById('radial-tip'),
  document.getElementById('radial-legend'));
let radialKey = '';
syncControls(state);
updateElement(app);

function advance() {
  state.frame++;
  state.seed = state.debugSeed ? state.frame : (Math.random() * 4294967296) >>> 0;
}

function resize() {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = Math.max(1, Math.round(glCanvas.clientWidth * dpr));
  const h = Math.max(1, Math.round(glCanvas.clientHeight * dpr));
  for (const c of [glCanvas, overlay]) {
    if (c.width !== w || c.height !== h) { c.width = w; c.height = h; dirty = true; }
  }
}

let fps = 0, lastT = performance.now();

function frame(t) {
  resize();
  if (!state.paused) advance();
  if (!state.paused || dirty) {
    if (state.paused && state.debugSeed) state.seed = state.frame; // keep reproducible
    lastGeo = renderer.render(state);
    drawOverlay(overlay, state, lastGeo);
    dirty = false;
    // Radial inset only changes with the element or the view scale.
    const key = `${state.element.z}|${state.fixedScale}|${lastGeo.halfExtent[1].toFixed(4)}`;
    if (key !== radialKey) { radialKey = key; radial.draw(state, lastGeo); }
  }
  const dt = t - lastT;
  lastT = t;
  if (dt > 0) fps = fps * 0.9 + (1000 / dt) * 0.1;
  if (lastGeo) updateStatus(state, lastGeo, state.paused ? 0 : fps);
  requestAnimationFrame(frame);
}

if (state.paused) advance();   // have a valid first frame/seed
requestAnimationFrame(frame);

// Expose for the console / tests.
window.orbitalsApp = app;
