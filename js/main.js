// App bootstrap and animation loop.
//
// URL parameters (handy for debugging / screenshots):
//   ?z=6          element (1..10)
//   ?debug=1      seed = frame number (reproducible frames)
//   ?mode=white   background mode ('dark' | 'white')
//   ?paused=1     start paused
//   ?density=0.2  exposure (mean detections per pixel; peak/20 if equalized)
//   ?equalize=1   equalize orbital brightness (default: true relative densities)
//   ?scale=fixed  same view scale for every element (compare sizes)
//   ?spin=0       start without the slow default rotation

import { elementByZ, ELEMENTS, ORBITALS } from './elements.js';
import { isValence } from './orbitals.js';
import { qAxisAngle, qMul, qNormalize } from './math.js';
import { createRenderer, baseHalf, viewGeometry } from './renderer.js';
import { createScaleAnimator, prefersReducedMotion } from './view.js';
import { loadTables } from './tables.js';
import { attachControls } from './controls.js';
import { buildUI, syncControls, updateElement, updateTable, updateStatus, drawOverlay, defaultExpandedShells } from './ui.js';
import { createRadialPlot, radialMode } from './radial.js';

const params = new URLSearchParams(location.search);

// Default orientation: tilted so x, y and z lobes are all distinguishable.
// Gentle spin at start (about view +y), unless the viewer prefers reduced
// motion or ?spin=0. A fling replaces it; grabbing the atom stops it.
function defaultSpin() {
  if (params.get('spin') === '0' || prefersReducedMotion()) return null;
  return { axis: [0, 1, 0], rate: 0.00025 };            // rad/ms (~14 deg/s)
}

function defaultOrientation() {
  return qMul(qAxisAngle([1, 0, 0], 0.45), qAxisAngle([0, 1, 0], -0.6));
}

const state = {
  element: elementByZ(clampInt(params.get('z'), 1, ELEMENTS.length, 6)),
  orientation: defaultOrientation(),
  zoom: 1,
  density: +params.get('density') || 1,
  equalize: params.get('equalize') === '1',
  fixedScale: params.get('scale') === 'fixed',
  radialMode: 'auto',            // 'auto' | 'shell' | 'subshell'
  spin: defaultSpin(),           // { axis, rate } in view space, or null
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

function syncRadialButtons() {
  const mode = radialMode(state, state.element);
  for (const [id, m] of [['radial-shells', 'shell'], ['radial-subshells', 'subshell']]) {
    const b = document.getElementById(id);
    b.classList.toggle('on', mode === m);
    b.setAttribute('aria-pressed', String(mode === m));
  }
}

function loadCollapsed() {
  try {
    const v = localStorage.getItem('orbitals.tableCollapsed');
    if (v === '1' || v === '0') return v === '1';
  } catch { /* storage unavailable */ }
  return matchMedia('(max-width: 640px), (orientation: portrait) and (max-width: 1024px)').matches;
}

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
    const el = z && elementByZ(z);
    if (!el || !el.supported) return;          // greyed-out elements can't be selected
    state.element = el;
    state.visible = ORBITALS.map(() => true);
    app.expandedShells = defaultExpandedShells(el);
    updateElement(app);
    dirty = true;
  },
  setDensity(d) { state.density = d; syncControls(state); dirty = true; },
  setKMax(k) { state.kMax = k; syncControls(state); dirty = true; },
  setGamma(g) { state.gamma = g; syncControls(state); dirty = true; },
  setEqualize(b) { state.equalize = b; syncControls(state); dirty = true; },
  setFixedScale(b) { state.fixedScale = b; syncControls(state); dirty = true; },
  setRadialMode(m) { state.radialMode = m; radialKey = ''; syncRadialButtons(); dirty = true; },
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
  // Periodic-table strip: collapsed shows only the current row. Remembered
  // per viewer; starts collapsed on phones.
  tableCollapsed: loadCollapsed(),
  setTableCollapsed(b) {
    app.tableCollapsed = b;
    try { localStorage.setItem('orbitals.tableCollapsed', b ? '1' : '0'); } catch { /* storage unavailable */ }
    updateTable(app);
  },
  expanded: new Set(),
  expandedShells: defaultExpandedShells(state.element),
  toggleShell(n) {
    if (app.expandedShells.has(n)) app.expandedShells.delete(n); else app.expandedShells.add(n);
    updateElement(app);
  },
  toggleExpanded(id) {
    if (app.expanded.has(id)) app.expanded.delete(id); else app.expanded.add(id);
    updateElement(app);
  },
  resetView() { state.orientation = defaultOrientation(); state.zoom = 1; state.spin = defaultSpin(); dirty = true; },
  step() { advance(); dirty = true; },
  // Download what is on screen: the WebGL frame (kept by
  // preserveDrawingBuffer) plus the 2D overlay when it is shown.
  savePng() {
    const out = document.createElement('canvas');
    out.width = glCanvas.width;
    out.height = glCanvas.height;
    const ctx = out.getContext('2d');
    ctx.drawImage(glCanvas, 0, 0);
    if (state.overlay) ctx.drawImage(overlay, 0, 0);
    out.toBlob(blob => {
      const a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = `orbitals-${state.element.symbol}-${state.mode}.png`;
      a.click();
      setTimeout(() => URL.revokeObjectURL(a.href), 1000);
    }, 'image/png');
  },
};

attachControls(glCanvas, {
  getOrientation: () => state.orientation,
  setOrientation: q => { state.orientation = q; dirty = true; },
  getZoom: () => state.zoom,
  setZoom: z => { state.zoom = z; dirty = true; },
  setSpin: s => { state.spin = s; },
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
const scaleAnim = createScaleAnimator({ duration: prefersReducedMotion() ? 0 : 450 });

let lastFrameT = null;

function frame(t) {
  resize();
  if (!state.paused) advance();
  // Keep spinning (paused stops it too; dt capped so a background tab
  // doesn't jump on return).
  const dtSpin = lastFrameT == null ? 0 : Math.min(100, t - lastFrameT);
  lastFrameT = t;
  if (state.spin && !state.paused && dtSpin > 0) {
    const dq = qAxisAngle(state.spin.axis, state.spin.rate * dtSpin);
    state.orientation = qNormalize(qMul(dq, state.orientation));
    dirty = true;
  }
  // Ease the view scale toward the element's (or fixed) scale.
  state.baseHalf = scaleAnim.value(baseHalf(state.element, state.fixedScale), t);
  if (scaleAnim.animating) dirty = true;
  if (!state.paused || dirty) {
    if (state.paused && state.debugSeed) state.seed = state.frame; // keep reproducible
    lastGeo = renderer.render(state);
    drawOverlay(overlay, state, lastGeo);
    dirty = false;
    // Radial inset only changes with the element or the *target* view scale,
    // so it updates once per switch rather than every animation frame.
    const target = viewGeometry(state.element, glCanvas.width, glCanvas.height, state.zoom, state.fixedScale);
    const radialW = document.getElementById('radial').clientWidth;   // redraw after layout changes
    const key = `${state.element.z}|${state.fixedScale}|${target.halfExtent[1].toFixed(4)}|${radialW}|${state.radialMode}`;
    if (key !== radialKey) { radialKey = key; radial.draw(state, target); syncRadialButtons(); }
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
