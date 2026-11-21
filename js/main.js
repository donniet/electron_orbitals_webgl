// App bootstrap and animation loop.
//
// URL parameters (handy for debugging / screenshots):
//   ?z=6          element (1..86)
//   ?z2=8         second element, shown side by side
//   ?debug=1      seed = frame number (reproducible frames)
//   ?mode=white   background mode ('dark' | 'white')
//   ?paused=1     start paused
//   ?density=0.2  exposure (mean detections per pixel; peak/20 if equalized)
//   ?equalize=1   equalize orbital brightness (default: true relative densities)
//   ?scale=fixed  same view scale for every element (compare sizes)
//   ?cut=1        cut-away: show only the half of the cloud behind the nucleus
//   ?colour=blend    mix every orbital's colours (default: nearest, the orbital
//                    whose detection is nearest the viewer)
//   ?spin=0       start without the slow default rotation
//   ?res=full     noise canvas resolution: auto (default) | full | half
//   ?kmax=4       counts for full brightness
//   ?gamma=0.5    tone gamma
//   ?record=1     open the video recording controls

import { elementByZ, ELEMENTS, ORBITALS } from './elements.js';
import { qAxisAngle, qMul, qNormalize } from './math.js';
import { createRenderer, baseHalf, viewGeometry, splitViews } from './renderer.js';
import { createScaleAnimator, prefersReducedMotion, defaultOrientation } from './view.js';
import { loadTables } from './tables.js';
import { attachControls } from './controls.js';
import { buildUI, syncControls, updateElement, updateTable, updatePerf, drawOverlay } from './ui.js';
import { atomsOf, valenceUnion, defaultExpandedShells } from './selection.js';
import { createRadialPlot, radialMode } from './radial.js';
import { createRecorderUI } from './recordui.js';

const params = new URLSearchParams(location.search);

// Default orientation: see view.js. Gentle spin at start (about view +y),
// unless the viewer prefers reduced motion or ?spin=0. A fling replaces it;
// grabbing the atom stops it.
function defaultSpin() {
  if (params.get('spin') === '0' || prefersReducedMotion()) return null;
  return { axis: [0, 1, 0], rate: 0.0025 };            // rad/ms (~14 deg/s)
}

const state = {
  element: elementByZ(clampInt(params.get('z'), 1, ELEMENTS.length, 6)),
  second: null,                  // optional second atom, shown side by side (shift-click)
  orientation: defaultOrientation(),
  zoom: 1,
  density: +params.get('density') || 1,
  equalize: params.get('equalize') === '1',
  fixedScale: params.get('scale') === 'fixed',
  cut: params.get('cut') === '1',  // cut-away: only the cloud behind the nucleus
  colour: params.get('colour') === 'blend' ? 'blend' : 'nearest',
  radialMode: 'auto',            // 'auto' | 'shell' | 'subshell'
  spin: defaultSpin(),           // { axis, rate } in view space, or null
  kMax: +params.get('kmax') || 4,
  gamma: +params.get('gamma') || 0.5,
  mode: params.get('mode') === 'white' ? 'white' : 'dark',
  visible: ORBITALS.map(() => true),
  paused: params.get('paused') === '1',
  resMode: ['full', 'half'].includes(params.get('res')) ? params.get('res') : 'auto',
  debugSeed: params.get('debug') === '1',
  overlay: true,
  frame: 0,
  seed: 0,
};

// After A or B changes: expand the selection's shells, re-try full
// resolution for the new workload, refresh the panel and redraw.
function selectionChanged() {
  app.expandedShells = defaultExpandedShells(atomsOf(state));
  if (state.resMode === 'auto') { renderK = 1; slowSince = null; }
  updateElement(app);
  dirty = true;
}

function syncRadialButtons() {
  const mode = radialMode(state, state.element);
  document.getElementById('radial-compare').hidden = !state.second;
  for (const [id, m] of [['radial-shells', 'shell'], ['radial-subshells', 'subshell'], ['radial-compare', 'compare']]) {
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
const view = document.getElementById('view');

// Exact device-pixel size of the view. Per-pixel noise must map 1:1 onto
// device pixels: any non-integer browser rescale (e.g. a 2x canvas on a 3x
// screen) turns white noise into a periodic moire.
let deviceSize = null;
try {
  const ro = new ResizeObserver(entries => {
    const box = entries[0].devicePixelContentBoxSize?.[0];
    deviceSize = box ? [box.inlineSize, box.blockSize] : null;
  });
  try { ro.observe(view, { box: 'device-pixel-content-box' }); } catch { ro.observe(view); }
} catch { /* no ResizeObserver: fall back to clientWidth * devicePixelRatio */ }

// Render scale k (1 = native): the noise is rendered at 1/k of the device
// resolution and blitted up by exactly k with nearest-neighbour sampling, so
// each sample becomes a clean k x k block (never a moire).
let renderK = state.resMode === 'half' ? 2 : 1;
let slowSince = null;
let lastDevice = '';                      // Auto: when frames started being slow
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
let lastViews = null;          // [{ element, geo }] of the last rendered frame

const app = {
  state,
  setResMode(m) { state.resMode = m; renderK = m === 'half' ? 2 : 1; slowSince = null; syncControls(state); dirty = true; },
  // Side-by-side comparison (undocumented): shift-click another element to
  // show it next to the primary; shift-click it again or press Esc to clear.
  // Shift-click / long-press an element, or tap one while "Compare" is
  // armed: it becomes atom B (or stops being B if it already is).
  toggleSecond(z) {
    const el = z && elementByZ(z);
    app.pickingSecond = false;
    if (el && el.supported && el !== state.element) {
      state.second = state.second === el ? null : el;
      selectionChanged();
    } else {
      updateTable(app);
      updateElement(app);
    }
  },
  clearSecond() {
    const had = !!state.second;
    state.second = null;
    if (had) selectionChanged();
  },
  swapAtoms() {
    if (!state.second) return;
    [state.element, state.second] = [state.second, state.element];
    selectionChanged();
  },
  pickingSecond: false,
  togglePicking() {
    app.pickingSecond = !app.pickingSecond;
    updateElement(app);
  },
  setElement(z) {
    const el = z && elementByZ(z);
    if (!el || !el.supported) return;          // greyed-out elements can't be selected
    if (state.second === el) state.second = null;   // stepping onto the second atom
    state.element = el;
    state.visible = ORBITALS.map(() => true);
    selectionChanged();
  },
  setDensity(d) { state.density = d; syncControls(state); dirty = true; },
  setKMax(k) { state.kMax = k; syncControls(state); dirty = true; },
  setGamma(g) { state.gamma = g; syncControls(state); dirty = true; },
  setEqualize(b) { state.equalize = b; syncControls(state); dirty = true; },
  setFixedScale(b) { state.fixedScale = b; syncControls(state); dirty = true; },
  setCut(b) { state.cut = b; syncControls(state); dirty = true; },
  setColour(m) { state.colour = m; syncControls(state); dirty = true; },
  setRadialMode(m) { state.radialMode = m; radialKey = ''; syncRadialButtons(); dirty = true; },
  setMode(m) { state.mode = m; syncControls(state); dirty = true; },
  setPaused(p) { state.paused = p; syncControls(state); },
  setDebugSeed(b) { state.debugSeed = b; syncControls(state); dirty = true; },
  setOverlay(b) { state.overlay = b; syncControls(state); dirty = true; },
  setVisible(i, b) { state.visible[i] = b; updateElement(app); dirty = true; },
  setVisibleMany(idx, b) { for (const i of idx) state.visible[i] = b; updateElement(app); dirty = true; },
  showAll() { app.setVisibleMany(ORBITALS.map((o, i) => i), true); },
  showValence() {
    const valence = valenceUnion(atomsOf(state));
    state.visible = ORBITALS.map((o, i) => valence.has(i));
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
  expandedShells: defaultExpandedShells([state.element]),
  toggleShell(n) {
    if (app.expandedShells.has(n)) app.expandedShells.delete(n); else app.expandedShells.add(n);
    updateElement(app);
  },
  toggleExpanded(id) {
    if (app.expanded.has(id)) app.expanded.delete(id); else app.expanded.add(id);
    updateElement(app);
  },
  resetView() { state.orientation = defaultOrientation(); state.zoom = 1; state.spin = defaultSpin(); dirty = true; recorder.refresh(); },
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
  // Printable page (print.html) in a new tab: this view on white, with
  // the element details. Everything it needs is in the URL (see
  // printParams), so the page can be reloaded or bookmarked.
  openPrintPage() {
    window.open(`print.html?${printParams(state)}`, '_blank', 'noopener');
  },
};

// URL parameters for print.html: the view (atoms, orientation, zoom), the
// exposure and tone settings and which orbitals are hidden. The background
// is always white there; the seed is not passed (the page draws a fresh
// sample, re-rollable; ?seed= fixes it).
function printParams(s) {
  const p = new URLSearchParams({ z: s.element.z });
  if (s.second) p.set('z2', s.second.z);
  p.set('q', s.orientation.map(c => +c.toFixed(6)).join(','));
  if (s.zoom !== 1) p.set('zoom', +s.zoom.toFixed(4));
  p.set('density', +s.density.toFixed(4));
  p.set('kmax', s.kMax);
  p.set('gamma', s.gamma);
  if (s.equalize) p.set('equalize', '1');
  if (s.fixedScale) p.set('scale', 'fixed');
  if (s.cut) p.set('cut', '1');
  if (s.colour !== 'nearest') p.set('colour', s.colour);
  if (!s.overlay) p.set('overlay', '0');
  if (s.radialMode !== 'auto') p.set('radial', s.radialMode);
  const hidden = s.visible.flatMap((v, i) => (v ? [] : [i]));
  if (hidden.length) p.set('hide', hidden.join(','));
  return p.toString();
}

attachControls(glCanvas, {
  getOrientation: () => state.orientation,
  setOrientation: q => { state.orientation = q; dirty = true; },
  getZoom: () => state.zoom,
  setZoom: z => { state.zoom = z; dirty = true; },
  setSpin: s => { state.spin = s; recorder?.refresh(); },
});

buildUI(app);
app.tables = tables;
// Video recording (red dot over the view). While it renders, the live view
// is frozen (see frame()).
const recorder = createRecorderUI({
  app,
  viewInfo: () => ({
    width: glCanvas.width, height: glCanvas.height,
    bohrPerPixel: lastViews?.[0]?.geo.bohrPerPixel ??
      viewGeometry(state.element, glCanvas.width, glCanvas.height, state.zoom, state.fixedScale).bohrPerPixel,
    dpr: glCanvas.width / view.clientWidth || 1,
    panes: lastViews?.map(v => v.geo.rect) ?? [{ x: 0, y: 0, w: glCanvas.width, h: glCanvas.height }],
  }),
  views: () => (state.second ? [state.element, state.second] : [state.element])
    .map((element, i) => ({ element, visible: i === 0 ? state.visible : undefined })),
  changed: () => { dirty = true; },
  // Lay the strip out collapsed for a moment (no paint happens in between)
  // to measure it, whatever its current state.
  // While a video renders, the view shows each finished frame, scaled to
  // fit and letterboxed, on the overlay canvas (the live view is paused
  // underneath; it is redrawn when rendering ends). Smoothing averages the
  // noise when scaling down, which avoids moiré in this preview.
  preview: frame => {
    const g = overlay.getContext('2d');
    const W = overlay.width, H = overlay.height;
    const s = Math.min(W / frame.width, H / frame.height);
    const w = Math.round(frame.width * s), h = Math.round(frame.height * s);
    g.save();
    g.fillStyle = '#000';
    g.fillRect(0, 0, W, H);
    g.imageSmoothingEnabled = true;
    g.imageSmoothingQuality = 'high';
    g.drawImage(frame, Math.round((W - w) / 2), Math.round((H - h) / 2), w, h);
    g.restore();
  },
  tableHeight: () => {
    const pt = document.getElementById('ptable');
    const hidden = pt.hidden, collapsed = app.tableCollapsed;
    pt.hidden = false;
    if (!collapsed) { app.tableCollapsed = true; updateTable(app); }
    const h = pt.getBoundingClientRect().height;
    if (!collapsed) { app.tableCollapsed = false; updateTable(app); }
    pt.hidden = hidden;
    return h;
  },
});
const radial = createRadialPlot(
  document.getElementById('radial'),
  document.getElementById('radial-tip'),
  document.getElementById('radial-legend'));
let radialKey = '';
syncControls(state);
updateElement(app);
if (params.get('z2')) app.toggleSecond(clampInt(params.get('z2'), 1, ELEMENTS.length, 0));
if (params.get('record') === '1') recorder.open();

function advance() {
  state.frame++;
  state.seed = state.debugSeed ? state.frame : (Math.random() * 4294967296) >>> 0;
}

function resize() {
  const dpr = window.devicePixelRatio || 1;
  const [dw, dh] = deviceSize ||
    [Math.round(view.clientWidth * dpr), Math.round(view.clientHeight * dpr)];
  if (`${dw}x${dh}` !== lastDevice) {             // new size: Auto re-tries full resolution
    lastDevice = `${dw}x${dh}`;
    if (state.resMode === 'auto') { renderK = 1; slowSince = null; }
    recorder.refresh();
  }
  if (overlay.width !== dw || overlay.height !== dh) { overlay.width = dw; overlay.height = dh; dirty = true; }
  // The canvas always matches device pixels exactly; render scale is
  // applied inside the renderer (low-res pass + exact integer blit).
  if (glCanvas.width !== dw || glCanvas.height !== dh) { glCanvas.width = dw; glCanvas.height = dh; dirty = true; }
}

// Auto resolution: if frames stay under ~45 fps for 1.5 s, halve the
// noise-canvas resolution (k = 1 -> 2 -> 3). It is reset to full on
// element changes and resizes rather than probing upward, so it never
// oscillates.
function autoResolution(t) {
  if (state.resMode !== 'auto' || state.paused) { slowSince = null; return; }
  if (fps > 0 && fps < 45) {
    if (slowSince == null) slowSince = t;
    else if (t - slowSince > 1500 && renderK < 3) { renderK++; slowSince = null; fps = 50; dirty = true; }
  } else {
    slowSince = null;
  }
}

let fps = 0, lastT = performance.now();
const scaleAnim = createScaleAnimator({ duration: prefersReducedMotion() ? 0 : 450 });

let lastFrameT = null;

function frame(t) {
  if (recorder.phase() === 'rendering') {      // frozen while a video renders
    lastFrameT = null;
    requestAnimationFrame(frame);
    return;
  }
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
  // Ease the view scale toward the element's (or fixed) scale; with two
  // atoms, the larger of the two so their sizes compare directly.
  const atoms = state.second ? [state.element, state.second] : [state.element];
  state.baseHalf = scaleAnim.value(Math.max(...atoms.map(el => baseHalf(el, state.fixedScale))), t);
  if (scaleAnim.animating) dirty = true;
  if (!state.paused || dirty) {
    if (state.paused && state.debugSeed) state.seed = state.frame; // keep reproducible
    const rects = splitViews(glCanvas.width, glCanvas.height, renderK, atoms.length);
    const views = atoms.map((element, i) => ({
      element, rect: rects[i], visible: i === 0 ? state.visible : undefined,   // second atom: all orbitals
    }));
    const geos = renderer.renderViews(state, renderK, views);
    lastViews = views.map((v, i) => ({ element: v.element, geo: geos[i] }));
    drawOverlay(overlay, state, lastViews);
    recorder.drawCrop(overlay);
    dirty = false;
    // Radial inset only changes with the element or the *target* view scale,
    // so it updates once per switch rather than every animation frame.
    const target = viewGeometry(state.element, glCanvas.width, glCanvas.height, state.zoom, state.fixedScale);
    const radialW = document.getElementById('radial').clientWidth;   // redraw after layout changes
    const key = `${state.element.z}|${state.second?.z ?? 0}|${state.fixedScale}|${target.halfExtent[1].toFixed(4)}|${radialW}|${state.radialMode}`;
    if (key !== radialKey) { radialKey = key; radial.draw(state, target); syncRadialButtons(); }
  }
  const dt = t - lastT;
  lastT = t;
  if (dt > 0) fps = fps * 0.9 + (1000 / dt) * 0.1;
  autoResolution(t);
  updatePerf(fps, renderK, state.paused);
  requestAnimationFrame(frame);
}

if (state.paused) advance();   // have a valid first frame/seed
requestAnimationFrame(frame);

// Expose for the console / tests.
window.orbitalsApp = app;
