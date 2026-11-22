// Printable page (print.html): one frame of the current view rendered on
// white at print resolution, in a box marking the bounds of the rendering,
// with ruled space for notes above and below and the side panel's details
// (element card, orbital table, radial density, model notes) underneath.
//
// The app's Print button passes the view in the URL (see printParams in
// main.js):
//   z, z2               atom A and optional atom B (side by side)
//   q=w,x,y,z | zoom    orientation quaternion and zoom (default: the app's)
//   density, kmax, gamma, equalize=1, scale=fixed, cut=1, color=blend,
//   streaks=local|axis|global,
//   overlay=0, radial=shell|subshell|compare   as in the app
//   hide=i,j,...        orbital indices hidden (in both panes, like the app)
// Print-page options (also in its toolbar):
//   paper=letter|a4     grain=1|2|4 (image pixels per noise sample at 300/in)
//   seed=n              fixed random seed (default: a fresh sample)

import { elementByZ, ELEMENTS, ORBITALS } from './elements.js';
import { SUBSHELLS, orbitalsOf } from './orbitals.js';
import { loadTables } from './tables.js';
import { createRenderer, baseHalf, viewGeometry, splitViews, streakParam } from './renderer.js';
import { defaultOrientation } from './view.js';
import { qNormalize } from './math.js';
import { drawOverlay, orbitalCss, periodGroup, en, rmsPm, modelNotes, CAPACITY, SHELL_LETTER } from './ui.js';
import { drawElementTile, TILE_W, TILE_H } from './tile.js';
import { createRadialPlot, radialMode } from './radial.js';
import { paneViews } from './selection.js';

const params = new URLSearchParams(location.search);
const $ = sel => document.querySelector(sel);

// Paper sizes and margins (inches). The image is 300 pixels per inch. A
// pair fills the content width at 4:3 (two side-by-side panes); a single
// atom gets a box the size of one of those panes, centred, so it prints
// at the same size and with the same overlay as in a pair. Sizes are
// multiples of 4 pixels so every grain divides them.
const PAPERS = {
  letter: { size: 'letter', w: 8.5, h: 11, margin: 0.5 },
  a4: { size: 'A4', w: 210 / 25.4, h: 297 / 25.4, margin: 12 / 25.4 },
};
const IMAGE_DPI = 300;
const CSS_DPI = 96;             // CSS px per inch
const ASPECT = 3 / 4;           // box height / width
const TILE_SCALE = 3;           // tile canvas pixels per CSS px (print-sharp)

function num(name, dflt) {
  const v = parseFloat(params.get(name));
  return Number.isFinite(v) ? v : dflt;
}

function atomParam(name) {
  const z = parseInt(params.get(name), 10);
  const el = Number.isFinite(z) && z >= 1 && z <= ELEMENTS.length ? elementByZ(z) : null;
  return el && el.supported ? el : null;
}

function orientationParam() {
  const q = (params.get('q') || '').split(',').map(Number);
  return q.length === 4 && q.every(Number.isFinite) && Math.hypot(...q) > 0
    ? qNormalize(q) : defaultOrientation();
}

const A = atomParam('z') || elementByZ(6);
const Bel = atomParam('z2');
const B = Bel && Bel !== A ? Bel : null;
const atoms = B ? [A, B] : [A];
const hidden = new Set((params.get('hide') || '').split(',').filter(Boolean).map(Number));

const state = {
  element: A,
  second: B,
  orientation: orientationParam(),
  zoom: num('zoom', 1) > 0 ? num('zoom', 1) : 1,
  density: num('density', 1) > 0 ? num('density', 1) : 1,
  kMax: num('kmax', 4),
  gamma: num('gamma', 0.5),
  equalize: params.get('equalize') === '1',
  fixedScale: params.get('scale') === 'fixed',
  cut: params.get('cut') === '1',
  color: params.get('color') === 'blend' ? 'blend' : 'nearest',
  streaks: streakParam(params.get('streaks')),
  overlay: params.get('overlay') !== '0',
  radialMode: ['shell', 'subshell', 'compare'].includes(params.get('radial')) ? params.get('radial') : 'auto',
  mode: 'white',
  visible: ORBITALS.map((o, i) => !hidden.has(i)),
  seed: 0,
};
// A pair shares one scale (the larger atom's), as in the app.
state.baseHalf = Math.max(...atoms.map(el => baseHalf(el, state.fixedScale)));

const opts = {
  paper: PAPERS[params.get('paper')] ? params.get('paper') : 'letter',
  grain: [1, 2, 4].includes(+params.get('grain')) ? +params.get('grain') : 2,
};
const fixedSeed = params.has('seed') ? parseInt(params.get('seed'), 10) >>> 0 : null;

// ------------------------------------------------------------ page size
// Image size (px) and its width on paper (CSS px) for this paper.
function imageSize(paper) {
  const p = PAPERS[paper];
  const contentIn = p.w - 2 * p.margin;
  const pairW = Math.floor((contentIn * IMAGE_DPI) / 4) * 4;
  const H = Math.floor((pairW * ASPECT) / 4) * 4;
  if (B) return { W: pairW, H, cssW: contentIn * CSS_DPI };
  const W = splitViews(pairW, H, 4, 2)[0].w;           // one pane of a pair
  return { W, H, cssW: (W / pairW) * contentIn * CSS_DPI };
}

function applyPaper() {
  const p = PAPERS[opts.paper];
  $('.sim').style.width = `${imageSize(opts.paper).cssW}px`;
  $('#page-size').textContent =
    `@page { size: ${p.size}; margin: ${p.margin}in; }\n` +
    `:root { --paper-w: ${p.w}in; --paper-h: ${p.h}in; --margin: ${p.margin}in; }`;
}

// ------------------------------------------------------------ rendering
let renderer = null;
const glCanvas = document.createElement('canvas');
const overlay = document.createElement('canvas');
const compose = document.createElement('canvas');
let imageUrl = null;

// Render one frame into the box image. Returns a promise that settles once
// the image has loaded (so printing right after shows it).
async function renderImage() {
  const { W, H, cssW } = imageSize(opts.paper);
  for (const c of [glCanvas, overlay, compose]) { c.width = W; c.height = H; }
  state.seed = fixedSeed ?? ((Math.random() * 4294967296) >>> 0);
  const k = opts.grain;
  const rects = splitViews(W, H, k, atoms.length);
  const views = paneViews(state).map((v, i) => ({ ...v, rect: rects[i] }));   // visibility shared, as in the app
  const geos = renderer.renderViews(state, k, views);
  // Overlay at the proportions it has in a view of the box's CSS width.
  drawOverlay(overlay, state, views.map((v, i) => ({ element: v.element, geo: geos[i] })), W / cssW);
  const ctx = compose.getContext('2d');
  ctx.drawImage(glCanvas, 0, 0);
  ctx.drawImage(overlay, 0, 0);

  const blob = await new Promise(r => compose.toBlob(r, 'image/png'));
  if (imageUrl) URL.revokeObjectURL(imageUrl);
  imageUrl = URL.createObjectURL(blob);
  const img = $('#sim-img');
  await new Promise((resolve, reject) => {
    img.onload = resolve;
    img.onerror = reject;
    img.src = imageUrl;
  });
  img.alt = `Electron orbitals of ${atoms.map(el => el.name).join(' and ')}: one frame of sampled detections`;
  $('#caption').textContent = caption(W, H);
}

function caption(W, H) {
  const k = opts.grain;
  const parts = [
    'One frame',
    `exposure ${state.density.toFixed(state.density < 0.1 ? 3 : 2)}` + (state.equalize ? ' (equalized, not to scale)' : ''),
    `full brightness at ${state.kMax} counts`,
    `gamma ${state.gamma.toFixed(2)}`,
    state.color === 'nearest' ? 'nearest color' : 'blended color',
  ];
  if (state.cut) parts.push('cut away: half behind the nucleus');
  if (state.streaks && state.color === 'nearest') {
    parts.push({ local: 'momentum lines: local spread', axis: 'momentum lines: local main axis',
      global: 'momentum lines: whole-orbital sample' }[state.streaks]);
  }
  if (state.fixedScale) parts.push('fixed scale');
  if (state.zoom !== 1) parts.push(`zoom ${+state.zoom.toFixed(2)}×`);
  parts.push(`${IMAGE_DPI / k} samples/in (${W / k} × ${H / k})`, `seed ${state.seed}`);
  return parts.join(' · ');
}

// ------------------------------------------------------------ details
function h(tag, attrs = {}, ...children) {
  const el = document.createElement(tag);
  for (const [key, v] of Object.entries(attrs)) {
    if (key === 'class') el.className = v;
    else if (key === 'html') el.innerHTML = v;
    else el.setAttribute(key, v);
  }
  for (const c of children) el.append(c);
  return el;
}

function tile(el, scale = 1) {
  const c = h('canvas', { class: 'tile', role: 'img', 'aria-label': `${el.name}, atomic number ${el.z}, mass ${el.mass}` });
  c.width = Math.round(TILE_W * scale * TILE_SCALE);
  c.height = Math.round(TILE_H * scale * TILE_SCALE);
  c.style.width = `${TILE_W * scale}px`;
  c.style.height = `${TILE_H * scale}px`;
  drawElementTile(c.getContext('2d'), el, 0, 0, scale * TILE_SCALE);
  return c;
}

const ion = el => `${el.ionizationEnergy.toFixed(3)} eV`;
const radius = el => `${el.atomicRadius.toFixed(0)} pm`;

function renderCard() {
  const card = $('#element-card');
  if (!B) {
    const row = (k, v) => h('div', { class: 'info-row' }, h('span', { class: 'k' }, k), h('span', { class: 'v' }, v));
    card.replaceChildren(tile(A), h('div', { class: 'card-info' },
      h('div', { class: 'card-name' }, A.name),
      h('div', { class: 'card-cat' }, A.category),
      row('Configuration', A.config),
      row('Period / Group', periodGroup(A)),
      row('Electronegativity', en(A)),
      row('1st ionization', ion(A)),
      row('Atomic radius (Rahm)', radius(A))));
    return;
  }
  card.classList.add('cmp');
  const side = (el, badge) => h('div', { class: 'cmp-atom' },
    h('span', { class: 'badge' }, badge), tile(el, 0.6),
    h('div', {}, h('div', { class: 'card-name' }, el.name), h('div', { class: 'card-cat' }, el.category)));
  const row = (k, f) => h('tr', {}, h('th', {}, k), h('td', {}, f(A)), h('td', {}, f(B)));
  card.replaceChildren(
    h('div', { class: 'cmp-head' }, side(A, 'A'), side(B, 'B')),
    h('table', { class: 'cmp-table' },
      h('thead', {}, h('tr', {}, h('th', {}), h('th', {}, `A · ${A.symbol}`), h('th', {}, `B · ${B.symbol}`))),
      h('tbody', {},
        row('Configuration', el => el.configShort),
        row('Period / Group', periodGroup),
        row('Electronegativity', en),
        row('1st ionization', ion),
        row('Atomic radius (Rahm)', radius))));
}

// Occupied subshells (of either atom), grouped by shell n, as in the app's
// legend fully expanded to subshells. Subshells hidden from the image are
// marked; a note lists them (and single hidden orbitals).
function renderLegend() {
  const cmp = !!B;
  const occOf = (el, idx) => idx.reduce((a, i) => a + el.occupancy[i], 0);
  const owner = i => atoms.find(el => el.occupancy[i] > 0) || A;
  const th = (text, cls = '') => h('th', { class: cls, html: text });
  const head = cmp
    ? [th(''), th('Subshell'), th(`e⁻ ${A.symbol}`, 'num'), th(`e⁻ ${B.symbol}`, 'num'),
      th(`r<sub>rms</sub> ${A.symbol}`, 'num'), th(`r<sub>rms</sub> ${B.symbol}`, 'num')]
    : [th(''), th('Subshell'), th('e⁻', 'num'), th('Z<sub>eff</sub>', 'num'), th('r<sub>rms</sub> pm', 'num')];
  const cols = head.length;

  const rows = [];
  const hiddenNames = [];
  const occupied = SUBSHELLS.filter(sh => atoms.some(el => el.subshells[sh.id]));
  for (const n of [...new Set(occupied.map(sh => sh.n))].sort((a, b) => a - b)) {
    const shells = occupied.filter(sh => sh.n === n);
    const idxAll = shells.flatMap(sh => orbitalsOf(sh.id));
    const totals = atoms.map(el => occOf(el, idxAll));
    rows.push(h('tr', { class: 'shell' },
      h('td', { colspan: '2', html: `n = ${n} <span class="shell-letter">${SHELL_LETTER[n]}</span>` }),
      ...totals.map(t => h('td', { class: 'num' }, t ? String(t) : '—')),
      ...Array.from({ length: cols - 2 - totals.length }, () => h('td'))));
    for (const sh of shells) {
      const idx = orbitalsOf(sh.id);
      const occ = atoms.map(el => occOf(el, idx));
      // Hidden (in both panes): of the orbitals either atom occupies.
      const occAny = idx.filter(i => atoms.some(el => el.occupancy[i] > 0));
      const hid = occAny.filter(i => hidden.has(i));
      const allHidden = occAny.length > 0 && hid.length === occAny.length;
      if (allHidden) hiddenNames.push(sh.id);
      else hiddenNames.push(...hid.map(i => ORBITALS[i].label));
      const label = h('td', { class: 'orb-label' }, sh.id);
      if (allHidden) label.append(h('span', { class: 'tag' }, ' hidden'));
      else if (hid.length) label.append(h('span', { class: 'tag' }, ' partly hidden'));
      const numbers = cmp
        ? [occ[0] ? String(occ[0]) : '—', occ[1] ? String(occ[1]) : '—', rmsPm(A, sh), rmsPm(B, sh)]
        : [occ[0] ? `${occ[0]}/${CAPACITY[sh.l]}` : '—', occ[0] ? A.zeffBySubshell[sh.id].toFixed(3) : '', rmsPm(A, sh)];
      rows.push(h('tr', { class: 'sub' + (allHidden ? ' hidden' : '') },
        h('td', { class: 'swatches' }, ...idx.map(i =>
          h('span', { class: 'swatch', style: `background:${orbitalCss(ORBITALS[i], owner(i))}` }))),
        label, ...numbers.map(v => h('td', { class: 'num' }, v))));
    }
  }
  $('#legend').replaceChildren(h('thead', {}, h('tr', {}, ...head)), h('tbody', {}, ...rows));
  const note = $('#hidden-note');
  note.hidden = hiddenNames.length === 0;
  note.innerHTML = hiddenNames.length
    ? `Not shown in the image: ${hiddenNames.join(', ')}.` : '';
}

function renderRadial() {
  const { W, H } = imageSize(opts.paper);
  const plot = createRadialPlot($('#radial'), null, $('#radial-legend'), { theme: 'light', dpr: 4 });
  // x range: A's view radius (the shorter half-dimension of a view at this
  // zoom), as in the app's inset, whatever the box's shape.
  const side = Math.min(W, H);
  plot.draw(state, viewGeometry(A, side, side, state.zoom, state.fixedScale));
  if (radialMode(state, A) === 'compare') $('#radial').setAttribute('aria-label',
    `Total radial electron density of ${A.name} and ${B.name} versus distance from the nucleus`);
}

function renderText() {
  const name = el => `${el.name} (${el.symbol}) · Z = ${el.z}`;
  for (const t of document.querySelectorAll('.title')) t.textContent = B ? `${name(A)}  vs  ${name(B)}` : name(A);
  document.title = `${atoms.map(el => el.symbol).join(' vs ')} · Electron Orbitals · Print`;
  $('#date').textContent = new Date().toLocaleDateString(undefined, { year: 'numeric', month: 'long', day: 'numeric' });
  const notes = modelNotes(atoms);
  $('#model-notes').hidden = notes.length === 0;
  $('#model-notes').textContent = notes.join(' ');
}

// ------------------------------------------------------------ toolbar
function setStatus(text) { $('#status').textContent = text; }

// Reflect the options in the URL so a reload keeps them.
function syncUrl() {
  const p = new URLSearchParams(location.search);
  p.set('paper', opts.paper);
  p.set('grain', String(opts.grain));
  history.replaceState(null, '', `${location.pathname}?${p}`);
}

async function rerender(what) {
  setStatus('Rendering…');
  try {
    await renderImage();
    if (what === 'paper') renderRadial();
    setStatus('');
  } catch (err) {
    showError(err);
  }
}

function showError(err) {
  $('#error').textContent = err.message || String(err);
  $('#error').hidden = false;
  setStatus('');
}

$('#paper').value = opts.paper;
$('#grain').value = String(opts.grain);
$('#paper').addEventListener('change', e => { opts.paper = e.target.value; applyPaper(); syncUrl(); rerender('paper'); });
$('#grain').addEventListener('change', e => { opts.grain = +e.target.value; syncUrl(); rerender(); });
$('#resample').addEventListener('click', () => rerender());
$('#do-print').addEventListener('click', () => window.print());

applyPaper();
renderText();
renderCard();
renderLegend();
renderRadial();
try {
  renderer = createRenderer(glCanvas, await loadTables());
  await rerender();
  document.body.dataset.ready = '1';           // for tests / headless capture
} catch (err) {
  showError(err);
}

window.printPage = { state, opts };
