// Side panel (periodic table, element card, orbital legend, controls) and the
// 2D overlay drawn over the WebGL canvas (horizon ring, axes, scale bar).

import { ELEMENTS, ORBITALS } from './elements.js';
import { SUBSHELLS, orbitalsOf } from './orbitals.js';
import { qToMat3, meanR2, hslCss } from './math.js';
import { categoryClass, drawElementTile, renderTileCanvas, TILE_W } from './tile.js';

const BOHR_PM = 52.9177;
const CAPACITY = { 0: 2, 1: 6, 2: 10 };

// Swatch colour for an orbital: paler for low-chroma (outer-shell) orbitals,
// approximating how they read in the render.
function orbitalCss(o) {
  return hslCss(o.hue, Math.round(90 * o.chroma), Math.round(55 + 25 * (1 - o.chroma)));
}

const $ = sel => document.querySelector(sel);

function h(tag, attrs = {}, ...children) {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === 'class') el.className = v;
    else if (k === 'html') el.innerHTML = v;
    else if (k.startsWith('on')) el.addEventListener(k.slice(2), v);
    else el.setAttribute(k, v);
  }
  for (const c of children) el.append(c);
  return el;
}


export function buildUI(app) {
  // ---------------------------------------------------- periodic table
  const pt = $('#ptable');
  for (const el of ELEMENTS) {
    const cell = h('button', {
      class: `pt-cell ${categoryClass(el.category)}`,
      title: el.name,
      'data-z': el.z,
      style: `grid-column: ${el.group}; grid-row: ${el.period}`,
      onclick: () => app.setElement(el.z),
    },
    h('span', { class: 'pt-z' }, String(el.z)),
    h('span', { class: 'pt-sym' }, el.symbol));
    pt.append(cell);
  }

  // ---------------------------------------------------- controls
  const exposure = $('#exposure');
  exposure.addEventListener('input', () => app.setDensity(sliderToDensity(+exposure.value)));
  const kmax = $('#kmax');
  kmax.addEventListener('input', () => app.setKMax(+kmax.value));
  const gamma = $('#gamma');
  gamma.addEventListener('input', () => app.setGamma(+gamma.value));
  $('#equalize').addEventListener('change', e => app.setEqualize(e.target.checked));
  $('#fixed-scale').addEventListener('change', e => app.setFixedScale(e.target.checked));
  for (const r of document.querySelectorAll('input[name=bgmode]')) {
    r.addEventListener('change', () => r.checked && app.setMode(r.value));
  }
  $('#pause').addEventListener('click', () => app.setPaused(!app.state.paused));
  $('#step').addEventListener('click', () => app.step());
  $('#debug-seed').addEventListener('change', e => app.setDebugSeed(e.target.checked));
  $('#show-overlay').addEventListener('change', e => app.setOverlay(e.target.checked));
  $('#reset-view').addEventListener('click', () => app.resetView());
  $('#save-png').addEventListener('click', () => app.savePng());
  $('#show-all').addEventListener('click', () => app.showAll());
  $('#show-valence').addEventListener('click', () => app.showValence());

  window.addEventListener('keydown', e => {
    if (e.target.tagName === 'INPUT') return;
    if (e.key === 'ArrowRight') app.setElement(Math.min(ELEMENTS.length, app.state.element.z + 1));
    else if (e.key === 'ArrowLeft') app.setElement(Math.max(1, app.state.element.z - 1));
    else if (e.key === ' ') { e.preventDefault(); app.setPaused(!app.state.paused); }
    else if (e.key === 'r') app.resetView();
    else if (e.key === 's') app.savePng();
  });
}

// Exposure slider is logarithmic: value in [0, 100] -> density in [0.01, 3].
export function sliderToDensity(v) { return 0.01 * Math.pow(300, v / 100); }
export function densityToSlider(d) { return 100 * Math.log(d / 0.01) / Math.log(300); }

export function syncControls(state) {
  $('#exposure').value = densityToSlider(state.density);
  $('#exposure-val').textContent = state.density.toFixed(state.density < 0.1 ? 3 : 2);
  $('#kmax').value = state.kMax;
  $('#kmax-val').textContent = state.kMax.toFixed(1);
  $('#gamma').value = state.gamma;
  $('#gamma-val').textContent = state.gamma.toFixed(2);
  $('#equalize').checked = state.equalize;
  $('#fixed-scale').checked = state.fixedScale;
  $('#exposure').title = state.equalize
    ? 'Peak expected detections per pixel = exposure x 20'
    : 'Mean detections per pixel per frame (if the whole atom were in view)';
  for (const r of document.querySelectorAll('input[name=bgmode]')) r.checked = r.value === state.mode;
  $('#pause').textContent = state.paused ? 'Play' : 'Pause';
  $('#step').disabled = !state.paused;
  $('#debug-seed').checked = state.debugSeed;
  $('#show-overlay').checked = state.overlay;
  document.body.classList.toggle('mode-white', state.mode === 'white');
}

export function updateElement(app) {
  const { state } = app;
  const el = state.element;
  for (const c of document.querySelectorAll('.pt-cell')) {
    c.classList.toggle('active', +c.dataset.z === el.z);
  }

  const tile = h('canvas', { class: 'card-symbol', role: 'img',
    'aria-label': `${el.name}, atomic number ${el.z}, mass ${el.mass}` });
  renderTileCanvas(tile, el);
  $('#element-card').replaceChildren(
    tile,
    h('div', { class: 'card-info' },
      h('div', { class: 'card-name' }, el.name),
      h('div', { class: 'card-cat' }, el.category),
      infoRow('Configuration', el.configShort, el.config),
      infoRow('Period / Group', `${el.period} / ${el.group}`),
      infoRow('Electronegativity', el.electronegativity == null ? '—' : el.electronegativity.toFixed(2)),
      infoRow('1st ionization', `${el.ionizationEnergy.toFixed(3)} eV`)),
  );

  // Subshells of every period up to this element's (at least 2), grouped;
  // p and d rows expand to their individual orbitals.
  const maxPeriod = Math.max(2, el.period);
  const rows = [];
  for (const sh of SUBSHELLS.filter(x => x.period <= maxPeriod)) {
    const idx = orbitalsOf(sh.id);
    const occ = idx.reduce((a, i) => a + el.occupancy[i], 0);
    const z = el.zeffBySubshell[sh.id];
    const r = occ ? Math.sqrt(meanR2(sh.n, sh.l, z)) : null;
    const shown = idx.filter(i => state.visible[i]).length;

    const cb = h('input', { type: 'checkbox', title: 'Show / hide this subshell' });
    cb.checked = shown === idx.length;
    cb.indeterminate = shown > 0 && shown < idx.length;
    cb.disabled = occ === 0;
    cb.addEventListener('change', () => app.setVisibleMany(idx, cb.checked));

    const expandable = sh.l > 0 && occ > 0;
    const open = expandable && app.expanded.has(sh.id);
    const caret = expandable
      ? h('button', { class: 'caret', title: open ? 'Collapse' : 'Show individual orbitals',
          'aria-expanded': String(open), onclick: () => app.toggleExpanded(sh.id) }, open ? '▾' : '▸')
      : h('span', { class: 'caret' });

    rows.push(h('tr', { class: occ ? '' : 'empty' },
      h('td', {}, cb),
      h('td', { class: 'swatches' }, ...idx.map(i =>
        h('span', { class: 'swatch', style: `background:${orbitalCss(ORBITALS[i])}` }))),
      h('td', { class: 'orb-label' }, caret, sh.id),
      h('td', { class: 'num' }, occ ? `${occ}/${CAPACITY[sh.l]}` : '—'),
      h('td', { class: 'num' }, occ ? z.toFixed(3) : ''),
      h('td', { class: 'num' }, occ ? (r * BOHR_PM).toFixed(0) : '')));

    if (open) {
      for (const i of idx) {
        const o = ORBITALS[i];
        const n = el.occupancy[i];
        const ocb = h('input', { type: 'checkbox', title: 'Show / hide this orbital' });
        ocb.checked = state.visible[i];
        ocb.disabled = n === 0;
        ocb.addEventListener('change', () => app.setVisible(i, ocb.checked));
        rows.push(h('tr', { class: 'child' + (n ? '' : ' empty') },
          h('td', {}, ocb),
          h('td', { class: 'swatches' }, h('span', { class: 'swatch', style: `background:${orbitalCss(o)}` })),
          h('td', { class: 'orb-label', html: o.label }),
          h('td', { class: 'pips' }, n === 0 ? '—' : (n === 2 ? '↑↓' : '↑')),
          h('td', {}), h('td', {})));
      }
    }
  }
  $('#legend tbody').replaceChildren(...rows);

  // Partly filled p/d subshells: which orbitals hold the electrons is one of
  // several equivalent (degenerate) choices; say which convention is shown.
  const partial = SUBSHELLS.filter(sh => sh.l > 0 &&
    (el.subshells[sh.id] || 0) > 0 && el.subshells[sh.id] < CAPACITY[sh.l]);
  const note = $('#fill-note');
  note.hidden = partial.length === 0;
  if (partial.length) {
    const order = partial.map(sh => (sh.l === 1 ? 'p: x, y, z' : 'd: xy, yz, xz, x²−y², z²'));
    note.textContent = `Partly filled ${partial.map(sh => sh.id).join(', ')}: shown is one of several ` +
      `equivalent fillings (Hund's rule in the order ${[...new Set(order)].join('; ')}).`;
  }
}

function infoRow(k, v, title) {
  const val = h('span', { class: 'v' }, v);
  if (title) val.title = title;
  return h('div', { class: 'info-row' }, h('span', { class: 'k' }, k), val);
}

export function updateStatus(state, geo, fps) {
  $('#status').textContent =
    `horizon radius = 2σ = ${(geo.horizon * BOHR_PM).toFixed(0)} pm (${geo.horizon.toFixed(2)} a₀)` +
    (fps ? `  ·  ${fps.toFixed(0)} fps` : '');
}

// ------------------------------------------------------------ overlay
export function drawOverlay(canvas, state, geo) {
  const ctx = canvas.getContext('2d');
  const w = canvas.width, hgt = canvas.height;
  ctx.clearRect(0, 0, w, hgt);
  if (!state.overlay) return;
  const dpr = w / canvas.clientWidth || 1;
  const fg = state.mode === 'white' ? 'rgba(0,0,0,0.55)' : 'rgba(255,255,255,0.55)';
  const pxPerBohr = hgt / (2 * geo.halfExtent[1]);
  const cx = w / 2, cy = hgt / 2;

  // Canvas text/line state persists between frames, so isolate this draw
  // and start from known settings.
  ctx.save();
  ctx.font = `${12 * dpr}px system-ui, sans-serif`;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'alphabetic';
  ctx.lineWidth = 1 * dpr;

  // Horizon ring (2 sigma).
  const r = geo.horizon * pxPerBohr;
  ctx.strokeStyle = fg;
  ctx.setLineDash([6 * dpr, 6 * dpr]);
  ctx.beginPath();
  ctx.arc(cx, cy, r, 0, 2 * Math.PI);
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = fg;
  ctx.fillText(`${Math.round(geo.horizon * BOHR_PM)} pm`, cx + r * 0.707 + 4 * dpr, cy - r * 0.707 - 4 * dpr);

  // Axes gizmo (atom frame axes projected to the screen), bottom-left.
  const m = qToMat3(state.orientation);   // atom -> view
  const gx = 50 * dpr, gy = hgt - 50 * dpr, gl = 34 * dpr;
  const axes = [0, 1, 2]
    .map(a => ({ a, v: [m[0][a], m[1][a], m[2][a]] }))
    .sort((p, q) => p.v[2] - q.v[2]);      // draw back-to-front
  for (const { a, v } of axes) {
    const o = ORBITALS.find(orb => orb.subshell === '2p' && orb.axis === a);
    ctx.strokeStyle = ctx.fillStyle = hslCss(o.hue, 90, state.mode === 'white' ? 40 : 60);
    ctx.lineWidth = 2 * dpr;
    ctx.beginPath();
    ctx.moveTo(gx, gy);
    ctx.lineTo(gx + v[0] * gl, gy - v[1] * gl);
    ctx.stroke();
    ctx.fillText('xyz'[a], gx + v[0] * (gl + 9 * dpr) - 3 * dpr, gy - v[1] * (gl + 9 * dpr) + 4 * dpr);
  }

  // Scale bar: 1 Å, or 1 a0 if that is too short / long.
  let len = 100 / BOHR_PM, label = '1 Å';
  if (len * pxPerBohr > w * 0.3) { len = 1; label = '1 a₀'; }
  const bx = w - 30 * dpr - len * pxPerBohr, by = hgt - 30 * dpr;
  ctx.strokeStyle = ctx.fillStyle = fg;
  ctx.lineWidth = 2 * dpr;
  ctx.beginPath();
  ctx.moveTo(bx, by); ctx.lineTo(bx + len * pxPerBohr, by);
  ctx.moveTo(bx, by - 4 * dpr); ctx.lineTo(bx, by + 4 * dpr);
  ctx.moveTo(bx + len * pxPerBohr, by - 4 * dpr); ctx.lineTo(bx + len * pxPerBohr, by + 4 * dpr);
  ctx.stroke();
  ctx.textAlign = 'center';
  ctx.textBaseline = 'bottom';
  ctx.fillText(label, bx + (len * pxPerBohr) / 2, by - 7 * dpr);   // centred above the bar

  // Element tile (same rendering as the side-panel card) with name,
  // category and configuration beside it, top-left.
  const el = state.element;
  const tx = 16 * dpr, ty = 16 * dpr;
  drawElementTile(ctx, el, tx, ty, dpr);
  const lx = tx + (TILE_W + 12) * dpr;
  ctx.fillStyle = state.mode === 'white' ? 'rgba(0,0,0,0.8)' : 'rgba(255,255,255,0.85)';
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.font = `600 ${18 * dpr}px system-ui, sans-serif`;
  ctx.fillText(el.name, lx, ty + 4 * dpr);
  ctx.fillStyle = fg;
  ctx.font = `${12 * dpr}px system-ui, sans-serif`;
  ctx.fillText(el.category, lx, ty + 28 * dpr);
  ctx.fillText(el.configShort, lx, ty + 46 * dpr);
  ctx.restore();
}
