// Side panel (periodic table, element card, orbital legend, controls) and the
// 2D overlay drawn over the WebGL canvas (horizon ring, axes, scale bar).

import { ELEMENTS, ORBITALS } from './elements.js';
import { SUBSHELL } from './orbitals.js';
import { qToMat3, meanR2, hslCss } from './math.js';

const BOHR_PM = 52.9177;
// Compact periodic table columns: groups 1, 2, 13..18.
const GROUP_COL = { 1: 1, 2: 2, 13: 3, 14: 4, 15: 5, 16: 6, 17: 7, 18: 8 };

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

function categoryClass(cat) {
  return 'cat-' + cat.toLowerCase().replace(/[^a-z]+/g, '-');
}

export function buildUI(app) {
  // ---------------------------------------------------- periodic table
  const pt = $('#ptable');
  for (const el of ELEMENTS) {
    const cell = h('button', {
      class: `pt-cell ${categoryClass(el.category)}`,
      title: el.name,
      'data-z': el.z,
      style: `grid-column: ${GROUP_COL[el.group]}; grid-row: ${el.period}`,
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

  window.addEventListener('keydown', e => {
    if (e.target.tagName === 'INPUT') return;
    if (e.key === 'ArrowRight') app.setElement(Math.min(10, app.state.element.z + 1));
    else if (e.key === 'ArrowLeft') app.setElement(Math.max(1, app.state.element.z - 1));
    else if (e.key === ' ') { e.preventDefault(); app.setPaused(!app.state.paused); }
    else if (e.key === 'r') app.resetView();
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

  $('#element-card').replaceChildren(
    h('div', { class: `card-symbol ${categoryClass(el.category)}` },
      h('span', { class: 'card-z' }, String(el.z)),
      h('span', { class: 'card-sym' }, el.symbol),
      h('span', { class: 'card-mass' }, el.mass.toFixed(3))),
    h('div', { class: 'card-info' },
      h('div', { class: 'card-name' }, el.name),
      h('div', { class: 'card-cat' }, el.category),
      infoRow('Configuration', el.config),
      infoRow('Period / Group', `${el.period} / ${el.group}`),
      infoRow('Electronegativity', el.electronegativity == null ? '—' : el.electronegativity.toFixed(2)),
      infoRow('1st ionization', `${el.ionizationEnergy.toFixed(3)} eV`)),
  );

  // Show the subshells of every period up to this element's (at least 2).
  const maxPeriod = Math.max(2, el.period);
  const rows = ORBITALS.map((o, i) => [o, i])
    .filter(([o]) => SUBSHELL[o.subshell].period <= maxPeriod)
    .map(([o, i]) => {
    const occ = el.occupancy[i];
    const z = el.zeff[i];
    const pips = occ === 0 ? '—' : (occ === 2 ? '↑↓' : '↑');
    const r = occ ? Math.sqrt(meanR2(o.n, o.l, z)) : null;
    const cb = h('input', { type: 'checkbox', title: 'Show / hide this orbital' });
    cb.checked = state.visible[i];
    cb.disabled = occ === 0;
    cb.addEventListener('change', () => app.setVisible(i, cb.checked));
    return h('tr', { class: occ ? '' : 'empty' },
      h('td', {}, cb),
      h('td', {}, h('span', { class: 'swatch', style: `background:${hslCss(o.hue)}` })),
      h('td', { class: 'orb-label', html: o.label }),
      h('td', { class: 'pips' }, pips),
      h('td', { class: 'num' }, occ ? z.toFixed(3) : ''),
      h('td', { class: 'num' }, occ ? (r * BOHR_PM).toFixed(0) : ''));
  });
  $('#legend tbody').replaceChildren(...rows);
}

function infoRow(k, v) {
  return h('div', { class: 'info-row' }, h('span', { class: 'k' }, k), h('span', { class: 'v' }, v));
}

export function updateStatus(state, geo, fps) {
  $('#status').textContent =
    `frame ${state.frame}  ·  seed ${state.seed >>> 0}  ·  horizon 2σ = ${geo.horizon.toFixed(2)} a₀ ` +
    `(${(geo.horizon * BOHR_PM).toFixed(0)} pm)` + (fps ? `  ·  ${fps.toFixed(0)} fps` : '');
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

  ctx.font = `${12 * dpr}px system-ui, sans-serif`;
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
  ctx.fillText('2σ', cx + r * 0.707 + 4 * dpr, cy - r * 0.707 - 4 * dpr);

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
  ctx.fillText(label, bx, by - 8 * dpr);

  // Element label, top-left.
  ctx.font = `600 ${28 * dpr}px system-ui, sans-serif`;
  ctx.fillText(state.element.symbol, 20 * dpr, 44 * dpr);
  ctx.font = `${12 * dpr}px system-ui, sans-serif`;
  ctx.fillText(`${state.element.name} · ${state.element.config}`, 20 * dpr, 64 * dpr);
}
