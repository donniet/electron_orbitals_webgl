// Side panel (periodic table, element card, orbital legend, controls) and the
// 2D overlay drawn over the WebGL canvas (horizon ring, axes, scale bar).

import { ELEMENTS, ORBITALS, stepElement } from './elements.js';
import { SUBSHELLS, orbitalsOf, orbitalChroma, valenceSubshells } from './orbitals.js';
import { qToMat3, meanR2, hslCss } from './math.js';
import { categoryClass, drawElementTile, renderTileCanvas, TILE_W } from './tile.js';

const BOHR_PM = 52.9177;
const CAPACITY = { 0: 2, 1: 6, 2: 10, 3: 14 };

// Swatch colour for an orbital: paler for low-chroma (outer-shell) orbitals,
// approximating how they read in the render.
function orbitalCss(o, el) {
  const c = orbitalChroma(o, el);
  return hslCss(o.hue, Math.round(90 * c), Math.round(55 + 25 * (1 - c)));
}

const SHELL_LETTER = { 1: 'K', 2: 'L', 3: 'M', 4: 'N', 5: 'O', 6: 'P' };

// Shells expanded by default: those holding valence subshells, or every
// shell for light atoms (few subshells).
export function defaultExpandedShells(el) {
  const occupied = Object.keys(el.subshells);
  if (occupied.length <= 6) return new Set(occupied.map(id => +id[0]));
  return new Set([...valenceSubshells(el)].map(id => +id[0]));
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
  // Standard 18-column layout; the f-block (La-Lu) sits in its own row below,
  // behind a "57–71" marker in group 3. Every cell records its row so the
  // strip can collapse to the current element's row.
  const pt = $('#ptable');
  const F_ROW = 8;                                   // row 7 is a visual gap
  for (const el of ELEMENTS) {
    const fBlock = el.group == null;
    const row = fBlock ? F_ROW : el.period;
    const col = fBlock ? 3 + (el.z - 57) : el.group;
    const cell = h('button', {
      class: `pt-cell ${categoryClass(el.category)}` + (el.supported ? '' : ' unsupported') +
        (fBlock ? ' f-row' : ''),
      title: el.supported ? el.name : `${el.name} — not available: ${el.unsupportedReason}`,
      'data-z': el.z, 'data-row': row, 'data-col': col,
      // A collapsed table shows the neighbouring rows as thin "peek" slivers;
      // clicking one expands the table rather than picking a tiny target.
      onclick: e => (e.currentTarget.classList.contains('peek')
        ? app.setTableCollapsed(false) : app.setElement(el.z)),
    },
    h('span', { class: 'pt-z' }, String(el.z)),
    h('span', { class: 'pt-sym' }, el.symbol));
    cell.dataset.title = cell.title;
    pt.append(cell);
  }
  if (ELEMENTS.some(el => el.group == null)) {
    pt.append(h('div', { class: 'pt-cell pt-marker cat-lanthanide', 'data-row': 6, 'data-col': 3,
      title: 'Lanthanides (57–71), shown in the row below' }, h('span', { class: 'pt-sym' }, '57–71')));
  }
  const toggle = h('button', { id: 'ptable-toggle', onclick: () => app.setTableCollapsed(!app.tableCollapsed) });
  pt.append(toggle);

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
  for (const r of document.querySelectorAll('input[name=resmode]')) {
    r.addEventListener('change', () => r.checked && app.setResMode(r.value));
  }
  $('#pause').addEventListener('click', () => app.setPaused(!app.state.paused));
  $('#step').addEventListener('click', () => app.step());
  $('#debug-seed').addEventListener('change', e => app.setDebugSeed(e.target.checked));
  $('#show-overlay').addEventListener('change', e => app.setOverlay(e.target.checked));
  $('#reset-view').addEventListener('click', () => app.resetView());
  $('#save-png').addEventListener('click', () => app.savePng());
  $('#show-all').addEventListener('click', () => app.showAll());
  $('#show-valence').addEventListener('click', () => app.showValence());
  $('#radial-shells').addEventListener('click', () => app.setRadialMode('shell'));
  $('#radial-subshells').addEventListener('click', () => app.setRadialMode('subshell'));

  window.addEventListener('keydown', e => {
    if (e.target.tagName === 'INPUT') return;
    if (e.key === 'ArrowRight') app.setElement(stepElement(app.state.element.z, +1));
    else if (e.key === 'ArrowLeft') app.setElement(stepElement(app.state.element.z, -1));
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
  for (const r of document.querySelectorAll('input[name=resmode]')) r.checked = r.value === state.resMode;
  $('#pause').textContent = state.paused ? 'Play' : 'Pause';
  $('#step').disabled = !state.paused;
  $('#debug-seed').checked = state.debugSeed;
  $('#show-overlay').checked = state.overlay;
  document.body.classList.toggle('mode-white', state.mode === 'white');
}

export function updateTable(app) {
  const el = app.state.element;
  const current = el.group == null ? 8 : el.period;
  const cells = [...document.querySelectorAll('#ptable .pt-cell')];
  // Collapsed: the current row in full, plus thin "peek" slivers of the rows
  // directly above and below (in table order: periods, then the f-block row)
  // hinting that there is a whole table to expand.
  const rows = [...new Set(cells.map(c => +c.dataset.row))].sort((a, b) => a - b);
  const at = rows.indexOf(current);
  const above = rows[at - 1], below = rows[at + 1];
  for (const c of cells) {
    c.classList.toggle('active', +c.dataset.z === el.z);
    const row = +c.dataset.row;
    const peek = app.tableCollapsed && (row === above || row === below);
    c.classList.toggle('peek', peek);
    c.hidden = app.tableCollapsed && row !== current && !peek;
    if (c.tagName === 'BUTTON') c.tabIndex = peek ? -1 : 0;
    if (peek) c.title = 'Show the whole periodic table';
    else if (c.dataset.title) c.title = c.dataset.title;
    c.style.gridColumn = c.dataset.col;
    c.style.gridRow = !app.tableCollapsed ? row : row === above ? 1 : row === current ? 2 : 3;
  }
  const t = $('#ptable-toggle');
  t.textContent = app.tableCollapsed ? '▾' : '▴';
  t.title = app.tableCollapsed ? 'Show the whole periodic table' : 'Collapse to the current row';
  t.setAttribute('aria-expanded', String(!app.tableCollapsed));
  $('#ptable').classList.toggle('collapsed', app.tableCollapsed);
}

export function updateElement(app) {
  const { state } = app;
  const el = state.element;
  updateTable(app);

  const tile = h('canvas', { class: 'card-symbol', role: 'img',
    'aria-label': `${el.name}, atomic number ${el.z}, mass ${el.mass}` });
  renderTileCanvas(tile, el);
  $('#element-card').replaceChildren(
    tile,
    h('div', { class: 'card-info' },
      h('div', { class: 'card-name-row' },
        h('div', { class: 'card-name' }, el.name),
        h('span', { class: 'stepper' },
          stepButton('‹', 'Previous element (←)', stepElement(el.z, -1), app),
          stepButton('›', 'Next element (→)', stepElement(el.z, +1), app))),
      h('div', { class: 'card-cat' }, el.category),
      infoRow('Configuration', el.configShort, el.config),
      infoRow('Period / Group', `${el.period} / ${el.group ?? 'f-block'}`),
      infoRow('Electronegativity', el.electronegativity == null ? '—' : el.electronegativity.toFixed(2)),
      infoRow('1st ionization', `${el.ionizationEnergy.toFixed(3)} eV`)),
  );

  // Subshells of every period up to this element's (at least 2), grouped by
  // shell (n). Shell rows collapse; p/d/f rows expand to single orbitals.
  const maxPeriod = Math.max(2, el.period);
  const rows = [];
  const shown = SUBSHELLS.filter(x => x.period <= maxPeriod);
  for (const n of [...new Set(shown.map(x => x.n))].sort((a, b) => a - b)) {
    const shells = shown.filter(x => x.n === n);
    const idxAll = shells.flatMap(sh => orbitalsOf(sh.id));
    const occN = idxAll.reduce((a, i) => a + el.occupancy[i], 0);
    const visN = idxAll.filter(i => state.visible[i] && el.occupancy[i] > 0).length;
    const occupiedN = idxAll.filter(i => el.occupancy[i] > 0).length;
    const openShell = app.expandedShells.has(n);

    const scb = h('input', { type: 'checkbox', title: `Show / hide shell n = ${n}` });
    scb.checked = occupiedN > 0 && visN === occupiedN;
    scb.indeterminate = visN > 0 && visN < occupiedN;
    scb.disabled = occN === 0;
    scb.addEventListener('change', () => app.setVisibleMany(idxAll, scb.checked));
    rows.push(h('tr', { class: 'shell' + (occN ? '' : ' empty') },
      h('td', {}, scb),
      h('td', { class: 'orb-label', colspan: '2' },
        h('button', { class: 'caret', title: openShell ? 'Collapse shell' : 'Show subshells',
          'aria-expanded': String(openShell), onclick: () => app.toggleShell(n) }, openShell ? '▾' : '▸'),
        `n = ${n}`, h('span', { class: 'shell-letter' }, ` ${SHELL_LETTER[n]}`)),
      h('td', { class: 'num' }, occN ? String(occN) : '—'),
      h('td', {}), h('td', {})));
    if (!openShell) continue;

    for (const sh of shells) {
      const idx = orbitalsOf(sh.id);
      const occ = idx.reduce((a, i) => a + el.occupancy[i], 0);
      const z = el.zeffBySubshell[sh.id];
      const r = occ ? Math.sqrt(meanR2(sh.n, sh.l, z)) : null;
      const vis = idx.filter(i => state.visible[i]).length;

      const cb = h('input', { type: 'checkbox', title: 'Show / hide this subshell' });
      cb.checked = vis === idx.length;
      cb.indeterminate = vis > 0 && vis < idx.length;
      cb.disabled = occ === 0;
      cb.addEventListener('change', () => app.setVisibleMany(idx, cb.checked));

      const expandable = sh.l > 0 && occ > 0;
      const open = expandable && app.expanded.has(sh.id);
      const caret = expandable
        ? h('button', { class: 'caret', title: open ? 'Collapse' : 'Show individual orbitals',
            'aria-expanded': String(open), onclick: () => app.toggleExpanded(sh.id) }, open ? '▾' : '▸')
        : h('span', { class: 'caret' });

      rows.push(h('tr', { class: 'sub' + (occ ? '' : ' empty') },
        h('td', {}, cb),
        h('td', { class: 'swatches' }, ...idx.map(i =>
          h('span', { class: 'swatch', style: `background:${orbitalCss(ORBITALS[i], el)}` }))),
        h('td', { class: 'orb-label' }, caret, sh.id),
        h('td', { class: 'num' }, occ ? `${occ}/${CAPACITY[sh.l]}` : '—'),
        h('td', { class: 'num' }, occ ? z.toFixed(3) : ''),
        h('td', { class: 'num' }, occ ? (r * BOHR_PM).toFixed(0) : '')));

      if (open) {
        for (const i of idx) {
          const o = ORBITALS[i];
          const k = el.occupancy[i];
          const ocb = h('input', { type: 'checkbox', title: 'Show / hide this orbital' });
          ocb.checked = state.visible[i];
          ocb.disabled = k === 0;
          ocb.addEventListener('change', () => app.setVisible(i, ocb.checked));
          rows.push(h('tr', { class: 'child' + (k ? '' : ' empty') },
            h('td', {}, ocb),
            h('td', { class: 'swatches' }, h('span', { class: 'swatch', style: `background:${orbitalCss(o, el)}` })),
            h('td', { class: 'orb-label', html: o.label }),
            h('td', { class: 'pips' }, k === 0 ? '—' : (k === 2 ? '↑↓' : '↑')),
            h('td', {}), h('td', {})));
        }
      }
    }
  }
  $('#legend tbody').replaceChildren(...rows);

  // Partly filled p/d subshells: which orbitals hold the electrons is one of
  // several equivalent (degenerate) choices; say which convention is shown.
  const partial = SUBSHELLS.filter(sh => sh.l > 0 &&
    (el.subshells[sh.id] || 0) > 0 && el.subshells[sh.id] < CAPACITY[sh.l]);
  const HUND_ORDER = { 1: 'p: x, y, z', 2: 'd: xy, yz, xz, x²−y², z²',
    3: 'f: z³, xz², yz², xyz, z(x²−y²), x(x²−3y²), y(3x²−y²)' };
  const notes = [];
  if (partial.length) {
    const order = [...new Set(partial.map(sh => HUND_ORDER[sh.l]))].join('; ');
    notes.push(`Partly filled ${partial.map(sh => sh.id).join(', ')}: shown is one of several ` +
      `equivalent fillings (Hund's rule in the order ${order}).`);
  }
  if (el.z >= 55) {
    notes.push('Non-relativistic hydrogen-like model: relativistic contraction / expansion of ' +
      'heavy-atom orbitals (e.g. 6s, 5d) is not included.');
  }
  const note = $('#fill-note');
  note.hidden = notes.length === 0;
  note.textContent = notes.join(' ');
}

function stepButton(label, title, z, app) {
  const b = h('button', { class: 'step', title, 'aria-label': title }, label);
  if (!z) b.disabled = true;
  else b.addEventListener('click', () => app.setElement(z));
  return b;
}

function infoRow(k, v, title) {
  const val = h('span', { class: 'v' }, v);
  if (title) val.title = title;
  return h('div', { class: 'info-row' }, h('span', { class: 'k' }, k), val);
}

export function updateStatus(state, geo, fps, renderK = 1) {
  $('#status').textContent =
    `horizon radius = 2σ = ${(geo.horizon * BOHR_PM).toFixed(0)} pm (${geo.horizon.toFixed(2)} a₀)` +
    (fps ? `  ·  ${fps.toFixed(0)} fps` : '') +
    (renderK > 1 ? `  ·  1/${renderK} res` : '');
}

// ------------------------------------------------------------ overlay
// Fill text with a halo in the background colour so labels stay legible
// over the horizon ring and the speckle.
function haloText(ctx, text, x, y, halo, width) {
  ctx.save();
  ctx.lineJoin = 'round';
  ctx.strokeStyle = halo;
  ctx.lineWidth = width;
  ctx.strokeText(text, x, y);
  ctx.restore();
  ctx.fillText(text, x, y);
}

export function drawOverlay(canvas, state, geo) {
  const ctx = canvas.getContext('2d');
  const w = canvas.width, hgt = canvas.height;
  ctx.clearRect(0, 0, w, hgt);
  if (!state.overlay) return;
  const dpr = w / canvas.clientWidth || 1;
  const fg = state.mode === 'white' ? 'rgba(0,0,0,0.55)' : 'rgba(255,255,255,0.55)';
  const pxPerBohr = hgt / (2 * geo.halfExtent[1]);
  const cx = w / 2, cy = hgt / 2;
  // Compact overlay on narrow canvases (phones in portrait).
  const compact = canvas.clientWidth < 520;
  const k = compact ? 0.7 : 1;            // tile, gizmo and margins
  const kt = compact ? 0.8 : 1;           // label text
  const halo = state.mode === 'white' ? 'rgba(255,255,255,0.9)' : 'rgba(0,0,0,0.9)';
  const hw = 4 * dpr;

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
  haloText(ctx, `${Math.round(geo.horizon * BOHR_PM)} pm`, cx + r * 0.707 + 4 * dpr, cy - r * 0.707 - 4 * dpr, halo, hw);

  // Axes gizmo (atom frame axes projected to the screen), bottom-left.
  const m = qToMat3(state.orientation);   // atom -> view
  const gx = 50 * k * dpr, gy = hgt - 50 * k * dpr, gl = 34 * k * dpr;
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
  const bx = w - 30 * k * dpr - len * pxPerBohr, by = hgt - 30 * k * dpr;
  ctx.strokeStyle = ctx.fillStyle = fg;
  ctx.lineWidth = 2 * dpr;
  ctx.beginPath();
  ctx.moveTo(bx, by); ctx.lineTo(bx + len * pxPerBohr, by);
  ctx.moveTo(bx, by - 4 * dpr); ctx.lineTo(bx, by + 4 * dpr);
  ctx.moveTo(bx + len * pxPerBohr, by - 4 * dpr); ctx.lineTo(bx + len * pxPerBohr, by + 4 * dpr);
  ctx.stroke();
  ctx.textAlign = 'center';
  ctx.textBaseline = 'bottom';
  haloText(ctx, label, bx + (len * pxPerBohr) / 2, by - 7 * dpr, halo, hw);   // centred above the bar

  // Element tile (same rendering as the side-panel card) with name,
  // category and configuration beside it, top-left.
  const el = state.element;
  const tx = 16 * k * dpr, ty = 16 * k * dpr;
  drawElementTile(ctx, el, tx, ty, k * dpr);
  const lx = tx + (TILE_W * k + 10) * dpr;
  ctx.fillStyle = state.mode === 'white' ? 'rgba(0,0,0,0.8)' : 'rgba(255,255,255,0.85)';
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.font = `600 ${18 * kt * dpr}px system-ui, sans-serif`;
  haloText(ctx, el.name, lx, ty + 4 * kt * dpr, halo, hw);
  ctx.fillStyle = fg;
  ctx.font = `${12 * kt * dpr}px system-ui, sans-serif`;
  haloText(ctx, el.category, lx, ty + 28 * kt * dpr, halo, hw);
  haloText(ctx, el.configShort, lx, ty + 46 * kt * dpr, halo, hw);
  ctx.restore();
}
