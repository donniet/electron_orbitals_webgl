// Side panel (periodic table, element card, orbital legend, controls) and the
// 2D overlay drawn over the WebGL canvas (atomic-radius ring, axes, scale bar).

import { ELEMENTS, ORBITALS, stepElement } from './elements.js';
import { SUBSHELLS, orbitalsOf, orbitalChroma } from './orbitals.js';
import { atomsOf, occupiedAny, maxPeriod, collapsedRows, tableRow } from './selection.js';
import { qToMat3, meanR2, hslCss } from './math.js';
import { categoryClass, drawElementTile, renderTileCanvas, TILE_W } from './tile.js';

const BOHR_PM = 52.9177;
export const CAPACITY = { 0: 2, 1: 6, 2: 10, 3: 14 };

// Swatch color for an orbital: paler for low-chroma (outer-shell) orbitals,
// approximating how they read in the render.
export function orbitalCss(o, el) {
  const c = orbitalChroma(o, el);
  return hslCss(o.hue, Math.round(90 * c), Math.round(55 + 25 * (1 - c)));
}

export const SHELL_LETTER = { 1: 'K', 2: 'L', 3: 'M', 4: 'N', 5: 'O', 6: 'P' };


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
      // Shift-click, long-press, or a tap while "Compare" is armed shows a
      // second atom side by side; a plain click goes back to a single atom.
      onclick: e => {
        const c = e.currentTarget;
        if (c.dataset.longpressed) { delete c.dataset.longpressed; return; }   // handled on press
        if (c.classList.contains('peek')) app.setTableCollapsed(false);
        else if (e.shiftKey || app.pickingSecond) app.toggleSecond(el.z);
        else { app.clearSecond(); app.setElement(el.z); }
      },
    },
    h('span', { class: 'pt-z' }, String(el.z)),
    h('span', { class: 'pt-sym' }, el.symbol));
    cell.dataset.title = cell.title;
    onLongPress(cell, () => { if (!cell.classList.contains('peek')) app.toggleSecond(el.z); });
    pt.append(cell);
  }
  if (ELEMENTS.some(el => el.group == null)) {
    pt.append(h('div', { class: 'pt-cell pt-marker cat-lanthanide', 'data-row': 6, 'data-col': 3,
      title: 'Lanthanides (57–71), shown in the row below' }, h('span', { class: 'pt-sym' }, '57–71')));
  }
  const toggle = h('button', { id: 'ptable-toggle', onclick: () => app.setTableCollapsed(!app.tableCollapsed) });
  pt.append(toggle);
  pt.append(h('div', { id: 'pick-hint', class: 'pick-hint', hidden: '' }, 'Tap an element to compare'));

  // ---------------------------------------------------- controls
  const exposure = $('#exposure');
  exposure.addEventListener('input', () => app.setDensity(sliderToDensity(+exposure.value)));
  const kmax = $('#kmax');
  kmax.addEventListener('input', () => app.setKMax(+kmax.value));
  const gamma = $('#gamma');
  gamma.addEventListener('input', () => app.setGamma(+gamma.value));
  $('#equalize').addEventListener('change', e => app.setEqualize(e.target.checked));
  $('#fixed-scale').addEventListener('change', e => app.setFixedScale(e.target.checked));
  $('#cut').addEventListener('change', e => app.setCut(e.target.checked));
  for (const r of document.querySelectorAll('input[name=colormode]')) {
    r.addEventListener('change', e => { if (e.target.checked) app.setColor(e.target.value); });
  }
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
  $('#print-page').addEventListener('click', () => app.openPrintPage());
  $('#show-all').addEventListener('click', () => app.showAll());
  $('#show-valence').addEventListener('click', () => app.showValence());
  $('#radial-shells').addEventListener('click', () => app.setRadialMode('shell'));
  $('#radial-subshells').addEventListener('click', () => app.setRadialMode('subshell'));
  $('#radial-compare').addEventListener('click', () => app.setRadialMode('compare'));

  window.addEventListener('keydown', e => {
    if (e.target.tagName === 'INPUT') return;
    if (e.key === 'ArrowRight') app.setElement(stepElement(app.state.element.z, +1));
    else if (e.key === 'ArrowLeft') app.setElement(stepElement(app.state.element.z, -1));
    else if (e.key === ' ') { e.preventDefault(); app.setPaused(!app.state.paused); }
    else if (e.key === 'r') app.resetView();
    else if (e.key === 's') app.savePng();
    else if (e.key === 'p') app.openPrintPage();
    else if (e.key === 'c') app.setCut(!app.state.cut);
    else if (e.key === 'n') app.setColor(app.state.color === 'nearest' ? 'blend' : 'nearest');
    else if (e.key === 'Escape') { if (app.pickingSecond) app.togglePicking(); else app.clearSecond(); }
  });
}

// Long-press (touch: the equivalent of shift-click). Fires after `ms` unless
// the pointer lifts or moves > 10 px first; the click that follows is
// swallowed by the cell's click handler, and the browser's long-press menu
// is suppressed.
function onLongPress(el, fn, ms = 500) {
  let timer = null, x0 = 0, y0 = 0;
  const cancel = () => { if (timer) { clearTimeout(timer); timer = null; } };
  el.addEventListener('pointerdown', e => {
    if (e.button !== 0) return;
    delete el.dataset.longpressed;
    x0 = e.clientX; y0 = e.clientY;
    timer = setTimeout(() => { timer = null; el.dataset.longpressed = '1'; fn(); }, ms);
  });
  el.addEventListener('pointermove', e => { if (timer && Math.hypot(e.clientX - x0, e.clientY - y0) > 10) cancel(); });
  for (const ev of ['pointerup', 'pointercancel', 'pointerleave']) el.addEventListener(ev, cancel);
  el.addEventListener('contextmenu', e => e.preventDefault());
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
  $('#cut').checked = state.cut;
  $('#exposure').title = state.equalize
    ? 'Peak expected detections per pixel = exposure x 20'
    : 'Mean detections per pixel per frame (if the whole atom were in view)';
  for (const r of document.querySelectorAll('input[name=colormode]')) r.checked = r.value === state.color;
  for (const r of document.querySelectorAll('input[name=bgmode]')) r.checked = r.value === state.mode;
  for (const r of document.querySelectorAll('input[name=resmode]')) r.checked = r.value === state.resMode;
  $('#pause').textContent = state.paused ? 'Play' : 'Pause';
  $('#step').disabled = !state.paused;
  $('#debug-seed').checked = state.debugSeed;
  $('#show-overlay').checked = state.overlay;
  document.body.classList.toggle('mode-white', state.mode === 'white');
}

export function updateTable(app) {
  const atoms = atomsOf(app.state);
  const [A, B] = atoms;
  const cells = [...document.querySelectorAll('#ptable .pt-cell')];
  // Collapsed: the rows of the selected atoms in full, plus thin "peek"
  // slivers of their neighbouring rows (in table order: periods, then the
  // f-block row), hinting that there is a whole table to expand.
  const layout = app.tableCollapsed
    ? new Map(collapsedRows(cells.map(c => +c.dataset.row), atoms.map(tableRow)).map(r => [r.row, r]))
    : null;
  for (const c of cells) {
    const z = +c.dataset.z;
    c.classList.toggle('active', z === A.z);
    c.classList.toggle('second', !!B && z === B.z);
    if (B && (z === A.z || z === B.z)) c.dataset.badge = z === A.z ? 'A' : 'B';
    else delete c.dataset.badge;
    const row = +c.dataset.row;
    const info = layout?.get(row);
    const peek = !!info?.peek;
    c.classList.toggle('peek', peek);
    c.hidden = !!layout && !info;
    if (c.tagName === 'BUTTON') c.tabIndex = peek ? -1 : 0;
    if (peek) c.title = 'Show the whole periodic table';
    else if (c.dataset.title) c.title = c.dataset.title;
    c.style.gridColumn = c.dataset.col;
    c.style.gridRow = layout ? (info ? info.gridRow : 1) : row;
  }
  const t = $('#ptable-toggle');
  t.textContent = app.tableCollapsed ? '▾' : '▴';
  t.title = app.tableCollapsed ? 'Show the whole periodic table' : 'Collapse to the current row';
  t.setAttribute('aria-expanded', String(!app.tableCollapsed));
  $('#ptable').classList.toggle('collapsed', app.tableCollapsed);
  $('#ptable').classList.toggle('picking', !!app.pickingSecond);
  $('#pick-hint').hidden = !app.pickingSecond;
}

export function updateElement(app) {
  const { state } = app;
  const atoms = atomsOf(state);
  const [A, B] = atoms;
  updateTable(app);
  renderCard(app, A, B);
  renderLegend(app, atoms);
  renderNotes(atoms);
}

// ------------------------------------------------------------ element card
// One atom: the tile and its data, with prev/next and Compare (⇄) buttons.
// Two atoms: both tiles (A, B) with swap / remove, and a two-column table.
function renderCard(app, A, B) {
  const card = $('#element-card');
  card.classList.toggle('cmp', !!B);
  const steppers = [
    stepButton('‹', 'Previous element (←)', stepElement(A.z, -1), app),
    stepButton('›', 'Next element (→)', stepElement(A.z, +1), app),
  ];
  if (!B) {
    const cmp = h('button', { class: 'step' + (app.pickingSecond ? ' on' : ''),
      title: 'Compare with another element: then tap it in the periodic table ' +
        '(or shift-click / long-press an element)', 'aria-pressed': String(!!app.pickingSecond),
      onclick: () => app.togglePicking() }, '⇄');
    card.replaceChildren(
      tileCanvas(A),
      h('div', { class: 'card-info' },
        h('div', { class: 'card-name-row' },
          h('div', { class: 'card-name' }, A.name),
          h('span', { class: 'stepper' }, ...steppers, cmp)),
        h('div', { class: 'card-cat' }, A.category),
        infoRow('Configuration', A.configShort, A.config),
        infoRow('Period / Group', periodGroup(A)),
        infoRow('Electronegativity', en(A)),
        infoRow('1st ionization', `${A.ionizationEnergy.toFixed(3)} eV`),
        infoRow('Atomic radius', `${A.atomicRadius.toFixed(0)} pm`, RADIUS_TITLE)),
    );
    return;
  }
  const side = (el, badge) => h('div', { class: 'cmp-atom' },
    h('span', { class: 'badge' }, badge), tileCanvas(el, 0.6),
    h('div', { class: 'cmp-name' }, el.name), h('div', { class: 'card-cat' }, el.category));
  const row = (k, a, b, titleA, titleB) => h('tr', {},
    h('th', {}, k), h('td', { title: titleA || '' }, a), h('td', { title: titleB || '' }, b));
  card.replaceChildren(
    h('div', { class: 'cmp-head' },
      side(A, 'A'), side(B, 'B'),
      h('span', { class: 'stepper cmp-actions' }, ...steppers,
        h('button', { class: 'step', title: 'Swap A and B', onclick: () => app.swapAtoms() }, '⇆'),
        h('button', { class: 'step', title: 'Stop comparing (Esc)', onclick: () => app.clearSecond() }, '×'))),
    h('table', { class: 'cmp-table' },
      h('thead', {}, h('tr', {}, h('th', {}), h('th', {}, `A · ${A.symbol}`), h('th', {}, `B · ${B.symbol}`))),
      h('tbody', {},
        row('Configuration', A.configShort, B.configShort, A.config, B.config),
        row('Period / Group', periodGroup(A), periodGroup(B)),
        row('Electronegativity', en(A), en(B)),
        row('1st ionization', `${A.ionizationEnergy.toFixed(3)} eV`, `${B.ionizationEnergy.toFixed(3)} eV`),
        row('Atomic radius', `${A.atomicRadius.toFixed(0)} pm`, `${B.atomicRadius.toFixed(0)} pm`, RADIUS_TITLE, RADIUS_TITLE))),
  );
}

function tileCanvas(el, scale = 1) {
  const c = h('canvas', { class: 'card-symbol', role: 'img',
    'aria-label': `${el.name}, atomic number ${el.z}, mass ${el.mass}` });
  renderTileCanvas(c, el, scale);
  return c;
}
export const periodGroup = el => `${el.period} / ${el.group ?? 'f-block'}`;
export const RADIUS_TITLE = 'Rahm atomic radius: where the free atom\'s electron density falls to ' +
  '0.001 e/bohr³ (Rahm, Hoffmann & Ashcroft 2016)';
export const en = el => (el.electronegativity == null ? '—' : el.electronegativity.toFixed(2));

// ------------------------------------------------------------ orbital legend
// The union of the selected atoms' subshells, grouped by shell (n). Shell
// rows collapse; p/d/f rows expand to single orbitals. Visibility is shared:
// hiding an orbital hides it in every pane. With two atoms the numeric
// columns are per atom (electrons and r_rms for A and for B).
function renderLegend(app, atoms) {
  const { state } = app;
  const cmp = atoms.length > 1;
  const [A, B] = atoms;
  $('#legend thead').replaceChildren(cmp
    ? h('tr', {}, h('th', {}), h('th', {}), h('th', {}, 'Subshell'),
      h('th', { class: 'num' }, `e⁻ ${A.symbol}`), h('th', { class: 'num' }, `e⁻ ${B.symbol}`),
      h('th', { class: 'num', html: `r<sub>rms</sub> ${A.symbol}` }), h('th', { class: 'num', html: `r<sub>rms</sub> ${B.symbol}` }))
    : h('tr', {}, h('th', {}), h('th', {}), h('th', {}, 'Subshell'), h('th', {}, 'e⁻'),
      h('th', { class: 'num', html: 'Z<sub>eff</sub>' }),
      h('th', { class: 'num', title: 'RMS radius √⟨r²⟩', html: 'r<sub>rms</sub> pm' })));

  const occOf = (el, idx) => idx.reduce((a, i) => a + el.occupancy[i], 0);
  // Swatch color (depth-based chroma) from the first atom that has it.
  const owner = i => atoms.find(el => el.occupancy[i] > 0) || A;
  const checkbox = (title, idx, onChange) => {
    const occupied = idx.filter(i => occupiedAny(atoms, i));
    const vis = occupied.filter(i => state.visible[i]).length;
    const cb = h('input', { type: 'checkbox', title });
    cb.checked = occupied.length > 0 && vis === occupied.length;
    cb.indeterminate = vis > 0 && vis < occupied.length;
    cb.disabled = occupied.length === 0;
    cb.addEventListener('change', () => onChange(cb.checked));
    return cb;
  };

  const rows = [];
  const shown = SUBSHELLS.filter(x => x.period <= Math.max(2, maxPeriod(atoms)));
  for (const n of [...new Set(shown.map(x => x.n))].sort((a, b) => a - b)) {
    const shells = shown.filter(x => x.n === n);
    const idxAll = shells.flatMap(sh => orbitalsOf(sh.id));
    const totals = atoms.map(el => occOf(el, idxAll));
    const any = totals.some(t => t > 0);
    const openShell = app.expandedShells.has(n);
    const scb = checkbox(`Show / hide shell n = ${n}`, idxAll, on => app.setVisibleMany(idxAll, on));
    const label = h('td', { class: 'orb-label', colspan: '2' },
      h('button', { class: 'caret', title: openShell ? 'Collapse shell' : 'Show subshells',
        'aria-expanded': String(openShell), onclick: () => app.toggleShell(n) }, openShell ? '▾' : '▸'),
      `n = ${n}`, h('span', { class: 'shell-letter' }, ` ${SHELL_LETTER[n]}`));
    rows.push(h('tr', { class: 'shell' + (any ? '' : ' empty') },
      h('td', {}, scb), label,
      ...(cmp
        ? [h('td', { class: 'num' }, totals[0] ? String(totals[0]) : '—'),
          h('td', { class: 'num' }, totals[1] ? String(totals[1]) : '—'), h('td', {}), h('td', {})]
        : [h('td', { class: 'num' }, totals[0] ? String(totals[0]) : '—'), h('td', {}), h('td', {})])));
    if (!openShell) continue;

    for (const sh of shells) {
      const idx = orbitalsOf(sh.id);
      const occs = atoms.map(el => occOf(el, idx));
      const any = occs.some(o => o > 0);
      const cb = checkbox('Show / hide this subshell', idx, on => app.setVisibleMany(idx, on));
      const expandable = sh.l > 0 && any;
      const open = expandable && app.expanded.has(sh.id);
      const caret = expandable
        ? h('button', { class: 'caret', title: open ? 'Collapse' : 'Show individual orbitals',
            'aria-expanded': String(open), onclick: () => app.toggleExpanded(sh.id) }, open ? '▾' : '▸')
        : h('span', { class: 'caret' });
      const cap = CAPACITY[sh.l];
      const numbers = cmp
        ? [h('td', { class: 'num' }, occs[0] ? String(occs[0]) : '—'),
          h('td', { class: 'num' }, occs[1] ? String(occs[1]) : '—'),
          h('td', { class: 'num' }, rmsPm(A, sh)), h('td', { class: 'num' }, rmsPm(B, sh))]
        : [h('td', { class: 'num' }, occs[0] ? `${occs[0]}/${cap}` : '—'),
          h('td', { class: 'num' }, occs[0] ? A.zeffBySubshell[sh.id].toFixed(3) : ''),
          h('td', { class: 'num' }, rmsPm(A, sh))];
      rows.push(h('tr', { class: 'sub' + (any ? '' : ' empty') },
        h('td', {}, cb),
        h('td', { class: 'swatches' }, ...idx.map(i =>
          h('span', { class: 'swatch', style: `background:${orbitalCss(ORBITALS[i], owner(i))}` }))),
        h('td', { class: 'orb-label' }, caret, sh.id), ...numbers));

      if (open) {
        for (const i of idx) {
          const o = ORBITALS[i];
          const pips = el => (el.occupancy[i] === 0 ? '—' : (el.occupancy[i] === 2 ? '↑↓' : '↑'));
          const ocb = checkbox('Show / hide this orbital', [i], on => app.setVisible(i, on));
          rows.push(h('tr', { class: 'child' + (occupiedAny(atoms, i) ? '' : ' empty') },
            h('td', {}, ocb),
            h('td', { class: 'swatches' }, h('span', { class: 'swatch', style: `background:${orbitalCss(o, owner(i))}` })),
            h('td', { class: 'orb-label', html: o.label }),
            ...(cmp
              ? [h('td', { class: 'pips num' }, pips(A)), h('td', { class: 'pips num' }, pips(B)), h('td', {}), h('td', {})]
              : [h('td', { class: 'pips' }, pips(A)), h('td', {}), h('td', {})])));
        }
      }
    }
  }
  $('#legend tbody').replaceChildren(...rows);
}

// RMS radius of a subshell in pm ('' if the atom has no such electrons).
export function rmsPm(el, sh) {
  const z = el.zeffBySubshell[sh.id];
  return el.subshells[sh.id] ? (Math.sqrt(meanR2(sh.n, sh.l, z)) * BOHR_PM).toFixed(0) : '';
}

function renderNotes(atoms) {
  const notes = modelNotes(atoms);
  const note = $('#fill-note');
  note.hidden = notes.length === 0;
  note.textContent = notes.join(' ');
}

// Model caveats for these atoms, as sentences. Partly filled p/d/f
// subshells: which orbitals hold the electrons is one of several equivalent
// (degenerate) choices; say which convention is shown.
export function modelNotes(atoms) {
  const HUND_ORDER = { 1: 'p: x, y, z', 2: 'd: xy, yz, xz, x²−y², z²',
    3: 'f: z³, xz², yz², xyz, z(x²−y²), x(x²−3y²), y(3x²−y²)' };
  const notes = [];
  const partial = [];
  for (const el of atoms) {
    for (const sh of SUBSHELLS) {
      const c = el.subshells[sh.id] || 0;
      if (sh.l > 0 && c > 0 && c < CAPACITY[sh.l]) partial.push({ el, sh });
    }
  }
  if (partial.length) {
    const names = partial.map(({ el, sh }) => (atoms.length > 1 ? `${el.symbol} ${sh.id}` : sh.id));
    const order = [...new Set(partial.map(({ sh }) => HUND_ORDER[sh.l]))].join('; ');
    notes.push(`Partly filled ${names.join(', ')}: shown is one of several ` +
      `equivalent fillings (Hund's rule in the order ${order}).`);
  }
  if (atoms.some(el => el.z >= 55)) {
    notes.push('Non-relativistic hydrogen-like model: relativistic contraction / expansion of ' +
      'heavy-atom orbitals (e.g. 6s, 5d) is not included.');
  }
  return notes;
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

// Frame rate and render resolution, under the Resolution selector.
export function updatePerf(fps, renderK = 1, paused = false) {
  const res = renderK > 1 ? `1/${renderK} resolution` : 'full resolution';
  $('#perf').textContent = (paused ? 'paused' : `${fps.toFixed(0)} fps`) + ` · ${res}`;
}

// ------------------------------------------------------------ overlay
// Fill text with a halo in the background color so labels stay legible
// over the atomic-radius ring and the speckle.
function haloText(ctx, text, x, y, halo, width) {
  ctx.save();
  ctx.lineJoin = 'round';
  ctx.strokeStyle = halo;
  ctx.lineWidth = width;
  ctx.strokeText(text, x, y);
  ctx.restore();
  ctx.fillText(text, x, y);
}

// views: [{ element, geo }] with geo.rect in device px, GL origin
// (bottom-left), from renderer.renderViews. Each pane gets its own atomic-radius
// ring, element tile and scale bar; the axes marker (shared orientation) is
// drawn once, and two panes are separated by a thin divider.
// `scale`: device pixels per CSS pixel. Defaults to the canvas' on-screen
// ratio; give it for canvases that are not in the page (video frames).
export function drawOverlay(canvas, state, views, scale = null) {
  const ctx = canvas.getContext('2d');
  const w = canvas.width, hgt = canvas.height;
  ctx.clearRect(0, 0, w, hgt);
  if (!state.overlay) return;
  const dpr = scale ?? (w / canvas.clientWidth || 1);
  const cssWidth = w / dpr;
  const fg = state.mode === 'white' ? 'rgba(0,0,0,0.55)' : 'rgba(255,255,255,0.55)';
  const halo = state.mode === 'white' ? 'rgba(255,255,255,0.9)' : 'rgba(0,0,0,0.9)';
  const hw = 4 * dpr;

  // Canvas text/line state persists between frames, so isolate this draw
  // and start from known settings.
  ctx.save();
  // The scale is shared between panes, so only the last pane (bottom-right
  // of the canvas) draws a scale bar.
  // Each pane's ring, labels and tile are clipped to the pane, so a ring
  // larger than its pane (zoomed in, or sized for the bigger atom) never
  // spills into the neighbouring pane.
  views.forEach(({ element, geo }, i) => {
    const r = geo.rect;
    const box = { left: r.x, top: hgt - r.y - r.h, w: r.w, h: r.h };
    ctx.save();
    ctx.beginPath();
    ctx.rect(box.left, box.top, box.w, box.h);
    ctx.clip();
    drawPane(ctx, state, element, geo, box,
      { dpr, fg, halo, hw, scaleBar: i === views.length - 1, badge: views.length > 1 ? 'AB'[i] : null });
    ctx.restore();
  });

  // Divider between two panes.
  if (views.length > 1) {
    const r = views[1].geo.rect;
    ctx.strokeStyle = fg;
    ctx.lineWidth = 1 * dpr;
    ctx.beginPath();
    if (r.x > 0) { ctx.moveTo(r.x + 0.5, 0); ctx.lineTo(r.x + 0.5, hgt); }
    else { const yTop = hgt - r.h; ctx.moveTo(0, yTop + 0.5); ctx.lineTo(w, yTop + 0.5); }
    ctx.stroke();
  }

  // Axes gizmo (atom frame axes projected to the screen), bottom-left.
  const k = cssWidth < 520 ? 0.7 : 1;
  const m = qToMat3(state.orientation);   // atom -> view
  const gx = 50 * k * dpr, gy = hgt - 50 * k * dpr, gl = 34 * k * dpr;
  ctx.font = `${12 * dpr}px system-ui, sans-serif`;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'alphabetic';
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
  ctx.restore();
}

// One pane: atomic-radius ring, scale bar and element tile within `box`
// (device px, top-left origin).
function drawPane(ctx, state, el, geo, box, { dpr, fg, halo, hw, scaleBar = true, badge = null }) {
  const pxPerBohr = box.h / (2 * geo.halfExtent[1]);
  const cx = box.left + box.w / 2, cy = box.top + box.h / 2;
  // Compact overlay on narrow panes (phones in portrait, split views).
  const compact = box.w / dpr < 520;
  const k = compact ? 0.7 : 1;            // tile and margins
  const kt = compact ? 0.8 : 1;           // label text

  ctx.font = `${12 * dpr}px system-ui, sans-serif`;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'alphabetic';
  ctx.lineWidth = 1 * dpr;

  // Ring at the element's Rahm atomic radius (free-atom density 0.001
  // e/bohr^3); the view itself is still sized to fit the model's cloud.
  const r = (el.atomicRadius / BOHR_PM) * pxPerBohr;
  ctx.strokeStyle = fg;
  ctx.setLineDash([6 * dpr, 6 * dpr]);
  ctx.beginPath();
  ctx.arc(cx, cy, r, 0, 2 * Math.PI);
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = fg;
  haloText(ctx, `${Math.round(el.atomicRadius)} pm`, cx + r * 0.707 + 4 * dpr, cy - r * 0.707 - 4 * dpr, halo, hw);

  // Scale bar: 1 Å, or 1 a0 if that is too long; bottom-right of the pane.
  if (scaleBar) drawScaleBar(ctx, box, pxPerBohr, k, { dpr, fg, halo, hw });
  drawPaneLabel(ctx, state, el, box, k, kt, { dpr, fg, halo, hw, badge, compact });
}

function drawScaleBar(ctx, box, pxPerBohr, k, { dpr, fg, halo, hw }) {
  let len = 100 / BOHR_PM, label = '1 Å';
  if (len * pxPerBohr > box.w * 0.3) { len = 1; label = '1 a₀'; }
  const bx = box.left + box.w - 30 * k * dpr - len * pxPerBohr, by = box.top + box.h - 30 * k * dpr;
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
}

// Element tile (same rendering as the side-panel card) with name, category
// and configuration beside it, top-left of the pane.
function drawPaneLabel(ctx, state, el, box, k, kt, { dpr, fg, halo, hw, badge, compact }) {
  const tx = box.left + 16 * k * dpr, ty = box.top + 16 * k * dpr;
  drawElementTile(ctx, el, tx, ty, k * dpr);
  if (badge) {                              // A / B, matching the periodic table
    const r = 9 * k * dpr, bx = tx + TILE_W * k * dpr, by = ty;
    ctx.beginPath();
    ctx.arc(bx, by, r, 0, 2 * Math.PI);
    ctx.fillStyle = '#fff';
    ctx.fill();
    ctx.fillStyle = '#111';
    ctx.font = `700 ${11 * k * dpr}px system-ui, sans-serif`;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(badge, bx, by + 0.5 * dpr);
  }
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

  // Rahm atomic radius of this pane's atom: top-right of the pane, or (on
  // compact panes, where it would collide with the name) a fourth line
  // under the configuration.
  const pm = Math.round(el.atomicRadius), a0 = (el.atomicRadius / BOHR_PM).toFixed(2);
  ctx.fillStyle = state.mode === 'white' ? 'rgba(0,0,0,0.6)' : 'rgba(255,255,255,0.6)';
  if (compact) {
    ctx.font = `${10 * dpr}px ui-monospace, monospace`;
    haloText(ctx, `Rahm radius ${pm} pm`, lx, ty + 64 * kt * dpr, halo, hw);
  } else {
    ctx.font = `${12 * dpr}px ui-monospace, monospace`;
    ctx.textAlign = 'right';
    haloText(ctx, `Rahm atomic radius = ${pm} pm (${a0} a₀)`, box.left + box.w - 12 * dpr, box.top + 10 * dpr, halo, hw);
  }
}
