// Radial electron density inset: occ * P(r) per subshell plus the total,
// in electrons per Bohr radius, always physical (never equalized).

import { ELEMENTS } from './elements.js';
import { radialDensity } from './math.js';
import { SUBSHELLS, orbitalsOf, SHELL_COLORS } from './orbitals.js';

const BOHR_PM = 52.9177;

// Color = shell n, validated on the dark panel (#151922) in shell order,
// i.e. neighbouring shells, which overlap in r (n5/n6 sit in the CVD
// 6-8 band, so curves are also direct-labelled). Two modes:
//  - subshells: one series per subshell, line style = l
//    (s solid, p dashed, d dotted, f dash-dot);
//  - shells: one solid series per shell n (for heavy atoms).
// SHELL_COLORS live in orbitals.js: the 'shell' color scheme uses them too.
const L_DASH = { 0: [], 1: [6, 3], 2: [1.5, 2.5], 3: [8, 3, 1.5, 3] };
const part = sh => ({ n: sh.n, l: sh.l, orbitals: orbitalsOf(sh.id) });
export const SERIES = SUBSHELLS.map(sh => ({
  key: sh.id, label: sh.id, color: SHELL_COLORS[sh.n], dash: L_DASH[sh.l], parts: [part(sh)],
}));
export const SHELL_SERIES = [...new Set(SUBSHELLS.map(sh => sh.n))].sort((a, b) => a - b).map(n => ({
  key: `n${n}`, label: `n=${n}`, color: SHELL_COLORS[n], dash: [],
  parts: SUBSHELLS.filter(sh => sh.n === n).map(part),
}));

// 'shell' | 'subshell' for this element: explicit choice, or shells when
// more than 8 subshells are occupied.
export function radialMode(state, el) {
  if (state.second && (state.radialMode === 'compare' || state.radialMode === 'auto')) return 'compare';
  if (state.radialMode === 'shell' || state.radialMode === 'subshell') return state.radialMode;
  return Object.keys(el.subshells).length > 8 ? 'shell' : 'subshell';
}

// A vs B: each atom's total radial density, A solid and B dashed, in two
// colors validated as a pair on the dark panel (also the n=3 / n=2 shell
// colors).
const COMPARE_COLORS = ['#0f91bd', '#b88b05'];
const COMPARE_DASH = [[], [6, 3]];

// Small inline line sample for legends/tooltips (same color + dash as the plot).
function lineSample(color, dash = []) {
  const da = dash.length ? ` stroke-dasharray="${dash.join(' ')}"` : '';
  return `<svg width="18" height="6" aria-hidden="true"><line x1="0" y1="3" x2="18" y2="3" ` +
    `stroke="${color}" stroke-width="2"${da}/></svg>`;
}

// Neutral inks (axes, grid, labels, the total curve) per surface: the dark
// panel, or white paper on the print page. The series colors are the same
// on both; they were also validated on #ffffff (same results: contrast,
// lightness and normal-vision PASS, CVD n5/n6 in the 6-8 band, which the
// direct labels cover).
const INKS = {
  dark: { muted: '#8b93a7', faint: '#2a3142', text: '#e6e9f0', total: '#c5cad6', surface: '#151922' },
  light: { muted: '#5c6375', faint: '#dfe3ea', text: '#1d2230', total: '#6b7280', surface: '#ffffff' },
};
const M = { l: 34, r: 8, t: 10, b: 22 };   // margins, CSS px
const SAMPLES = 240;

// Electrons per a0 for one series (a subshell or a whole shell) at r (Bohr).
export function seriesDensity(el, s, r) {
  let sum = 0;
  for (const p of s.parts) {
    const occ = p.orbitals.reduce((a, i) => a + el.occupancy[i], 0);
    if (occ) sum += occ * radialDensity(p.n, p.l, el.zeff[p.orbitals[0]], r);
  }
  return sum;
}

const occupiedSeries = (list, el) =>
  list.filter(s => s.parts.some(p => p.orbitals.some(i => el.occupancy[i] > 0)));

function totalDensity(el, r) {
  return SERIES.reduce((a, s) => a + seriesDensity(el, s, r), 0);
}

function peakOver(el, rMax) {
  let m = 0;
  for (let i = 1; i <= 600; i++) m = Math.max(m, totalDensity(el, (i / 600) * rMax));
  return m;
}

function niceStep(range, target) {
  const raw = range / target;
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  const f = raw / p;
  return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p;
}

// theme: 'dark' | 'light'; dpr: canvas pixels per CSS pixel (default: the
// screen's, capped at 2; the print page asks for more). Without a tooltip
// element there is no hover.
export function createRadialPlot(canvas, tooltip, legendEl, { theme = 'dark', dpr: fixedDpr = null } = {}) {
  const INK = INKS[theme];
  let model = null;   // last drawn data, for hover

  let legendFor = null;
  function updateLegend(active, showTotal) {
    const key = active.map(s => s.key + s.label).join() + showTotal;
    if (key === legendFor) return;
    legendFor = key;
    const items = showTotal ? [...active, { label: 'total', color: INK.total }] : active;
    legendEl.replaceChildren(...items.map(s => {
      const item = document.createElement('span');
      item.className = 'rl-item';
      item.innerHTML = `${lineSample(s.color, s.dash)}${s.label}`;
      return item;
    }));
  }

  function draw(state, geo, hoverX = null) {
    const el = state.element;
    const dpr = fixedDpr ?? Math.min(window.devicePixelRatio || 1, 2);
    const cw = canvas.clientWidth, ch = canvas.clientHeight;
    if (canvas.width !== Math.round(cw * dpr)) canvas.width = Math.round(cw * dpr);
    if (canvas.height !== Math.round(ch * dpr)) canvas.height = Math.round(ch * dpr);
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cw, ch);

    const rMax = geo.halfExtent[1];                       // Bohr, matches view
    // Fixed scale: share the y range across the element's period so e.g. He
    // vs H heights compare directly. Auto: fit this element.
    const mode = radialMode(state, el);
    const atoms = mode === 'compare' ? [el, state.second] : [el];
    const yMax = 1.08 * (state.fixedScale
      ? Math.max(...ELEMENTS.filter(e => atoms.some(a => a.period === e.period) && e.supported)
        .map(e => peakOver(e, rMax)))
      : Math.max(...atoms.map(a => peakOver(a, rMax))));
    const pw = cw - M.l - M.r, ph = ch - M.t - M.b;
    const X = r => M.l + (r / rMax) * pw;
    const Y = v => M.t + ph - (v / yMax) * ph;

    // grid + axes (recessive)
    ctx.font = '10px system-ui, sans-serif';
    ctx.lineWidth = 1;
    const yStep = niceStep(yMax, 3);
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (let v = 0; v <= yMax + 1e-9; v += yStep) {
      ctx.strokeStyle = INK.faint;
      ctx.beginPath(); ctx.moveTo(M.l, Y(v) + 0.5); ctx.lineTo(M.l + pw, Y(v) + 0.5); ctx.stroke();
      ctx.fillStyle = INK.muted;
      ctx.fillText(+v.toPrecision(3) + '', M.l - 4, Y(v));
    }
    const xStepPm = niceStep(rMax * BOHR_PM, 4);
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    for (let pm = 0; pm <= rMax * BOHR_PM + 1e-9; pm += xStepPm) {
      ctx.fillStyle = INK.muted;
      ctx.fillText(String(Math.round(pm)), X(pm / BOHR_PM), M.t + ph + 4);
    }

    // Atomic-radius marker(s): Rahm et al. radius of each plotted atom.
    atoms.forEach((a, i) => {
      const rr = a.atomicRadius / BOHR_PM;
      if (!(rr <= rMax)) return;
      const hx = X(rr);
      ctx.strokeStyle = INK.muted;
      ctx.setLineDash([3, 3]);
      ctx.beginPath(); ctx.moveTo(hx + 0.5, M.t); ctx.lineTo(hx + 0.5, M.t + ph); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = INK.muted;
      ctx.textBaseline = 'top';
      const label = (atoms.length > 1 ? `${a.symbol} ` : '') + `${Math.round(a.atomicRadius)} pm`;
      // Right of the line, or left of it when it would run off the plot.
      const fits = hx + 3 + ctx.measureText(label).width <= M.l + pw;
      ctx.textAlign = fits ? 'left' : 'right';
      ctx.fillText(label, fits ? hx + 3 : hx - 3, M.t + 12 * i);
    });

    // curves
    const rs = Array.from({ length: SAMPLES + 1 }, (_, i) => (i / SAMPLES) * rMax);
    const curve = (fn, color, width, dash = []) => {
      ctx.strokeStyle = color;
      ctx.lineWidth = width;
      ctx.lineJoin = 'round';
      ctx.setLineDash(dash);
      ctx.beginPath();
      rs.forEach((r, i) => { const y = Y(Math.min(fn(r), yMax * 1.2)); i ? ctx.lineTo(X(r), y) : ctx.moveTo(X(r), y); });
      ctx.stroke();
    };
    // Plotted series, each with its own density function r -> e- per a0.
    const active = mode === 'compare'
      ? atoms.map((a, i) => ({ key: `atom${i}`, label: a.symbol, color: COMPARE_COLORS[i],
          dash: COMPARE_DASH[i], f: r => totalDensity(a, r) }))
      : occupiedSeries(mode === 'shell' ? SHELL_SERIES : SERIES, el)
          .map(sr => ({ ...sr, f: r => seriesDensity(el, sr, r) }));
    const showTotal = mode !== 'compare' && active.length > 1;
    updateLegend(active, showTotal);
    ctx.save();
    ctx.beginPath(); ctx.rect(M.l, M.t - 2, pw, ph + 4); ctx.clip();
    if (showTotal) curve(r => totalDensity(el, r), INK.total, 1.25);
    for (const s of active) curve(s.f, s.color, 2, s.dash);
    ctx.setLineDash([]);
    ctx.restore();

    // selective direct labels at each series peak (skip if they would collide)
    ctx.font = '11px system-ui, sans-serif';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'bottom';
    const placed = [];
    for (const s of active) {
      let best = 0, br = 0;
      for (const r of rs) { const v = s.f(r); if (v > best) { best = v; br = r; } }
      const lx = X(br) + 4, ly = Math.max(M.t + 12, Y(best) - 2);
      if (placed.some(p => Math.abs(p[0] - lx) < 22 && Math.abs(p[1] - ly) < 12)) continue;
      placed.push([lx, ly]);
      ctx.fillStyle = INK.text;
      ctx.fillText(s.label, lx, ly);
    }

    // hover crosshair
    model = { state, geo, rMax, X, Y, pw, ph, active, showTotal };
    if (hoverX != null && hoverX >= M.l && hoverX <= M.l + pw) {
      const r = ((hoverX - M.l) / pw) * rMax;
      ctx.strokeStyle = INK.muted;
      ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(hoverX + 0.5, M.t); ctx.lineTo(hoverX + 0.5, M.t + ph); ctx.stroke();
      const rows = [];
      for (const s of active) {
        const v = s.f(r);
        ctx.fillStyle = s.color;
        ctx.strokeStyle = INK.surface;
        ctx.lineWidth = 2;
        ctx.beginPath(); ctx.arc(X(r), Y(Math.min(v, yMax)), 4, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
        rows.push(`<div>${lineSample(s.color, s.dash)}${s.label}<b>${v.toFixed(3)}</b></div>`);
      }
      if (showTotal) rows.push(`<div>${lineSample(INK.total)}total<b>${totalDensity(el, r).toFixed(3)}</b></div>`);
      tooltip.innerHTML = `<div class="tt-head">r = ${(r * BOHR_PM).toFixed(0)} pm · e⁻/a₀</div>` + rows.join('');
      tooltip.hidden = false;
      const left = hoverX + 12 + tooltip.offsetWidth > cw ? hoverX - 12 - tooltip.offsetWidth : hoverX + 12;
      tooltip.style.left = `${left}px`;
      tooltip.style.top = `${M.t}px`;
    } else if (tooltip) {
      tooltip.hidden = true;
    }
  }

  if (!tooltip) return { draw };
  // Hit target is the whole plot area, wider than the 2px lines.
  canvas.addEventListener('pointermove', e => {
    if (!model) return;
    const rect = canvas.getBoundingClientRect();
    draw(model.state, model.geo, e.clientX - rect.left);
  });
  canvas.addEventListener('pointerleave', () => { if (model) draw(model.state, model.geo); });

  return { draw };
}
