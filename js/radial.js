// Radial electron density inset: occ * P(r) per subshell plus the total,
// in electrons per Bohr radius, always physical (never equalized).

import { ELEMENTS } from './elements.js';
import { radialDensity } from './math.js';
import { SUBSHELLS, orbitalsOf } from './orbitals.js';

const BOHR_PM = 52.9177;

// Same hues as the render, stepped to chart-safe lightness for the dark
// panel (validated: lightness band, CVD separation, contrast vs #151922).
// One series per subshell. Colours for shells 3-4 are provisional until the
// period 3/4 colour pass re-validates the palette.
const CHART_COLORS = {
  '1s': '#ca1616', '2s': '#b88b05', '2p': '#0f91bd',
  '3s': '#c2410c', '3p': '#0e7490', '3d': '#a21caf', '4s': '#be123c', '4p': '#155e75',
};
export const SERIES = SUBSHELLS.map(sh => ({
  key: sh.id, label: sh.id, color: CHART_COLORS[sh.id], n: sh.n, l: sh.l,
  orbitals: orbitalsOf(sh.id),
}));
const INK = { muted: '#8b93a7', faint: '#2a3142', text: '#e6e9f0', total: '#c5cad6' };
const M = { l: 34, r: 8, t: 10, b: 22 };   // margins, CSS px
const SAMPLES = 240;

// Electrons per a0 for one subshell series of `el` at radius r (Bohr).
export function seriesDensity(el, s, r) {
  const occ = s.orbitals.reduce((a, i) => a + el.occupancy[i], 0);
  if (!occ) return 0;
  return occ * radialDensity(s.n, s.l, el.zeff[s.orbitals[0]], r);
}

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

export function createRadialPlot(canvas, tooltip, legendEl) {
  let model = null;   // last drawn data, for hover

  let legendFor = null;
  function updateLegend(active) {
    const key = active.map(s => s.key).join();
    if (key === legendFor) return;
    legendFor = key;
    const items = active.length > 1 ? [...active, { label: 'total', color: INK.total }] : active;
    legendEl.replaceChildren(...items.map(s => {
      const item = document.createElement('span');
      item.className = 'rl-item';
      item.innerHTML = `<i style="background:${s.color}"></i>${s.label}`;
      return item;
    }));
  }

  function draw(state, geo, hoverX = null) {
    const el = state.element;
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const cw = canvas.clientWidth, ch = canvas.clientHeight;
    if (canvas.width !== Math.round(cw * dpr)) canvas.width = Math.round(cw * dpr);
    if (canvas.height !== Math.round(ch * dpr)) canvas.height = Math.round(ch * dpr);
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cw, ch);

    const rMax = geo.halfExtent[1];                       // Bohr, matches view
    // Fixed scale: share the y range across the element's period so e.g. He
    // vs H heights compare directly. Auto: fit this element.
    const yMax = 1.08 * (state.fixedScale
      ? Math.max(...ELEMENTS.filter(e => e.period === el.period).map(e => peakOver(e, rMax)))
      : peakOver(el, rMax));
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

    // horizon marker
    const hx = X(geo.horizon);
    if (geo.horizon <= rMax) {
      ctx.strokeStyle = INK.muted;
      ctx.setLineDash([3, 3]);
      ctx.beginPath(); ctx.moveTo(hx + 0.5, M.t); ctx.lineTo(hx + 0.5, M.t + ph); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = INK.muted;
      ctx.textAlign = 'left';
      ctx.textBaseline = 'top';
      ctx.fillText('2σ', hx + 3, M.t);
    }

    // curves
    const rs = Array.from({ length: SAMPLES + 1 }, (_, i) => (i / SAMPLES) * rMax);
    const curve = (fn, color, width) => {
      ctx.strokeStyle = color;
      ctx.lineWidth = width;
      ctx.lineJoin = 'round';
      ctx.beginPath();
      rs.forEach((r, i) => { const y = Y(Math.min(fn(r), yMax * 1.2)); i ? ctx.lineTo(X(r), y) : ctx.moveTo(X(r), y); });
      ctx.stroke();
    };
    const active = SERIES.filter(s => s.orbitals.some(i => el.occupancy[i] > 0));
    updateLegend(active);
    ctx.save();
    ctx.beginPath(); ctx.rect(M.l, M.t - 2, pw, ph + 4); ctx.clip();
    if (active.length > 1) curve(r => totalDensity(el, r), INK.total, 1.25);
    for (const s of active) curve(r => seriesDensity(el, s, r), s.color, 2);
    ctx.restore();

    // selective direct labels at each series peak (skip if they would collide)
    ctx.font = '11px system-ui, sans-serif';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'bottom';
    const placed = [];
    for (const s of active) {
      let best = 0, br = 0;
      for (const r of rs) { const v = seriesDensity(el, s, r); if (v > best) { best = v; br = r; } }
      const lx = X(br) + 4, ly = Math.max(M.t + 12, Y(best) - 2);
      if (placed.some(p => Math.abs(p[0] - lx) < 22 && Math.abs(p[1] - ly) < 12)) continue;
      placed.push([lx, ly]);
      ctx.fillStyle = INK.text;
      ctx.fillText(s.label, lx, ly);
    }

    // hover crosshair
    model = { state, geo, rMax, X, Y, pw, ph, active };
    if (hoverX != null && hoverX >= M.l && hoverX <= M.l + pw) {
      const r = ((hoverX - M.l) / pw) * rMax;
      ctx.strokeStyle = INK.muted;
      ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(hoverX + 0.5, M.t); ctx.lineTo(hoverX + 0.5, M.t + ph); ctx.stroke();
      const rows = [];
      for (const s of active) {
        const v = seriesDensity(el, s, r);
        ctx.fillStyle = s.color;
        ctx.strokeStyle = '#151922';
        ctx.lineWidth = 2;
        ctx.beginPath(); ctx.arc(X(r), Y(Math.min(v, yMax)), 4, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
        rows.push(`<div><i style="background:${s.color}"></i>${s.label}<b>${v.toFixed(3)}</b></div>`);
      }
      if (active.length > 1) rows.push(`<div><i style="background:${INK.total}"></i>total<b>${totalDensity(el, r).toFixed(3)}</b></div>`);
      tooltip.innerHTML = `<div class="tt-head">r = ${(r * BOHR_PM).toFixed(0)} pm · e⁻/a₀</div>` + rows.join('');
      tooltip.hidden = false;
      const left = hoverX + 12 + tooltip.offsetWidth > cw ? hoverX - 12 - tooltip.offsetWidth : hoverX + 12;
      tooltip.style.left = `${left}px`;
      tooltip.style.top = `${M.t}px`;
    } else {
      tooltip.hidden = true;
    }
  }

  // Hit target is the whole plot area, wider than the 2px lines.
  canvas.addEventListener('pointermove', e => {
    if (!model) return;
    const rect = canvas.getBoundingClientRect();
    draw(model.state, model.geo, e.clientX - rect.left);
  });
  canvas.addEventListener('pointerleave', () => { if (model) draw(model.state, model.geo); });

  return { draw };
}
