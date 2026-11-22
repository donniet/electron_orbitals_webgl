// Orbital definitions (data-driven) and hydrogen-like radial functions.
// Mirrors tools/orbitals.py; see that file for the math.
//
// Scaling: each radial function is tabulated for reference charge Z0 = n,
// in u = Z r / n. For effective charge Z, kappa = Z / n and
//   col_Z(p, d) = kappa^2 col_ref(kappa p, d).

// Subshells in written (configuration) order. `row` is the moment-table row
// (must match data/orbital_tables.json; checked in tests).
export const SUBSHELLS = [
  { id: '1s', n: 1, l: 0, period: 1, row: 0 },
  { id: '2s', n: 2, l: 0, period: 2, row: 1 },
  { id: '2p', n: 2, l: 1, period: 2, row: 4 },
  { id: '3s', n: 3, l: 0, period: 3, row: 2 },
  { id: '3p', n: 3, l: 1, period: 3, row: 5 },
  { id: '3d', n: 3, l: 2, period: 4, row: 7 },
  { id: '4s', n: 4, l: 0, period: 4, row: 3 },
  { id: '4p', n: 4, l: 1, period: 4, row: 6 },
  { id: '4d', n: 4, l: 2, period: 5, row: 12 },
  { id: '4f', n: 4, l: 3, period: 6, row: 14 },
  { id: '5s', n: 5, l: 0, period: 5, row: 8 },
  { id: '5p', n: 5, l: 1, period: 5, row: 10 },
  { id: '5d', n: 5, l: 2, period: 6, row: 13 },
  { id: '6s', n: 6, l: 0, period: 6, row: 9 },
  { id: '6p', n: 6, l: 1, period: 6, row: 11 },
];
export const SUBSHELL = Object.fromEntries(SUBSHELLS.map(s => [s.id, s]));

// Aufbau filling order.
export const AUFBAU = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s', '4d', '5p', '6s', '4f', '5d', '6p'];

const P_NORM = Math.sqrt(3 / (4 * Math.PI));
const S15 = Math.sqrt(15 / (4 * Math.PI)) / 2;
const Z2 = Math.sqrt(5 / (16 * Math.PI));

// p: e = sqrt(3/4pi) * axis.  d: symmetric Q with (x.Qx)/r^2 = real Y_2m.
const P_SHAPES = [
  { m: 'x', label: 'x', e: [P_NORM, 0, 0] },
  { m: 'y', label: 'y', e: [0, P_NORM, 0] },
  { m: 'z', label: 'z', e: [0, 0, P_NORM] },
];
const D_SHAPES = [
  { m: 'xy',    label: 'xy',       Q: [[0, S15, 0], [S15, 0, 0], [0, 0, 0]] },
  { m: 'yz',    label: 'yz',       Q: [[0, 0, 0], [0, 0, S15], [0, S15, 0]] },
  { m: 'xz',    label: 'xz',       Q: [[0, 0, S15], [0, 0, 0], [S15, 0, 0]] },
  { m: 'x2-y2', label: 'x²−y²',    Q: [[S15, 0, 0], [0, -S15, 0], [0, 0, 0]] },
  { m: 'z2',    label: 'z²',       Q: [[-Z2, 0, 0], [0, -Z2, 0], [0, 0, 2 * Z2]] },
];

// f: real l = 3 harmonics as symmetric cubic tensors C (C(x,x,x)/r^3 = Y_3m),
// built from monomials {[ex, ey, ez]: coefficient}. Hund fill order follows
// this list. Mirrors tools/orbitals.py.
function cubicTensor(monomials) {
  const C = [0, 1, 2].map(() => [0, 1, 2].map(() => [0, 0, 0]));
  for (const [exps, coef] of monomials) {
    const idx = [];
    exps.forEach((e, ax) => { for (let k = 0; k < e; k++) idx.push(ax); });
    const perms = new Set();
    for (const [a, b, c] of [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]) {
      perms.add(`${idx[a]}${idx[b]}${idx[c]}`);
    }
    for (const key of perms) C[+key[0]][+key[1]][+key[2]] += coef / perms.size;
  }
  return C;
}
const FA = Math.sqrt(7 / (16 * Math.PI)), FB = Math.sqrt(21 / (32 * Math.PI));
const FC = Math.sqrt(105 / (4 * Math.PI)), FD = Math.sqrt(105 / (16 * Math.PI));
const FE = Math.sqrt(35 / (32 * Math.PI));
const F_SHAPES = [
  { m: 'z3',        label: 'z³',        C: cubicTensor([[[0, 0, 3], 2 * FA], [[2, 0, 1], -3 * FA], [[0, 2, 1], -3 * FA]]) },
  { m: 'xz2',       label: 'xz²',       C: cubicTensor([[[1, 0, 2], 4 * FB], [[3, 0, 0], -FB], [[1, 2, 0], -FB]]) },
  { m: 'yz2',       label: 'yz²',       C: cubicTensor([[[0, 1, 2], 4 * FB], [[2, 1, 0], -FB], [[0, 3, 0], -FB]]) },
  { m: 'xyz',       label: 'xyz',       C: cubicTensor([[[1, 1, 1], FC]]) },
  { m: 'z(x2-y2)',  label: 'z(x²−y²)',  C: cubicTensor([[[2, 0, 1], FD], [[0, 2, 1], -FD]]) },
  { m: 'x(x2-3y2)', label: 'x(x²−3y²)', C: cubicTensor([[[3, 0, 0], FE], [[1, 2, 0], -3 * FE]]) },
  { m: 'y(3x2-y2)', label: 'y(3x²−y²)', C: cubicTensor([[[2, 1, 0], 3 * FE], [[0, 3, 0], -FE]]) },
];

// Render colors. Hue encodes the angular shape (chroma, the weight of the
// hue, depends on the element: see orbitalChroma below).
//  s: warm hues by shell (concentric shells stay distinguishable)
//  p: one hue per axis in every shell
//  d, f: hues of their own
const S_HUE = { 1: 0, 2: 45, 3: 340, 4: 20, 5: 45, 6: 340 };
const P_HUE = { x: 195, y: 275, z: 125 };
const D_HUE = { xy: 300, yz: 235, xz: 160, 'x2-y2': 85, z2: 55 };
const F_HUE = { z3: 15, xz2: 65, yz2: 105, xyz: 145, 'z(x2-y2)': 215, 'x(x2-3y2)': 255, 'y(3x2-y2)': 320 };

function makeOrbitals(subshellId) {
  const s = SUBSHELL[subshellId];
  const base = { subshell: s.id, n: s.n, l: s.l, row: s.row };
  if (s.l === 0) return [{ ...base, id: s.id, m: '', label: s.id, hue: S_HUE[s.n] }];
  if (s.l === 1) {
    return P_SHAPES.map((p, a) => ({
      ...base, id: s.id + p.m, m: p.m, axis: a, e: p.e,
      label: `${s.id}<sub>${p.label}</sub>`, hue: P_HUE[p.m],
    }));
  }
  if (s.l === 2) {
    return D_SHAPES.map(d => ({
      ...base, id: s.id + d.m, m: d.m, Q: d.Q,
      label: `${s.id}<sub>${d.label}</sub>`, hue: D_HUE[d.m],
    }));
  }
  return F_SHAPES.map(f => ({
    ...base, id: s.id + f.m, m: f.m, C: f.C,
    label: `${s.id}<sub>${f.label}</sub>`, hue: F_HUE[f.m],
  }));
}

// The 12 angular coefficients the shader needs for an orbital:
//  p: e (3); d: unique Q entries q00 q01 q02 q11 q12 q22 (6);
//  f: unique C entries c000 c001 c002 c011 c012 c022 c111 c112 c122 c222 (10).
export function angularCoefs(o) {
  const out = new Array(12).fill(0);
  if (o.l === 1) o.e.forEach((v, i) => { out[i] = v; });
  if (o.l === 2) {
    const Q = o.Q;
    [Q[0][0], Q[0][1], Q[0][2], Q[1][1], Q[1][2], Q[2][2]].forEach((v, i) => { out[i] = v; });
  }
  if (o.l === 3) {
    const C = o.C;
    [C[0][0][0], C[0][0][1], C[0][0][2], C[0][1][1], C[0][1][2], C[0][2][2],
      C[1][1][1], C[1][1][2], C[1][2][2], C[2][2][2]].forEach((v, i) => { out[i] = v; });
  }
  return out;
}

// DT: Why is this different from Aufbau filling order (AUFBAU constant above)?
//     Do we need both?
// Canonical (shader) order. Periods 1-4 keep their original order and new
// subshells are appended, so frames for H-Kr consume random numbers in the
// same order as before.
export const ORBITALS = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p',
  '5s', '4d', '5p', '6s', '4f', '5d', '6p'].flatMap(makeOrbitals);
export const NUM_ORBITALS = ORBITALS.length;   // 43

export function orbitalIndex(id) {
  return ORBITALS.findIndex(o => o.id === id);
}

export function orbitalsOf(subshellId) {
  return ORBITALS.map((o, i) => i).filter(i => ORBITALS[i].subshell === subshellId);
}

// ------------------------------------------------------------ radial functions
function laguerre(k, alpha) {
  // ascending coefficients of the generalised Laguerre L_k^alpha(x)
  let prev = [1], cur = [1 + alpha, -1];
  if (k === 0) return prev;
  for (let j = 1; j < k; j++) {
    const nxt = new Array(cur.length + 1).fill(0);
    cur.forEach((c, i) => { nxt[i] += (2 * j + 1 + alpha) * c; nxt[i + 1] -= c; });
    prev.forEach((c, i) => { nxt[i] -= (j + alpha) * c; });
    prev = cur;
    cur = nxt.map(c => c / (j + 1));
  }
  return cur;
}

function factorial(n) { let f = 1; for (let i = 2; i <= n; i++) f *= i; return f; }

const RADIAL_CACHE = new Map();
function radialRef(n, l) {
  const key = n * 10 + l;
  if (!RADIAL_CACHE.has(key)) {
    const norm = Math.sqrt(8 * factorial(n - l - 1) / (2 * n * factorial(n + l)));
    const poly = laguerre(n - l - 1, 2 * l + 1).map((c, k) => c * 2 ** k * 2 ** l);
    RADIAL_CACHE.set(key, { norm, poly });
  }
  return RADIAL_CACHE.get(key);
}

// R_ref(u) = norm exp(-u) P(u) u^l: { norm, poly } with poly the ascending
// coefficients of P (a scaled Laguerre polynomial, degree n - l - 1).
export function radialPolyRef(n, l) { return radialRef(n, l); }

// R_ref(u) for reference charge Z0 = n.
export function radialRefR(n, l, u) {
  const { norm, poly } = radialRef(n, l);
  let v = 0;
  for (let i = poly.length - 1; i >= 0; i--) v = v * u + poly[i];
  return norm * Math.exp(-u) * v * u ** l;
}

// Normalised radial probability density P(r) = r^2 R_nl(r)^2 (per a0) for
// effective charge z. Integrates to 1 over r.
export function radialDensity(n, l, z, r) {
  const k = z / n;
  const u = k * r;
  const R = radialRefR(n, l, u);
  return k * u * u * R * R;
}

// Render chroma (weight of the orbital's hue) by depth relative to the
// element: valence orbitals are vivid, the shell just inside paler, deeper
// shells palest. s orbitals stay vivid (concentric s shells are told apart
// by their hues), so H-Ne look as before.
export function orbitalChroma(orbital, element) {
  return orbital.l === 0 ? 1 : depthChroma(orbital, element);
}

function depthChroma(orbital, element) {
  if (isValence(orbital, element)) return 1;
  return element.period - orbital.n <= 1 ? 0.6 : 0.35;
}

// Color schemes: what an orbital's hue encodes. The renderer composites
// (hue, chroma) per orbital as always, so every scheme works in blend and
// nearest color and keeps the white theme the exact dual of the dark one.
//  shape     hue = angular shape (S/P/D/F_HUE above), chroma by depth
//  gray      no hue: brightness only
//  subshell  hue = l (s, p, d, f), chroma by depth (s included)
//  shell     hue = n, the radial chart's shell colors, full chroma
//  valence   core (paler blue) vs valence (amber), see isValence
// Per detection, nearest color only (the renderer colors the nearest
// detection itself; per orbital these have no color):
//  phase     sign of psi at the detection: + red, - blue
//  momentum  local RMS momentum at the detection: slow blue -> fast amber
export const COLOR_SCHEMES = ['shape', 'gray', 'subshell', 'shell', 'valence', 'phase', 'momentum'];
export const SCHEME_LABELS = {
  shape: 'Shape', gray: 'Gray', subshell: 'Subshell (s p d f)', shell: 'Shell (n)',
  valence: 'Core / valence', phase: 'Phase (sign of ψ)', momentum: 'Momentum (local)',
};
export const NEAREST_ONLY_SCHEMES = new Set(['phase', 'momentum']);

// Two-color palettes [hue (degrees), chroma]. As drawn by the renderer (full
// value; white theme the dual) each pair passes the dataviz palette
// validator's CVD check in both themes (worst dE 26 core/valence, 31 phase,
// 31 momentum). Phase follows the textbook red +, blue -. The momentum
// hues are opposite, so mixing them passes through neutral gray halfway.
export const PALETTE_COLORS = {
  valence: { core: [205, 0.6], valence: [35, 1] },
  phase: { plus: [0, 1], minus: [210, 1] },
  momentum: { slow: [220, 1], fast: [40, 1] },
};

// URL value -> scheme ('shape' for anything unknown).
export function schemeParam(v) { return COLOR_SCHEMES.includes(v) ? v : 'shape'; }

// s red, p sky blue, d green, f violet: drawn as the renderer draws them
// (full value, both themes) all pairs pass the dataviz palette validator's
// CVD check (worst: d/f deutan dE 17).
export const L_HUE = { 0: 350, 1: 200, 2: 110, 3: 270 };

// Shell colors, shared with the radial chart (js/radial.js), validated there
// on the dark panel in shell order; at the renderer's full value their hues
// also pass the CVD check for neighbouring shells (worst n1/n2 dE 20).
export const SHELL_COLORS = { 1: '#ca1616', 2: '#b88b05', 3: '#0f91bd', 4: '#8b6cf0', 5: '#3f9a6a', 6: '#c2569a' };

function hexHue(hex) {
  const [r, g, b] = [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16) / 255);
  const max = Math.max(r, g, b), d = max - Math.min(r, g, b);
  if (d === 0) return 0;
  const h = max === r ? (g - b) / d : max === g ? 2 + (b - r) / d : 4 + (r - g) / d;
  return (h * 60 + 360) % 360;
}
const SHELL_HUE = Object.fromEntries(Object.entries(SHELL_COLORS).map(([n, c]) => [n, hexHue(c)]));

// { hue (degrees), chroma in [0, 1] } of an orbital of `element` in `scheme`.
export function orbitalColor(orbital, element, scheme = 'shape') {
  switch (scheme) {
    case 'gray': return { hue: 0, chroma: 0 };
    case 'subshell': return { hue: L_HUE[orbital.l], chroma: depthChroma(orbital, element) };
    case 'shell': return { hue: SHELL_HUE[orbital.n], chroma: 1 };
    case 'valence': {
      const [hue, chroma] = PALETTE_COLORS.valence[isValence(orbital, element) ? 'valence' : 'core'];
      return { hue, chroma };
    }
    case 'phase': case 'momentum': return { hue: 0, chroma: 0 };
    default: return { hue: orbital.hue, chroma: orbitalChroma(orbital, element) };
  }
}

// Valence orbitals: the outermost shell plus an open (n-1)d subshell.
export function isValence(orbital, element) {
  const subs = element.subshells || {};
  const count = subs[orbital.subshell] || 0;
  if (!count) return false;
  return valenceSubshells(element).has(orbital.subshell);
}

// Valence subshells: the period's ns/np, plus an open (n-1)d or (n-2)f.
// If none of those is occupied (e.g. Pd, [Kr] 4d10), fall back to the
// last occupied subshell in Aufbau order.
export function valenceSubshells(element) {
  const subs = element.subshells || {};
  const p = element.period;
  const out = new Set();
  for (const id of [`${p}s`, `${p}p`]) if (subs[id]) out.add(id);
  const d = `${p - 1}d`, f = `${p - 2}f`;
  if (subs[d] > 0 && subs[d] < 10) out.add(d);
  if (subs[f] > 0 && subs[f] < 14) out.add(f);
  if (out.size === 0) {
    const last = AUFBAU.filter(id => subs[id]).pop();
    if (last) out.add(last);
  }
  return out;
}
