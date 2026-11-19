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
];
export const SUBSHELL = Object.fromEntries(SUBSHELLS.map(s => [s.id, s]));

// Aufbau filling order.
export const AUFBAU = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p'];

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

// Hues (degrees). s: warm sequence by shell; p: one hue per axis in every
// shell; d: five hues of their own. Shells beyond 2 are provisional until
// the period 3/4 colour pass.
const S_HUE = { 1: 0, 2: 45, 3: 20, 4: 345 };
const P_HUE = { x: 195, y: 275, z: 125 };
const D_HUE = { xy: 310, yz: 230, xz: 165, 'x2-y2': 85, z2: 60 };

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
  return D_SHAPES.map(d => ({
    ...base, id: s.id + d.m, m: d.m, Q: d.Q,
    label: `${s.id}<sub>${d.label}</sub>`, hue: D_HUE[d.m],
  }));
}

// Shader order. The first five match the original H-Ne layout so frames for
// those elements consume random numbers in the same order.
export const ORBITALS = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p'].flatMap(makeOrbitals);
export const NUM_ORBITALS = ORBITALS.length;   // 18

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
