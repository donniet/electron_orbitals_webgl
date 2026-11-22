// Math helpers + CPU reference implementations that mirror the shader.

import { radialPolyRef } from './orbitals.js';

// ------------------------------------------------------------ quaternions
// q = [x, y, z, w]; represents the atom's orientation in view space:
// v_view = rotate(q, v_atom).

export function qIdentity() { return [0, 0, 0, 1]; }

export function qAxisAngle(axis, angle) {
  const n = Math.hypot(axis[0], axis[1], axis[2]) || 1;
  const s = Math.sin(angle / 2) / n;
  return [axis[0] * s, axis[1] * s, axis[2] * s, Math.cos(angle / 2)];
}

export function qMul(a, b) {
  const [ax, ay, az, aw] = a, [bx, by, bz, bw] = b;
  return [
    aw * bx + ax * bw + ay * bz - az * by,
    aw * by - ax * bz + ay * bw + az * bx,
    aw * bz + ax * by - ay * bx + az * bw,
    aw * bw - ax * bx - ay * by - az * bz,
  ];
}

export function qNormalize(q) {
  const n = Math.hypot(q[0], q[1], q[2], q[3]) || 1;
  return [q[0] / n, q[1] / n, q[2] / n, q[3] / n];
}

// Row-major 3x3 rotation matrix of q (as nested arrays): M[r][c].
export function qToMat3(q) {
  const [x, y, z, w] = q;
  return [
    [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
    [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
    [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
  ];
}

export function mat3MulVec(m, v) {
  return [
    m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
    m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
    m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
  ];
}

export function mat3Transpose(m) {
  return [0, 1, 2].map(r => [0, 1, 2].map(c => m[c][r]));
}

// Column-major Float32Array for gl.uniformMatrix3fv.
export function mat3ToGL(m) {
  return new Float32Array([
    m[0][0], m[1][0], m[2][0],
    m[0][1], m[1][1], m[2][1],
    m[0][2], m[1][2], m[2][2],
  ]);
}

// View->atom matrix (the shader's uRot) for orientation q.
export function viewToAtom(q) {
  return mat3Transpose(qToMat3(q));
}

// ------------------------------------------------------------ orbitals
// Hydrogen-like <r^2> = n^2 / (2 Z^2) * (5 n^2 + 1 - 3 l (l + 1)) (Bohr^2).
export function meanR2(n, l, z) {
  return (n * n) / (2 * z * z) * (5 * n * n + 1 - 3 * l * (l + 1));
}

// Radial densities live in orbitals.js; re-exported for convenience.
export { radialDensity } from './orbitals.js';

// "Horizon": 2 standard deviations of the electron position
// (sigma = sqrt(<r^2>)) for the most diffuse occupied orbital.
export function horizon(element, orbitals) {
  let h = 0;
  orbitals.forEach((o, i) => {
    if (element.occupancy[i] > 0) {
      h = Math.max(h, 2 * Math.sqrt(meanR2(o.n, o.l, element.zeff[i])));
    }
  });
  return h;
}

// ------------------------------------------------------------ tables
// Linear interpolation of moment I_k (k = 0, 2, 4, 6) of table row `row` at
// scaled impact parameter u, or of the odd half-line moment J_k (k = 1, 3, 5).
// Mirrors momentAt() in the fragment shader.
export function momentLookup(tables, row, k, u) {
  const { N, uMax } = tables;
  const data = k % 2 === 0 ? tables.data : tables.odd;
  const x = (u / uMax) * (N - 1);
  if (x >= N - 1) return 0;
  const i = Math.floor(x);
  const f = x - i;
  const base = (row * N + i) * 4 + (k >> 1);
  return data[base] * (1 - f) + data[base + 4] * f;
}

function dot3(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
function mulQ(Q, v) {
  return [dot3(Q[0], v), dot3(Q[1], v), dot3(Q[2], v)];
}

// Column density in reference units (Z0 = n) for scaled closest point P
// and unit direction d. Mirrors the shader. cut = 0: whole ray; cut = -1 / +1:
// twice the half of the ray behind / in front of the plane through the
// nucleus (t <= 0 / t >= 0), i.e. even terms + 2 cut * odd terms.
export function columnRef(tables, orbital, P, d, cut = 0) {
  const u = Math.hypot(P[0], P[1], P[2]);
  const I = k => momentLookup(tables, orbital.row, k, u);   // odd k: J_k
  if (orbital.l === 0) return I(0) / (4 * Math.PI);
  let even, odd;
  if (orbital.l === 1) {
    const ep = dot3(orbital.e, P), ed = dot3(orbital.e, d);
    even = ep * ep * I(0) + ed * ed * I(2);
    odd = () => 2 * ep * ed * I(1);
  } else if (orbital.l === 2) {
    const QP = mulQ(orbital.Q, P), Qd = mulQ(orbital.Q, d);
    const a = dot3(P, QP), beta = dot3(P, Qd), c = dot3(d, Qd);
    even = a * a * I(0) + (4 * beta * beta + 2 * a * c) * I(2) + c * c * I(4);
    odd = () => 4 * a * beta * I(1) + 4 * beta * c * I(3);
  } else {
    // f: C(P + t d)^3 = a0 + a1 t + a2 t^2 + a3 t^3
    const a0 = cform(orbital.C, P, P, P), a1 = 3 * cform(orbital.C, P, P, d);
    const a2 = 3 * cform(orbital.C, P, d, d), a3 = cform(orbital.C, d, d, d);
    even = a0 * a0 * I(0) + (a1 * a1 + 2 * a0 * a2) * I(2) +
           (a2 * a2 + 2 * a1 * a3) * I(4) + a3 * a3 * I(6);
    odd = () => 2 * a0 * a1 * I(1) + 2 * (a0 * a3 + a1 * a2) * I(3) + 2 * a2 * a3 * I(5);
  }
  return cut ? even + 2 * cut * odd() : even;
}

// Angular amplitude along the scaled ray P + t d as a cubic in t:
// A(t) = c[0] + c[1] t + c[2] t^2 + c[3] t^3, |psi|^2 = q2(r) A(t)^2.
// Mirrors angularPoly() in the shader.
export function angularPoly(orbital, P, d) {
  if (orbital.l === 0) return [1 / Math.sqrt(4 * Math.PI), 0, 0, 0];
  if (orbital.l === 1) return [dot3(orbital.e, P), dot3(orbital.e, d), 0, 0];
  if (orbital.l === 2) {
    const QP = mulQ(orbital.Q, P), Qd = mulQ(orbital.Q, d);
    return [dot3(P, QP), 2 * dot3(P, Qd), dot3(d, Qd), 0];
  }
  return [cform(orbital.C, P, P, P), 3 * cform(orbital.C, P, P, d),
          3 * cform(orbital.C, P, d, d), cform(orbital.C, d, d, d)];
}

// Radial factor q2(u) = R^2 / u^(2l) of table row `row` (linear interpolation,
// 0 beyond the table). Mirrors radialAt() in the shader.
export function radialLookup(tables, row, u) {
  const { N, uMax, radial } = tables;
  const x = (u / uMax) * (N - 1);
  if (x >= N - 1) return 0;
  const i = Math.floor(x);
  const f = x - i;
  return radial[row * N + i] * (1 - f) + radial[row * N + i + 1] * f;
}

export const NEAREST_STEPS = 48;
export const NO_DEPTH = -1e30;

// Depth (reference units along d, larger = nearer the viewer) of the nearest
// of k detections of `orbital` on the scaled ray P + t d, by marching the
// density from the front in s, t = H s^3; see nearestDepthRef() in the
// shader, which this mirrors. massRef: column of the drawn part of the ray;
// v uniform in [0, 1); cut: only t <= 0 is drawn; returns NO_DEPTH once
// behind tStop.
export function nearestDepthRef(tables, orbital, P, d, k, massRef, v, cut = false,
                                tStop = NO_DEPTH, steps = NEAREST_STEPS) {
  const u2 = dot3(P, P);
  const ext = tables.rows[orbital.row].extent;
  if (u2 >= ext * ext) return 0;
  const H = Math.sqrt(ext * ext - u2);
  const sHi = cut ? 0 : 1;
  const hs = (sHi + 1) / steps;
  const A = angularPoly(orbital, P, d);
  const dens = s => {                       // q2(r) A(t)^2 dt/ds
    const t = H * s * s * s;
    const a = A[0] + t * (A[1] + t * (A[2] + t * A[3]));
    return radialLookup(tables, orbital.row, Math.sqrt(u2 + t * t)) * a * a * 3 * H * s * s;
  };
  const target = (1 - Math.pow(v, 1 / k)) * massRef;
  let fPrev = dens(sHi), acc = 0;
  for (let j = 1; j <= steps; j++) {
    const s = sHi - j * hs;
    const f = dens(s);
    const seg = 0.5 * (fPrev + f) * hs;
    if (acc + seg >= target) {
      const sc = s + hs * (1 - (target - acc) / Math.max(seg, 1e-30));
      return H * sc * sc * sc;
    }
    acc += seg;
    fPrev = f;
    if (H * s * s * s < tStop) return NO_DEPTH;
  }
  return -H;
}

// ------------------------------------------------------------ momentum (streaks)
// A streak is one sample of the orbital's momentum distribution |phi(p)|^2
// (phi the Fourier transform of psi): direction with density A(p)^2 (psi's
// own angular part), scaled magnitude q = p / kappa from the inverse-CDF
// table. Independent of the detection's position; see the shader's
// "Momentum (streaks)" section, which these mirror.
export const MOMENTUM_TRIES = 48;
export const STREAM_MOMENTUM_DIR = 2 << 20;
export const STREAM_MOMENTUM_Q = 3 << 20;

const AMAX_CACHE = new WeakMap();

// Bound on A(x)^2 over unit vectors x: dense Fibonacci sphere + 2% margin
// (exact maxima are 3/4pi for p, max |eigenvalue|^2 of Q for d). The shader
// gets it as c2.w; s needs none (uniform directions).
export function angularMax2(orbital) {
  if (orbital.l === 0) return 1 / (4 * Math.PI);
  if (AMAX_CACHE.has(orbital)) return AMAX_CACHE.get(orbital);
  const n = 20000, zero = [0, 0, 0];
  let best = 0;
  for (let i = 0; i < n; i++) {
    const z = 1 - (2 * (i + 0.5)) / n, r = Math.sqrt(Math.max(0, 1 - z * z));
    const ph = i * 2.399963229728653;
    const a = angularPoly(orbital, [r * Math.cos(ph), r * Math.sin(ph), z], zero)[0];
    best = Math.max(best, a * a);
  }
  AMAX_CACHE.set(orbital, best * 1.02);
  return best * 1.02;
}

const u01 = v => (v >>> 8) / 16777216;

// Unit momentum direction for pixel key [x, y, seed]. Mirrors momentumDirection().
export function momentumDirection(orbital, key, amax2 = angularMax2(orbital)) {
  let dir = [0, 0, 1];
  for (let j = 0; j < MOMENTUM_TRIES; j++) {
    const r = pcg4d(key[0] >>> 0, key[1] >>> 0, key[2] >>> 0, (STREAM_MOMENTUM_DIR + j) >>> 0);
    const z = 2 * u01(r[0]) - 1, s = Math.sqrt(Math.max(0, 1 - z * z));
    const ph = 2 * Math.PI * u01(r[1]);
    dir = [s * Math.cos(ph), s * Math.sin(ph), z];
    if (orbital.l === 0) return dir;
    const a = angularPoly(orbital, dir, [0, 0, 0])[0];
    if (u01(r[2]) * amax2 <= a * a) return dir;
  }
  return dir;
}

// Scaled momentum q of table row `row` for uniform F. Mirrors momentumQ().
export function momentumQ(tables, row, F) {
  const { N, momentum } = tables;
  const x = Math.min(Math.max(F * N - 0.5, 0), N - 1);
  const i = Math.min(Math.floor(x), N - 2);
  const f = x - i;
  return momentum[row * N + i] * (1 - f) + momentum[row * N + i + 1] * f;
}

export const STREAM_MOMENTUM_LOCAL = 4 << 20;

// Hessian (row-major 3x3) of ln|psi_ref| at scaled point X; see
// hessLogPsi() in the shader, which this mirrors. psi_ref ~ exp(-u) P(u) A(X).
export function hessLogPsi(orbital, X) {
  const u = Math.max(Math.hypot(X[0], X[1], X[2]), 1e-4);
  const r = X.map(c => c / u);
  const cs = radialPolyRef(orbital.n, orbital.l).poly;
  let p0 = 0, p1 = 0, p2 = 0;
  for (let k = cs.length - 1; k >= 0; k--) {
    p2 = p2 * u + 2 * p1;
    p1 = p1 * u + p0;
    p0 = p0 * u + cs[k];
  }
  p0 = p0 < 0 ? Math.min(p0, -1e-6) : Math.max(p0, 1e-6);
  const lp = p1 / p0, g1 = -1 + lp, g2 = p2 / p0 - lp * lp;
  const H = [0, 1, 2].map(i => [0, 1, 2].map(j =>
    g2 * r[i] * r[j] + (g1 / u) * ((i === j ? 1 : 0) - r[i] * r[j])));
  if (orbital.l === 0) return H;
  let A, g, HA;
  if (orbital.l === 1) {
    A = dot3(orbital.e, X); g = orbital.e; HA = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  } else if (orbital.l === 2) {
    const QX = mulQ(orbital.Q, X);
    A = dot3(X, QX); g = QX.map(c => 2 * c); HA = orbital.Q.map(row => row.map(c => 2 * c));
  } else {
    const C = orbital.C;
    const M = [0, 1, 2].map(i => [0, 1, 2].map(j => C[0][i][j] * X[0] + C[1][i][j] * X[1] + C[2][i][j] * X[2]));
    const MX = mulQ(M, X);
    A = dot3(X, MX); g = MX.map(c => 3 * c); HA = M.map(row => row.map(c => 6 * c));
  }
  const eps = 1e-6 * u ** orbital.l;
  A = A < 0 ? Math.min(A, -eps) : Math.max(A, eps);
  return H.map((row, i) => row.map((h, j) => h + HA[i][j] / A - g[i] * g[j] / (A * A)));
}

// Local momentum second moment T = -Hess ln|psi| / 2 (units of kappa^2).
export function localMomentumTensor(orbital, X) {
  return hessLogPsi(orbital, X).map(row => row.map(c => -0.5 * c));
}

// Sign of psi_ref at scaled point X: sign P(u) * sign A(X). Mirrors
// psiSign() in the shader (phase palette).
export function psiSign(orbital, X) {
  const u = Math.hypot(X[0], X[1], X[2]);
  const cs = radialPolyRef(orbital.n, orbital.l).poly;
  let P = 0;
  for (let k = cs.length - 1; k >= 0; k--) P = P * u + cs[k];
  let A = 1;
  if (orbital.l === 1) A = dot3(orbital.e, X);
  else if (orbital.l === 2) A = dot3(X, mulQ(orbital.Q, X));
  else if (orbital.l === 3) {
    const C = orbital.C;
    const M = [0, 1, 2].map(i => [0, 1, 2].map(j => C[0][i][j] * X[0] + C[1][i][j] * X[1] + C[2][i][j] * X[2]));
    A = dot3(X, mulQ(M, X));
  }
  return (P < 0) !== (A < 0) ? -1 : 1;
}

// Local RMS momentum (atomic units) at scaled point X of an orbital with
// kappa = Z_eff / n: kappa sqrt(trace of T, negative part dropped). Mirrors
// localMomentum() in the shader (momentum palette).
export function localMomentum(orbital, X, kappa) {
  const { lam } = jacobiEigen(localMomentumTensor(orbital, X));
  const p = kappa * Math.sqrt(lam.reduce((a, l) => a + Math.max(l, 0), 0));
  return p < 1e30 ? Math.max(p, 1e-30) : 1e30;
}

// Cyclic Jacobi eigen-decomposition of symmetric 3x3 A (row-major):
// { lam: eigenvalues, V: row-major matrix whose columns are the
// eigenvectors }. Mirrors jacobiEigen() in the shader.
export function jacobiEigen(A0) {
  let A = A0.map(r => r.slice());
  let V = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  const mul = (X, Y) => [0, 1, 2].map(i => [0, 1, 2].map(j => X[i][0] * Y[0][j] + X[i][1] * Y[1][j] + X[i][2] * Y[2][j]));
  for (let sweep = 0; sweep < 6; sweep++) {
    for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
      const apq = A[p][q];
      if (Math.abs(apq) <= 1e-12 * (Math.abs(A[p][p]) + Math.abs(A[q][q]))) continue;
      const theta = (A[q][q] - A[p][p]) / (2 * apq);
      const t = (theta >= 0 ? 1 : -1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      const J = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
      J[p][p] = c; J[q][q] = c; J[p][q] = s; J[q][p] = -s;
      A = mul(mul(mat3Transpose(J), A), J);
      V = mul(V, J);
    }
  }
  return { lam: [A[0][0], A[1][1], A[2][2]], V };
}

// Momentum (units of kappa, atom frame) of a streak for pixel key
// [x, y, seed] and the detection at scaled point X; mode 'global' | 'local'
// | 'axis'. Mirrors the momentum part of streakHalf().
export function streakMomentum(tables, orbital, key, mode = 'global', X = null) {
  if (mode === 'global') {
    const dir = momentumDirection(orbital, key);
    const q = momentumQ(tables, orbital.row,
      u01(pcg4d(key[0] >>> 0, key[1] >>> 0, key[2] >>> 0, STREAM_MOMENTUM_Q >>> 0)[0]));
    return dir.map(c => q * c);
  }
  const { lam, V } = jacobiEigen(localMomentumTensor(orbital, X));
  const l = lam.map(v => Math.max(v, 0));
  const r = pcg4d(key[0] >>> 0, key[1] >>> 0, key[2] >>> 0, STREAM_MOMENTUM_LOCAL >>> 0).map(u01);
  const m0 = Math.sqrt(-2 * Math.log(1 - r[0])), m1 = Math.sqrt(-2 * Math.log(1 - r[2]));
  const z = [m0 * Math.cos(2 * Math.PI * r[1]), m0 * Math.sin(2 * Math.PI * r[1]), m1 * Math.cos(2 * Math.PI * r[3])];
  const zV = [0, 1, 2].map(k => V[0][k] * z[0] + V[1][k] * z[1] + V[2][k] * z[2]);
  const fromBasis = w => [0, 1, 2].map(i => V[i][0] * w[0] + V[i][1] * w[1] + V[i][2] * w[2]);
  if (mode === 'axis') {
    // principal axis; random within the eigenspace of coinciding largest
    // eigenvalues (spectral weight (lam / max)^8 on z)
    const lmax = Math.max(...l);
    const a = fromBasis(l.map((v, k) => (v / Math.max(lmax, 1e-30)) ** 8 * zV[k]));
    const al = Math.hypot(...a);
    return al > 0 ? a.map(c => Math.sqrt(lmax) * c / al) : [0, 0, 0];
  }
  // symmetric square root V sqrt(L) V^T z (independent of the eigenbasis)
  return fromBasis(l.map((v, k) => Math.sqrt(v) * zV[k]));
}

// Streak half-vector (render-buffer pixels): the streak momentum rotated by
// atomToView (row-major, = transpose of the shader's uRot), projected on the
// screen, times len, clamped to rMax. Mirrors streakHalf().
export function streakHalf(tables, orbital, key, atomToView, len, rMax, mode = 'global', X = null) {
  const v = mat3MulVec(atomToView, streakMomentum(tables, orbital, key, mode, X));
  let h = [len * v[0], len * v[1]];
  const l = Math.hypot(h[0], h[1]);
  if (!(l < 1e30)) return [0, 0];
  if (l > rMax) h = [h[0] * rMax / l, h[1] * rMax / l];
  return h;
}

function cform(C, a, b, c) {
  let s = 0;
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) for (let k = 0; k < 3; k++) {
    s += C[i][j][k] * a[i] * b[j] * c[k];
  }
  return s;
}

// Column density (per electron in the orbital, 1 / Bohr^2) along the ray
// with closest point p (atom frame, p.d = 0) and unit direction d.
export function columnDensity(tables, orbital, zeff, p, d, cut = 0) {
  const k = zeff / orbital.n;
  return k * k * columnRef(tables, orbital, [k * p[0], k * p[1], k * p[2]], d, cut);
}

// ------------------------------------------------------------ random
// PCG hash (same as shader). Returns uint32.
export function pcg(v) {
  const state = (Math.imul(v >>> 0, 747796405) + 2891336453) >>> 0;
  const shift = ((state >>> 28) + 4) >>> 0;
  const word = Math.imul(((state >>> shift) ^ state) >>> 0, 277803737) >>> 0;
  return ((word >>> 22) ^ word) >>> 0;
}

// pcg4d (Jarzynski & Olano 2020): 128-bit bijective hash. Mirrors the shader.
export function pcg4d(x, y, z, w) {
  x = (Math.imul(x, 1664525) + 1013904223) >>> 0;
  y = (Math.imul(y, 1664525) + 1013904223) >>> 0;
  z = (Math.imul(z, 1664525) + 1013904223) >>> 0;
  w = (Math.imul(w, 1664525) + 1013904223) >>> 0;
  x = (x + Math.imul(y, w)) >>> 0; y = (y + Math.imul(z, x)) >>> 0;
  z = (z + Math.imul(x, y)) >>> 0; w = (w + Math.imul(y, z)) >>> 0;
  x = (x ^ (x >>> 16)) >>> 0; y = (y ^ (y >>> 16)) >>> 0;
  z = (z ^ (z >>> 16)) >>> 0; w = (w ^ (w >>> 16)) >>> 0;
  x = (x + Math.imul(y, w)) >>> 0; y = (y + Math.imul(z, x)) >>> 0;
  z = (z + Math.imul(x, y)) >>> 0; w = (w + Math.imul(y, z)) >>> 0;
  return [x, y, z, w];
}

// Counter-based RNG for stream (x, y, seed), uniform in [0, 1) from the top
// 24 bits, consecutive pcg4d outputs. The shader uses the same hash and
// uniform conversion, keyed per orbital sample (x, y, seed, slot) so each
// sample's uniforms come from one hash; the uniforms are i.i.d. either way,
// so the CPU sampler below draws from the same distribution.
export function makeRng(seed, x = 0, y = 0) {
  let ctr = 0, buf = null, idx = 4;
  return () => {
    if (idx >= 4) { buf = pcg4d(x >>> 0, y >>> 0, seed >>> 0, ctr >>> 0); ctr++; idx = 0; }
    return (buf[idx++] >>> 8) / 16777216;
  };
}

// log(k!): exact table for small k, Stirling series above. Mirrors the shader.
const LOG_FACT = [0, 0, 0.69314718, 1.79175947, 3.17805383,
  4.78749174, 6.57925121, 8.52516136, 10.6046029, 12.8018275];
export function logFactorial(k) {
  if (k < 10) return LOG_FACT[k];
  const ik = 1 / k, ik2 = ik * ik;
  return (k + 0.5) * Math.log(k) - k + 0.91893853320467 +
    ik * (1 / 12 - ik2 * (1 / 360 - ik2 / 1260));
}

// Exact Poisson sample, same algorithm as the shader:
//  lambda < 10 : inversion with one uniform;
//  lambda >= 10: PTRS transformed rejection (Hormann 1993).
export function poisson(lambda, rand) {
  if (!(lambda > 0)) return 0;
  if (lambda < 10) {
    const u = rand();
    let p = Math.exp(-lambda), F = p, k = 0;
    for (let i = 0; i < 64; i++) {
      if (u <= F) return k;
      k++;
      p *= lambda / k;
      F += p;
    }
    return Math.floor(lambda);       // beyond the cumulative limit: mode, as in the shader
  }
  const slam = Math.sqrt(lambda), loglam = Math.log(lambda);
  const b = 0.931 + 2.53 * slam;
  const a = -0.059 + 0.02483 * b;
  const invalpha = 1.1239 + 1.1328 / (b - 3.4);
  const vr = 0.9277 - 3.6224 / (b - 2);
  for (let i = 0; i < 16; i++) {
    const U = rand() - 0.5;
    const V = rand();
    const us = 0.5 - Math.abs(U);
    const k = Math.floor((2 * a / us + b) * U + lambda + 0.43);
    if (us >= 0.07 && V <= vr) return k;
    if (k < 0 || (us < 0.013 && V > us)) continue;
    if (Math.log(V) + Math.log(invalpha) - Math.log(a / (us * us) + b) <=
        -lambda + k * loglam - logFactorial(k)) return k;
  }
  return Math.floor(lambda + 0.5);
}

// ------------------------------------------------------------ color
export function hsvToRgb(h, s, v) {
  // h in degrees
  const f = n => {
    const k = (n + h / 60) % 6;
    return v - v * s * Math.max(0, Math.min(k, 4 - k, 1));
  };
  return [f(5), f(3), f(1)];
}

export function hslCss(hue, s = 90, l = 55) {
  return `hsl(${hue} ${s}% ${l}%)`;
}
