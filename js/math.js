// Math helpers + CPU reference implementations that mirror the shader.

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

// Normalised radial probability density P(r) = r^2 R_nl(r)^2 (per a0) for a
// hydrogen-like orbital with effective charge z. Integrates to 1 over r.
export function radialDensity(n, l, z, r) {
  const x = z * r;
  if (n === 1) return 4 * z * x * x * Math.exp(-2 * x);
  if (l === 0) return (z / 8) * x * x * (2 - x) * (2 - x) * Math.exp(-x);
  return (z / 24) * x * x * x * x * Math.exp(-x);
}

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
// Linear interpolation of channel `ch` at scaled impact parameter b.
// Mirrors tableAt() in the fragment shader.
export function tableLookup(tables, ch, b) {
  const { N, bMax, data } = tables;
  const x = (b / bMax) * (N - 1);
  if (x >= N - 1) return 0;
  const i = Math.floor(x);
  const f = x - i;
  return data[4 * i + ch] * (1 - f) + data[4 * (i + 1) + ch] * f;
}

// Column density (electrons / Bohr^2 per electron in the orbital) along the
// ray with closest point p (atom frame, p.d = 0) and unit direction d.
export function columnDensity(tables, orbital, zeff, p, d) {
  const ps = [p[0] * zeff, p[1] * zeff, p[2] * zeff];
  const b = Math.hypot(ps[0], ps[1], ps[2]);
  const z2 = zeff * zeff;
  if (orbital.kind === 's1') return z2 * tableLookup(tables, 2, b);
  if (orbital.kind === 's2') return z2 * tableLookup(tables, 3, b);
  const a = orbital.axis;
  return z2 * (ps[a] * ps[a] * tableLookup(tables, 0, b) +
               d[a] * d[a] * tableLookup(tables, 1, b));
}

// ------------------------------------------------------------ random
// PCG hash (same as shader). Returns uint32.
export function pcg(v) {
  const state = (Math.imul(v >>> 0, 747796405) + 2891336453) >>> 0;
  const shift = ((state >>> 28) + 4) >>> 0;
  const word = Math.imul(((state >>> shift) ^ state) >>> 0, 277803737) >>> 0;
  return ((word >>> 22) ^ word) >>> 0;
}

// Stateful RNG over pcg, uniform in [0, 1). Mirrors rand() in the shader.
export function makeRng(seed) {
  let s = seed >>> 0;
  return () => {
    s = pcg(s);
    return s / 4294967296;
  };
}

// Poisson sample, same algorithm as the shader.
export function poisson(lambda, rand) {
  if (!(lambda > 0)) return 0;
  if (lambda < 30) {
    const L = Math.exp(-lambda);
    let p = 1, k = 0;
    for (let i = 0; i < 64; i++) {
      p *= rand();
      if (p <= L) break;
      k++;
    }
    return k;
  }
  const u1 = Math.max(rand(), 1e-7), u2 = rand();
  const n = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  return Math.max(0, Math.floor(lambda + Math.sqrt(lambda) * n + 0.5));
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
