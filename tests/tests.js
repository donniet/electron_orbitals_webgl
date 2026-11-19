// In-browser tests. Open http://localhost:8000/tests/test.html
// Results are written into the page; document.title becomes
// "PASS n/n" or "FAIL k/n" so headless runs can read it.

import { ELEMENTS, ORBITALS, aufbau, hundFill } from '../js/elements.js';
import { SUBSHELLS, orbitalIndex, orbitalsOf, isValence } from '../js/orbitals.js';
import {
  qAxisAngle, qMul, qToMat3, mat3MulVec, viewToAtom, mat3ToGL,
  meanR2, horizon, momentLookup, radialDensity, columnDensity, makeRng, poisson, pcg, hsvToRgb,
} from '../js/math.js';
import { createRenderer, viewGeometry, orbitalPeak, orbitalGains, MAX_HORIZON } from '../js/renderer.js';
import { loadTables } from '../js/tables.js';

const TABLES = await loadTables();
const P = id => orbitalIndex(id);

const results = [];
function test(name, fn) {
  try {
    fn();
    results.push({ name, ok: true });
  } catch (e) {
    results.push({ name, ok: false, msg: e.message });
  }
}
function assert(cond, msg = 'assertion failed') { if (!cond) throw new Error(msg); }
function near(a, b, tol, msg = '') {
  if (!(Math.abs(a - b) <= tol)) throw new Error(`${msg} expected ${b} got ${a} (tol ${tol})`);
}

// ------------------------------------------------------------ math
test('quaternion rotates x to y about z by 90 deg', () => {
  const m = qToMat3(qAxisAngle([0, 0, 1], Math.PI / 2));
  const v = mat3MulVec(m, [1, 0, 0]);
  near(v[0], 0, 1e-12); near(v[1], 1, 1e-12); near(v[2], 0, 1e-12);
});

test('quaternion composition matches matrix product', () => {
  const a = qAxisAngle([1, 2, 3], 0.7), b = qAxisAngle([-1, 0, 2], 1.3);
  const v = [0.3, -1.2, 0.8];
  const lhs = mat3MulVec(qToMat3(qMul(a, b)), v);
  const rhs = mat3MulVec(qToMat3(a), mat3MulVec(qToMat3(b), v));
  for (let i = 0; i < 3; i++) near(lhs[i], rhs[i], 1e-12);
});

test('viewToAtom is the inverse rotation', () => {
  const q = qAxisAngle([0.2, 1, -0.4], 2.1);
  const v = [0.5, 0.25, -2];
  const back = mat3MulVec(viewToAtom(q), mat3MulVec(qToMat3(q), v));
  for (let i = 0; i < 3; i++) near(back[i], v[i], 1e-12);
});

test('mat3ToGL is column-major', () => {
  const m = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
  assert(mat3ToGL(m).join() === '1,4,7,2,5,8,3,6,9');
});

test('<r^2> hydrogen 1s = 3 a0^2, 2s = 42, 2p = 30', () => {
  near(meanR2(1, 0, 1), 3, 1e-12);
  near(meanR2(2, 0, 1), 42, 1e-12);
  near(meanR2(2, 1, 1), 30, 1e-12);
  near(meanR2(1, 0, 2), 0.75, 1e-12);
});

test('horizon uses the most diffuse occupied orbital', () => {
  near(horizon(ELEMENTS[0], ORBITALS), 2 * Math.sqrt(3), 1e-12);
  const li = ELEMENTS[2];
  near(horizon(li, ORBITALS), 2 * Math.sqrt(42) / li.zeff[P('2s')], 1e-12);
});

function integrate(f, a, b, n = 4000) {
  const h = (b - a) / n;
  let s = f(a) + f(b);
  for (let i = 1; i < n; i++) s += (i & 1 ? 4 : 2) * f(a + i * h);
  return (s * h) / 3;
}

test('radialDensity is normalised and matches <r^2>', () => {
  for (const [n, l] of [[1, 0], [2, 0], [2, 1]]) {
    for (const z of [1, 1.688, 5.758]) {
      near(integrate(r => radialDensity(n, l, z, r), 0, 60 / z), 1, 1e-6, `norm ${n}${l} z=${z}`);
      near(integrate(r => r * r * radialDensity(n, l, z, r), 0, 60 / z), meanR2(n, l, z),
        1e-6 * meanR2(n, l, z), `<r2> ${n}${l} z=${z}`);
    }
  }
});

test('fixed scale: same view for every element; auto: He ~0.59x H', () => {
  const [H, He] = ELEMENTS;
  const fH = viewGeometry(H, 400, 300, 1, true), fHe = viewGeometry(He, 400, 300, 1, true);
  near(fH.halfExtent[1], fHe.halfExtent[1], 1e-12);
  near(fH.halfExtent[1], MAX_HORIZON / 0.9, 1e-12);
  const largest = ELEMENTS.reduce((a, e) => (horizon(e, ORBITALS) > horizon(a, ORBITALS) ? e : a));
  near(MAX_HORIZON, horizon(largest, ORBITALS), 1e-12);
  assert(ELEMENTS.length > 18 || largest.symbol === 'Na', `largest is ${largest.symbol}`);
  const aH = viewGeometry(H, 400, 300, 1, false), aHe = viewGeometry(He, 400, 300, 1, false);
  near(aHe.halfExtent[1] / aH.halfExtent[1], 1 / 1.688, 1e-9);
});

// ------------------------------------------------------------ elements
test('element data: sequential Z, occupancy sums to Z', () => {
  assert(ELEMENTS.length >= 18, `only ${ELEMENTS.length} elements`);
  ELEMENTS.forEach((el, i) => assert(el.z === i + 1, `${el.symbol} out of order`));
  for (const el of ELEMENTS) {
    const sum = el.occupancy.reduce((a, b) => a + b, 0);
    assert(sum === el.z, `${el.symbol}: ${sum} != ${el.z}`);
    el.occupancy.forEach((o, i) => assert(o === 0 || el.zeff[i] > 0, `${el.symbol} missing Z_eff`));
  }
});

test("p orbitals fill by Hund's rule", () => {
  const p = z => orbitalsOf('2p').map(i => ELEMENTS[z - 1].occupancy[i]).join('');
  assert(p(5) === '100' && p(6) === '110' && p(7) === '111');
  assert(p(8) === '211' && p(9) === '221' && p(10) === '222');
});

test("3p fills by Hund's rule (Al..Ar)", () => {
  const p = z => orbitalsOf('3p').map(i => ELEMENTS[z - 1].occupancy[i]).join('');
  assert(['100', '110', '111', '211', '221', '222'].every((v, k) => p(13 + k) === v));
});

test('shorthand configurations', () => {
  assert(ELEMENTS[10].configShort === '[Ne] 3s¹', ELEMENTS[10].configShort);
  assert(ELEMENTS[17].configShort === '[Ne] 3s² 3p⁶', ELEMENTS[17].configShort);
  assert(ELEMENTS[5].configShort === '1s² 2s² 2p²', 'period 2 stays in full');
});

test('valence orbitals: outer shell only', () => {
  const ids = z => ORBITALS.filter(o => isValence(o, ELEMENTS[z - 1])).map(o => o.subshell);
  assert([...new Set(ids(11))].join() === '3s', ids(11).join());
  assert([...new Set(ids(17))].join() === '3s,3p', ids(17).join());
});

test('outer p shells are paler (chroma) with the same axis hues', () => {
  const o = id => ORBITALS[orbitalIndex(id)];
  assert(o('2px').chroma === 1 && o('3px').chroma < 1 && o('4px').chroma < o('3px').chroma);
  assert(o('2px').hue === o('3px').hue && o('2pz').hue === o('4pz').hue);
});

test('electron configurations', () => {
  assert(ELEMENTS[0].config === '1s¹');
  assert(ELEMENTS[5].config === '1s² 2s² 2p²');
  assert(ELEMENTS[9].config === '1s² 2s² 2p⁶');
});

test('aufbau and Hund helpers', () => {
  assert(JSON.stringify(aufbau(7)) === JSON.stringify({ '1s': 2, '2s': 2, '2p': 3 }));
  assert(JSON.stringify(aufbau(21)) === JSON.stringify({ '1s': 2, '2s': 2, '2p': 6, '3s': 2, '3p': 6, '4s': 2, '3d': 1 }));
  assert(hundFill(7, 5).join('') === '22111');
});

// ------------------------------------------------------------ tables
test('table rows match orbitals.js subshell rows', () => {
  for (const sh of SUBSHELLS) {
    const r = TABLES.rows[sh.row];
    assert(r.name === sh.id && r.n === sh.n && r.l === sh.l, `${sh.id} -> row ${sh.row} is ${r.name}`);
  }
});

test('moment tables at u = 0 match closed forms', () => {
  // 1s: R^2 = 4 e^{-2r}; 2p (Z = 2): R^2 / r^2 = (4/3) e^{-2r}
  const r1s = TABLES.rows.findIndex(r => r.name === '1s');
  const r2p = TABLES.rows.findIndex(r => r.name === '2p');
  [4, 2, 6].forEach((v, j) => near(momentLookup(TABLES, r1s, 2 * j, 0), v, 1e-5));
  [4 / 3, 2 / 3, 2].forEach((v, j) => near(momentLookup(TABLES, r2p, 2 * j, 0), v, 1e-5));
  assert(momentLookup(TABLES, r1s, 0, 1e6) === 0, 'lookup past uMax must be 0');
});

test('radialDensity (general n, l) normalised with correct <r^2>', () => {
  for (const sh of SUBSHELLS) {
    for (const z of [1, 3.3]) {
      const hi = 60 * sh.n / z;
      near(integrate(r => radialDensity(sh.n, sh.l, z, r), 0, hi), 1, 1e-6, `norm ${sh.id}`);
      near(integrate(r => r * r * radialDensity(sh.n, sh.l, z, r), 0, hi), meanR2(sh.n, sh.l, z),
        1e-5 * meanR2(sh.n, sh.l, z), `<r2> ${sh.id}`);
    }
  }
});

test('column density integrates to 1 electron over the image plane', () => {
  // integrate over a grid in the view plane for a tilted ray direction
  const q = qAxisAngle([0.3, 1, 0.2], 0.9);
  const R = viewToAtom(q);
  const d = mat3MulVec(R, [0, 0, 1]);
  ORBITALS.forEach((o, i) => {
    const z = 1.3 + 0.4 * i;
    const ext = 30 * o.n / z, n = 260, hstep = 2 * ext / n;
    let acc = 0;
    for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) {
      const x = -ext + (a + 0.5) * hstep, y = -ext + (b + 0.5) * hstep;
      acc += columnDensity(TABLES, ORBITALS[i], z, mat3MulVec(R, [x, y, 0]), d);
    }
    near(acc * hstep * hstep, 1, 0.01, o.id);
  });
});

test('full p and d subshells are isotropic (Unsold)', () => {
  for (const sh of SUBSHELLS.filter(s => s.l > 0)) {
    const vals = [];
    for (const ang of [0, 0.4, 1.1, 2.0]) {
      const R = viewToAtom(qAxisAngle([1, 0.5, -0.3], ang));
      const p = mat3MulVec(R, [0.8, 0.3, 0]), d = mat3MulVec(R, [0, 0, 1]);
      vals.push(orbitalsOf(sh.id).reduce((s, i) => s + columnDensity(TABLES, ORBITALS[i], 2, p, d), 0));
    }
    for (const v of vals) near(v, vals[0], 1e-6 + 1e-5 * vals[0], sh.id);
  }
});

// ------------------------------------------------------------ random
test('pcg hash matches reference values', () => {
  // Reference: GLSL/C implementation of the same pcg_hash.
  assert(pcg(0) === 129708002, `pcg(0) = ${pcg(0)}`);
  assert(pcg(1) !== pcg(2));
});

test('Poisson sampler mean and variance ~ lambda', () => {
  const rand = makeRng(12345);
  for (const lam of [0.05, 0.7, 3, 12, 50]) {
    const n = 40000;
    let s = 0, s2 = 0;
    for (let i = 0; i < n; i++) { const k = poisson(lam, rand); s += k; s2 += k * k; }
    const mean = s / n, v = s2 / n - mean * mean;
    near(mean, lam, 5 * Math.sqrt(lam / n) + 0.01 * lam, `mean@${lam}`);
    near(v, lam, 0.06 * lam + 0.01, `var@${lam}`);
  }
});

test('hsvToRgb primaries', () => {
  assert(hsvToRgb(0, 1, 1).map(Math.round).join() === '1,0,0');
  assert(hsvToRgb(120, 1, 1).map(Math.round).join() === '0,1,0');
  assert(hsvToRgb(240, 0, 1).map(Math.round).join() === '1,1,1');
});

// ------------------------------------------------------------ equalize
test('orbital peaks: 1s/2s at b=0, 2p end-on (closed forms)', () => {
  const [H, , Li, , B] = ELEMENTS;
  near(orbitalPeak(TABLES, H, P('1s')), 1 / Math.PI, 1e-6);
  const z2s = Li.zeff[P('2s')];
  near(orbitalPeak(TABLES, Li, P('2s')), z2s * z2s / (8 * Math.PI), 1e-6);
  const z2p = B.zeff[P('2px')];
  near(orbitalPeak(TABLES, B, P('2px')), z2p * z2p / (8 * Math.PI), 1e-6);
});

test('equalized gains give every occupied orbital the same peak lambda', () => {
  const n = ELEMENTS[6];
  const g = orbitalGains(TABLES, n, 0.15, 50, true);
  n.occupancy.forEach((o, i) => {
    if (o > 0) near(orbitalPeak(TABLES, n, i) * g[i], 0.15 * 20, 1e-9, ORBITALS[i].id);
  });
  const phys = orbitalGains(TABLES, n, 0.15, 50, false);
  assert(phys.every(v => v === phys[0]), 'physical gains must be shared');
});

// ------------------------------------------------------------ GPU
function renderFrame(z, opts = {}) {
  const canvas = document.createElement('canvas');
  canvas.width = 256; canvas.height = 256;
  const r = createRenderer(canvas, TABLES);
  const state = {
    element: ELEMENTS[z - 1], orientation: opts.q || [0, 0, 0, 1], zoom: 1,
    seed: opts.seed ?? 1, density: opts.density ?? 0.3,
    equalize: opts.equalize ?? false, kMax: opts.kMax ?? 4, gamma: opts.gamma ?? 0.5,
    mode: opts.mode || 'dark', visible: ORBITALS.map(() => true), fixedScale: opts.fixedScale ?? false,
  };
  r.render(state);
  const px = new Uint8Array(256 * 256 * 4);
  r.gl.readPixels(0, 0, 256, 256, r.gl.RGBA, r.gl.UNSIGNED_BYTE, px);
  r.gl.getExtension('WEBGL_lose_context')?.loseContext();
  return px;
}

function brightness(px, x, y) { const i = 4 * (y * 256 + x); return px[i] + px[i + 1] + px[i + 2]; }

test('GPU: same seed reproduces, new seed differs', () => {
  const a = renderFrame(1, { seed: 7 }), b = renderFrame(1, { seed: 7 }), c = renderFrame(1, { seed: 8 });
  let same = true, diff = 0;
  for (let i = 0; i < a.length; i++) { if (a[i] !== b[i]) same = false; if (a[i] !== c[i]) diff++; }
  assert(same, 'same seed differs');
  assert(diff > 1000, 'different seed too similar');
});

test('GPU: hydrogen image is centred and radially falling', () => {
  // Average over several seeds; compare an inner disc to an outer ring.
  let inner = 0, outer = 0, nIn = 0, nOut = 0, lit = 0;
  for (let s = 1; s <= 4; s++) {
    const px = renderFrame(1, { seed: s });
    for (let y = 0; y < 256; y++) for (let x = 0; x < 256; x++) {
      const r = Math.hypot(x - 127.5, y - 127.5);
      const v = brightness(px, x, y);
      if (v > 0) lit++;
      if (r < 20) { inner += v; nIn++; } else if (r > 100 && r < 120) { outer += v; nOut++; }
    }
  }
  assert(lit > 1000, `too few lit pixels: ${lit}`);
  assert(inner / nIn > 3 * (outer / nOut), `inner ${inner / nIn} outer ${outer / nOut}`);
});

test('GPU: boron 2px lobes lie along view x when unrotated', () => {
  // Only the p electron is visible far from the centre along x.
  let alongX = 0, alongY = 0;
  for (let s = 1; s <= 4; s++) {
    const px = renderFrame(5, { seed: s, density: 1 });
    const g = viewGeometry(ELEMENTS[4], 256, 256, 1);
    const r = Math.round(2.0 / ELEMENTS[4].zeff[P('2px')] / g.bohrPerPixel); // ~ lobe peak
    for (let o = -6; o <= 6; o++) for (let t = -6; t <= 6; t++) {
      alongX += brightness(px, 128 + r + o, 128 + t) + brightness(px, 128 - r + o, 128 + t);
      alongY += brightness(px, 128 + t, 128 + r + o) + brightness(px, 128 + t, 128 - r + o);
    }
  }
  assert(alongX > 1.3 * alongY, `alongX ${alongX} alongY ${alongY}`);
});

test('GPU: equalize brightens nitrogen 2p relative to the core', () => {
  const ratio = eq => {
    let core = 0, lobe = 0;
    for (let s = 1; s <= 3; s++) {
      const px = renderFrame(7, { seed: s, equalize: eq, density: 0.15 });
      for (let y = 0; y < 256; y++) for (let x = 0; x < 256; x++) {
        const r = Math.hypot(x - 127.5, y - 127.5);
        const v = brightness(px, x, y);
        if (r < 6) core += v / 113; else if (r > 40 && r < 70) lobe += v / 9896;
      }
    }
    return lobe / core;
  };
  const phys = ratio(false), eq = ratio(true);
  assert(eq > 1.5 * phys, `equalized ${eq} vs physical ${phys}`);
});

test('GPU: fixed scale shows He smaller than H', () => {
  // Mean brightness in a ring at ~2.5 a0 from the nucleus.
  const ring = z => {
    const g = viewGeometry(ELEMENTS[z - 1], 256, 256, 1, true);
    const r0 = 2.1 / g.bohrPerPixel, r1 = 2.9 / g.bohrPerPixel;
    let sum = 0, n = 0;
    for (let s = 1; s <= 4; s++) {
      const px = renderFrame(z, { seed: s, fixedScale: true, density: 1, kMax: 1 });
      for (let y = 0; y < 256; y++) for (let x = 0; x < 256; x++) {
        const r = Math.hypot(x - 127.5, y - 127.5);
        if (r >= r0 && r < r1) { sum += brightness(px, x, y); n++; }
      }
    }
    return sum / n;
  };
  const h = ring(1), he = ring(2);
  assert(h > 2 * he, `H ring ${h} vs He ring ${he}`);
});

test('GPU: white mode background is white', () => {
  const px = renderFrame(1, { mode: 'white' });
  assert(brightness(px, 2, 2) === 765, 'corner not white');
});

// ------------------------------------------------------------ report
const pass = results.filter(r => r.ok).length;
const out = document.getElementById('results');
for (const r of results) {
  const li = document.createElement('li');
  li.className = r.ok ? 'ok' : 'fail';
  li.textContent = (r.ok ? '✔ ' : '✘ ') + r.name + (r.ok ? '' : ' — ' + r.msg);
  out.append(li);
}
const summary = `${pass === results.length ? 'PASS' : 'FAIL'} ${pass}/${results.length}`;
document.getElementById('summary').textContent = summary;
document.title = summary;
