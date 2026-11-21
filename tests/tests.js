// In-browser tests. Serve the repo root (python3 -m http.server 8000) and open
// http://localhost:8000/tests/test.html
// Results are written into the page; document.title becomes
// "PASS n/n" or "FAIL k/n" so headless runs can read it.

import { ELEMENTS, ORBITALS, aufbau, hundFill, stepElement } from '../js/elements.js';
import { SUBSHELLS, orbitalIndex, orbitalsOf, isValence, angularCoefs, orbitalChroma } from '../js/orbitals.js';
import {
  qAxisAngle, qMul, qToMat3, mat3MulVec, viewToAtom, mat3ToGL,
  meanR2, horizon, momentLookup, radialDensity, columnDensity, makeRng, poisson, pcg, hsvToRgb,
  columnRef, angularPoly, radialLookup, nearestDepthRef, NO_DEPTH,
} from '../js/math.js';
import { createRenderer, viewGeometry, orbitalPeak, orbitalGains, MAX_VIEW_RADIUS, shapePeak, splitViews } from '../js/renderer.js';
import { loadTables } from '../js/tables.js';
import { atomsOf, paneViews, subshellUnion, occupiedAny, valenceUnion, defaultExpandedShells, collapsedRows, tableRow } from '../js/selection.js';
import { createScaleAnimator } from '../js/view.js';
import { pinchZoom, ZOOM_MIN, ZOOM_MAX, flingSpin, FLING_MIN_SPEED, SPIN_MAX_RATE } from '../js/controls.js';

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

test('view sized to the Rahm radius; fixed scale uses the largest one', () => {
  const [H, He] = ELEMENTS;
  const fH = viewGeometry(H, 400, 300, 1, true), fHe = viewGeometry(He, 400, 300, 1, true);
  near(fH.halfExtent[1], fHe.halfExtent[1], 1e-12);
  near(fH.halfExtent[1], MAX_VIEW_RADIUS / 0.9, 1e-12);
  const largest = Math.max(...ELEMENTS.filter(e => e.supported).map(e => e.atomicRadius));
  near(MAX_VIEW_RADIUS, largest / 52.9177, 1e-12);
  // auto: the ring (Rahm radius) fills 90% of the shorter half-dimension
  const aH = viewGeometry(H, 400, 300, 1, false), aHe = viewGeometry(He, 400, 300, 1, false);
  near(aH.radius, 154 / 52.9177, 1e-12);
  near(aH.halfExtent[1], aH.radius / 0.9, 1e-12);
  near(aHe.halfExtent[1] / aH.halfExtent[1], 134 / 154, 1e-12);
});

// ------------------------------------------------------------ elements
test('element data: sequential Z, occupancy sums to Z', () => {
  assert(ELEMENTS.length >= 18, `only ${ELEMENTS.length} elements`);
  ELEMENTS.forEach((el, i) => assert(el.z === i + 1, `${el.symbol} out of order`));
  for (const el of ELEMENTS) {
    const sum = el.occupancy.reduce((a, b) => a + b, 0);
    assert(sum === el.z, `${el.symbol}: ${sum} != ${el.z}`);
    if (el.supported) el.occupancy.forEach((o, i) => assert(o === 0 || el.zeff[i] > 0, `${el.symbol} missing Z_eff`));
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

test('chroma by depth: valence vivid, inner shells paler; same axis hues', () => {
  const o = id => ORBITALS[orbitalIndex(id)];
  const kr = ELEMENTS[35], ne = ELEMENTS[9];
  assert(orbitalChroma(o('4px'), kr) === 1 && orbitalChroma(o('3px'), kr) === 0.6 &&
    orbitalChroma(o('2px'), kr) === 0.35, 'Kr p shells');
  assert(orbitalChroma(o('2px'), ne) === 1 && orbitalChroma(o('1s'), kr) === 1, 's stays vivid; Ne 2p valence');
  assert(o('2px').hue === o('3px').hue && o('2pz').hue === o('4pz').hue);
});

test('period 4: Cr and Cu exceptions, shorthand', () => {
  assert(ELEMENTS.length >= 36, `${ELEMENTS.length} elements`);
  const cfg = z => ELEMENTS[z - 1].configShort;
  assert(cfg(19) === '[Ar] 4s¹', cfg(19));
  assert(cfg(24) === '[Ar] 3d⁵ 4s¹', cfg(24));
  assert(cfg(26) === '[Ar] 3d⁶ 4s²', cfg(26));
  assert(cfg(29) === '[Ar] 3d¹⁰ 4s¹', cfg(29));
  assert(cfg(36) === '[Ar] 3d¹⁰ 4s² 4p⁶', cfg(36));
  const d = z => orbitalsOf('3d').map(i => ELEMENTS[z - 1].occupancy[i]).join('');
  assert(d(21) === '10000' && d(24) === '11111' && d(26) === '21111' && d(30) === '22222');
  const ids = z => [...new Set(ORBITALS.filter(o => isValence(o, ELEMENTS[z - 1])).map(o => o.subshell))].sort().join();
  assert(ids(26) === '3d,4s' && ids(30) === '4s' && ids(31) === '4p,4s', `${ids(26)} | ${ids(30)} | ${ids(31)}`);
});

test('d peak search is not beaten by random rays', () => {
  const rand = makeRng(77);
  for (const id of ['3dxy', '3dz2', '4fxyz', '4fz3']) {
    const i = orbitalIndex(id);
    let best = 0;
    for (let k = 0; k < 3000; k++) {
      // random direction, random closest point perpendicular to it
      const z = 2 * rand() - 1, ph = 2 * Math.PI * rand(), r = Math.sqrt(1 - z * z);
      const d = [r * Math.cos(ph), r * Math.sin(ph), z];
      let q = [rand() - 0.5, rand() - 0.5, rand() - 0.5];
      const t = q[0] * d[0] + q[1] * d[1] + q[2] * d[2];
      q = q.map((c, j) => c - t * d[j]);
      const n = Math.hypot(...q), u = 6 * rand();
      const P = q.map(c => (c / n) * u);
      best = Math.max(best, columnDensity(TABLES, ORBITALS[i], ORBITALS[i].n, P, d));   // kappa = 1
    }
    const peak = shapePeak(TABLES, i);
    assert(peak >= 0.97 * best, `${id}: search peak ${peak} < random best ${best}`);
    assert(peak <= 1.3 * best, `${id}: search peak ${peak} >> random best ${best}`);
  }
});

test('period 5: 4d anomalies and [Kr] shorthand from the data', () => {
  const cfg = sym => ELEMENTS.find(e => e.symbol === sym).configShort;
  assert(cfg('Rb') === '[Kr] 5s¹', cfg('Rb'));
  assert(cfg('Nb') === '[Kr] 4d⁴ 5s¹', cfg('Nb'));
  assert(cfg('Mo') === '[Kr] 4d⁵ 5s¹', cfg('Mo'));
  assert(cfg('Pd') === '[Kr] 4d¹⁰', cfg('Pd'));
  assert(cfg('Ag') === '[Kr] 4d¹⁰ 5s¹', cfg('Ag'));
  assert(cfg('Xe') === '[Kr] 4d¹⁰ 5s² 5p⁶', cfg('Xe'));
  const pd = ELEMENTS.find(e => e.symbol === 'Pd');
  const ids = [...new Set(ORBITALS.filter(o => isValence(o, pd)).map(o => o.subshell))].join();
  assert(ids === '4d', `Pd valence ${ids}`);
});

test('data: Clementi Z_eff for every occupied subshell, Kr 4p corrected', () => {
  for (const el of ELEMENTS.filter(e => e.supported)) {
    for (const id of Object.keys(el.subshells)) assert(el.zeffBySubshell[id] > 0, `${el.symbol} ${id}`);
  }
  near(ELEMENTS[35].zeffBySubshell['4p'], 9.769, 1e-9);
});

test('stepping skips unsupported elements', () => {
  const firstUnsupported = ELEMENTS.find(e => !e.supported);
  if (firstUnsupported) {
    const z = firstUnsupported.z;
    assert(stepElement(z - 1, +1) !== z, 'stepped onto an unsupported element');
  }
  assert(stepElement(1, -1) === null && stepElement(ELEMENTS.length, +1) === null);
});

test('period 6: lanthanides, Pt/Au, [Xe] shorthand, f-block valence', () => {
  const E = sym => ELEMENTS.find(e => e.symbol === sym);
  assert(ELEMENTS.length === 86, `${ELEMENTS.length} elements`);
  assert(E('Cs').configShort === '[Xe] 6s¹', E('Cs').configShort);
  assert(E('Eu').configShort === '[Xe] 4f⁷ 6s²', E('Eu').configShort);
  assert(E('Pt').configShort === '[Xe] 4f¹⁴ 5d⁹ 6s¹', E('Pt').configShort);
  assert(E('Au').configShort === '[Xe] 4f¹⁴ 5d¹⁰ 6s¹', E('Au').configShort);
  assert(E('Rn').configShort === '[Xe] 4f¹⁴ 5d¹⁰ 6s² 6p⁶', E('Rn').configShort);
  assert(['La', 'Ce', 'Gd'].every(s => !E(s).supported), 'La/Ce/Gd greyed');
  assert(E('Nd').group === null && E('Nd').period === 6);
  const val = sym => [...new Set(ORBITALS.filter(o => isValence(o, E(sym))).map(o => o.subshell))].sort().join();
  assert(val('Nd') === '4f,6s', `Nd ${val('Nd')}`);
  assert(val('Lu') === '5d,6s', `Lu ${val('Lu')}`);
  assert(val('Hg') === '6s', `Hg ${val('Hg')}`);
  // f Hund order: Nd 4f4 fills z3, xz2, yz2, xyz
  const f = orbitalsOf('4f').map(i => E('Nd').occupancy[i]).join('');
  assert(f === '1111000', `Nd 4f ${f}`);
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

test('odd half-line moments at u = 0 match closed forms', () => {
  // J_k(0) = int_0^inf R^2 t^k / t^(2l) dt
  const r1s = TABLES.rows.findIndex(r => r.name === '1s');
  const r2p = TABLES.rows.findIndex(r => r.name === '2p');
  [1, 1.5, 7.5].forEach((v, j) => near(momentLookup(TABLES, r1s, 2 * j + 1, 0), v, 1e-5));
  [1 / 3, 0.5, 2.5].forEach((v, j) => near(momentLookup(TABLES, r2p, 2 * j + 1, 0), v, 1e-5));
});

test('cut-away columns: back(p) = front(-p), back + front = 2 whole, s unchanged', () => {
  const R = viewToAtom(qMul(qAxisAngle([1, 0, 0], 0.7), qAxisAngle([0, 1, 0], -0.9)));
  const d = mat3MulVec(R, [0, 0, 1]);
  let asym = 0;
  ['1s', '2px', '3dxz', '3dz2', '4fxyz', '4fz3', '5py'].forEach(id => {
    const o = ORBITALS[P(id)];
    for (const [x, y] of [[0.4, 0.1], [-0.3, 0.8], [1.1, -0.6]]) {
      const p = mat3MulVec(R, [x, y, 0]), mp = p.map(c => -c);
      const whole = columnDensity(TABLES, o, o.n, p, d);
      const back = columnDensity(TABLES, o, o.n, p, d, -1);
      const front = columnDensity(TABLES, o, o.n, p, d, 1);
      const tol = 1e-5 * Math.max(1, whole);
      near(back, columnDensity(TABLES, o, o.n, mp, d, 1), tol, `${id} back(p) vs front(-p)`);
      near(back + front, 2 * whole, tol, `${id} back + front`);
      if (o.l === 0) near(back, whole, tol, `${id} s orbital unchanged by the cut`);
      asym = Math.max(asym, Math.abs(back - front) / whole);
    }
  });
  assert(asym > 0.1, `cut should change tilted p/d/f columns (max rel. change ${asym})`);
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

test('f shapes: 7 orbitals, packed coefficients rebuild the cubic tensor', () => {
  const fs = ORBITALS.filter(o => o.l === 3);
  assert(fs.length === 7, `${fs.length} f orbitals`);
  for (const o of fs) {
    const c = angularCoefs(o);
    // unique entries in the shader's packing order
    const at = (i, j, k) => {
      const key = [i, j, k].sort().join('');
      return c[['000', '001', '002', '011', '012', '022', '111', '112', '122', '222'].indexOf(key)];
    };
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) for (let k = 0; k < 3; k++) {
      near(at(i, j, k), o.C[i][j][k], 1e-12, `${o.id} C[${i}${j}${k}]`);
    }
  }
});

test('full p, d and f subshells are isotropic (Unsold)', () => {
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

// ------------------------------------------------------------ nearest detection
test('radial table: q2 at u = 0 matches closed forms; extents positive', () => {
  const r1s = TABLES.rows.findIndex(r => r.name === '1s');
  const r2p = TABLES.rows.findIndex(r => r.name === '2p');
  near(radialLookup(TABLES, r1s, 0), 4, 1e-5);
  near(radialLookup(TABLES, r2p, 0), 4 / 3, 1e-5);
  assert(TABLES.rows.every(r => r.extent > 0 && r.extent <= TABLES.uMax), 'extents');
});

// Rays (scaled, reference units) for the march tests.
function marchCases() {
  const R = viewToAtom(qMul(qAxisAngle([1, 0, 0], 0.7), qAxisAngle([0, 1, 0], -0.9)));
  const d = mat3MulVec(R, [0, 0, 1]);
  const out = [];
  ['1s', '4s', '2px', '5py', '3dxz', '3dz2', '5dxy', '4fxyz', '4fz3'].forEach(id => {
    const o = ORBITALS[P(id)];
    for (const [x, y] of [[0.05, 0.02], [0.9, -0.4], [-2.5, 1.8]]) {
      out.push({ id, o, d, P: mat3MulVec(R, [x * o.n, y * o.n, 0]) });   // ~ shell-sized offsets
    }
  });
  return out;
}

test('march: density along the ray integrates to the table column (whole and cut)', () => {
  // Fine trapezoid of q2(r) A(t)^2 over the orbital's extent vs the moment tables.
  for (const { id, o, d, P } of marchCases()) {
    const ext = TABLES.rows[o.row].extent, u2 = P[0] ** 2 + P[1] ** 2 + P[2] ** 2;
    if (u2 >= ext * ext) continue;
    const H = Math.sqrt(ext * ext - u2), A = angularPoly(o, P, d), n = 4000;
    const f = t => { const a = A[0] + t * (A[1] + t * (A[2] + t * A[3]));
      return radialLookup(TABLES, o.row, Math.sqrt(u2 + t * t)) * a * a; };
    let whole = 0, back = 0;
    for (let j = 0; j < n; j++) {
      const t0 = -H + (2 * H * j) / n, t1 = t0 + (2 * H) / n, seg = 0.5 * (f(t0) + f(t1)) * (t1 - t0);
      whole += seg; if (t1 <= 1e-12) back += seg;
    }
    const col = columnRef(TABLES, o, P, d), colCut = columnRef(TABLES, o, P, d, -1) / 2;
    near(whole, col, 2e-3 * col + 1e-7, `${id} whole`);
    near(back, colCut, 2e-3 * col + 1e-7, `${id} back half`);
  }
});

test('march: nearest depth with 48 steps matches a fine march', () => {
  let worst = 0;
  for (const { id, o, d, P } of marchCases()) {
    const ext = TABLES.rows[o.row].extent;
    for (const cut of [false, true]) {
      const mass = cut ? columnRef(TABLES, o, P, d, -1) / 2 : columnRef(TABLES, o, P, d);
      for (const k of [1, 4]) for (const v of [0.1, 0.5, 0.9]) {
        const coarse = nearestDepthRef(TABLES, o, P, d, k, mass, v, cut);
        const fine = nearestDepthRef(TABLES, o, P, d, k, mass, v, cut, NO_DEPTH, 20000);
        worst = Math.max(worst, Math.abs(coarse - fine) / ext);
        assert(!cut || coarse <= 1e-9, `${id}: cut depth must be behind the plane`);
      }
    }
  }
  // Measured worst 0.0115 of the extent (front tail, where the t = H s^3 steps
  // are widest; ~0.2 reference units); the median is far smaller.
  assert(worst < 0.02, `worst depth error ${worst.toFixed(4)} of the extent`);
});

test('march: nearest of k detections follows 1 - (1 - F)^k; stopping early is exact', () => {
  // v = 1 - q^k puts the nearest at front fraction q... check the order statistic
  // through the inverse: depth(k, v) == depth(1, v^(1/k)).
  const { o, d, P } = marchCases()[7];
  const mass = columnRef(TABLES, o, P, d);
  for (const v of [0.2, 0.7]) {
    near(nearestDepthRef(TABLES, o, P, d, 3, mass, v), nearestDepthRef(TABLES, o, P, d, 1, mass, Math.pow(v, 1 / 3)), 1e-9);
    const t = nearestDepthRef(TABLES, o, P, d, 1, mass, v);
    assert(nearestDepthRef(TABLES, o, P, d, 1, mass, v, false, t - 1e-6) === t, 'stop behind: same depth');
    assert(nearestDepthRef(TABLES, o, P, d, 1, mass, v, false, t + 0.5) === NO_DEPTH, 'stop in front: no depth');
  }
});

// ------------------------------------------------------------ GPU
function renderFrame(z, opts = {}) {
  const canvas = document.createElement('canvas');
  canvas.width = 256; canvas.height = 256;
  const r = createRenderer(canvas, TABLES);
  const state = {
    element: ELEMENTS[z - 1], orientation: opts.q || [0, 0, 0, 1], zoom: opts.zoom ?? 1,
    seed: opts.seed ?? 1, density: opts.density ?? 0.3,
    equalize: opts.equalize ?? false, kMax: opts.kMax ?? 4, gamma: opts.gamma ?? 0.5,
    mode: opts.mode || 'dark', visible: opts.visible || ORBITALS.map(() => true), fixedScale: opts.fixedScale ?? false,
    cut: opts.cut ?? false, color: opts.color || 'blend',
  };
  r.render(state);
  const px = new Uint8Array(256 * 256 * 4);
  r.gl.readPixels(0, 0, 256, 256, r.gl.RGBA, r.gl.UNSIGNED_BYTE, px);
  r.gl.getExtension('WEBGL_lose_context')?.loseContext();
  return px;
}

function brightness(px, x, y) { const i = 4 * (y * 256 + x); return px[i] + px[i + 1] + px[i + 2]; }

test('scale animator: snap, log-space ease, retarget without jumps', () => {
  const a = createScaleAnimator({ duration: 400 });
  near(a.value(2, 0), 2, 0, 'first call snaps');
  assert(!a.animating);
  let prev = a.value(16, 1000);                    // new target at t = 1000
  near(prev, 2, 1e-12, 'starts from the shown value');
  assert(a.animating);
  for (let t = 1010; t < 1400; t += 10) {
    const v = a.value(16, t);
    assert(v >= prev && v <= 16, `not monotonic at ${t}: ${v}`);
    prev = v;
  }
  near(a.value(16, 1200), Math.sqrt(2 * 16), 1e-9, 'log-space midpoint is the geometric mean');
  near(a.value(16, 1400), 16, 0, 'reaches the target');
  assert(!a.animating);
  // retarget mid-tween: continues from the current value
  a.value(4, 2000);
  const mid = a.value(4, 2200);
  near(a.value(32, 2200), mid, 1e-12, 'retarget does not jump');
  near(a.value(32, 2600), 32, 0);
  const instant = createScaleAnimator({ duration: 0 });
  instant.value(1, 0);
  near(instant.value(9, 1), 9, 0, 'duration 0 snaps');
});

test('pinch zoom follows the finger-distance ratio and clamps', () => {
  near(pinchZoom(2, 100, 150), 3, 1e-12);
  near(pinchZoom(2, 100, 50), 1, 1e-12);
  near(pinchZoom(50, 100, 1000), ZOOM_MAX, 0);
  near(pinchZoom(1, 100, 1), ZOOM_MIN, 0);
  near(pinchZoom(3, 0, 10), 3, 0, 'degenerate start distance keeps the zoom');
});

test('fling: spin only above the speed threshold, axis like a drag, rate capped', () => {
  assert(flingSpin(0.5 * FLING_MIN_SPEED, 0) === null, 'slow release keeps no spin');
  assert(flingSpin(0, 0) === null);
  const s = flingSpin(0.5, 0);                         // flick right
  near(s.axis[1], 1, 1e-12, 'horizontal flick spins about view +y');
  near(s.rate, 0.5 * Math.PI / 400, 1e-12, 'rate = drag mapping');
  const v = flingSpin(0, -0.6);                        // flick up
  near(v.axis[0], -1, 1e-12, 'vertical flick spins about view x');
  near(flingSpin(50, 0).rate, SPIN_MAX_RATE, 0, 'capped');
});

test('viewGeometry baseOverride scales the view but keeps the atom\'s own radius', () => {
  const ne = ELEMENTS[9];
  const g = viewGeometry(ne, 300, 300, 2, false, 10);
  near(g.halfExtent[1], 5, 1e-12);
  near(g.radius, ne.atomicRadius / 52.9177, 1e-12);
});

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
      const px = renderFrame(z, { seed: s, fixedScale: true, density: 0.2, kMax: 1 });
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

test('GPU: Sc 3d_xy viewed down z is a four-leaf clover', () => {
  const i = orbitalIndex('3dxy');
  const visible = ORBITALS.map((o, j) => j === i);
  const sc = ELEMENTS[20];
  // Zoom so the lobes (~0.9 a0 out) land ~64 px from the centre, whatever
  // the view's base scale.
  const rr = 0.9;
  const zoom = (64 * viewGeometry(sc, 256, 256, 1).bohrPerPixel) / rr;
  const g = viewGeometry(sc, 256, 256, zoom);
  const r = Math.round(rr / g.bohrPerPixel);
  let diag = 0, axes = 0;
  for (let s = 1; s <= 4; s++) {
    const px = renderFrame(21, { seed: s, equalize: true, density: 0.3, visible, zoom });
    const c = 128, k = Math.round(r / Math.SQRT2);
    for (let o = -4; o <= 4; o++) for (let t = -4; t <= 4; t++) {
      for (const [sx, sy] of [[1, 1], [1, -1], [-1, 1], [-1, -1]]) diag += brightness(px, c + sx * k + o, c + sy * k + t);
      for (const [ax, ay] of [[r, 0], [-r, 0], [0, r], [0, -r]]) axes += brightness(px, c + ax + o, c + ay + t);
    }
  }
  assert(diag > 3 * axes, `diagonals ${diag} vs axes ${axes}`);
});

test('GPU: white mode is the exact dual of dark mode', () => {
  // white = 1 - dark(hue + 180deg). For an HSV color, rotating the hue by
  // 180deg maps each channel x -> max + min - x, so per pixel:
  //   white = 255 - (max + min - dark)
  const opts = { seed: 3, equalize: true, density: 0.3 };
  const dark = renderFrame(8, { ...opts, mode: 'dark' });
  const white = renderFrame(8, { ...opts, mode: 'white' });
  let lit = 0, worst = 0, overlapDarker = 0, overlaps = 0;
  for (let i = 0; i < dark.length; i += 4) {
    const d = [dark[i], dark[i + 1], dark[i + 2]];
    const M = Math.max(...d), m = Math.min(...d);
    if (M > 0) lit++;
    for (let c = 0; c < 3; c++) worst = Math.max(worst, Math.abs(white[i + c] - (255 - (M + m - d[c]))));
    // low-saturation dark pixels (hue mixes) must be darker than paper in white
    if (M > 100 && M - m < 0.3 * M) {
      overlaps++;
      if (white[i] + white[i + 1] + white[i + 2] < 3 * 255 - 100) overlapDarker++;
    }
  }
  assert(lit > 2000, `too few lit pixels ${lit}`);
  assert(worst <= 3, `max channel deviation ${worst}`);
  assert(overlaps === 0 || overlapDarker === overlaps, `${overlapDarker}/${overlaps} mixes darker`);
});

test('GPU: Nd 4f_xyz viewed down z is a four-leaf clover', () => {
  const i = orbitalIndex('4fxyz');
  const nd = ELEMENTS.find(e => e.symbol === 'Nd');
  const visible = ORBITALS.map((o, j) => j === i);
  const rr = 0.8 * Math.sqrt(meanR2(4, 3, nd.zeff[i]));          // ~ lobe distance (a0)
  const zoom = (64 * viewGeometry(nd, 256, 256, 1).bohrPerPixel) / rr;   // lobes ~64 px out
  const g = viewGeometry(nd, 256, 256, zoom);
  const r = Math.round(rr / g.bohrPerPixel);
  assert(r > 15 && r < 110, `lobe radius ${r}px`);
  let diag = 0, axes = 0;
  for (let s = 1; s <= 4; s++) {
    const px = renderFrame(nd.z, { seed: s, equalize: true, density: 0.3, visible, zoom });
    const c = 128, k = Math.round(r / Math.SQRT2);
    for (let o = -3; o <= 3; o++) for (let t = -3; t <= 3; t++) {
      for (const [sx, sy] of [[1, 1], [1, -1], [-1, 1], [-1, -1]]) diag += brightness(px, c + sx * k + o, c + sy * k + t);
      for (const [ax, ay] of [[r, 0], [-r, 0], [0, r], [0, -r]]) axes += brightness(px, c + ax + o, c + ay + t);
    }
  }
  assert(diag > 3 * axes, `diagonals ${diag} vs axes ${axes}`);
});

test('selection: unions over A and B', () => {
  const E = sym => ELEMENTS.find(e => e.symbol === sym);
  const fe = E('Fe'), kr = E('Kr'), h = E('H');
  assert(atomsOf({ element: fe, second: null }).length === 1 && atomsOf({ element: fe, second: kr }).length === 2);
  const u = subshellUnion([fe, kr]);
  assert(u.has('3d') && u.has('4p') && !u.has('4d'), [...u].join());
  assert(occupiedAny([fe, kr], orbitalIndex('4pz')) && !occupiedAny([fe], orbitalIndex('4pz')));
  const v = valenceUnion([fe, kr]);
  assert(v.has(orbitalIndex('3dxy')) && v.has(orbitalIndex('4px')) && v.has(orbitalIndex('4s')), 'Fe 3d + Kr 4p + 4s');
  assert(!v.has(orbitalIndex('2px')), 'core 2p is not valence');
  const sh = defaultExpandedShells([h, fe]);
  assert(sh.has(1) && sh.has(3) && sh.has(4), [...sh].join());     // H valence 1s, Fe valence 3d + 4s
});

test('collapsedRows: full rows for A and B, peeks next to them', () => {
  const rows = [1, 2, 3, 4, 5, 6, 8];
  const show = full => collapsedRows(rows, full).map(r => (r.peek ? `(${r.row})` : String(r.row))).join(' ');
  assert(show([4]) === '(3) 4 (5)', show([4]));
  assert(show([4, 4]) === '(3) 4 (5)', 'same period: one row');
  assert(show([2, 5]) === '(1) 2 (3) (4) 5 (6)', show([2, 5]));
  assert(show([1, 8]) === '1 (2) (6) 8', show([1, 8]));
  assert(collapsedRows(rows, [3, 4]).map(r => r.gridRow).join() === '1,2,3,4', 'consecutive grid rows');
  assert(tableRow(ELEMENTS.find(e => e.symbol === 'Nd')) === 8 && tableRow(ELEMENTS[25]) === 4);
});

test('splitViews: panes tile the canvas on multiples of k, orientation by aspect', () => {
  for (const [W, H, k] of [[1000, 600, 1], [1001, 601, 2], [997, 400, 3], [390, 844, 1], [391, 845, 3]]) {
    const v = splitViews(W, H, k, 2);
    const area = v.reduce((a, r) => a + r.w * r.h, 0);
    assert(area === W * H, `area ${W}x${H}`);
    for (const r of v) assert(r.x % k === 0 && r.y % k === 0, `origin not on k grid ${JSON.stringify(r)}`);
    if (W >= H) assert(v[0].x === 0 && v[1].x === v[0].w && v[0].h === H && v[1].h === H, 'left|right');
    else assert(v[1].y === 0 && v[0].y === v[1].h && v[0].w === W, 'top over bottom (primary on top)');
  }
  const one = splitViews(800, 600, 2, 1);
  assert(one.length === 1 && one[0].w === 800 && one[0].h === 600);
});

test('GPU: two atoms side by side are each centred in their own pane', () => {
  const W = 320, H = 160;
  for (const k of [1, 2]) {
    const canvas = document.createElement('canvas');
    canvas.width = W; canvas.height = H;
    const r = createRenderer(canvas, TABLES);
    const rects = splitViews(W, H, k, 2);
    const state = { element: ELEMENTS[0], orientation: [0, 0, 0, 1], zoom: 1, seed: 3, density: 0.6,
      equalize: false, kMax: 1, gamma: 1, mode: 'dark', fixedScale: false,
      baseHalf: Math.max(...[ELEMENTS[0], ELEMENTS[1]].map(e => e.atomicRadius / 52.9177 / 0.9)) };
    r.renderViews(state, k, [{ element: ELEMENTS[0], rect: rects[0] }, { element: ELEMENTS[1], rect: rects[1] }]);
    const px = new Uint8Array(W * H * 4);
    r.gl.readPixels(0, 0, W, H, r.gl.RGBA, r.gl.UNSIGNED_BYTE, px);
    r.gl.getExtension('WEBGL_lose_context')?.loseContext();
    rects.forEach((rc, j) => {
      let sx = 0, sy = 0, sw = 0;
      for (let y = rc.y; y < rc.y + rc.h; y++) for (let x = rc.x; x < rc.x + rc.w; x++) {
        const i = 4 * (y * W + x), v = px[i] + px[i + 1] + px[i + 2];
        sx += v * (x + 0.5); sy += v * (y + 0.5); sw += v;
      }
      assert(sw > 0, `pane ${j} is dark (k=${k})`);
      near(sx / sw, rc.x + rc.w / 2, 3, `pane ${j} centroid x (k=${k})`);
      near(sy / sw, rc.y + rc.h / 2, 3, `pane ${j} centroid y (k=${k})`);
    });
  }
});

test('pane views: both atoms share the orbital visibility', () => {
  const visible = ORBITALS.map((o, i) => i !== P('2px'));
  const views = paneViews({ element: ELEMENTS[5], second: ELEMENTS[7], visible });
  assert(views.length === 2 && views[0].element.z === 6 && views[1].element.z === 8, 'A, B');
  assert(views.every(v => v.visible === visible), 'both panes get state.visible');
  assert(paneViews({ element: ELEMENTS[5], second: null, visible }).length === 1, 'single atom');
});

test('GPU: hiding orbitals hides them in both panes (C | O)', () => {
  const W = 320, H = 160;
  const canvas = document.createElement('canvas');
  canvas.width = W; canvas.height = H;
  const r = createRenderer(canvas, TABLES);
  const rects = splitViews(W, H, 1, 2);
  // Lit pixels per pane with the given orbitals shown.
  const lit = shown => {
    const state = { element: ELEMENTS[5], second: ELEMENTS[7], orientation: [0, 0, 0, 1], zoom: 1, seed: 3,
      density: 2, equalize: false, kMax: 1, gamma: 1, mode: 'dark', fixedScale: false,
      visible: ORBITALS.map((o, i) => shown.includes(i)) };
    r.renderViews(state, 1, paneViews(state).map((v, j) => ({ ...v, rect: rects[j] })));
    const px = new Uint8Array(W * H * 4);
    r.gl.readPixels(0, 0, W, H, r.gl.RGBA, r.gl.UNSIGNED_BYTE, px);
    return rects.map(rc => {
      let n = 0;
      for (let y = rc.y; y < rc.y + rc.h; y++) for (let x = rc.x; x < rc.x + rc.w; x++) {
        const i = 4 * (y * W + x);
        if (px[i] + px[i + 1] + px[i + 2] > 0) n++;
      }
      return n;
    });
  };
  const all = lit(ORBITALS.map((o, i) => i));
  const none = lit([]);
  const only1s = lit([P('1s')]);
  r.gl.getExtension('WEBGL_lose_context')?.loseContext();
  assert(all[0] > 0 && all[1] > 0, `all shown: ${all}`);
  assert(none[0] === 0 && none[1] === 0, `all hidden: A ${none[0]}, B ${none[1]} lit pixels`);
  assert(only1s[0] > 0 && only1s[1] > 0 && only1s[1] < all[1], `1s only: ${only1s} (all ${all})`);
});

test('GPU: cut-away keeps the lobe behind the nucleus (boron 2px tilted away on the right)', () => {
  // Rotating +0.8 rad about view y takes atom +x to view (cos, 0, -sin):
  // right of centre and behind the plane (view z points at the viewer).
  const q = qAxisAngle([0, 1, 0], 0.8);
  const visible = ORBITALS.map((o, i) => i === P('2px'));
  const g = viewGeometry(ELEMENTS[4], 256, 256, 1);
  const r = Math.round(Math.cos(0.8) * 2.0 / ELEMENTS[4].zeff[P('2px')] / g.bohrPerPixel);
  const sides = cut => {
    let right = 0, left = 0;
    for (let s = 1; s <= 4; s++) {
      const px = renderFrame(5, { seed: s, density: 1, q, visible, cut });
      for (let o = -8; o <= 8; o++) for (let t = -8; t <= 8; t++) {
        right += brightness(px, 128 + r + o, 128 + t);
        left += brightness(px, 128 - r + o, 128 + t);
      }
    }
    return [right, left];
  };
  const [r0, l0] = sides(false), [r1, l1] = sides(true);
  assert(r0 < 1.3 * l0 && l0 < 1.3 * r0, `uncut lobes should match: right ${r0} left ${l0}`);
  assert(r1 > 3 * l1, `cut: right (behind) ${r1} should dominate left (in front) ${l1}`);
});

test('GPU: cut-away applies to both panes', () => {
  // Boron tilted as above: with the cut each pane's brightness centroid moves
  // right of its centre (the 2px lobe behind the nucleus); 1s/2s stay centred.
  const W = 320, H = 160;
  const canvas = document.createElement('canvas');
  canvas.width = W; canvas.height = H;
  const rend = createRenderer(canvas, TABLES);
  const rects = splitViews(W, H, 1, 2);
  const state = { element: ELEMENTS[4], orientation: qAxisAngle([0, 1, 0], 0.8), zoom: 1, seed: 5,
    density: 3, equalize: true, kMax: 1, gamma: 1, mode: 'dark', fixedScale: false, cut: true };
  const shift = [0, 0];
  for (let seed = 1; seed <= 3; seed++) {
    state.seed = seed;
    rend.renderViews(state, 1, rects.map(rect => ({ element: ELEMENTS[4], rect })));
    const px = new Uint8Array(W * H * 4);
    rend.gl.readPixels(0, 0, W, H, rend.gl.RGBA, rend.gl.UNSIGNED_BYTE, px);
    rects.forEach((rc, j) => {
      let sx = 0, sw = 0;
      for (let y = rc.y; y < rc.y + rc.h; y++) for (let x = rc.x; x < rc.x + rc.w; x++) {
        const i = 4 * (y * W + x), v = px[i] + px[i + 1] + px[i + 2];
        sx += v * (x + 0.5); sw += v;
      }
      shift[j] += (sx / sw - (rc.x + rc.w / 2)) / 3;
    });
  }
  rend.gl.getExtension('WEBGL_lose_context')?.loseContext();
  shift.forEach((s, j) => assert(s > 2, `pane ${j} centroid shift ${s.toFixed(2)} px should be > 2`));
});

test('GPU: nearest color shows the front lobe (neon 2px vs 2pz, tilted)', () => {
  // Rotating +0.8 rad about view y: atom +x goes right and behind, atom +z
  // right and in front, so on the right 2pz is nearer and on the left 2px.
  // Their projections overlap, so blend mode mixes both on each side.
  const q = qAxisAngle([0, 1, 0], 0.8);
  const ix = P('2px'), iz = P('2pz');
  const visible = ORBITALS.map((o, i) => i === ix || i === iz);
  const hx = ORBITALS[ix].hue / 360, hz = ORBITALS[iz].hue / 360;
  const hueOf = (r, g, b) => {
    const mx = Math.max(r, g, b), mn = Math.min(r, g, b);
    if (mx - mn < 40) return null;                       // grey / mixed
    let h = mx === r ? (g - b) / (mx - mn) : mx === g ? 2 + (b - r) / (mx - mn) : 4 + (r - g) / (mx - mn);
    return ((h / 6) % 1 + 1) % 1;
  };
  const hd = (a, b) => Math.min(Math.abs(a - b), 1 - Math.abs(a - b));
  const share = color => {     // fraction of 2pz-colored pixels among classified ones, [left, right]
    const n = [[0, 0], [0, 0]];
    for (let s = 1; s <= 3; s++) {
      const px = renderFrame(10, { seed: s, density: 2, q, visible, color, kMax: 1, gamma: 1, equalize: true });
      for (let y = 64; y < 192; y++) for (let x = 0; x < 256; x++) {
        if (Math.abs(x - 128) < 12) continue;
        const i = 4 * (y * 256 + x), h = hueOf(px[i], px[i + 1], px[i + 2]);
        if (h === null) continue;
        const side = x < 128 ? 0 : 1;
        if (hd(h, hz) < 0.06) n[side][1]++; else if (hd(h, hx) < 0.06) n[side][0]++;
      }
    }
    return n.map(([a, b]) => b / Math.max(1, a + b));
  };
  // Not 0 / 1: only orbitals detected in a pixel compete, and in the faint
  // edges often just one is, whichever side it is on. Measured at this
  // exposure: left 0.27, right 0.74 (0.21 / 0.80 at density 30).
  const [l1, r1] = share('nearest');
  assert(r1 > 0.65 && l1 < 0.35 && r1 - l1 > 0.35,
    `nearest: 2pz share left ${l1.toFixed(2)} right ${r1.toFixed(2)}`);
  const [l0, r0] = share('blend');
  assert(Math.abs(l0 - r0) < 0.2, `blend should be symmetric: left ${l0.toFixed(2)} right ${r0.toFixed(2)}`);
});

test('GPU: white mode background is white', () => {
  const px = renderFrame(1, { mode: 'white' });
  assert(brightness(px, 2, 2) === 765, 'corner not white');
});

// ------------------------------------------------------------ print page
async function testAsync(name, fn) {
  try {
    await fn();
    results.push({ name, ok: true });
  } catch (e) {
    results.push({ name, ok: false, msg: e.message });
  }
}

// Load a page in a hidden iframe (letter-sized so its layout is the real
// one) and wait until it says it is ready.
async function loadFrame(url, timeout = 60000) {
  const f = document.createElement('iframe');
  f.style.cssText = 'position:absolute;left:-10000px;width:900px;height:1200px';
  f.src = url;
  document.body.append(f);
  const t0 = performance.now();
  while (f.contentDocument?.body?.dataset.ready !== '1') {
    if (performance.now() - t0 > timeout) throw new Error(`${url} did not finish loading`);
    await new Promise(r => setTimeout(r, 50));
  }
  return f;
}

await testAsync('print page: boxed 300 dpi image on white, details for both atoms, hidden orbitals listed', async () => {
  // C (A) with 2px (orbital index P('2px')) hidden, next to O (B).
  const f = await loadFrame(`../print.html?z=6&z2=8&hide=${P('2px')}&seed=7&paper=letter&grain=2`);
  const doc = f.contentDocument;
  const img = doc.getElementById('sim-img');
  // Letter: 7.5 in content width at 300 px/in, 4:3, multiples of 4.
  assert(img.naturalWidth === 2248 && img.naturalHeight === 1684,
    `image ${img.naturalWidth}x${img.naturalHeight}`);
  const box = doc.getElementById('sim-box').getBoundingClientRect();
  const ir = img.getBoundingClientRect();
  near(ir.left, box.left, 0.5, 'box left'); near(ir.right, box.right, 0.5, 'box right');
  near(ir.top, box.top, 0.5, 'box top'); near(ir.bottom, box.bottom, 0.5, 'box bottom');
  assert(f.contentWindow.getComputedStyle(doc.getElementById('sim-box')).outlineStyle === 'solid', 'box outline');
  // White background: the image's corner pixel.
  const c = document.createElement('canvas');
  c.width = c.height = 1;
  c.getContext('2d').drawImage(img, 0, img.naturalHeight - 1, 1, 1, 0, 0, 1, 1);
  const px = c.getContext('2d').getImageData(0, 0, 1, 1).data;
  assert(px[0] === 255 && px[1] === 255 && px[2] === 255, `corner ${[...px]}`);
  // Notes above and below the image on page 1.
  const notes = [...doc.querySelectorAll('.page1 .notes')].map(n => n.getBoundingClientRect());
  assert(notes.length === 2 && notes[0].bottom <= box.top && notes[1].top >= box.bottom, 'notes above and below');
  assert(notes.every(n => n.height >= 0.9 * 96), 'notes at least 0.9 in tall');
  // Details: both atoms, the occupied subshells, and the hidden 2px.
  assert(doc.querySelector('.title').textContent.includes('Carbon') &&
    doc.querySelector('.title').textContent.includes('Oxygen'), 'title');
  const subs = [...doc.querySelectorAll('#legend tr.sub .orb-label')].map(td => td.firstChild.textContent);
  assert(subs.join() === '1s,2s,2p', `subshells ${subs}`);
  const hiddenNote = doc.getElementById('hidden-note');
  assert(!hiddenNote.hidden && hiddenNote.textContent.includes('2px'), `hidden note: ${hiddenNote.textContent}`);
  assert(doc.getElementById('radial').width > 0 && doc.querySelectorAll('#radial-legend .rl-item').length === 2,
    'radial A vs B');
  f.remove();
});

await testAsync('print page: a single atom gets a centred box the size of one pane of a pair', async () => {
  const f = await loadFrame('../print.html?z=26&seed=7&paper=letter&grain=2');
  const doc = f.contentDocument;
  const img = doc.getElementById('sim-img');
  // Letter pair 2248 x 1684 -> one pane 1124 x 1684.
  assert(img.naturalWidth === 1124 && img.naturalHeight === 1684, `image ${img.naturalWidth}x${img.naturalHeight}`);
  const box = doc.getElementById('sim-box').getBoundingClientRect();
  const page = doc.querySelector('.page1 .head').getBoundingClientRect();
  near(box.width, page.width / 2, 1, 'half the content width');
  near(box.left + box.width / 2, page.left + page.width / 2, 1, 'centred');
  const ir = img.getBoundingClientRect();
  near(ir.left, box.left, 0.5, 'box left'); near(ir.right, box.right, 0.5, 'box right');
  near(ir.top, box.top, 0.5, 'box top'); near(ir.bottom, box.bottom, 0.5, 'box bottom');
  f.remove();
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
// Also on the console (headless: --enable-logging=stderr).
console.log(`TESTS ${summary}` + results.filter(r => !r.ok).map(r => ` | FAIL ${r.name}: ${r.msg}`).join(''));
