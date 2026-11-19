// Draws one frame of the probabilistic orbital image.

import { createContext, createProgram, createTableTexture, createEmptyVAO } from './gl.js';
import { VERT, FRAG } from './shaders.js';
import { ORBITALS, NUM_ORBITALS } from './orbitals.js';
import { ELEMENTS } from './elements.js';
import { horizon, viewToAtom, mat3ToGL, momentLookup, columnRef } from './math.js';

// ------------------------------------------------------------ peaks
// Peak column density per electron of each orbital *shape* in reference
// units (Z0 = n), maximised over all rays. Cached per tables object.
//  s: max_u I0 / 4pi
//  p: (e.P)^2 + (e.d)^2 <= |e|^2 when |P| = u, so max = |e|^2 max_u max(u^2 I0, I2)
//  d: numerical search over directions and image-plane points.
const PEAK_CACHE = new WeakMap();

function searchPeakD(tables, orbital) {
  let best = 0;
  const nDir = 120, nR = 48, nPhi = 24, uHi = 12;
  for (let k = 0; k < nDir; k++) {
    // Fibonacci hemisphere (column density is symmetric under d -> -d)
    const z = 1 - (k + 0.5) / nDir;
    const r = Math.sqrt(1 - z * z), phi = k * 2.399963229728653;
    const d = [r * Math.cos(phi), r * Math.sin(phi), z];
    const a = Math.abs(d[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
    let e1 = [a[1] * d[2] - a[2] * d[1], a[2] * d[0] - a[0] * d[2], a[0] * d[1] - a[1] * d[0]];
    const n1 = Math.hypot(...e1);
    e1 = e1.map(c => c / n1);
    const e2 = [d[1] * e1[2] - d[2] * e1[1], d[2] * e1[0] - d[0] * e1[2], d[0] * e1[1] - d[1] * e1[0]];
    for (let i = 0; i < nR; i++) {
      const u = (i / (nR - 1)) * uHi;
      for (let j = 0; j < nPhi; j++) {
        const t = (2 * Math.PI * j) / nPhi;
        const P = [0, 1, 2].map(c => u * (Math.cos(t) * e1[c] + Math.sin(t) * e2[c]));
        best = Math.max(best, columnRef(tables, orbital, P, d));
      }
    }
  }
  return best;
}

export function shapePeaks(tables) {
  if (PEAK_CACHE.has(tables)) return PEAK_CACHE.get(tables);
  const du = tables.uMax / (tables.N - 1);
  const peaks = ORBITALS.map(o => {
    if (o.l === 2) return searchPeakD(tables, o);
    let best = 0;
    for (let i = 0; i < tables.N; i++) {
      const u = i * du;
      const I0 = momentLookup(tables, o.row, 0, u);
      if (o.l === 0) best = Math.max(best, I0 / (4 * Math.PI));
      else best = Math.max(best, u * u * I0, momentLookup(tables, o.row, 2, u));
    }
    return o.l === 1 ? best * (3 / (4 * Math.PI)) : best;
  });
  PEAK_CACHE.set(tables, peaks);
  return peaks;
}

// Peak column density (electrons / Bohr^2) of orbital i of `element`.
export function orbitalPeak(tables, element, i) {
  const k = element.zeff[i] / ORBITALS[i].n;
  return element.occupancy[i] * k * k * shapePeaks(tables)[i];
}

// In equalized mode the exposure value `density` maps to this peak count.
export const EQUALIZE_PEAK_PER_DENSITY = 20;

// Per-orbital gain: lambda_i = occ * col_i * gain_i.
//  - physical: one shared gain; sum over pixels of lambda ~= density * pixels
//  - equalized: each orbital's peak lambda = density * EQUALIZE_PEAK_PER_DENSITY
export function orbitalGains(tables, element, density, viewArea, equalize) {
  return ORBITALS.map((o, i) => {
    if (!equalize) return (density * viewArea) / element.z;
    const peak = orbitalPeak(tables, element, i);
    return peak > 0 ? (density * EQUALIZE_PEAK_PER_DENSITY) / peak : 0;
  });
}

// ------------------------------------------------------------ view
// Fraction of the shorter screen half-dimension covered by the horizon.
export const HORIZON_FILL = 0.9;

// Largest horizon over all supported elements. Fixed-scale mode uses it for
// every element so that sizes can be compared directly.
export const MAX_HORIZON = Math.max(...ELEMENTS.map(el => horizon(el, ORBITALS)));

// View geometry (Bohr) for a canvas of the given pixel size.
export function viewGeometry(element, width, height, zoom = 1, fixedScale = false) {
  const h = horizon(element, ORBITALS);
  const half = (fixedScale ? MAX_HORIZON : h) / HORIZON_FILL / zoom;  // shorter half-dimension
  const aspect = width / height;
  const halfExtent = aspect >= 1 ? [half * aspect, half] : [half, half / aspect];
  return { horizon: h, halfExtent, bohrPerPixel: (2 * halfExtent[1]) / height };
}

// ------------------------------------------------------------ renderer
export function createRenderer(canvas, tables) {
  const gl = createContext(canvas);
  const { prog, uniforms: u } = createProgram(gl, VERT, FRAG);
  const tex = createTableTexture(gl, tables.data, tables.N, tables.rows.length);
  createEmptyVAO(gl);

  gl.useProgram(prog);
  gl.activeTexture(gl.TEXTURE0);
  gl.bindTexture(gl.TEXTURE_2D, tex);
  gl.uniform1i(u.uTables, 0);
  gl.uniform1i(u.uTableN, tables.N);
  gl.uniform1f(u.uTableUMax, tables.uMax);

  // Static per-orbital shape uniforms.
  gl.uniform2iv(u.uKind, new Int32Array(ORBITALS.flatMap(o => [o.l, o.row])));
  gl.uniform3fv(u.uAxis, new Float32Array(ORBITALS.flatMap(o => o.e || [0, 0, 0])));
  gl.uniformMatrix3fv(u.uQ, false, new Float32Array(ORBITALS.flatMap(o => {
    const Q = o.Q || [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    return Array.from(mat3ToGL(Q));
  })));

  const orbBuf = new Float32Array(4 * NUM_ORBITALS);

  /**
   * state: {
   *   element, orientation (quaternion), zoom, seed (uint32),
   *   density   - exposure: mean detections per pixel if the whole atom
   *               were in view (or peak counts / 20 when equalized),
   *   equalize  - give every orbital the same peak expected count,
   *   kMax      - detections giving full brightness,
   *   gamma     - tone curve exponent: amount = (k / kMax)^gamma,
   *   fixedScale - same view scale for every element (compare sizes),
   *   mode      - 'dark' | 'white',
   *   visible   - boolean per orbital
   * }
   */
  function render(state) {
    const w = canvas.width, h = canvas.height;
    gl.viewport(0, 0, w, h);
    const el = state.element;
    const geo = viewGeometry(el, w, h, state.zoom, state.fixedScale);
    const viewArea = 4 * geo.halfExtent[0] * geo.halfExtent[1];
    const gains = orbitalGains(tables, el, state.density, viewArea, state.equalize);

    ORBITALS.forEach((o, i) => {
      orbBuf[4 * i] = state.visible[i] ? el.occupancy[i] : 0;
      orbBuf[4 * i + 1] = el.zeff[i] / o.n;
      orbBuf[4 * i + 2] = o.hue / 360;
      orbBuf[4 * i + 3] = gains[i];
    });

    gl.uniform2f(u.uResolution, w, h);
    gl.uniform2f(u.uHalfExtent, geo.halfExtent[0], geo.halfExtent[1]);
    gl.uniformMatrix3fv(u.uRot, false, mat3ToGL(viewToAtom(state.orientation)));
    gl.uniform1ui(u.uSeed, state.seed >>> 0);
    gl.uniform4fv(u.uOrb, orbBuf);
    gl.uniform1f(u.uKMax, state.kMax);
    gl.uniform1f(u.uGamma, state.gamma);
    gl.uniform1i(u.uMode, state.mode === 'white' ? 1 : 0);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
    geo.gains = gains;
    return geo;
  }

  return { gl, render };
}
