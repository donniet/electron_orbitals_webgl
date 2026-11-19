// Draws one frame of the probabilistic orbital image.

import { createContext, createProgram, createTableTexture, createEmptyVAO } from './gl.js';
import { VERT, FRAG } from './shaders.js';
import { ORBITALS, NUM_ORBITALS, angularCoefs } from './orbitals.js';
import { ELEMENTS } from './elements.js';
import { horizon, viewToAtom, mat3ToGL, momentLookup, columnRef } from './math.js';

// ------------------------------------------------------------ peaks
// Peak column density per electron of each orbital *shape* in reference
// units (Z0 = n), maximised over all rays. Cached per tables object.
//  s: max_u I0 / 4pi
//  p: (e.P)^2 + (e.d)^2 <= |e|^2 when |P| = u, so max = |e|^2 max_u max(u^2 I0, I2)
//  d, f: numerical search over directions and image-plane points.
const PEAK_CACHE = new WeakMap();

function searchPeak(tables, orbital) {
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

function computePeak(tables, o) {
  const du = tables.uMax / (tables.N - 1);
  if (o.l >= 2) return searchPeak(tables, o);
  {
    let best = 0;
    for (let i = 0; i < tables.N; i++) {
      const u = i * du;
      const I0 = momentLookup(tables, o.row, 0, u);
      if (o.l === 0) best = Math.max(best, I0 / (4 * Math.PI));
      else best = Math.max(best, u * u * I0, momentLookup(tables, o.row, 2, u));
    }
    return o.l === 1 ? best * (3 / (4 * Math.PI)) : best;
  }
}

// Reference-unit peak of orbital shape i, computed lazily (the d search
// takes a moment, so atoms without d electrons never pay for it).
export function shapePeak(tables, i) {
  if (!PEAK_CACHE.has(tables)) PEAK_CACHE.set(tables, new Array(ORBITALS.length).fill(null));
  const peaks = PEAK_CACHE.get(tables);
  if (peaks[i] === null) peaks[i] = computePeak(tables, ORBITALS[i]);
  return peaks[i];
}

// Peak column density (electrons / Bohr^2) of orbital i of `element`.
export function orbitalPeak(tables, element, i) {
  if (!element.occupancy[i]) return 0;
  const k = element.zeff[i] / ORBITALS[i].n;
  return element.occupancy[i] * k * k * shapePeak(tables, i);
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
export const MAX_HORIZON = Math.max(...ELEMENTS.filter(el => el.supported).map(el => horizon(el, ORBITALS)));

// Target half-extent (Bohr, shorter screen half-dimension) before zoom.
export function baseHalf(element, fixedScale = false) {
  return (fixedScale ? MAX_HORIZON : horizon(element, ORBITALS)) / HORIZON_FILL;
}

// View geometry (Bohr) for a canvas of the given pixel size. `baseOverride`
// replaces baseHalf() (used to animate scale changes); `horizon` is always
// the element's own, so the ring stays at its true physical size.
export function viewGeometry(element, width, height, zoom = 1, fixedScale = false, baseOverride = null) {
  const h = horizon(element, ORBITALS);
  const half = (baseOverride ?? baseHalf(element, fixedScale)) / zoom;  // shorter half-dimension
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

  // Per-orbital data lives in a std140 uniform block: 5 vec4 per orbital
  // (orb, kind, c0, c1, c2). Static shape data is prepared once; each frame
  // only the occupied, visible orbitals are packed (in canonical order) into
  // the first uCount slots.
  const STRIDE = 20;                                   // floats per Orbital
  const blockIndex = gl.getUniformBlockIndex(prog, 'Orbitals');
  gl.uniformBlockBinding(prog, blockIndex, 0);
  const ubo = gl.createBuffer();
  gl.bindBuffer(gl.UNIFORM_BUFFER, ubo);
  gl.bufferData(gl.UNIFORM_BUFFER, STRIDE * 4 * NUM_ORBITALS, gl.DYNAMIC_DRAW);
  gl.bindBufferBase(gl.UNIFORM_BUFFER, 0, ubo);
  const blockData = new Float32Array(STRIDE * NUM_ORBITALS);
  const shapeData = ORBITALS.map(o => [o.l, o.row, ...angularCoefs(o)]);   // 14 floats

  /**
   * state: {
   *   element, orientation (quaternion), zoom, seed (uint32),
   *   density   - exposure: mean detections per pixel if the whole atom
   *               were in view (or peak counts / 20 when equalized),
   *   equalize  - give every orbital the same peak expected count,
   *   kMax      - detections giving full brightness,
   *   gamma     - tone curve exponent: amount = (k / kMax)^gamma,
   *   fixedScale - same view scale for every element (compare sizes),
   *   baseHalf  - optional animated base half-extent (see view.js),
   *   mode      - 'dark' | 'white',
   *   visible   - boolean per orbital
   * }
   */
  function render(state) {
    const w = canvas.width, h = canvas.height;
    gl.viewport(0, 0, w, h);
    const el = state.element;
    const geo = viewGeometry(el, w, h, state.zoom, state.fixedScale, state.baseHalf ?? null);
    const viewArea = 4 * geo.halfExtent[0] * geo.halfExtent[1];
    const gains = orbitalGains(tables, el, state.density, viewArea, state.equalize);

    let count = 0;
    ORBITALS.forEach((o, i) => {
      if (!(el.occupancy[i] > 0) || !state.visible[i]) return;
      const b = STRIDE * count++;
      const [l, row, ...coefs] = shapeData[i];
      blockData[b] = el.occupancy[i];
      blockData[b + 1] = el.zeff[i] / o.n;
      blockData[b + 2] = o.hue / 360;
      blockData[b + 3] = gains[i];
      blockData[b + 4] = l;
      blockData[b + 5] = row;
      blockData[b + 6] = o.chroma;
      blockData[b + 7] = 0;
      blockData.set(coefs, b + 8);
    });
    gl.bindBuffer(gl.UNIFORM_BUFFER, ubo);
    gl.bufferSubData(gl.UNIFORM_BUFFER, 0, blockData, 0, Math.max(1, count) * STRIDE);

    gl.uniform2f(u.uResolution, w, h);
    gl.uniform2f(u.uHalfExtent, geo.halfExtent[0], geo.halfExtent[1]);
    gl.uniformMatrix3fv(u.uRot, false, mat3ToGL(viewToAtom(state.orientation)));
    gl.uniform1ui(u.uSeed, state.seed >>> 0);
    gl.uniform1i(u.uCount, count);
    gl.uniform1f(u.uKMax, state.kMax);
    gl.uniform1f(u.uGamma, state.gamma);
    gl.uniform1i(u.uMode, state.mode === 'white' ? 1 : 0);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
    geo.gains = gains;
    return geo;
  }

  return { gl, render };
}
