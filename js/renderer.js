// Draws one frame of the probabilistic orbital image.

import { createContext, createProgram, createTableTexture, createEmptyVAO } from './gl.js';
import { VERT, FRAG } from './shaders.js';
import { TABLE_N, TABLE_B_MAX, TABLE_DATA } from './orbital_tables.js';
import { ORBITALS, ELEMENTS } from './elements.js';
import { horizon, viewToAtom, mat3ToGL, tableLookup } from './math.js';

export const TABLES = { N: TABLE_N, bMax: TABLE_B_MAX, data: TABLE_DATA };

// Peak column density per electron in scaled units (Z = 1), over all rays.
// For p orbitals col = p_a^2 A(b) + d_a^2 B(b) with p_a^2/b^2 + d_a^2 <= 1,
// so the maximum over orientations is max_b max(b^2 A(b), B(b)).
export const SCALED_PEAK = (() => {
  const peak = { s1: 0, s2: 0, p: 0 };
  for (let i = 0; i < TABLE_N; i++) {
    const b = (i / (TABLE_N - 1)) * TABLE_B_MAX;
    peak.s1 = Math.max(peak.s1, tableLookup(TABLES, 2, b));
    peak.s2 = Math.max(peak.s2, tableLookup(TABLES, 3, b));
    peak.p = Math.max(peak.p, b * b * tableLookup(TABLES, 0, b), tableLookup(TABLES, 1, b));
  }
  return peak;
})();

// Peak column density (electrons / Bohr^2) of orbital i of `element`.
export function orbitalPeak(element, i) {
  const z = element.zeff[i];
  return element.occupancy[i] * z * z * SCALED_PEAK[ORBITALS[i].kind];
}

// In equalized mode the exposure value `density` maps to this peak count.
export const EQUALIZE_PEAK_PER_DENSITY = 20;

// Per-orbital gain: lambda_i = occ * Z^2 * col_scaled * gain_i.
//  - physical: one shared gain; sum over pixels of lambda ~= density * pixels
//  - equalized: each orbital's peak lambda = density * EQUALIZE_PEAK_PER_DENSITY
export function orbitalGains(element, density, viewArea, equalize) {
  return ORBITALS.map((o, i) => {
    if (!equalize) return (density * viewArea) / element.z;
    const peak = orbitalPeak(element, i);
    return peak > 0 ? (density * EQUALIZE_PEAK_PER_DENSITY) / peak : 0;
  });
}

// Fraction of the shorter screen half-dimension covered by the horizon.
export const HORIZON_FILL = 0.9;

// Largest horizon over all supported elements (Li). Fixed-scale mode uses it
// for every element so that sizes can be compared directly.
export const MAX_HORIZON = Math.max(...ELEMENTS.map(el => horizon(el, ORBITALS)));

// View geometry (Bohr) for a canvas of the given pixel size.
export function viewGeometry(element, width, height, zoom = 1, fixedScale = false) {
  const h = horizon(element, ORBITALS);
  const half = (fixedScale ? MAX_HORIZON : h) / HORIZON_FILL / zoom;  // shorter half-dimension
  const aspect = width / height;
  const halfExtent = aspect >= 1 ? [half * aspect, half] : [half, half / aspect];
  return { horizon: h, halfExtent, bohrPerPixel: (2 * halfExtent[1]) / height };
}

export function createRenderer(canvas) {
  const gl = createContext(canvas);
  const { prog, uniforms: u } = createProgram(gl, VERT, FRAG);
  const tex = createTableTexture(gl, TABLE_DATA, TABLE_N);
  createEmptyVAO(gl);

  gl.useProgram(prog);
  gl.activeTexture(gl.TEXTURE0);
  gl.bindTexture(gl.TEXTURE_2D, tex);
  gl.uniform1i(u.uTables, 0);
  gl.uniform1i(u.uTableN, TABLE_N);
  gl.uniform1f(u.uTableBMax, TABLE_B_MAX);
  gl.uniform1fv(u.uHue, ORBITALS.map(o => o.hue / 360));

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
    const gains = orbitalGains(el, state.density, viewArea, state.equalize);

    gl.uniform2f(u.uResolution, w, h);
    gl.uniform2f(u.uHalfExtent, geo.halfExtent[0], geo.halfExtent[1]);
    gl.uniformMatrix3fv(u.uRot, false, mat3ToGL(viewToAtom(state.orientation)));
    gl.uniform1ui(u.uSeed, state.seed >>> 0);
    gl.uniform1fv(u.uGain, gains);
    gl.uniform1f(u.uKMax, state.kMax);
    gl.uniform1f(u.uGamma, state.gamma);
    gl.uniform1i(u.uMode, state.mode === 'white' ? 1 : 0);
    gl.uniform1fv(u.uOcc, el.occupancy.map((o, i) => (state.visible[i] ? o : 0)));
    gl.uniform1fv(u.uZeff, el.zeff);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
    geo.gains = gains;
    return geo;
  }

  return { gl, render };
}
