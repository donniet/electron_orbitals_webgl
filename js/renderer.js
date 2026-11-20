// Draws one frame of the probabilistic orbital image.

import { createContext, createProgram, createTableTexture, createEmptyVAO } from './gl.js';
import { VERT, FRAG, SHADER_DEFINES } from './shaders.js';
import { ORBITALS, NUM_ORBITALS, angularCoefs, orbitalChroma } from './orbitals.js';
import { ELEMENTS } from './elements.js';
import { viewToAtom, mat3ToGL, momentLookup, columnRef } from './math.js';

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
const BOHR_PM = 52.9177;

// The view is sized to the element's Rahm atomic radius (free-atom density
// 0.001 e/bohr^3, in elements.json): that ring fills VIEW_FILL of the
// shorter half-dimension. (The model's own cloud can extend beyond it,
// e.g. the single-Z_eff 4s of the 3d metals.)
export const VIEW_FILL = 0.9;

// Rahm atomic radius in Bohr.
export function viewRadius(element) {
  return element.atomicRadius / BOHR_PM;
}

// Largest atomic radius over all supported elements. Fixed-scale mode uses
// it for every element so that sizes can be compared directly.
export const MAX_VIEW_RADIUS = Math.max(...ELEMENTS.filter(el => el.supported).map(viewRadius));

// Target half-extent (Bohr, shorter screen half-dimension) before zoom.
export function baseHalf(element, fixedScale = false) {
  return (fixedScale ? MAX_VIEW_RADIUS : viewRadius(element)) / VIEW_FILL;
}

// View geometry (Bohr) for a canvas of the given pixel size. `baseOverride`
// replaces baseHalf() (used to animate scale changes and to share one scale
// between two panes); `radius` is always the element's own atomic radius.
export function viewGeometry(element, width, height, zoom = 1, fixedScale = false, baseOverride = null) {
  const half = (baseOverride ?? baseHalf(element, fixedScale)) / zoom;  // shorter half-dimension
  const aspect = width / height;
  const halfExtent = aspect >= 1 ? [half * aspect, half] : [half, half / aspect];
  return { radius: viewRadius(element), halfExtent, bohrPerPixel: (2 * halfExtent[1]) / height };
}

// Split a W x H canvas (device px) into n panes for side-by-side atoms:
// left | right when wider than tall, otherwise top over bottom. Rects use
// the GL origin (bottom-left) and every split falls on a multiple of the
// render scale k, so each pane maps exactly onto k x k blocks. Pane 0 is
// the primary atom (left / top).
export function splitViews(W, H, k = 1, n = 1) {
  if (n < 2) return [{ x: 0, y: 0, w: W, h: H }];
  if (W >= H) {
    const wl = Math.max(k, Math.floor(W / (2 * k)) * k);
    return [{ x: 0, y: 0, w: wl, h: H }, { x: wl, y: 0, w: W - wl, h: H }];
  }
  const hb = Math.max(k, Math.floor(H / (2 * k)) * k);       // bottom pane height
  return [{ x: 0, y: hb, w: W, h: H - hb }, { x: 0, y: 0, w: W, h: hb }];
}

// ------------------------------------------------------------ renderer
export function createRenderer(canvas, tables) {
  const gl = createContext(canvas);
  const { prog, uniforms: u } = createProgram(gl, VERT, FRAG, { defines: SHADER_DEFINES });
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
  // Low-resolution target for render scale k > 1. The result is copied to
  // the canvas with an exact integer nearest-neighbour blit (k x k blocks),
  // independent of CSS / browser scaling.
  let fbo = null, rbo = null, fboW = 0, fboH = 0;
  function lowResTarget(bw, bh) {
    if (!fbo) { fbo = gl.createFramebuffer(); rbo = gl.createRenderbuffer(); }
    if (bw !== fboW || bh !== fboH) {
      gl.bindRenderbuffer(gl.RENDERBUFFER, rbo);
      gl.renderbufferStorage(gl.RENDERBUFFER, gl.RGBA8, bw, bh);
      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
      gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.RENDERBUFFER, rbo);
      fboW = bw; fboH = bh;
    }
    return fbo;
  }

  const ALL_VISIBLE = ORBITALS.map(() => true);

  // Pack the occupied, visible orbitals of `el` into the uniform block.
  function uploadOrbitals(el, visible, gains) {
    let count = 0;
    ORBITALS.forEach((o, i) => {
      if (!(el.occupancy[i] > 0) || !visible[i]) return;
      const b = STRIDE * count++;
      const [l, row, ...coefs] = shapeData[i];
      blockData[b] = el.occupancy[i];
      blockData[b + 1] = el.zeff[i] / o.n;
      blockData[b + 2] = o.hue / 360;
      blockData[b + 3] = gains[i];
      blockData[b + 4] = l;
      blockData[b + 5] = row;
      blockData[b + 6] = orbitalChroma(o, el);
      blockData[b + 7] = 0;
      blockData.set(coefs, b + 8);
    });
    gl.bindBuffer(gl.UNIFORM_BUFFER, ubo);
    gl.bufferSubData(gl.UNIFORM_BUFFER, 0, blockData, 0, Math.max(1, count) * STRIDE);
    return count;
  }

  /**
   * Draw one or more panes into the canvas. Each view is
   *   { element, visible (per orbital, optional), rect: { x, y, w, h } }
   * with rect in device pixels, GL origin (bottom-left); see splitViews().
   * All panes share orientation, seed and tone settings. At render scale
   * k > 1 every pane draws into the one 1/k buffer at rect / k (rect origins
   * are multiples of k), followed by a single exact blit.
   * Returns one geometry per view (radius, halfExtent, gains, rect).
   */
  function renderViews(state, k, views) {
    const W = canvas.width, H = canvas.height;
    const BW = Math.ceil(W / k), BH = Math.ceil(H / k);
    gl.bindFramebuffer(gl.FRAMEBUFFER, k > 1 ? lowResTarget(BW, BH) : null);

    gl.uniformMatrix3fv(u.uRot, false, mat3ToGL(viewToAtom(state.orientation)));
    gl.uniform1ui(u.uSeed, state.seed >>> 0);
    gl.uniform1f(u.uKMax, state.kMax);
    gl.uniform1f(u.uGamma, state.gamma);
    gl.uniform1i(u.uMode, state.mode === 'white' ? 1 : 0);

    const geos = views.map(({ element: el, visible = ALL_VISIBLE, rect }) => {
      const bw = Math.ceil(rect.w / k), bh = Math.ceil(rect.h / k);
      gl.viewport(rect.x / k, rect.y / k, bw, bh);
      const geo = viewGeometry(el, rect.w, rect.h, state.zoom, state.fixedScale, state.baseHalf ?? null);
      const viewArea = 4 * geo.halfExtent[0] * geo.halfExtent[1];
      const gains = orbitalGains(tables, el, state.density, viewArea, state.equalize);
      const count = uploadOrbitals(el, visible, gains);

      // The pane's 1/k buffer region covers bw*k x bh*k device pixels
      // (>= rect); the vertex shader scales its NDC so geometry matches
      // the native view.
      gl.uniform2f(u.uNdcScale, (bw * k) / rect.w, (bh * k) / rect.h);
      gl.uniform2f(u.uHalfExtent, geo.halfExtent[0], geo.halfExtent[1]);
      gl.uniform1i(u.uCount, count);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
      geo.gains = gains;
      geo.rect = rect;
      return geo;
    });

    if (k > 1) {
      gl.bindFramebuffer(gl.READ_FRAMEBUFFER, fbo);
      gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, null);
      gl.blitFramebuffer(0, 0, BW, BH, 0, 0, BW * k, BH * k, gl.COLOR_BUFFER_BIT, gl.NEAREST);
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    }
    return geos;
  }

  // Single atom filling the canvas.
  function render(state, k = 1) {
    const rect = { x: 0, y: 0, w: canvas.width, h: canvas.height };
    return renderViews(state, k, [{ element: state.element, visible: state.visible, rect }])[0];
  }

  return { gl, render, renderViews };
}
