#version 300 es
// Fragment shader: per-pixel Poisson detections from orbital column densities.
//
// Layout (read top to bottom in order of importance):
//   Inputs                 uniforms (incl. the cut-away), the per-orbital uniform block, the pixel's ray
//   main()                 physics per pixel, then rendering to RGB
//   1. Physics             orbital column densities from the moment tables
//   2. Rendering helpers   colour conversion
//   3. Numerical methods   3a. uniform random numbers, 3b. Poisson sampling
//
// Host-provided defines (injected by js/gl.js createProgram after #version):
//   NUM_ORBITALS  size of the Orbitals uniform block (js/orbitals.js)
#ifndef NUM_ORBITALS
#error NUM_ORBITALS must be defined by the host (createProgram defines)
#endif

// Steps of the march along a ray that places an orbital's nearest detection
// (colour mode "nearest" only; see nearestDepthRef).
#define NEAREST_STEPS 48

#define TWO_PI (6.28318530718)

precision highp float;
precision highp int;
precision highp sampler2D;


// ==========================================================================
// Inputs
// ==========================================================================
// Moment tables: row r, texel i = (I0, I2, I4, I6) of radial function r at
// scaled impact parameter u = i * uTableUMax / (uTableN - 1), u = Z b / n.
uniform sampler2D uTables;
uniform int   uTableN;
uniform float uTableUMax;
// Odd half-line moments for the cut-away view, same layout:
// texel = (J1, J3, J5, 0), J_k = integral over t >= 0 only.
uniform sampler2D uOddTables;
// Radial factor q2(u) = R^2 / u^(2l) (reference units), same rows and u grid;
// single channel. Used to march along a ray.
uniform sampler2D uRadialTables;

// Cut-away: 0 = whole ray; -1 = only the half behind the plane through the
// nucleus parallel to the screen (t <= 0; view z points at the viewer),
// at double gain so front/back-symmetric parts keep their brightness.
uniform float uCut;

uniform uint  uSeed;         // per-frame random seed
uniform float uKMax;         // detections giving full brightness
uniform float uGamma;        // tone curve: amount = (k / uKMax)^uGamma
uniform int   uMode;         // 0 = dark background, 1 = white background
uniform int   uColourMode;   // 0 = blend every orbital's detections;
                             // 1 = only the orbital with the nearest detection

// Occupied, visible orbitals, packed into the first uCount slots of a
// uniform block (std140) in canonical order.
struct Orbital {
  vec4 orb;    // (occupancy, kappa = Z/n, hue in [0,1), gain)
  vec4 kind;   // (l, table row, chroma, extent: radius (reference units)
               //  beyond which the radial density is negligible)
  vec4 c0;     // angular coefficients (see js/orbitals.js angularCoefs):
  vec4 c1;     //  p: e.xyz; d: q00 q01 q02 q11 | q12 q22;
  vec4 c2;     //  f: c000 c001 c002 c011 | c012 c022 c111 c112 | c122 c222
};
layout(std140) uniform Orbitals {
  Orbital uOrbitals[NUM_ORBITALS];
};
uniform int uCount;

// This pixel's ray, from the vertex shader (atom frame, Bohr):
in vec3 vP;                  // closest point of the ray to the nucleus
flat in vec3 vD;             // unit direction

out vec4 outColor;

// ==========================================================================
// Forward declarations (definitions follow main)
// ==========================================================================
// 1. Physics
vec4  momentAt(sampler2D tables, int row, float u);
vec4  angularPoly(Orbital o, vec3 P, vec3 d);
float columnRef(Orbital o, vec3 P, vec3 d);
float radialAt(int row, float u);
float rayDensityS(int row, vec4 A, float u2, float H, float s);
float nearestDepthRef(Orbital o, vec3 P, vec3 d, int k, float massRef, float v, float tStop);
// 2. Rendering helpers
vec3  hsv2rgb(vec3 c);
// 3. Numerical methods
int   poisson(float lambda, uvec3 key, uint stream);
float logFactorial(float k);
uvec4 pcg4d(uvec4 v);
float u01(uint v);

void main() {
  // --- Physics: expected detections lambda for each orbital along this
  //     pixel's ray, sampled as Poisson counts.
  // The ray through this pixel: p + t d (p . d = 0), from the vertex shader.
  vec3 p = vP;
  vec3 d = vD;

  uvec3 key = uvec3(uvec2(gl_FragCoord.xy), uSeed);

  vec2 chroma = vec2(0.0);
  float total = 0.0;

  // Nearest detection so far (colour mode 1): depth t along d (Bohr; larger
  // is nearer the viewer) and that orbital's tone-mapped amount and colour.
  const float NONE = -1e30;
  float nearT = NONE;
  float nearSat = 0.0;
  vec2 nearChroma = vec2(0.0);

  for (int i = 0; i < NUM_ORBITALS; i++) {
    if (i >= uCount) break;
    Orbital o = uOrbitals[i];
    vec4 orb = o.orb;
    float kappa = orb.y;
    vec3 P = kappa * p;                         // scaled ray: P + (kappa t) d
    float colRef = columnRef(o, P, d);
    float col = kappa * kappa * colRef;
    float lambda = orb.x * col * orb.w;
    int k = poisson(lambda, key, uint(i));      // stream = orbital slot (< 64)
    float sat = k > 0 ? min(1.0, pow(float(k) / uKMax, uGamma)) : 0.0;
    float h = orb.z * TWO_PI;
    vec2 c = sat * o.kind.z * vec2(cos(h), sin(h));
    chroma += c;
    total += sat;

    if (uColourMode == 1 && k > 0) {
      // Where is this orbital's nearest detection? Its k detections are
      // independent draws from its density along the ray, so the nearest
      // sits at cumulative fraction 1 - v^(1/k) from the front (v uniform).
      // The march stops early once it is behind the nearest found so far.
      float v = u01(pcg4d(uvec4(key, (1u << 20) + uint(i))).x);  // own stream range
      float massRef = uCut != 0.0 ? 0.5 * colRef : colRef;       // the marched part
      float tRef = nearestDepthRef(o, P, d, k, massRef, v, nearT * kappa);
      if (tRef / kappa > nearT) {
        nearT = tRef / kappa;
        nearSat = sat;
        nearChroma = c;
      }
    }
  }
  if (uColourMode == 1) {
    // Only the nearest orbital's detections are shown, in its own colour.
    chroma = nearChroma;
    total = nearSat;
  }

  // --- Rendering: tone curve per orbital, hue compositing, then RGB.
  float hue = fract(atan(chroma.y, chroma.x) / TWO_PI + 1.0);
  float purity = total > 0.0 ? length(chroma) / total : 0.0;
  float amount = min(total, 1.0);

  vec3 rgb;
  if (uMode == 0) {
    // Dark (light on black): detections emit in the orbital hue; overlaps
    // mix toward white.
    rgb = hsv2rgb(vec3(hue, purity, amount));
  } else {
    // White (ink on paper): the exact dual of the dark image -- the negative
    // of the dark image computed with every hue rotated 180 degrees. Rotating
    // all orbital hues rotates the summed hue vector, so that is just
    // hue + 1/2. A single detection becomes its orbital colour as ink, empty
    // space is white, and overlaps mix toward black.
    rgb = vec3(1.0) - hsv2rgb(vec3(fract(hue + 0.5), purity, amount));
  }
  outColor = vec4(rgb, 1.0);
}

// ==========================================================================
// 1. Physics: orbital column densities
// ==========================================================================
// Manual linear interpolation (float textures need not be filterable).
vec4 momentAt(sampler2D tables, int row, float u) {
  float x = u / uTableUMax * float(uTableN - 1);
  if (x >= float(uTableN - 1)) return vec4(0.0);
  int i = int(x);
  float f = x - float(i);
  return mix(texelFetch(tables, ivec2(i, row), 0),
             texelFetch(tables, ivec2(i + 1, row), 0), f);
}

// Angular amplitude along the scaled ray P + t d as a polynomial in t:
// A(t) = c.x + c.y t + c.z t^2 + c.w t^3, with |psi|^2 = q2(r) A(t)^2.
//  s: 1 / sqrt(4 pi)
//  p: (e.P, e.d)
//  d: (a, 2 beta, c),  a = PQP, beta = PQd, c = dQd
//  f: (a0, a1, a2, a3), a0 = C(P,P,P), a1 = 3C(P,P,d), a2 = 3C(P,d,d), a3 = C(d,d,d)
vec4 angularPoly(Orbital o, vec3 P, vec3 d) {
  int l = int(o.kind.x + 0.5);
  if (l == 0) return vec4(0.28209479177, 0.0, 0.0, 0.0);
  if (l == 1) return vec4(dot(o.c0.xyz, P), dot(o.c0.xyz, d), 0.0, 0.0);
  if (l == 2) {
    mat3 Q = mat3(o.c0.x, o.c0.y, o.c0.z,
                  o.c0.y, o.c0.w, o.c1.x,
                  o.c0.z, o.c1.x, o.c1.y);
    vec3 Qd = Q * d;
    return vec4(dot(P, Q * P), 2.0 * dot(P, Qd), dot(d, Qd), 0.0);
  }
  // f: contract the cubic tensor once with P and once with d
  //   C(a, b, c) = b . M(a) c,  M(a) = a.x Cx + a.y Cy + a.z Cz
  mat3 Cx = mat3(o.c0.x, o.c0.y, o.c0.z,  o.c0.y, o.c0.w, o.c1.x,  o.c0.z, o.c1.x, o.c1.y);
  mat3 Cy = mat3(o.c0.y, o.c0.w, o.c1.x,  o.c0.w, o.c1.z, o.c1.w,  o.c1.x, o.c1.w, o.c2.x);
  mat3 Cz = mat3(o.c0.z, o.c1.x, o.c1.y,  o.c1.x, o.c1.w, o.c2.x,  o.c1.y, o.c2.x, o.c2.y);
  mat3 MP = P.x * Cx + P.y * Cy + P.z * Cz;
  mat3 Md = d.x * Cx + d.y * Cy + d.z * Cz;
  vec3 MPd = MP * d;
  return vec4(dot(P, MP * P), 3.0 * dot(P, MPd), 3.0 * dot(d, MPd), dot(d, Md * d));
}

// Column density in reference units (Z0 = n) for scaled closest point P.
// Whole ray (even moments I, odd powers of t cancel), with angularPoly's terms:
//  s: I0 / 4pi
//  p: (e.P)^2 I0 + (e.d)^2 I2
//  d: a^2 I0 + (4 beta^2 + 2 a c) I2 + c^2 I4
//  f: a0^2 I0 + (a1^2 + 2 a0 a2) I2 + (a2^2 + 2 a1 a3) I4 + a3^2 I6
// Cut-away adds 2 uCut times the odd terms (half-line moments J):
//  s: 0
//  p: 2 (e.P)(e.d) J1
//  d: 4 a beta J1 + 4 beta c J3
//  f: 2 a0 a1 J1 + 2 (a0 a3 + a1 a2) J3 + 2 a2 a3 J5
// (whole = E; half behind = E/2 - O; shown at double gain: E - 2 O).
float columnRef(Orbital o, vec3 P, vec3 d) {
  int l = int(o.kind.x + 0.5);
  int row = int(o.kind.y + 0.5);
  float u = length(P);
  vec4 m = momentAt(uTables, row, u);
  if (l == 0) return m.x * 0.0795774715459;   // 1 / (4 pi)
  vec4 A = angularPoly(o, P, d);
  if (l == 1) {
    float ep = A.x, ed = A.y;
    float col = ep * ep * m.x + ed * ed * m.y;
    if (uCut != 0.0) col += 2.0 * uCut * (2.0 * ep * ed * momentAt(uOddTables, row, u).x);
    return col;
  }
  if (l == 2) {
    float a = A.x, beta = 0.5 * A.y, c = A.z;
    float col = a * a * m.x + (4.0 * beta * beta + 2.0 * a * c) * m.y + c * c * m.z;
    if (uCut != 0.0) {
      vec4 j = momentAt(uOddTables, row, u);
      col += 2.0 * uCut * (4.0 * a * beta * j.x + 4.0 * beta * c * j.y);
    }
    return col;
  }
  float a0 = A.x, a1 = A.y, a2 = A.z, a3 = A.w;
  float col = a0 * a0 * m.x + (a1 * a1 + 2.0 * a0 * a2) * m.y +
              (a2 * a2 + 2.0 * a1 * a3) * m.z + a3 * a3 * m.w;
  if (uCut != 0.0) {
    vec4 j = momentAt(uOddTables, row, u);
    col += 2.0 * uCut * (2.0 * a0 * a1 * j.x + 2.0 * (a0 * a3 + a1 * a2) * j.y + 2.0 * a2 * a3 * j.z);
  }
  return col;
}

// Radial factor q2 at u (manual lerp, 0 beyond the table).
float radialAt(int row, float u) {
  float x = u / uTableUMax * float(uTableN - 1);
  if (x >= float(uTableN - 1)) return 0.0;
  int i = int(x);
  float f = x - float(i);
  return mix(texelFetch(uRadialTables, ivec2(i, row), 0).r,
             texelFetch(uRadialTables, ivec2(i + 1, row), 0).r, f);
}

// Density along the scaled ray per unit of the march variable s, where
// t = H s^3 (dt/ds = 3 H s^2): q2(r) A(t)^2 dt/ds.
float rayDensityS(int row, vec4 A, float u2, float H, float s) {
  float t = H * s * s * s;
  float a = A.x + t * (A.y + t * (A.z + t * A.w));
  return radialAt(row, sqrt(u2 + t * t)) * a * a * 3.0 * H * s * s;
}

// Depth t (reference units along d; larger = nearer the viewer) of the
// nearest of k detections of orbital o on the scaled ray P + t d. massRef is
// the column of the part of the ray that is drawn (half of it when cut);
// v is uniform in [0, 1).
// The density is marched from the front (s = 1, or the cut plane s = 0) to
// the back (s = -1) with the trapezoid rule in s, t = H s^3,
// H = sqrt(extent^2 - |P|^2): the cubic packs the steps near the closest
// point, where the narrow inner lobes of high-n orbitals lie. The march
// stops when the cumulative mass reaches (1 - v^(1/k)) massRef; the crossing
// is interpolated linearly in s. If the march gets behind tStop first, the
// result cannot be the nearest overall: returns -1e30. Rays outside the
// extent return 0 (the closest point).
float nearestDepthRef(Orbital o, vec3 P, vec3 d, int k, float massRef, float v, float tStop) {
  float u2 = dot(P, P);
  float ext = o.kind.w;
  if (u2 >= ext * ext) return 0.0;
  float H = sqrt(ext * ext - u2);
  float sHi = uCut != 0.0 ? 0.0 : 1.0;
  float hs = (sHi + 1.0) / float(NEAREST_STEPS);
  int row = int(o.kind.y + 0.5);
  vec4 A = angularPoly(o, P, d);
  float target = (1.0 - pow(v, 1.0 / float(k))) * massRef;

  float fPrev = rayDensityS(row, A, u2, H, sHi);
  float acc = 0.0;
  for (int j = 1; j <= NEAREST_STEPS; j++) {
    float s = sHi - float(j) * hs;
    float f = rayDensityS(row, A, u2, H, s);
    float seg = 0.5 * (fPrev + f) * hs;
    if (acc + seg >= target) {
      float sc = s + hs * (1.0 - (target - acc) / max(seg, 1e-30));
      return H * sc * sc * sc;
    }
    acc += seg;
    fPrev = f;
    if (H * s * s * s < tStop) return -1e30;
  }
  return -H;    // short of the target by the quadrature error: the back end
}

// ==========================================================================
// 2. Rendering helpers
// ==========================================================================
vec3 hsv2rgb(vec3 c) {
  vec3 p = abs(fract(c.xxx + vec3(1.0, 2.0 / 3.0, 1.0 / 3.0)) * 6.0 - 3.0);
  return c.z * mix(vec3(1.0), clamp(p - 1.0, 0.0, 1.0), c.y);
}

// ==========================================================================
// 3. Numerical methods
// ==========================================================================
// --- 3a. Uniform random numbers -------------------------------------------
// Counter-based randomness: each orbital sample in each pixel hashes a
// unique (x, y, frame seed, stream) with pcg4d (Jarzynski & Olano 2020), a
// bijection on 128 bits, so no two pixels or samples can share or merge a
// stream. One hash gives four uniforms (top 24 bits: exact in float32,
// in [0, 1)); PTRS re-hashes with an offset stream if it needs more.
uvec4 pcg4d(uvec4 v) {
  v = v * 1664525u + 1013904223u;
  v.x += v.y * v.w; v.y += v.z * v.x; v.z += v.x * v.y; v.w += v.y * v.z;
  v ^= v >> 16u;
  v.x += v.y * v.w; v.y += v.z * v.x; v.z += v.x * v.y; v.w += v.y * v.z;
  return v;
}

float u01(uint v) { return float(v >> 8u) * (1.0 / 16777216.0); }   // top 24 bits, [0, 1)

// --- 3b. Poisson sampling -------------------------------------------------
// log(k!) : exact table for small k, Stirling series above (float32-safe).
float logFactorial(float k) {
  if (k < 10.0) {
    const float LF[10] = float[10](0.0, 0.0, 0.69314718, 1.79175947, 3.17805383,
                                   4.78749174, 6.57925121, 8.52516136, 10.6046029, 12.8018275);
    return LF[int(k)];
  }
  float ik = 1.0 / k, ik2 = ik * ik;
  return (k + 0.5) * log(k) - k + 0.91893853320467        // 0.5 log(2 pi)
       + ik * (1.0 / 12.0 - ik2 * (1.0 / 360.0 - ik2 / 1260.0));
}

// Exact Poisson sampler with bounded cost.
//  lambda < 10 : inversion with one uniform (pmf by recurrence).
//  lambda >= 10: PTRS, transformed rejection with squeeze (Hormann 1993);
//                ~1.1-1.3 iterations of 2 uniforms on average.
int poisson(float lambda, uvec3 key, uint stream) {
  uvec4 h = pcg4d(uvec4(key, stream));
  if (lambda <= 0.0) return 0;
  if (lambda < 10.0) {
    float u = u01(h.x);
    float p = exp(-lambda), F = p;
    int k = 0;
    for (int i = 0; i < 64; i++) {
      if (u <= F) return k;
      k++;
      p *= lambda / float(k);
      F += p;
    }
    // u fell beyond the float32 limit of the cumulative sum (~1e-7 below 1);
    // return the mode rather than the loop cap, so there are no hot pixels.
    return int(lambda);
  }
  float slam = sqrt(lambda), loglam = log(lambda);
  float b = 0.931 + 2.53 * slam;
  float a = -0.059 + 0.02483 * b;
  float invalpha = 1.1239 + 1.1328 / (b - 3.4);
  float vr = 0.9277 - 3.6224 / (b - 2.0);
  for (int i = 0; i < 16; i++) {
    // two uniforms per iteration: (x, y) then (z, w) of each hash
    if (i > 0 && (i & 1) == 0) h = pcg4d(uvec4(key, stream + uint(i) * 1024u));
    bool first = (i & 1) == 0;
    float U = u01(first ? h.x : h.z) - 0.5;
    float V = u01(first ? h.y : h.w);
    float us = 0.5 - abs(U);
    float k = floor((2.0 * a / us + b) * U + lambda + 0.43);
    if (us >= 0.07 && V <= vr) return int(k);
    if (k < 0.0 || (us < 0.013 && V > us)) continue;
    if (log(V) + log(invalpha) - log(a / (us * us) + b) <= -lambda + k * loglam - logFactorial(k)) {
      return int(k);
    }
  }
  return int(floor(lambda + 0.5));     // practically unreachable (p < 1e-15)
}
