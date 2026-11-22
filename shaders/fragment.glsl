#version 300 es
// Fragment shader: per-pixel Poisson detections from orbital column densities.
//
// Layout (read top to bottom in order of importance):
//   Inputs                 uniforms (incl. the cut-away), the per-orbital uniform block, the pixel's ray
//   main()                 physics per pixel, then rendering to RGB
//   1. Physics             orbital column densities from the moment tables
//   2. Rendering helpers   color conversion
//   3. Numerical methods   3a. uniform random numbers, 3b. Poisson sampling
//
// Host-provided defines (injected by js/gl.js createProgram after #version):
//   NUM_ORBITALS  size of the Orbitals uniform block (js/orbitals.js)
//   STREAKS       (optional, prototype) nearest-color detections also get a
//                 momentum streak: instead of a color this pass writes a
//                 packed texel (color, streak half-vector, depth) that
//                 shaders/streak.glsl draws as short lines.
//   NUM_ROWS      (with STREAKS) number of radial table rows
#ifndef NUM_ORBITALS
#error NUM_ORBITALS must be defined by the host (createProgram defines)
#endif

// Steps of the march along a ray that places an orbital's nearest detection
// (color mode "nearest" only; see nearestDepthRef).
#define NEAREST_STEPS 48

#define TWO_PI (6.28318530718)

#ifdef STREAKS
// Longest streak half-length (render-buffer pixels); must match STREAK_R in
// shaders/streak.glsl (both come from js/shaders.js).
#ifndef STREAK_R
#error STREAK_R must be defined by the host
#endif
// Rejection-sampling attempts for a momentum direction (acceptance >= 1/7
// for f: P(no accept) = (6/7)^48 < 1e-3; then the last proposal is used).
#define MOMENTUM_TRIES 48
// Random streams of the momentum draws (disjoint from the Poisson streams
// < 2^20 and the nearest-depth streams 2^20 + slot).
#define STREAM_MOMENTUM_DIR (2u << 20)
#define STREAM_MOMENTUM_Q   (3u << 20)
#define STREAM_MOMENTUM_LOCAL (4u << 20)
#endif

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
uniform int   uColorMode;   // 0 = blend every orbital's detections;
                             // 1 = only the orbital with the nearest detection
// Depth fog (nearest color only): with d the distance behind the front of
// the atom (+uFogRadius, or the cut plane) as a fraction of the depth to
// its back (-uFogRadius), clamped to [0, 1], a detection's brightness is
// exp(-(FOG_SCALE uFog d)^2) (Gaussian, like GL's EXP2 fog): clear in
// front, falling off with distance squared, to exp(-(3 uFog)^2) at the back
// (0.04 at the default 0.6). A visual cue, not physics. uFog = 0: off.
uniform float uFog;
uniform float uFogRadius;    // Bohr: the element's atomic radius
const float FOG_SCALE = 3.0; // uFog = 1: e^-9 at the back

#ifdef STREAKS
// Momentum streaks: inverse CDF of each row's radial momentum density
// (scaled momentum q = p n / Z at F = (i + 0.5) / N, uTableN x rows, R32F),
// the view -> atom rotation (as in the vertex shader) and the streak
// half-length per unit of projected q (render-buffer pixels).
uniform sampler2D uMomentumTables;
uniform mat3  uRot;
uniform float uStreakLen;
// 0: global -- a sample of the orbital's whole momentum distribution;
// 1: local  -- a sample of the local momentum spread at the detection;
// 2: axis   -- the principal axis of that local spread.
uniform int   uStreakMode;
// Radial polynomial P of each row, R_ref = norm exp(-u) P(u) u^l: ascending
// coefficients c0..c7 in two vec4 per row (local modes).
uniform vec4  uRadialPoly[2 * NUM_ROWS];
#endif

// Occupied, visible orbitals, packed into the first uCount slots of a
// uniform block (std140) in canonical order.
struct Orbital {
  vec4 orb;    // (occupancy, kappa = Z/n, hue in [0,1), gain)
  vec4 kind;   // (l, table row, chroma, extent: radius (reference units)
               //  beyond which the radial density is negligible)
  vec4 c0;     // angular coefficients (see js/orbitals.js angularCoefs):
  vec4 c1;     //  p: e.xyz; d: q00 q01 q02 q11 | q12 q22;
  vec4 c2;     //  f: c000 c001 c002 c011 | c012 c022 c111 c112 | c122 c222
               // c2.w: max of A^2 over unit vectors (momentum direction bound)
};
layout(std140) uniform Orbitals {
  Orbital uOrbitals[NUM_ORBITALS];
};
uniform int uCount;

// This pixel's ray, from the vertex shader (atom frame, Bohr):
in vec3 vP;                  // closest point of the ray to the nucleus
flat in vec3 vD;             // unit direction

#ifdef STREAKS
// Packed detection (see main): x = RGB8 + has-detection bit 24, y = streak
// half-vector (pixels, half floats), z = depth bits, w = the nearest
// orbital's slot (for tests).
out uvec4 outData;
#else
out vec4 outColor;
#endif

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
#ifdef STREAKS
vec3  momentumDirection(Orbital o, uvec3 key);
float momentumQ(int row, float F);
vec2  streakHalf(Orbital o, uvec3 key, vec3 X);
mat3  hessLogPsi(Orbital o, vec3 X);
vec3  jacobiEigen(mat3 A, out mat3 V);
#endif
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

  // Nearest detection so far (color mode 1): depth t along d (Bohr; larger
  // is nearer the viewer) and that orbital's tone-mapped amount and color.
  const float NONE = -1e30;
  float nearT = NONE;
  float nearSat = 0.0;
  vec2 nearChroma = vec2(0.0);
  int nearSlot = -1;
#ifdef STREAKS
  vec3 nearX = vec3(0.0);     // the nearest detection, scaled by its orbital's kappa
#endif

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

    if (uColorMode == 1 && k > 0) {
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
        nearSlot = i;
#ifdef STREAKS
        nearX = P + tRef * d;
#endif
      }
    }
  }
  if (uColorMode == 1) {
    // Only the nearest orbital's detections are shown, in its own color.
    chroma = nearChroma;
    total = nearSat;
    if (uFog > 0.0 && nearSlot >= 0) {
      // Scaling chroma and total together dims the detection, keeping its
      // hue and purity; applied before the theme, so white stays the dual.
      float front = uCut != 0.0 ? 0.0 : uFogRadius;
      float dist = clamp((front - nearT) / (front + uFogRadius), 0.0, 1.0);
      float k = FOG_SCALE * uFog * dist;
      float fog = exp(-k * k);
      chroma *= fog;
      total *= fog;
    }
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
    // hue + 1/2. A single detection becomes its orbital color as ink, empty
    // space is white, and overlaps mix toward black.
    rgb = vec3(1.0) - hsv2rgb(vec3(fract(hue + 0.5), purity, amount));
  }
#ifdef STREAKS
  // The nearest detection's orbital also gives a momentum streak, drawn by
  // shaders/streak.glsl. Empty pixels keep their color, without a streak.
  vec2 h = vec2(0.0);
  uint has = 0u;
  if (uColorMode == 1 && nearSlot >= 0) {
    h = streakHalf(uOrbitals[nearSlot], key, nearX);
    has = 1u;
  }
  uvec3 c8 = uvec3(round(clamp(rgb, 0.0, 1.0) * 255.0));
  outData = uvec4(c8.r | (c8.g << 8) | (c8.b << 16) | (has << 24),
                  packHalf2x16(h), floatBitsToUint(nearT), uint(max(nearSlot, 0)));
#else
  outColor = vec4(rgb, 1.0);
#endif
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

#ifdef STREAKS
// --- Momentum (streaks) ---------------------------------------------------
// A detected position has no definite momentum, and the probability current
// of these real orbitals is zero, so there is no local velocity to draw.
// What is defined is the orbital's momentum distribution |phi(p)|^2, phi the
// Fourier transform of psi: it keeps psi's angular part A^2 and has its own
// radial part (tools/orbitals.py). A streak is one sample of that momentum,
// independent of where the detection lies: "what a momentum measurement of
// this electron would give instead". |phi(p)|^2 = |phi(-p)|^2, so streaks
// are undirected segments centred on the detection.

// Unit momentum direction with density A(p)^2: rejection from the uniform
// sphere against the bound c2.w (s: uniform). Mirrored in js/math.js.
vec3 momentumDirection(Orbital o, uvec3 key) {
  bool isS = o.kind.x < 0.5;
  vec3 dir = vec3(0.0, 0.0, 1.0);
  for (int j = 0; j < MOMENTUM_TRIES; j++) {
    uvec4 r = pcg4d(uvec4(key, STREAM_MOMENTUM_DIR + uint(j)));
    float z = 2.0 * u01(r.x) - 1.0;
    float s = sqrt(max(0.0, 1.0 - z * z));
    float ph = TWO_PI * u01(r.y);
    dir = vec3(s * cos(ph), s * sin(ph), z);
    if (isS) return dir;
    float a = angularPoly(o, dir, vec3(0.0)).x;       // A at the unit vector
    if (u01(r.z) * o.c2.w <= a * a) return dir;
  }
  return dir;
}

// Scaled momentum q for uniform F from the inverse-CDF table (lerp,
// clamped at the ends). Mirrored in js/math.js.
float momentumQ(int row, float F) {
  float x = clamp(F * float(uTableN) - 0.5, 0.0, float(uTableN - 1));
  int i = min(int(x), uTableN - 2);
  float f = x - float(i);
  return mix(texelFetch(uMomentumTables, ivec2(i, row), 0).r,
             texelFetch(uMomentumTables, ivec2(i + 1, row), 0).r, f);
}

// Local momentum spread at the detection (modes 1, 2). With x and p not
// commuting there is no joint distribution, but the Wigner function gives
// the second moment of momentum at a point x. For real psi it is
//     T(x) = <p p^T>_x = -(hbar^2 / 2) Hess ln|psi(x)|,
// and integrated over |psi|^2 it gives exactly the second moments of
// |phi(p)|^2. In the scaled frame of the orbital (X = kappa x) and in units
// of kappa (as q): T = -Hess_X ln|psi_ref(X)| / 2. It can be indefinite (the
// Wigner function is not a probability); negative eigenvalues are dropped.
// Mirrored in js/math.js.
//
// Hessian of ln|psi_ref| at scaled point X: psi_ref ~ exp(-u) P(u) A(X), so
// ln|psi| = -u + ln|P(u)| + ln|A(X)|, A homogeneous of degree l.
mat3 hessLogPsi(Orbital o, vec3 X) {
  int l = int(o.kind.x + 0.5);
  int row = int(o.kind.y + 0.5);
  float u = max(length(X), 1e-4);
  vec3 r = X / u;
  // P, P', P'' by Horner (unrolled: constant indices only)
  vec4 c0 = uRadialPoly[2 * row], c1 = uRadialPoly[2 * row + 1];
  float p0 = 0.0, p1 = 0.0, p2 = 0.0;
#define HORNER(c) p2 = p2 * u + 2.0 * p1; p1 = p1 * u + p0; p0 = p0 * u + (c);
  HORNER(c1.w) HORNER(c1.z) HORNER(c1.y) HORNER(c1.x)
  HORNER(c0.w) HORNER(c0.z) HORNER(c0.y) HORNER(c0.x)
#undef HORNER
  p0 = p0 < 0.0 ? min(p0, -1e-6) : max(p0, 1e-6);       // away from radial nodes
  float lp = p1 / p0;
  float g1 = -1.0 + lp;                // d/du ln|q|
  float g2 = p2 / p0 - lp * lp;        // d2/du2 ln|q|
  mat3 rr = outerProduct(r, r);
  mat3 H = g2 * rr + (g1 / u) * (mat3(1.0) - rr);
  if (l == 0) return H;
  // ln|A|: Hess = Hess A / A - grad A grad A^T / A^2
  float A;
  vec3 g;
  mat3 HA;
  if (l == 1) {
    A = dot(o.c0.xyz, X); g = o.c0.xyz; HA = mat3(0.0);
  } else if (l == 2) {
    mat3 Q = mat3(o.c0.x, o.c0.y, o.c0.z,  o.c0.y, o.c0.w, o.c1.x,  o.c0.z, o.c1.x, o.c1.y);
    A = dot(X, Q * X); g = 2.0 * (Q * X); HA = 2.0 * Q;
  } else {
    mat3 Cx = mat3(o.c0.x, o.c0.y, o.c0.z,  o.c0.y, o.c0.w, o.c1.x,  o.c0.z, o.c1.x, o.c1.y);
    mat3 Cy = mat3(o.c0.y, o.c0.w, o.c1.x,  o.c0.w, o.c1.z, o.c1.w,  o.c1.x, o.c1.w, o.c2.x);
    mat3 Cz = mat3(o.c0.z, o.c1.x, o.c1.y,  o.c1.x, o.c1.w, o.c2.x,  o.c1.y, o.c2.x, o.c2.y);
    mat3 M = X.x * Cx + X.y * Cy + X.z * Cz;
    A = dot(X, M * X); g = 3.0 * (M * X); HA = 6.0 * M;
  }
  float eps = 1e-6 * pow(u, float(l));                   // away from angular nodes
  A = A < 0.0 ? min(A, -eps) : max(A, eps);
  return H + HA / A - outerProduct(g, g) / (A * A);
}

// Eigen-decomposition of symmetric A by cyclic Jacobi rotations: returns
// the eigenvalues, V the eigenvectors (columns). Mirrored in js/math.js.
// One rotation zeroing A[p][q] (constant p < q; J[col][row]).
#define JACOBI_ROT(p, q) { \
    float apq = A[q][p]; \
    if (abs(apq) > 1e-12 * (abs(A[p][p]) + abs(A[q][q]))) { \
      float theta = (A[q][q] - A[p][p]) / (2.0 * apq); \
      float t = (theta >= 0.0 ? 1.0 : -1.0) / (abs(theta) + sqrt(theta * theta + 1.0)); \
      float c = inversesqrt(t * t + 1.0), s = t * c; \
      mat3 J = mat3(1.0); \
      J[p][p] = c; J[q][q] = c; J[q][p] = s; J[p][q] = -s; \
      A = transpose(J) * A * J; \
      V = V * J; \
    } }
vec3 jacobiEigen(mat3 A, out mat3 V) {
  V = mat3(1.0);
  for (int sweep = 0; sweep < 6; sweep++) {
    JACOBI_ROT(0, 1) JACOBI_ROT(0, 2) JACOBI_ROT(1, 2)
  }
  return vec3(A[0][0], A[1][1], A[2][2]);
}
#undef JACOBI_ROT

// Streak half-vector on screen (render-buffer pixels) for the detection at
// scaled point X: a momentum in units of kappa (so every shell shows its own
// pattern), rotated to view space, projected onto the screen, scaled by
// uStreakLen and clamped to STREAK_R.
//  global: q * direction from the whole momentum distribution;
//  local:  a Gaussian sample with covariance T (negative part dropped);
//  axis:   sqrt(largest eigenvalue of T) along its eigenvector (random
//          within the eigenspace where the largest eigenvalues coincide).
vec2 streakHalf(Orbital o, uvec3 key, vec3 X) {
  vec3 p;
  if (uStreakMode == 0) {
    vec3 dir = momentumDirection(o, key);
    p = momentumQ(int(o.kind.y + 0.5), u01(pcg4d(uvec4(key, STREAM_MOMENTUM_Q)).x)) * dir;
  } else {
    mat3 V;
    vec3 lam = max(jacobiEigen(-0.5 * hessLogPsi(o, X), V), 0.0);
    vec4 r = vec4(pcg4d(uvec4(key, STREAM_MOMENTUM_LOCAL)) >> 8u) * (1.0 / 16777216.0);
    vec2 m = sqrt(-2.0 * log(1.0 - r.xz));                 // Box-Muller, 1 - u in (0, 1]
    vec3 z = vec3(m.x * cos(TWO_PI * r.y), m.x * sin(TWO_PI * r.y), m.y * cos(TWO_PI * r.w));
    vec3 zV = transpose(V) * z;                            // z in the eigenbasis
    if (uStreakMode == 2) {
      // Principal axis, length sqrt(largest eigenvalue). Where the largest
      // eigenvalues coincide (the tangential plane of s) no axis is singled
      // out: the direction is random in that eigenspace, from the spectral
      // weight (lam / max)^8 applied to z -- a function of T alone.
      float lmax = max(lam.x, max(lam.y, lam.z));
      vec3 w = lam / max(lmax, 1e-30);
      w *= w; w *= w; w *= w;
      vec3 a = V * (w * zV);
      float al = length(a);
      p = al > 0.0 ? sqrt(lmax) * a / al : vec3(0.0);
    } else {
      // symmetric square root V sqrt(L) V^T: unique even where eigenvalues
      // coincide, so the sample does not depend on the eigenbasis Jacobi
      // happens to pick
      p = V * (sqrt(lam) * zV);
    }
  }
  vec2 h = uStreakLen * (transpose(uRot) * p).xy;         // atom -> view
  float len = length(h);
  if (!(len < 1e30)) return vec2(0.0);                    // non-finite at a node
  return len > float(STREAK_R) ? h * (float(STREAK_R) / len) : h;
}
#endif

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
