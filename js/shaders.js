// GLSL ES 3.00 sources.

export const VERT = `#version 300 es
// Full-screen triangle; no vertex buffers needed.
void main() {
  vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
  gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}
`;

import { NUM_ORBITALS } from './orbitals.js';
export { NUM_ORBITALS };

export const FRAG = `#version 300 es
precision highp float;
precision highp int;
precision highp sampler2D;

#define NUM_ORBITALS ${NUM_ORBITALS}

// Moment tables: row r, texel i = (I0, I2, I4, I6) of radial function r at
// scaled impact parameter u = i * uTableUMax / (uTableN - 1), u = Z b / n.
uniform sampler2D uTables;
uniform int   uTableN;
uniform float uTableUMax;

uniform vec2  uResolution;   // full canvas size (device px)
uniform float uScale;        // render scale k: this pass draws 1/k resolution
uniform vec2  uHalfExtent;   // half width/height of the view (Bohr)
uniform mat3  uRot;          // view -> atom frame
uniform uint  uSeed;         // per-frame random seed
uniform float uKMax;         // detections giving full brightness
uniform float uGamma;        // tone curve: amount = (k / uKMax)^uGamma
uniform int   uMode;         // 0 = dark background, 1 = white background

// Occupied, visible orbitals, packed into the first uCount slots of a
// uniform block (std140) in canonical order.
struct Orbital {
  vec4 orb;    // (occupancy, kappa = Z/n, hue in [0,1), gain)
  vec4 kind;   // (l, table row, chroma, unused)
  vec4 c0;     // angular coefficients (see js/orbitals.js angularCoefs):
  vec4 c1;     //  p: e.xyz; d: q00 q01 q02 q11 | q12 q22;
  vec4 c2;     //  f: c000 c001 c002 c011 | c012 c022 c111 c112 | c122 c222
};
layout(std140) uniform Orbitals {
  Orbital uOrbitals[NUM_ORBITALS];
};
uniform int uCount;

out vec4 outColor;

// ---------------------------------------------------------------- random
uint pcg(uint v) {
  uint state = v * 747796405u + 2891336453u;
  uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
  return (word >> 22u) ^ word;
}

float rand(inout uint s) {
  s = pcg(s);
  return float(s) * (1.0 / 4294967296.0);
}

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
int poisson(float lambda, inout uint s) {
  if (lambda <= 0.0) return 0;
  if (lambda < 10.0) {
    float u = rand(s);
    float p = exp(-lambda), F = p;
    int k = 0;
    for (int i = 0; i < 64; i++) {
      if (u <= F) break;
      k++;
      p *= lambda / float(k);
      F += p;
    }
    return k;
  }
  float slam = sqrt(lambda), loglam = log(lambda);
  float b = 0.931 + 2.53 * slam;
  float a = -0.059 + 0.02483 * b;
  float invalpha = 1.1239 + 1.1328 / (b - 3.4);
  float vr = 0.9277 - 3.6224 / (b - 2.0);
  for (int i = 0; i < 16; i++) {
    float U = rand(s) - 0.5;
    float V = rand(s);
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

// ---------------------------------------------------------------- tables
// Manual linear interpolation (float textures need not be filterable).
vec4 momentAt(int row, float u) {
  float x = u / uTableUMax * float(uTableN - 1);
  if (x >= float(uTableN - 1)) return vec4(0.0);
  int i = int(x);
  float f = x - float(i);
  return mix(texelFetch(uTables, ivec2(i, row), 0),
             texelFetch(uTables, ivec2(i + 1, row), 0), f);
}

// Column density in reference units (Z0 = n) for scaled closest point P.
//  s: I0 / 4pi
//  p: (e.P)^2 I0 + (e.d)^2 I2
//  d: a^2 I0 + (4 beta^2 + 2 a c) I2 + c^2 I4,  a = PQP, beta = PQd, c = dQd
//  f: a0^2 I0 + (a1^2 + 2 a0 a2) I2 + (a2^2 + 2 a1 a3) I4 + a3^2 I6,
//     a0 = C(P,P,P), a1 = 3C(P,P,d), a2 = 3C(P,d,d), a3 = C(d,d,d)
float columnRef(Orbital o, vec3 P, vec3 d) {
  int l = int(o.kind.x + 0.5);
  vec4 m = momentAt(int(o.kind.y + 0.5), length(P));
  if (l == 0) return m.x * 0.0795774715459;   // 1 / (4 pi)
  if (l == 1) {
    float ep = dot(o.c0.xyz, P), ed = dot(o.c0.xyz, d);
    return ep * ep * m.x + ed * ed * m.y;
  }
  if (l == 2) {
    mat3 Q = mat3(o.c0.x, o.c0.y, o.c0.z,
                  o.c0.y, o.c0.w, o.c1.x,
                  o.c0.z, o.c1.x, o.c1.y);
    vec3 Qd = Q * d;
    float a = dot(P, Q * P), beta = dot(P, Qd), c = dot(d, Qd);
    return a * a * m.x + (4.0 * beta * beta + 2.0 * a * c) * m.y + c * c * m.z;
  }
  // f: contract the cubic tensor once with P and once with d
  //   C(a, b, c) = b . M(a) c,  M(a) = a.x Cx + a.y Cy + a.z Cz
  mat3 Cx = mat3(o.c0.x, o.c0.y, o.c0.z,  o.c0.y, o.c0.w, o.c1.x,  o.c0.z, o.c1.x, o.c1.y);
  mat3 Cy = mat3(o.c0.y, o.c0.w, o.c1.x,  o.c0.w, o.c1.z, o.c1.w,  o.c1.x, o.c1.w, o.c2.x);
  mat3 Cz = mat3(o.c0.z, o.c1.x, o.c1.y,  o.c1.x, o.c1.w, o.c2.x,  o.c1.y, o.c2.x, o.c2.y);
  mat3 MP = P.x * Cx + P.y * Cy + P.z * Cz;
  mat3 Md = d.x * Cx + d.y * Cy + d.z * Cz;
  vec3 MPd = MP * d;
  float a0 = dot(P, MP * P), a1 = 3.0 * dot(P, MPd), a2 = 3.0 * dot(d, MPd), a3 = dot(d, Md * d);
  return a0 * a0 * m.x + (a1 * a1 + 2.0 * a0 * a2) * m.y +
         (a2 * a2 + 2.0 * a1 * a3) * m.z + a3 * a3 * m.w;
}

vec3 hsv2rgb(vec3 c) {
  vec3 p = abs(fract(c.xxx + vec3(1.0, 2.0 / 3.0, 1.0 / 3.0)) * 6.0 - 3.0);
  return c.z * mix(vec3(1.0), clamp(p - 1.0, 0.0, 1.0), c.y);
}

void main() {
  // Orthographic ray through this pixel. In view space the ray is
  // (x, y, 0) + t (0, 0, 1), so (x, y, 0) is already the closest point to
  // the nucleus and b = |(x, y)|.
  // Each (low-res) fragment stands for the centre of a k x k block of device
  // pixels, so geometry matches the native-resolution view exactly.
  vec2 ndc = gl_FragCoord.xy * uScale / uResolution * 2.0 - 1.0;
  vec2 xy = ndc * uHalfExtent;
  vec3 p = uRot * vec3(xy, 0.0);
  vec3 d = uRot * vec3(0.0, 0.0, 1.0);

  uint s = pcg(uint(gl_FragCoord.x) + pcg(uint(gl_FragCoord.y) + pcg(uSeed)));

  vec2 chroma = vec2(0.0);
  float total = 0.0;

  for (int i = 0; i < NUM_ORBITALS; i++) {
    if (i >= uCount) break;
    Orbital o = uOrbitals[i];
    vec4 orb = o.orb;
    float kappa = orb.y;
    float col = kappa * kappa * columnRef(o, kappa * p, d);
    float lambda = orb.x * col * orb.w;
    int k = poisson(lambda, s);
    float sat = k > 0 ? min(1.0, pow(float(k) / uKMax, uGamma)) : 0.0;
    float h = orb.z * 6.28318530718;
    chroma += sat * o.kind.z * vec2(cos(h), sin(h));
    total += sat;
  }

  float hue = fract(atan(chroma.y, chroma.x) / 6.28318530718 + 1.0);
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
`;
