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

// Moment tables: row r, texel i = (I0, I2, I4, 0) of radial function r at
// scaled impact parameter u = i * uTableUMax / (uTableN - 1), u = Z b / n.
uniform sampler2D uTables;
uniform int   uTableN;
uniform float uTableUMax;

uniform vec2  uResolution;   // framebuffer size (px)
uniform vec2  uHalfExtent;   // half width/height of the view (Bohr)
uniform mat3  uRot;          // view -> atom frame
uniform uint  uSeed;         // per-frame random seed
uniform float uKMax;         // detections giving full brightness
uniform float uGamma;        // tone curve: amount = (k / uKMax)^uGamma
uniform int   uMode;         // 0 = dark background, 1 = white background

// Per orbital (packed to stay well inside the uniform-vector budget):
uniform vec4  uOrb[NUM_ORBITALS];   // (occupancy (0 = skip), kappa = Z/n, hue in [0,1), gain)
uniform ivec2 uKind[NUM_ORBITALS];  // (l, table row)
uniform vec3  uAxis[NUM_ORBITALS];  // p: sqrt(3/4pi) * lobe axis
uniform mat3  uQ[NUM_ORBITALS];     // d: quadratic form of the real harmonic
uniform float uChroma[NUM_ORBITALS]; // weight of the orbital's hue (outer p shells paler)

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

int poisson(float lambda, inout uint s) {
  if (lambda <= 0.0) return 0;
  if (lambda < 30.0) {
    // Knuth: multiply uniforms until the product drops below e^-lambda.
    float L = exp(-lambda);
    float p = 1.0;
    int k = 0;
    for (int i = 0; i < 64; i++) {
      p *= rand(s);
      if (p <= L) break;
      k++;
    }
    return k;
  }
  // Large lambda: normal approximation.
  float u1 = max(rand(s), 1e-7);
  float u2 = rand(s);
  float n = sqrt(-2.0 * log(u1)) * cos(6.28318530718 * u2);
  return max(0, int(floor(lambda + sqrt(lambda) * n + 0.5)));
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
float columnRef(int i, vec3 P, vec3 d) {
  ivec2 kind = uKind[i];
  vec4 m = momentAt(kind.y, length(P));
  if (kind.x == 0) return m.x * 0.0795774715459;   // 1 / (4 pi)
  if (kind.x == 1) {
    float ep = dot(uAxis[i], P), ed = dot(uAxis[i], d);
    return ep * ep * m.x + ed * ed * m.y;
  }
  mat3 Q = uQ[i];
  vec3 Qd = Q * d;
  float a = dot(P, Q * P), beta = dot(P, Qd), c = dot(d, Qd);
  return a * a * m.x + (4.0 * beta * beta + 2.0 * a * c) * m.y + c * c * m.z;
}

vec3 hsv2rgb(vec3 c) {
  vec3 p = abs(fract(c.xxx + vec3(1.0, 2.0 / 3.0, 1.0 / 3.0)) * 6.0 - 3.0);
  return c.z * mix(vec3(1.0), clamp(p - 1.0, 0.0, 1.0), c.y);
}

void main() {
  // Orthographic ray through this pixel. In view space the ray is
  // (x, y, 0) + t (0, 0, 1), so (x, y, 0) is already the closest point to
  // the nucleus and b = |(x, y)|.
  vec2 ndc = gl_FragCoord.xy / uResolution * 2.0 - 1.0;
  vec2 xy = ndc * uHalfExtent;
  vec3 p = uRot * vec3(xy, 0.0);
  vec3 d = uRot * vec3(0.0, 0.0, 1.0);

  uint s = pcg(uint(gl_FragCoord.x) + pcg(uint(gl_FragCoord.y) + pcg(uSeed)));

  vec2 chroma = vec2(0.0);
  float total = 0.0;

  for (int i = 0; i < NUM_ORBITALS; i++) {
    vec4 orb = uOrb[i];
    if (orb.x <= 0.0) continue;
    float kappa = orb.y;
    float col = kappa * kappa * columnRef(i, kappa * p, d);
    float lambda = orb.x * col * orb.w;
    int k = poisson(lambda, s);
    float sat = k > 0 ? min(1.0, pow(float(k) / uKMax, uGamma)) : 0.0;
    float h = orb.z * 6.28318530718;
    chroma += sat * uChroma[i] * vec2(cos(h), sin(h));
    total += sat;
  }

  float hue = fract(atan(chroma.y, chroma.x) / 6.28318530718 + 1.0);
  float purity = total > 0.0 ? length(chroma) / total : 0.0;
  float amount = min(total, 1.0);

  vec3 rgb;
  if (uMode == 0) {
    // Dark: detections light up in the orbital hue.
    rgb = hsv2rgb(vec3(hue, purity, amount));
  } else {
    // White (literal): count drives saturation at full value. Hue mixes
    // that cancel (low purity) are darkened slightly so they stay visible.
    rgb = hsv2rgb(vec3(hue, amount * purity, 1.0 - 0.35 * amount * (1.0 - purity)));
  }
  outColor = vec4(rgb, 1.0);
}
`;
