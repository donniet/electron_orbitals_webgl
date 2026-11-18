// GLSL ES 3.00 sources.

export const VERT = `#version 300 es
// Full-screen triangle; no vertex buffers needed.
void main() {
  vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
  gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}
`;

export const NUM_ORBITALS = 5;

export const FRAG = `#version 300 es
precision highp float;
precision highp int;
precision highp sampler2D;

#define NUM_ORBITALS ${NUM_ORBITALS}

// Column-density tables, RGBA = (A_2p, B_2p, S_1s, S_2s), scaled units.
uniform sampler2D uTables;
uniform int   uTableN;
uniform float uTableBMax;

uniform vec2  uResolution;   // framebuffer size (px)
uniform vec2  uHalfExtent;   // half width/height of the view (Bohr)
uniform mat3  uRot;          // view -> atom frame
uniform uint  uSeed;         // per-frame random seed
uniform float uKMax;         // detections giving full brightness
uniform float uGamma;        // tone curve: amount = (k / uKMax)^uGamma
uniform int   uMode;         // 0 = dark background, 1 = white background

uniform float uOcc[NUM_ORBITALS];   // electrons in each orbital (0 = skip)
uniform float uZeff[NUM_ORBITALS];  // effective nuclear charge
uniform float uHue[NUM_ORBITALS];   // hue in [0, 1)
uniform float uGain[NUM_ORBITALS];  // lambda = occ * Z^2 * col * uGain[i]

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
vec4 tableAt(float b) {
  float x = b / uTableBMax * float(uTableN - 1);
  if (x >= float(uTableN - 1)) return vec4(0.0);
  int i = int(x);
  float f = x - float(i);
  return mix(texelFetch(uTables, ivec2(i, 0), 0),
             texelFetch(uTables, ivec2(i + 1, 0), 0), f);
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
  float b = length(xy);

  uint s = pcg(uint(gl_FragCoord.x) + pcg(uint(gl_FragCoord.y) + pcg(uSeed)));

  vec2 chroma = vec2(0.0);
  float total = 0.0;

  for (int i = 0; i < NUM_ORBITALS; i++) {
    float occ = uOcc[i];
    if (occ <= 0.0) continue;
    float z = uZeff[i];
    vec4 t = tableAt(b * z);
    float col;
    if (i == 0) {
      col = t.z;                       // 1s
    } else if (i == 1) {
      col = t.w;                       // 2s
    } else {
      int axis = i - 2;                // 2px, 2py, 2pz
      float pa = p[axis] * z;
      float da = d[axis];
      col = pa * pa * t.x + da * da * t.y;
    }
    float lambda = occ * z * z * col * uGain[i];
    int k = poisson(lambda, s);
    float sat = k > 0 ? min(1.0, pow(float(k) / uKMax, uGamma)) : 0.0;
    float h = uHue[i] * 6.28318530718;
    chroma += sat * vec2(cos(h), sin(h));
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
