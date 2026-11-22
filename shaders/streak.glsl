#version 300 es
// Streak pass (prototype; color mode "nearest" with streaks): draws every
// detection of the packed detection texture (fragment.glsl compiled with
// STREAKS) as a short line through its pixel, along its sampled momentum.
//
// Gather, not scatter: each output pixel looks at the (2R + 1)^2 detections
// around it and keeps the nearest (largest depth) of those whose segment
// covers it -- the same "nearest wins" rule as the color mode, with no
// blending, depth buffer or draw order. Pixels covered by no segment keep
// their own color (the background where there is no detection).
//
// Host-provided defines (js/gl.js createProgram):
//   STREAK_R  longest streak half-length in render-buffer pixels
#ifndef STREAK_R
#error STREAK_R must be defined by the host
#endif

precision highp float;
precision highp int;
precision highp usampler2D;

uniform usampler2D uDetections;  // x: RGB8 + has bit 24, y: half2 half-vector, z: depth bits
uniform ivec4 uPane;             // pane rect in render-buffer pixels (x, y, w, h)

out vec4 outColor;

vec3 unpackRGB(uint v) {
  return vec3(uvec3(v, v >> 8, v >> 16) & 255u) * (1.0 / 255.0);
}

void main() {
  ivec2 c = ivec2(gl_FragCoord.xy);
  vec3 rgb = unpackRGB(texelFetch(uDetections, c, 0).x);
  float best = -3.0e38;
  // Detections from other panes never reach into this one.
  ivec2 lo = max(c - STREAK_R, uPane.xy);
  ivec2 hi = min(c + STREAK_R, uPane.xy + uPane.zw - 1);

  for (int dy = -STREAK_R; dy <= STREAK_R; dy++) {
    for (int dx = -STREAK_R; dx <= STREAK_R; dx++) {
      ivec2 n = c + ivec2(dx, dy);
      if (any(lessThan(n, lo)) || any(greaterThan(n, hi))) continue;
      uvec4 t = texelFetch(uDetections, n, 0);
      if ((t.x >> 24) == 0u) continue;
      // Distance from this pixel's centre to the segment n - h .. n + h:
      // under half a pixel draws a connected one-pixel-wide line.
      vec2 h = unpackHalf2x16(t.y);
      vec2 o = vec2(-dx, -dy);
      float hh = dot(h, h);
      float s = hh > 0.0 ? clamp(dot(o, h) / hh, -1.0, 1.0) : 0.0;
      vec2 e = o - s * h;
      if (dot(e, e) >= 0.25) continue;
      float depth = uintBitsToFloat(t.z);
      if (depth > best) {
        best = depth;
        rgb = unpackRGB(t.x);
      }
    }
  }
  outColor = vec4(rgb, 1.0);
}
