#version 300 es
// Vertex shader: full-screen triangle and the geometry of each pixel's ray.
//
// Rays are orthographic along the view axis. In view space the ray through a
// pixel is (x, y, 0) + t (0, 0, 1), so (x, y, 0) is already its closest point
// to the nucleus. That point is linear in screen position, so it is computed
// here at the triangle's corners and interpolated exactly for every pixel.

uniform vec2 uNdcScale;    // (render-buffer size * k) / canvas size; (1, 1) at full resolution
uniform vec2 uHalfExtent;  // half width/height of the view (Bohr)
uniform mat3 uRot;         // view -> atom frame

out vec3 vP;               // closest point of the pixel's ray to the nucleus (atom frame, Bohr)
flat out vec3 vD;          // ray direction (atom frame); the same for every pixel

void main() {
  // Full-screen triangle; no vertex buffers needed. p is 0 or 2 per axis.
  vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
  gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);

  // At render scale k each fragment of the 1/k buffer stands for the centre
  // of a k x k block of device pixels; the buffer overhangs the canvas by
  // < k pixels, which uNdcScale accounts for, so the geometry matches the
  // native-resolution view exactly.
  vec2 ndc = p * 2.0 * uNdcScale - 1.0;
  vec2 xy = ndc * uHalfExtent;
  vP = uRot * vec3(xy, 0.0);
  vD = uRot * vec3(0.0, 0.0, 1.0);
}
