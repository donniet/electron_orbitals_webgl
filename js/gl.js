// Minimal WebGL2 helpers.

export function createContext(canvas) {
  const gl = canvas.getContext('webgl2', { antialias: false, preserveDrawingBuffer: true });
  if (!gl) throw new Error('WebGL2 is not available in this browser.');
  return gl;
}

function compile(gl, type, src) {
  const sh = gl.createShader(type);
  gl.shaderSource(sh, src);
  gl.compileShader(sh);
  if (!gl.getShaderParameter(sh, gl.COMPILE_STATUS)) {
    const log = gl.getShaderInfoLog(sh);
    gl.deleteShader(sh);
    throw new Error('Shader compile error:\n' + log);
  }
  return sh;
}

// Insert `#define NAME value` lines right after the first line (which must
// stay `#version ...`), then `#line 2` so compiler errors still report the
// line numbers of the .glsl file.
export function withDefines(src, defines = {}) {
  const entries = Object.entries(defines);
  if (entries.length === 0) return src;
  const nl = src.indexOf('\n');
  const head = src.slice(0, nl + 1), body = src.slice(nl + 1);
  if (!/^\s*#version/.test(head)) throw new Error('Shader source must start with #version');
  const lines = entries.map(([k, v]) => `#define ${k} ${v}`).join('\n');
  return `${head}${lines}\n#line 2\n${body}`;
}

export function createProgram(gl, vertSrc, fragSrc, { defines = {} } = {}) {
  const prog = gl.createProgram();
  gl.attachShader(prog, compile(gl, gl.VERTEX_SHADER, withDefines(vertSrc, defines)));
  gl.attachShader(prog, compile(gl, gl.FRAGMENT_SHADER, withDefines(fragSrc, defines)));
  gl.linkProgram(prog);
  if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) {
    throw new Error('Program link error:\n' + gl.getProgramInfoLog(prog));
  }
  // Cache uniform locations by name (arrays are looked up as "name[0]").
  const uniforms = {};
  const n = gl.getProgramParameter(prog, gl.ACTIVE_UNIFORMS);
  for (let i = 0; i < n; i++) {
    const info = gl.getActiveUniform(prog, i);
    const name = info.name.replace(/\[0\]$/, '');
    uniforms[name] = gl.getUniformLocation(prog, info.name);
  }
  return { prog, uniforms };
}

// Upload the RGBA float tables as an N x rows RGBA32F texture. Sampled with
// texelFetch + manual lerp, so NEAREST filtering is all that is required.
export function createTableTexture(gl, data, n, rows = 1) {
  const tex = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, tex);
  gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, n, rows, 0, gl.RGBA, gl.FLOAT, data);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  return tex;
}

// Full-screen triangle uses gl_VertexID; WebGL2 still wants a bound VAO.
export function createEmptyVAO(gl) {
  const vao = gl.createVertexArray();
  gl.bindVertexArray(vao);
  return vao;
}
