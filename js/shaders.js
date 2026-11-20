// Loads the GLSL ES 3.00 sources from shaders/*.glsl.
//
// The fragment shader expects host-provided defines (see SHADER_DEFINES),
// which js/gl.js createProgram() injects right after the #version line.

import { NUM_ORBITALS } from './orbitals.js';

export { NUM_ORBITALS };

export const SHADER_DEFINES = { NUM_ORBITALS };

async function loadText(rel) {
  const url = new URL(rel, import.meta.url);
  const res = await fetch(url);
  if (!res.ok) throw new Error(`Could not load shader ${url.pathname} (${res.status})`);
  return res.text();
}

export const [VERT, FRAG] = await Promise.all([
  loadText('../shaders/vertex.glsl'),
  loadText('../shaders/fragment.glsl'),
]);
