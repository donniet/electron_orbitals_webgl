# Electron Orbitals Rendering in WebGL

Renders atomic electron orbitals in WebGL using a probabalistic ray-tracing 
algorithm. Each frame a random seed is chosen (DEBUG: seed == frame_number). 
Then for each ray the probability of detecting one or more electrons along the 
ray is calculated. The random seed is used to choose a random number and the
number of electrons detected along that ray is generated. This value sets
the saturation of the "electron color" for the fragment of that ray.

We will support a limited number of elements.  For now just the s and p orbitals
up to Neon (Ne 10). 

# Implementation

We will use the proper electron orbital probabilities for the rendering, but
assume a fixed orientation. Helper tables for each orbital are pre-calculated
via simulation to speed up the computation.

# Rendering

Each orbital is assigned a hue and the saturation is set by the number of
electrons in that orbital detected.  Each orbital is calculated independently in
the shader and the results added together, then converted to RGB.  The horizon 
for 2 standard deviations of electron positions is calculated and used to scale
the visualization on the screen.

The user can rotate the orientation by clicking and draging.  Each time step is 
a new seed, and so the cloud should appear static-y. 

The nucleus for each atom can be skipped for now, but the UI should clearly show
the atom being rendered, the orbital hues, and the period table information for
that atom.

# Notes

This should be written in pure javascript/HTML/CSS with no libraries or
dependencies and a simple python server used for the dev server and any testing.

# Running

Serve the repository root with Python's built-in server:

    python3 -m http.server 8000

Then open http://localhost:8000/ for the app and
http://localhost:8000/tests/test.html for the in-browser tests. Python tests:
`python3 -m unittest discover tests`. Regenerate the orbital tables with
`python3 tools/gen_tables.py --verify`.

http.server does not send no-cache headers, so hard-reload (Ctrl+Shift+R)
after editing files (including the shaders) if the browser shows stale code.

The GLSL lives in `shaders/vertex.glsl` and `shaders/fragment.glsl` and is
loaded by `js/shaders.js`. Build-time constants such as `NUM_ORBITALS` are not
written in the .glsl files: `createProgram` in `js/gl.js` injects them as
`#define`s right after the `#version` line (followed by `#line 2`, so compiler
errors still point at the right line of the file).