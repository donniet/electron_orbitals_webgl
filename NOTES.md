# Electron Orbitals Rendering in WebGL

Renders atomic electron orbitals in WebGL using a probabilistic ray-tracing
algorithm. Each frame a new random seed is chosen (debug option: seed == frame
number). For each pixel's ray and each orbital, the expected number of electron
detections along the ray is computed from the orbital's probability density,
and a Poisson-distributed count is sampled from it. The counts set the colour
of that pixel, so the cloud looks like live, static-y noise.

Supported elements: hydrogen through radon (Z = 1–86), with s, p, d and f
orbitals. La, Ce and Gd are shown but greyed out: the published effective
nuclear charges lack a value their ground states need (see AGENTS.md).

# Implementation

Orbitals are hydrogen-like, using Clementi effective nuclear charges, in a
fixed orientation. Exact column densities come from pre-calculated 1-D moment
tables (numerical integration, verified against a Monte Carlo simulation), so
the shader stays cheap. Element data (configurations, charges, radii, and so
on) is generated from cited sources into `data/elements.json`.

# Rendering

Each orbital is assigned a hue: one per angular shape, paler for inner shells.
Each orbital is calculated independently in the shader; detection counts go
through a tone curve, the per-orbital colours are added together, and the
result is converted to RGB. On the dark background detections light up; the
white background is its exact "ink on paper" dual. The view is scaled so the
element's Rahm atomic radius (where the free atom's electron density falls to
0.001 e/bohr³) fills about 90% of the view, and that radius is drawn as a ring.

The user can rotate the orientation by clicking and dragging (or pinch-zoom on
touch). The atom spins slowly at start, and a flick keeps it spinning. Each
time step is a new seed, so the cloud appears static-y. Two atoms can be shown
side by side for comparison (shift-click, long-press or the ⇄ button).

The nucleus is not drawn. The UI shows the atom being rendered, its orbitals
and their hues, its periodic-table information, and a radial-density chart.

# Notes

This is written in pure JavaScript/HTML/CSS with no libraries or dependencies.
Python (standard library only) is used for the dev server, the data and table
generators, and tests. See AGENTS.md for the architecture, conventions and
verification steps.

# Running

Serve the repository root with Python's built-in server:

    python3 -m http.server 8000

Then open http://localhost:8000/ for the app and
http://localhost:8000/tests/test.html for the in-browser tests. Python tests:
`python3 -m unittest discover tests`. Regenerate the orbital tables with
`python3 tools/gen_tables.py --verify`, and the element data with
`python3 tools/build_elements.py`.

http.server does not send no-cache headers, so hard-reload (Ctrl+Shift+R)
after editing files (including the shaders) if the browser shows stale code.

The GLSL lives in `shaders/vertex.glsl` and `shaders/fragment.glsl` and is
loaded by `js/shaders.js`. Build-time constants such as `NUM_ORBITALS` are not
written in the .glsl files: `createProgram` in `js/gl.js` injects them as
`#define`s right after the `#version` line (followed by `#line 2`, so compiler
errors still point at the right line of the file).
