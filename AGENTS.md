# AGENTS.md

Guidance for coding agents (and new contributors) working in this repository.
`NOTES.md` is the owner's original design brief; this file records how the
project actually works today and how changes are expected to be made.

## What this is

A WebGL2 visualization of atomic electron orbitals for H–Rn (Z = 1–86).
Every frame, each pixel's ray samples a Poisson number of electron
"detections" from the exact column density of hydrogen-like orbitals
(Clementi effective nuclear charges), so the cloud looks like live,
static-y noise. Two atoms can be shown side by side (A/B comparison).

**Hard constraints**

- Pure JavaScript (ES modules) / HTML / CSS. **No libraries, no build step,
  no package manager.** Python tools use the **standard library only**.
- Serve the repo root with the standard server — there is no custom server:

  ```sh
  python3 -m http.server 8000
  # app:    http://localhost:8000/
  # print:  http://localhost:8000/print.html?z=26   (the app's Print… button passes the view)
  # tests:  http://localhost:8000/tests/test.html
  # video:  http://localhost:8000/tests/video.html   (real time; see below)
  # bench:  http://localhost:8000/tests/bench.html?z=36,86   (also &color=blend, &cut=1)
  ```

  `http.server` lets browsers cache. After editing JS or `.glsl` files,
  hard-reload; in headless Chrome use a **fresh `--user-data-dir`**, or a
  stale module can mix with new ones (this has caused false test failures).

## Layout

| Path | Role |
|---|---|
| `index.html`, `css/style.css` | Page and styles (responsive: stacked layout in portrait/narrow, compact strip on short screens). `css/categories.css`: element-category colors, shared with the print page. |
| `print.html`, `js/print.js`, `css/print.css` | Print page: one frame on white at 300 px/in in a box (Letter/A4, grain 1/2/4; a pair fills the width at 4:3, a single atom gets one pane's box, centred), ruled notes above and below; page 2 has the element card, orbital table, radial density (light theme) and model notes. State comes in the URL (`printParams` in `main.js`). |
| `js/main.js` | App state, frame loop, resize, render scale (Auto/Full/Half), spin, selection actions. |
| `js/renderer.js` | WebGL renderer: `renderViews` (one pane per atom), `splitViews`, view geometry (`viewRadius`, `baseHalf`), equalize gains, d/f peak search. |
| `js/gl.js` | WebGL helpers; `createProgram(..., { defines })` injects `#define`s after `#version` plus `#line 2`. |
| `js/shaders.js` | Loads `shaders/*.glsl` (top-level await) and exports `SHADER_DEFINES`. |
| `shaders/vertex.glsl` | Full-screen triangle **and** per-pixel ray geometry (`vP`, `vD`). |
| `shaders/fragment.glsl` | Ordered by importance: inputs → `main()` (physics, then rendering) → 1. physics → 2. rendering helpers → 3. numerical methods (3a. uniform RNG, 3b. Poisson). |
| `js/orbitals.js` | The 43 real orbitals (s/p/d/f), angular forms (p axis, d quadratic `Q`, f cubic tensor `C`), radial functions, chroma, valence rules. |
| `js/elements.js` | Loads `data/elements.json`; occupancy (Hund fill), configurations, `MAX_Z`. |
| `js/selection.js` | The selected atoms A (+ optional B) and everything derived from their union. |
| `js/math.js` | CPU mirrors of the shader (column density, `pcg4d`, Poisson) used by tests, plus quaternions. |
| `js/record.js`, `js/recordui.js`, `js/mux-mp4.js`, `js/mux-webm.js` | Video recording: frame plan, offscreen rendering at a fixed scale (centre crop of the view), WebCodecs encoding, codec levels; the red dot and recording controls; hand-written MP4 (H.264) and WebM (VP9/AV1/VP8) muxers. The UI calls it *Render*; the code says record. |
| `js/ui.js`, `js/radial.js`, `js/tile.js`, `js/controls.js`, `js/view.js`, `js/tables.js` | Panel/overlay, radial-density chart, element tile, input, scale animation, table loader. |
| `tools/gen_tables.py`, `tools/orbitals.py` | Generate the moment tables (`data/orbital_tables*`, incl. the odd cut-away moments). |
| `tools/build_elements.py` | Builds `data/elements.json` from cached, cited sources in `data/sources/`. |
| `tests/` | `test.html` + `tests.js` (browser, incl. GPU), `video.html` + `video.js` (encoders, muxers, recording; real time), `test_*.py` (stdlib unittest), `bench.*` (frame-time probe). |

## Physics and math (don't break these)

- **Column densities from 1-D moment tables.** For each radial function the
  tables hold `I_k(u) = ∫ R² t^k / r^(2l) dt`, k = 0, 2, 4, 6, in the scaled
  variable `u = Z_eff·b/n` (texture: 1024 × rows, RGBA = I0, I2, I4, I6).
  s/p/d/f column densities are exact combinations of these (see the header of
  `tools/orbitals.py`). New radial rows are **appended**; never reorder rows.
- **Cut-away** (front half removed): the odd half-line moments
  `J_k = ∫_0^∞ R² t^k / r^(2l) dt`, k = 1, 3, 5, live in
  `data/orbital_tables_odd.bin` (same layout, RGBA = J1, J3, J5, 0). The shader
  adds `2·uCut·(odd terms)` to the whole-ray column (uCut = −1: the half
  behind the plane, view z < 0, shown at double gain). It works only because
  the plane passes through the nucleus; a movable plane would need 2-D tables.
- **Nearest color** (`uColorMode` 1): per-orbital counts are unchanged; for
  each orbital with k > 0 the nearest of its k detections sits at cumulative
  fraction 1 − v^(1/k) from the front. `nearestDepthRef` finds it by marching
  q2(r)·A(t)² (`data/orbital_tables_radial.bin`, q2 = R²/u^(2l); `angularPoly`
  gives A) over the row's `extent`, with t = H·s³ to resolve the inner lobes
  near the nucleus, and stops once behind the nearest so far. The pixel shows
  only that orbital (its hue, chroma and tone-mapped count). Mirrored in
  `js/math.js`.
- **Ray geometry** lives in the vertex shader; the fragment shader starts at
  the physics. `gl_FragCoord` is used only for the random-number key.
- **Randomness:** counter-based `pcg4d(x, y, frameSeed, stream)`; uniforms
  from the top 24 bits (`[0,1)`, never 1.0). **Poisson:** inversion for
  λ < 10 (with a tail guard), PTRS (Hörmann 1993) for λ ≥ 10. Keep the shader
  and `js/math.js` mirrors in sync.
- **Display:** noise must map 1:1 onto device pixels (ResizeObserver
  `devicePixelContentBoxSize`). Reduced resolution renders 1/k into a
  framebuffer and blits by exactly k; pane splits fall on multiples of k.
  Non-integer browser rescaling of per-pixel noise causes moiré.
- **Color:** hue = angular shape; chroma by depth (valence vivid, inner
  shells paler; s stays vivid). White background is the exact dual of dark
  (negative of the 180°-hue-rotated dark image) — not a bespoke mapping.
- **View scale:** the Rahm atomic radius ring fills 90% of the shorter
  half-dimension; a pair shares one scale (the larger radius). The model's
  hydrogen-like cloud may extend past the ring — that is expected.
- Known model limits are surfaced in the UI, not hidden: non-relativistic
  hydrogen-like orbitals (note shown for Z ≥ 55); partly filled subshells use
  a fixed Hund order (note shown).

## Data and provenance

- Element data is **generated**: edit `tools/build_elements.py` or the cached
  sources, then run `python3 tools/build_elements.py`. Never hand-edit
  `data/elements.json`.
- Every number has a cited source (PubChem periodic table; Clementi & Raimondi
  1963 / Clementi, Raimondi & Reinhardt 1967 Z_eff; Rahm, Hoffmann & Ashcroft
  2016 radii; reviewed overrides for Z ≤ 36). Record new sources in
  `data/sources/` and in the generated JSON's `sources`.
- **Do not invent or estimate values** to fill gaps. If a source lacks a value
  the model needs, ask the owner. Precedent: La, Ce and Gd have no Clementi 5d
  Z_eff, so they are marked `supported: false` and greyed out.
- Check fetched data against independent values before trusting it (e.g. the
  Wikipedia Kr 4p Z_eff was wrong; the Clementi table value replaced it).
- Downloaded material goes in a scratch directory first and is read before use.
- Regenerating tables: `timeout 600 python3 tools/gen_tables.py --verify`
  (Monte Carlo cross-check must PASS). Rows that already existed must stay
  bit-identical.

## Verifying changes

Run all of these for anything non-trivial; report results honestly.

```sh
python3 -m unittest discover tests
python3 -m http.server 8000 &
google-chrome --headless=new --no-sandbox --use-angle=swiftshader \
  --enable-unsafe-swiftshader --user-data-dir="$(mktemp -d)" \
  --virtual-time-budget=200000 --dump-dom http://localhost:8000/tests/test.html \
  | grep -oE '<title>[^<]*|class="fail">[^<]*'      # expect: PASS n/n
# Video tests play media, which needs real time: no --virtual-time-budget.
# Chrome stays open, hence the timeout; the result is a console line.
timeout 240 google-chrome --headless=new --no-sandbox --use-angle=swiftshader \
  --enable-unsafe-swiftshader --user-data-dir="$(mktemp -d)" \
  --enable-logging=stderr --v=0 http://localhost:8000/tests/video.html 2>&1 \
  | grep -m1 -o 'TESTS.*'                          # expect: TESTS PASS n/n
```

- **Recordings:** check real files with `ffprobe -count_frames` (codec, size,
  frame count, duration) and `ffmpeg -v error -i FILE -f null -` (decodes
  cleanly). Quality presets were calibrated by comparing source frames
  (`record({ onFrame })`) with decoded ones; re-check by eye if they change.

- **Screenshots** with `--screenshot` and `--window-size` (desktop 1400×1000,
  1366×768, phone `--window-size=390,844 --force-device-scale-factor=2`).
  Look at them; check overlaps, clipping and both themes.
- **Reproducible frames:** URL options `?z=26&debug=1&paused=1&spin=0&res=full`
  (also `z2=` for a second atom, `mode=white`, `equalize=1`, `scale=fixed`, `cut=1`,
  `color=blend` (default is nearest), `record=1`,
  `density=`). Use these for **pixel-identical before/after comparisons** when a
  change should not alter output (refactors, file moves).
- For numeric/visual claims (noise quality, alignment, performance), measure:
  spectra/autocorrelation of screenshots, GPU readbacks in a test, the bench
  page. In headless Chrome, `performance.now()` stalls under
  `--virtual-time-budget`; capture bench output via `console.log` and
  `--enable-logging=stderr` instead.
- Things headless can't do (long-press, shift-click, real-device fps, touch):
  say so and ask the owner to check.
- When a test fails, find out whether the test or the code is wrong before
  changing either; tune tests only with a stated reason.
- Wrap long-running commands in `timeout`.

## Charts and color

Chart colors are validated with the dataviz palette validator (lightness band,
CVD separation, contrast against the dark panel `#151922`). The radial chart
uses color = shell n and dash = l (s solid, p dashed, d dotted, f dash-dot);
A vs B uses solid/dashed with a validated pair. Re-validate when adding series.

## Working agreements

- **Branch per change**, off the current `main` tip. Commit there; the owner
  reviews, then asks to fast-forward `main` and push. **Never push, merge or
  move `main` unless asked.** Keep history linear (`git merge --ff-only`).
- Tag revisions the owner wants to compare (e.g. `poisson-knuth`,
  `poisson-ptrs`, `display-1to1`, `rng-pcg4d`); tags stay local unless asked.
- Commit messages: a short summary line, then what changed and why, how it was
  verified, and the co-author trailer used in this repo.
- Don't commit generated junk (`__pycache__/`, `*.pyc` are ignored).
- Plan larger changes first (problem, approach, files, verification), ask about
  genuine design decisions (with a recommendation), then implement.
- Prefer clarity and separation of concerns over micro-optimization; keep
  comments accurate when moving code.
- Features may land undocumented at first (e.g. shift-click compare); user
  documentation is written once a feature settles.
