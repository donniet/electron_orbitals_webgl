# Electron Orbitals

A live WebGL2 picture of the electron cloud of any atom from hydrogen to radon
(Z = 1–86). Each frame, every pixel looks along its line of sight and counts
simulated electron "detections". The counts are drawn at random from the
orbital's probability density |ψ|². The cloud flickers like static, and its
average over time is the true density.

It is plain JavaScript, HTML and CSS: no libraries, no build step, nothing to
install.

## Running it

Serve the repository root with Python's built-in web server:

```sh
python3 -m http.server 8000
```

Then open <http://localhost:8000/>. You need a browser with WebGL2 (any
current Chrome, Firefox, Safari or Edge).

The server lets the browser cache files. If you edit the code, hard-reload
(Ctrl+Shift+R) to make sure the new version loads.

## Using it

- **Choose an element** in the periodic table at the top, or step through
  them with ← / →. La, Ce and Gd are greyed out (see *Limits* below).
- **Rotate** by dragging. A flick keeps the atom spinning. **Zoom** with the
  scroll wheel or a pinch.
- **Compare two atoms** side by side: press ⇄ on the element card and pick a
  second element, or shift-click (long-press on touch) an element in the
  table. ⇆ swaps them; × or Esc goes back to one atom.
- **Orbitals table:** show or hide whole shells, subshells, or single orbitals.
  Click ▸ to expand a row. *Valence* shows only the outer shell (plus an open
  d subshell). The table also lists each subshell's electron count,
  effective nuclear charge Z<sub>eff</sub> and RMS radius.
- **Radial density chart:** electrons per bohr radius against distance from
  the nucleus, per shell or per subshell. When two atoms are shown, *A vs B*
  compares their totals. Hover the chart to read values.
- **Rendering controls:** exposure (mean detections per pixel), *Equalize
  orbitals* (brightens faint outer orbitals; densities are then no longer to
  scale), *Fixed scale* (the same zoom for every element, so sizes can be
  compared), the tone curve, a dark or white background, and noise
  resolution.
- **Keys:** space pause · r reset view · s save PNG · ←/→ element ·
  Esc stop comparing.

### Reading the picture

- **Hue shows the orbital's shape** (s, each p, each d, each f orientation).
  Inner shells are drawn paler, so the outer (valence) electrons stand out.
- **The ring is the atom's size:** the Rahm atomic radius, where the free
  atom's electron density falls to 0.001 electrons per bohr³. The view is
  zoomed so this ring fills about 90% of the view. Some of the modelled cloud
  can reach past the ring; that is expected.
- By default, brightness is to scale, so the dense inner shells dominate.
  Turn on *Equalize orbitals* to see the faint outer orbitals clearly.
- The nucleus is not drawn.

### URL options

Settings can be given in the address, which is handy for sharing a view or
for reproducible screenshots:

| Option | Meaning |
|---|---|
| `z=26` | Element by atomic number. |
| `z2=8` | Second element, shown side by side. |
| `mode=white` | White background (default `dark`). |
| `equalize=1` | Equalize orbital brightness. |
| `scale=fixed` | Same view scale for every element. |
| `density=0.5` | Exposure. |
| `kmax=4`, `gamma=0.5` | Tone curve: counts for full brightness, and gamma. |
| `res=full` | Noise resolution: `auto` (default), `full` or `half`. |
| `paused=1` | Start paused. |
| `spin=0` | Start without the slow rotation. |
| `debug=1` | Seed = frame number, so frames are reproducible. |

Example: `http://localhost:8000/?z=26&z2=44&equalize=1`

## How it works

- **The model.** Each electron is in a hydrogen-like orbital whose nuclear
  charge is the effective charge Z<sub>eff</sub> (Clementi & Raimondi). This
  includes s, p, d and f orbitals in real (not complex) form. Electrons fill
  the subshells in the element's ground-state configuration.
- **Exact densities along each ray.** The amount of each orbital along a line
  of sight is an integral. It is computed exactly from precomputed 1-D tables
  (`data/orbital_tables.*`, made by `tools/gen_tables.py` and checked against
  a Monte Carlo simulation). That keeps the shader fast.
- **Random detections.** Every frame, each pixel draws a Poisson-distributed
  count from that expected value, using a counter-based random-number
  generator. The counts are coloured by orbital, added, and tone-mapped.
- **One sample per screen pixel.** The noise is drawn at the display's real
  pixel resolution (or an exact half), so it stays crisp without moiré.

`AGENTS.md` describes the architecture and conventions in detail. `NOTES.md`
is the original design brief.

## Data sources

All element data is generated from cited sources by
`tools/build_elements.py` into `data/elements.json`. The cached source files
are in `data/sources/`.

- Names, masses, configurations, electronegativity, ionization energy:
  [PubChem periodic table](https://pubchem.ncbi.nlm.nih.gov/periodic-table/).
- Effective nuclear charges: Clementi & Raimondi, *J. Chem. Phys.* **38**,
  2686 (1963); Clementi, Raimondi & Reinhardt, *J. Chem. Phys.* **47**, 1300
  (1967).
- Atomic radii: Rahm, Hoffmann & Ashcroft, *Chem. Eur. J.* **22**, 14625
  (2016); corrigendum **23**, 4017 (2017).
- Reviewed corrections for Z ≤ 36: `data/sources/overrides.json`.

## Limits

These are shown in the app as well:

- **Not relativistic.** For heavy atoms (Z ≥ 55), the real contraction or
  expansion of orbitals such as 6s and 5d is not included.
- **Partly filled subshells** use one fixed filling order (Hund's rule). Other
  equally valid fillings would look different.
- **La, Ce and Gd are unavailable.** The published Z<sub>eff</sub> tables have
  no value for their 5d electron, and the app does not estimate one.
- Hydrogen-like orbitals are a simple model. They show the shapes and the
  shell structure well, but are not exact many-electron densities.

## Development

```sh
python3 -m unittest discover tests        # Python tests (data, tables)
# with the server running:
#   http://localhost:8000/tests/test.html   browser and GPU tests
#   http://localhost:8000/tests/bench.html  frame-time benchmark
python3 tools/build_elements.py           # regenerate data/elements.json
timeout 600 python3 tools/gen_tables.py --verify   # regenerate orbital tables
```

The Python tools use only the standard library. See `AGENTS.md` for the
project layout, the rules for changing physics and data, and how to verify
changes.

## License

Copyright © 2026 Donnie Tolbert.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License, version 3 or (at your
option) any later version. See [`LICENSE`](LICENSE).

The AGPL also covers use over a network: if you host a modified version for
others to use, you must offer them its source code.
