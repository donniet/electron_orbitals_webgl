#!/usr/bin/env python3
"""Generate the orbital column-density moment tables (stdlib only).

For each radial function (1s 2s 3s 4s 2p 3p 4p 3d) the moments

    I_k(u) = int R(r)^2 t^k / r^(2l) dt,   r = sqrt(u^2 + t^2),  k = 0, 2, 4

are tabulated in the scaled variable u = Z b / n (reference charge Z0 = n),
1024 samples on [0, 30]. Output:

    data/orbital_tables.bin   float32 little-endian, [row][i][RGBA]
    data/orbital_tables.json  metadata (N, uMax, rows, channels)

Each texel holds (I0, I2, I4, 0). tools/orbitals.py has the math that turns
moments into column densities for s, p and d orbitals.

`--verify` runs a Monte Carlo simulation: electron positions are sampled
from |psi|^2, projected onto an image plane and histogrammed; the histograms
are compared to the tables.
"""
import argparse
import array
import json
import math
import multiprocessing
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orbitals as orb  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ================================================================ moment tables
M_N = 1024            # samples per row
M_U_MAX = 30.0        # scaled impact parameter range u = Z b / n
M_T_MAX = 40.0        # integrate t in [0, M_T_MAX] (exp(-80) beyond)
M_T_STEPS = 4000      # Simpson intervals (even)
BIN_PATH = os.path.join(ROOT, "data", "orbital_tables.bin")
META_PATH = os.path.join(ROOT, "data", "orbital_tables.json")


def moments_at(u):
    """[[I0, I2, I4] per radial row] at scaled impact parameter u."""
    h = M_T_MAX / M_T_STEPS
    rads = orb.RADIALS
    acc = [[0.0, 0.0, 0.0] for _ in rads]
    polys = [(r.norm, r.poly) for r in rads]
    for i in range(M_T_STEPS + 1):
        t = i * h
        w = 1.0 if i in (0, M_T_STEPS) else (4.0 if i & 1 else 2.0)
        r = math.sqrt(u * u + t * t)
        we = w * math.exp(-2.0 * r)
        t2 = t * t
        for row, (norm, poly) in enumerate(polys):
            v = 0.0
            for c in reversed(poly):
                v = v * r + c
            q2 = we * norm * norm * v * v
            a = acc[row]
            a[0] += q2
            a[1] += q2 * t2
            a[2] += q2 * t2 * t2
    f = 2.0 * h / 3.0     # Simpson weight x2 for the symmetric half-line
    return [[v * f for v in a] for a in acc]


def build_moment_tables(n=M_N, u_max=M_U_MAX, processes=None):
    """tables[row][i] = (I0, I2, I4) at u_i = i * u_max / (n - 1)."""
    us = [i * u_max / (n - 1) for i in range(n)]
    with multiprocessing.Pool(processes) as pool:
        per_u = pool.map(moments_at, us, chunksize=16)
    return [[per_u[i][row] for i in range(n)] for row in range(len(orb.RADIALS))]


def write_moment_tables(tables, bin_path=BIN_PATH, meta_path=META_PATH):
    n = len(tables[0])
    arr = array.array("f")
    for row in tables:
        for i0, i2, i4 in row:
            arr.extend((i0, i2, i4, 0.0))
    if sys.byteorder == "big":
        arr.byteswap()
    os.makedirs(os.path.dirname(bin_path), exist_ok=True)
    with open(bin_path, "wb") as fh:
        arr.tofile(fh)
    meta = {
        "generator": "tools/gen_tables.py",
        "format": "float32 little-endian, layout [row][i][channel], RGBA texel per sample",
        "N": n,
        "uMax": M_U_MAX,
        "channels": ["I0", "I2", "I4", "unused"],
        "rows": [{"name": nm, "n": nn, "l": ll}
                 for nm, (nn, ll) in zip(orb.ROW_NAMES, orb.RADIAL)],
    }
    with open(meta_path, "w") as fh:
        json.dump(meta, fh, indent=2)
        fh.write("\n")


def load_moment_tables(bin_path=BIN_PATH, meta_path=META_PATH):
    with open(meta_path) as fh:
        meta = json.load(fh)
    arr = array.array("f")
    with open(bin_path, "rb") as fh:
        arr.frombytes(fh.read())
    if sys.byteorder == "big":
        arr.byteswap()
    n, rows = meta["N"], len(meta["rows"])
    assert len(arr) == rows * n * 4, (len(arr), rows, n)
    tables = [[(arr[(r * n + i) * 4], arr[(r * n + i) * 4 + 1], arr[(r * n + i) * 4 + 2])
               for i in range(n)] for r in range(rows)]
    return tables, meta


def moment_lookup(tables, row, k, u, u_max=M_U_MAX):
    """Linear interpolation of I_k for `row` at u (0 beyond the table)."""
    t = tables[row]
    n = len(t)
    x = u / u_max * (n - 1)
    if x >= n - 1:
        return 0.0
    i = int(x)
    f = x - i
    ch = k // 2
    return t[i][ch] * (1.0 - f) + t[i + 1][ch] * f


def moments_fn(tables, row):
    return lambda k, b: moment_lookup(tables, row, k, b)


# ---------------------------------------------------------------- Monte Carlo
def _iso(rng, r):
    c = 2.0 * rng.random() - 1.0
    s = math.sqrt(1.0 - c * c)
    phi = 2.0 * math.pi * rng.random()
    return [r * s * math.cos(phi), r * s * math.sin(phi), r * c]


class RadialSampler:
    """Inverse-CDF sampler for the radial pdf u^2 R_ref(u)^2 (reference units)."""

    def __init__(self, radial, u_max=40.0, n=20000):
        self.us = [u_max * i / n for i in range(n + 1)]
        cdf = [0.0]
        prev = 0.0
        for i in range(1, n + 1):
            u = self.us[i]
            cur = u * u * radial.R(u) ** 2
            cdf.append(cdf[-1] + 0.5 * (prev + cur) * (u_max / n))
            prev = cur
        self.cdf = [c / cdf[-1] for c in cdf]

    def sample(self, rng):
        y = rng.random()
        lo, hi = 0, len(self.cdf) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if self.cdf[mid] < y:
                lo = mid
            else:
                hi = mid
        c0, c1 = self.cdf[lo], self.cdf[hi]
        f = 0.0 if c1 == c0 else (y - c0) / (c1 - c0)
        return self.us[lo] + f * (self.us[hi] - self.us[lo])


def sample_shape(shape, rng, radial_sampler, amax):
    """Electron position in reference units from |psi|^2 of `shape`."""
    r = radial_sampler.sample(rng)
    while True:
        v = _iso(rng, 1.0)
        if rng.random() * amax <= shape.angular2(v):
            return [c * r for c in v]


def _basis(d):
    a = [1.0, 0.0, 0.0] if abs(d[0]) < 0.9 else [0.0, 1.0, 0.0]
    e1 = [a[1] * d[2] - a[2] * d[1], a[2] * d[0] - a[0] * d[2], a[0] * d[1] - a[1] * d[0]]
    n1 = math.sqrt(orb.dot(e1, e1))
    e1 = [c / n1 for c in e1]
    e2 = [d[1] * e1[2] - d[2] * e1[1], d[2] * e1[0] - d[0] * e1[2], d[0] * e1[1] - d[1] * e1[0]]
    return e1, e2


def ring_masses_from_tables(tables, shape, d, bin_w, nbins, sub=20, nphi=48):
    e1, e2 = _basis(d)
    mom = moments_fn(tables, shape.row)
    out = []
    h = bin_w / sub
    for k in range(nbins):
        acc = 0.0
        for j in range(sub):
            b = k * bin_w + (j + 0.5) * h
            ring = 0.0
            for m in range(nphi):
                ph = 2 * math.pi * (m + 0.5) / nphi
                p = [b * (math.cos(ph) * e1[i] + math.sin(ph) * e2[i]) for i in range(3)]
                ring += shape.column_ref(mom, p, d)
            acc += ring * (2 * math.pi / nphi) * b
        out.append(acc * h)
    return out


def mc_shape_check(tables, shape, d, n_samples=40000, seed=1, bin_w=None, nbins=16):
    """Max |MC histogram - table| of projected impact parameter (ref units)."""
    rng = random.Random(seed)
    sampler = RadialSampler(shape.radial)
    amax = shape.angular2_max()
    rms = math.sqrt(orb.mean_r2(shape.n, shape.l, shape.n))
    bin_w = bin_w or rms / 6
    counts = [0] * nbins
    for _ in range(n_samples):
        x = sample_shape(shape, rng, sampler, amax)
        t = orb.dot(x, d)
        b = math.sqrt(max(0.0, orb.dot(x, x) - t * t))
        kk = int(b / bin_w)
        if kk < nbins:
            counts[kk] += 1
    mc = [c / n_samples for c in counts]
    tb = ring_masses_from_tables(tables, shape, d, bin_w, nbins)
    return max(abs(a - b) for a, b in zip(mc, tb)), sum(tb)


MC_CASES = [  # (shape label, view direction)
    ("1s", [0, 0, 1]), ("3s", [0, 0, 1]), ("4s", [0.6, 0, 0.8]),
    ("3pz", [0, 0, 1]), ("4px", [0, 0.6, 0.8]), ("2py", [1, 0, 0]),
    ("3dxy", [0, 0, 1]), ("3dz2", [0.6, 0, 0.8]), ("3dx2-y2", [0.36, 0.48, 0.8]),
]


def monte_carlo_verify_moments(tables, n_samples=60000, verbose=True):
    shapes = {s.label: s for s in orb.all_shapes()}
    worst = 0.0
    for i, (lab, d) in enumerate(MC_CASES):
        err, mass = mc_shape_check(tables, shapes[lab], d, n_samples, seed=10 + i)
        worst = max(worst, err)
        if verbose:
            print(f"  MC {lab:8s} d={d}  max |hist - table| = {err:.4f}  (mass in bins {mass:.3f})")
    return worst


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--verify", action="store_true",
                    help="run Monte Carlo cross-check after generating")
    args = ap.parse_args(argv)

    print(f"Building moment tables: {len(orb.RADIALS)} radial functions x {M_N} samples, "
          f"u in [0, {M_U_MAX}] ...")
    tables = build_moment_tables()
    write_moment_tables(tables)
    print(f"Wrote {BIN_PATH} and {META_PATH}")
    for name, row in zip(orb.ROW_NAMES, tables):
        print(f"  {name}: I0(0)={row[0][0]:.6g}  I2(0)={row[0][1]:.6g}  I4(0)={row[0][2]:.6g}")

    if args.verify:
        print("Monte Carlo verification ...")
        worst = monte_carlo_verify_moments(tables)
        ok = worst < 0.01
        print("PASS" if ok else "FAIL", f"(worst bin error {worst:.4f})")
        return 0 if ok else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
