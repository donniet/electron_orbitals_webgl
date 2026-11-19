"""Tests for the general (periods 1-4) moment tables and tools/orbitals.py.

Checks the checked-in data/orbital_tables.bin + .json, so stale or corrupted
output is caught. Run: python3 -m unittest discover tests
"""
import math
import os
import random
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import gen_tables as gt  # noqa: E402
import orbitals as orb  # noqa: E402

TABLES, META = gt.load_moment_tables()
SHAPES = orb.all_shapes()
BY_LABEL = {s.label: s for s in SHAPES}


def moments(row):
    return gt.moments_fn(TABLES, row)


def simpson(f, a, b, n=4000):
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4 if i & 1 else 2) * f(a + i * h)
    return s * h / 3


def random_ray(rng, b_hi):
    while True:
        d = [rng.gauss(0, 1) for _ in range(3)]
        n = math.sqrt(orb.dot(d, d))
        if n > 1e-6:
            d = [c / n for c in d]
            break
    while True:
        q = [rng.gauss(0, 1) for _ in range(3)]
        t = orb.dot(q, d)
        q = [a - t * b for a, b in zip(q, d)]
        n = math.sqrt(orb.dot(q, q))
        if n > 1e-6:
            break
    b = rng.uniform(0, b_hi)
    return [c / n * b for c in q], d


class TestMetadata(unittest.TestCase):
    def test_rows_match_orbitals_module(self):
        self.assertEqual([r["name"] for r in META["rows"]], orb.ROW_NAMES)
        self.assertEqual([(r["n"], r["l"]) for r in META["rows"]], orb.RADIAL)
        self.assertEqual(META["N"], len(TABLES[0]))
        self.assertEqual(META["uMax"], gt.M_U_MAX)

    def test_binary_size(self):
        size = os.path.getsize(gt.BIN_PATH)
        self.assertEqual(size, len(orb.RADIAL) * META["N"] * 4 * 4)

    def test_non_negative_and_tail(self):
        for row in TABLES:
            for v in row:
                for x in v:
                    self.assertGreaterEqual(x, 0.0)
            for x in row[-1]:
                self.assertLess(x, 1e-10)


class TestRadialFunctions(unittest.TestCase):
    def test_laguerre_low_orders(self):
        self.assertEqual(orb.laguerre(0, 3), [1.0])
        self.assertEqual(orb.laguerre(1, 3), [4.0, -1.0])
        # L_2^a(x) = ((a+1)(a+2) - 2(a+2) x + x^2) / 2
        for a, c in zip(orb.laguerre(2, 1), [3.0, -3.0, 0.5]):
            self.assertAlmostEqual(a, c)

    def test_normalised_and_mean_r2(self):
        for rad in orb.RADIALS:
            n, l = rad.n, rad.l
            self.assertAlmostEqual(simpson(lambda u: u * u * rad.R(u) ** 2, 0, 60), 1.0, places=6,
                                   msg=f"{n}{l}")
            self.assertAlmostEqual(simpson(lambda u: u ** 4 * rad.R(u) ** 2, 0, 60),
                                   orb.mean_r2(n, l, n), places=6, msg=f"<r2> {n}{l}")

    def test_radial_node_counts(self):
        for rad in orb.RADIALS:
            us = [0.01 + 0.01 * i for i in range(3000)]
            signs = [rad.q(u) > 0 for u in us]
            nodes = sum(1 for a, b in zip(signs, signs[1:]) if a != b)
            self.assertEqual(nodes, rad.n - rad.l - 1, f"{rad.n}{rad.l}")


class TestAngular(unittest.TestCase):
    def test_each_shape_normalised_on_sphere(self):
        nt, nphi = 200, 200
        for s in SHAPES:
            acc = 0.0
            for i in range(nt):
                th = math.pi * (i + 0.5) / nt
                for j in range(nphi):
                    ph = 2 * math.pi * (j + 0.5) / nphi
                    v = [math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)]
                    acc += s.angular2(v) * math.sin(th)
            acc *= (math.pi / nt) * (2 * math.pi / nphi)
            self.assertAlmostEqual(acc, 1.0, delta=2e-4, msg=s.label)  # midpoint-rule error ~5e-5

    def test_d_eigenvalues(self):
        s15 = math.sqrt(15 / (4 * math.pi))
        ev = sorted(orb.sym_eigvals(orb.D_SHAPES["xy"]))
        for a, b in zip(ev, [-s15 / 2, 0, s15 / 2]):
            self.assertAlmostEqual(a, b)
        z = math.sqrt(5 / (16 * math.pi))
        ev = sorted(orb.sym_eigvals(orb.D_SHAPES["z2"]))
        for a, b in zip(ev, [-z, -z, 2 * z]):
            self.assertAlmostEqual(a, b)


class TestClosedForms(unittest.TestCase):
    def test_b_zero(self):
        # 1s: R^2 = 4 e^{-2r}; 2p (Z=2): R^2 / r^2 = (4/3) e^{-2r}
        want = {"1s": (4, 2, 6), "2p": (4 / 3, 2 / 3, 2)}
        for name, vals in want.items():
            row = orb.ROW_NAMES.index(name)
            for got, w in zip(TABLES[row][0], vals):
                self.assertAlmostEqual(got, w, places=5, msg=name)


class TestPlaneNormalisation(unittest.TestCase):
    """Each orbital integrated over the image plane holds one electron."""

    def test_all_shapes_tilted_view(self):
        d = [0.36, 0.48, 0.8]
        e1, e2 = gt._basis(d)
        n = META["N"]
        h = gt.M_U_MAX / (n - 1)
        nphi = 48
        for s in SHAPES:
            mom = moments(s.row)
            acc = 0.0
            for i in range(n):
                b = i * h
                w = 0.5 if i in (0, n - 1) else 1.0
                ring = 0.0
                for m in range(nphi):
                    ph = 2 * math.pi * (m + 0.5) / nphi
                    p = [b * (math.cos(ph) * e1[k] + math.sin(ph) * e2[k]) for k in range(3)]
                    ring += s.column_ref(mom, p, d)
                acc += w * ring * b
            acc *= h * 2 * math.pi / nphi
            self.assertAlmostEqual(acc, 1.0, delta=2e-3, msg=s.label)


class TestBruteForceRays(unittest.TestCase):
    """Table formula == direct numerical line integral of |psi|^2."""

    def test_random_rays_every_shape(self):
        rng = random.Random(2024)
        for s in SHAPES:
            for z in (1.3, 6.7):
                k = z / s.n
                rms = math.sqrt(orb.mean_r2(s.n, s.l, z))
                for _ in range(2):
                    p, d = random_ray(rng, 1.5 * rms)
                    tmax = 40.0 / k
                    want = simpson(lambda t: s.density(z, [p[i] + t * d[i] for i in range(3)]),
                                   -tmax, tmax, 20000)
                    got = s.column(moments(s.row), z, p, d)
                    scale = k * k * max(TABLES[s.row][0])  # typical magnitude
                    self.assertAlmostEqual(got, want, delta=2e-3 * max(want, 1e-2 * scale) + 1e-9,
                                           msg=f"{s.label} z={z} p={p} d={d}")


class TestUnsold(unittest.TestCase):
    """Filled (or half-filled) p and d subshells are spherically symmetric."""

    def check(self, shapes, rng):
        for b in (0.0, 0.7, 2.0, 5.0):
            vals = []
            for _ in range(6):
                _, d = random_ray(rng, 0)
                e1, _ = gt._basis(d)          # any direction perpendicular to d
                p = [b * c for c in e1]
                vals.append(sum(sh.column_ref(moments(sh.row), p, d) for sh in shapes))
            spread = max(vals) - min(vals)
            self.assertLess(spread, 1e-6 + 1e-5 * max(vals), msg=f"{[s.label for s in shapes]} b={b}")

    def test_p_and_d_subshells(self):
        rng = random.Random(5)
        for n, l in orb.RADIAL:
            if l == 0:
                continue
            self.check([s for s in SHAPES if (s.n, s.l) == (n, l)], rng)


class TestLegacyConsistency(unittest.TestCase):
    """New 1s/2s/2p tables reproduce the legacy H-Ne tables used by the app."""

    def test_matches_legacy(self):
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from test_tables import col_from_tables  # legacy formula
        rng = random.Random(9)
        for lab, legacy in (("1s", "1s"), ("2s", "2s"), ("2px", "2px"), ("2pz", "2pz")):
            s = BY_LABEL[lab]
            for z in (1.0, 3.2, 5.758):
                for _ in range(3):
                    p, d = random_ray(rng, 4.0 / z)
                    got = s.column(moments(s.row), z, p, d)
                    want = col_from_tables(legacy, p, d, z)
                    self.assertAlmostEqual(got, want, delta=1e-3 * max(want, 1e-3 * z * z) + 1e-8,
                                           msg=f"{lab} z={z}")


class TestMonteCarlo(unittest.TestCase):
    def test_sampled_histograms_match_tables(self):
        for i, (lab, d) in enumerate([("3s", [0, 0, 1]), ("4px", [0, 0.6, 0.8]),
                                      ("3dxy", [0, 0, 1]), ("3dz2", [0.6, 0, 0.8])]):
            err, _ = gt.mc_shape_check(TABLES, BY_LABEL[lab], d, n_samples=20000, seed=100 + i)
            self.assertLess(err, 0.012, lab)


if __name__ == "__main__":
    unittest.main()
