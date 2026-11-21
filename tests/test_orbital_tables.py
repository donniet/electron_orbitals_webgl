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
ODD = gt.load_odd_tables()
SHAPES = orb.all_shapes()
BY_LABEL = {s.label: s for s in SHAPES}


def moments(row):
    return gt.moments_fn(TABLES, row, ODD)


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
            # negligible at u = uMax relative to each channel's peak
            for ch, x in enumerate(row[-1]):
                peak = max(v[ch] for v in row)
                self.assertLess(x, 1e-10 * max(peak, 1.0))


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


class TestFShapes(unittest.TestCase):
    def test_seven_f_shapes_with_symmetric_tensors(self):
        fs = [sh for sh in SHAPES if sh.l == 3]
        self.assertEqual(len(fs), 7)
        for sh in fs:
            C = sh.C
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        self.assertAlmostEqual(C[i][j][k], C[j][i][k])
                        self.assertAlmostEqual(C[i][j][k], C[i][k][j])
            # traceless (harmonic): sum_i C[i][i][k] = 0 for every k
            for k in range(3):
                self.assertAlmostEqual(sum(C[i][i][k] for i in range(3)), 0.0, places=12, msg=sh.label)

    def test_f_shapes_orthogonal(self):
        fs = [sh for sh in SHAPES if sh.l == 3]
        nt, nphi = 120, 120
        pts = []
        for i in range(nt):
            th = math.pi * (i + 0.5) / nt
            for j in range(nphi):
                ph = 2 * math.pi * (j + 0.5) / nphi
                pts.append(([math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)],
                            math.sin(th) * (math.pi / nt) * (2 * math.pi / nphi)))
        for a in range(len(fs)):
            for b in range(a + 1, len(fs)):
                acc = sum(w * orb.cform(fs[a].C, v, v, v) * orb.cform(fs[b].C, v, v, v) for v, w in pts)
                self.assertAlmostEqual(acc, 0.0, delta=1e-3, msg=f"{fs[a].label}.{fs[b].label}")


class TestClosedForms(unittest.TestCase):
    def test_b_zero(self):
        # 1s: R^2 = 4 e^{-2r}; 2p (Z=2): R^2 / r^2 = (4/3) e^{-2r}
        want = {"1s": (4, 2, 6, 45), "2p": (4 / 3, 2 / 3, 2, 15)}
        for name, vals in want.items():
            row = orb.ROW_NAMES.index(name)
            for got, w in zip(TABLES[row][0], vals):
                self.assertAlmostEqual(got, w, places=5, msg=name)


class TestOddClosedForms(unittest.TestCase):
    def test_b_zero(self):
        # J_k(0) = int_0^inf R^2 t^k / t^(2l) dt; 1s: 4 e^{-2t}, 2p: (4/3) e^{-2t}
        want = {"1s": (1, 1.5, 7.5), "2p": (1 / 3, 0.5, 2.5)}
        for name, vals in want.items():
            row = orb.ROW_NAMES.index(name)
            for got, w in zip(ODD[row][0], vals):
                self.assertAlmostEqual(got, w, places=5, msg=name)

    def test_metadata_and_size(self):
        self.assertEqual(META["odd"]["channels"], ["J1", "J3", "J5", "0"])
        path = os.path.join(ROOT, "data", META["odd"]["file"])
        self.assertEqual(os.path.getsize(path), len(META["rows"]) * META["N"] * 16)
        for row in ODD:
            self.assertTrue(all(t[3] == 0.0 for t in row))


class TestRadialTable(unittest.TestCase):
    """q2 = R^2 / u^(2l) table and row extents (nearest-detection coloring)."""

    def test_values_and_closed_forms(self):
        q2 = gt.load_radial_tables()
        n = META["N"]
        self.assertEqual(len(q2), len(META["rows"]))
        for r, rad in enumerate(orb.RADIALS):
            for i in (0, 17, 300, 900):
                u = i * gt.M_U_MAX / (n - 1)
                self.assertAlmostEqual(q2[r][i], rad.R2_over_r2l(u),
                                       delta=1e-6 * max(1.0, rad.R2_over_r2l(u)), msg=f"row {r} i {i}")
        self.assertAlmostEqual(q2[orb.ROW_NAMES.index("1s")][0], 4.0, places=5)
        self.assertAlmostEqual(q2[orb.ROW_NAMES.index("2p")][0], 4 / 3, places=5)

    def test_extents_hold_all_but_the_tail(self):
        for row, rad in zip(META["rows"], orb.RADIALS):
            ext = row["extent"]
            self.assertTrue(0 < ext <= gt.M_U_MAX, row)
            tail = simpson(lambda u: u * u * rad.R(u) ** 2, ext, gt.M_U_MAX + 10, 4000)
            self.assertLess(tail, 2 * META["radial"]["tail"], row["name"])
            self.assertGreater(tail, 0.2 * META["radial"]["tail"], row["name"])


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

    def test_cut_away_halves(self):
        # The density is symmetric through the nucleus, so each half holds
        # half an electron: the doubled cut-away column integrates to 1.
        d = [0.48, 0.6, 0.64]
        e1, e2 = gt._basis(d)
        n = META["N"]
        h = gt.M_U_MAX / (n - 1)
        nphi = 24
        for s in SHAPES:
            if s.l == 0:
                continue
            mom = moments(s.row)
            for cut in (-1, 1):
                acc = 0.0
                for i in range(0, n, 2):
                    b = i * h
                    w = 0.5 if i in (0, n - 1) else 1.0
                    ring = 0.0
                    for m in range(nphi):
                        ph = 2 * math.pi * (m + 0.5) / nphi
                        p = [b * (math.cos(ph) * e1[k] + math.sin(ph) * e2[k]) for k in range(3)]
                        ring += s.column_ref(mom, p, d, cut)
                    acc += w * ring * b
                acc *= 2 * h * 2 * math.pi / nphi
                self.assertAlmostEqual(acc, 1.0, delta=5e-3, msg=f"{s.label} cut={cut}")


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

    def test_random_half_rays_every_shape(self):
        # Cut-away: twice the integral over the back half (t <= 0) / front half.
        rng = random.Random(77)
        for s in SHAPES:
            if s.l == 0:
                continue
            z = 2.9
            k = z / s.n
            rms = math.sqrt(orb.mean_r2(s.n, s.l, z))
            for _ in range(3):
                p, d = random_ray(rng, 1.5 * rms)
                tmax = 40.0 / k
                f = lambda t: s.density(z, [p[i] + t * d[i] for i in range(3)])
                back, front = simpson(f, -tmax, 0, 10000), simpson(f, 0, tmax, 10000)
                scale = k * k * max(TABLES[s.row][0])
                for cut, half in ((-1, back), (1, front)):
                    want = 2 * half
                    got = s.column(moments(s.row), z, p, d, cut)
                    self.assertAlmostEqual(got, want, delta=2e-3 * max(want, 1e-2 * scale) + 1e-9,
                                           msg=f"{s.label} cut={cut} p={p} d={d}")
                self.assertGreater(abs(back - front), 0.0)


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


class TestMonteCarlo(unittest.TestCase):
    def test_sampled_histograms_match_tables(self):
        for i, (lab, d) in enumerate([("3s", [0, 0, 1]), ("4px", [0, 0.6, 0.8]),
                                      ("3dxy", [0, 0, 1]), ("3dz2", [0.6, 0, 0.8]),
                                      ("6s", [0.6, 0, 0.8]), ("4fxyz", [0, 0, 1]),
                                      ("4fz3", [0.36, 0.48, 0.8])]):
            err, _ = gt.mc_shape_check(TABLES, BY_LABEL[lab], d, n_samples=20000, seed=100 + i)
            self.assertLess(err, 0.012, lab)

    def test_cut_away_quadrant_histograms(self):
        for i, (lab, d, cut) in enumerate([("2px", [0.6, 0, 0.8], -1), ("3dz2", [0.36, 0.48, 0.8], 1),
                                           ("4fxz2", [0.36, 0.48, 0.8], -1)]):
            err, _ = gt.mc_shape_check(TABLES, BY_LABEL[lab], d, n_samples=20000, seed=200 + i,
                                       odd=ODD, cut=cut)
            self.assertLess(err, 0.012, f"{lab} cut={cut}")


if __name__ == "__main__":
    unittest.main()
