"""Tests for the generated orbital column-density tables (stdlib unittest).

These check the *checked-in* js/orbital_tables.js, parsed directly, so a
stale or corrupted table file is caught.

Run: python3 -m unittest discover tests
"""
import math
import os
import random
import re
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import gen_tables as gt  # noqa: E402

TABLE_JS = os.path.join(ROOT, "js", "orbital_tables.js")


def load_js_tables():
    src = open(TABLE_JS).read()
    n = int(re.search(r"TABLE_N\s*=\s*(\d+)", src).group(1))
    b_max = float(re.search(r"TABLE_B_MAX\s*=\s*([0-9.eE+-]+)", src).group(1))
    body = re.search(r"Float32Array\(\[(.*?)\]\)", src, re.S).group(1)
    vals = [float(v) for v in body.replace("\n", " ").split(",") if v.strip()]
    assert len(vals) == 4 * n, (len(vals), n)
    rows = [tuple(vals[4 * i:4 * i + 4]) for i in range(n)]
    return rows, b_max


ROWS, B_MAX = load_js_tables()


def lookup(ch, b):
    return gt.lookup(ROWS, ch, b, B_MAX)


def brute_line_integral(rho, p, d, t_max=60.0, steps=24000):
    """Simpson integral of rho(p + t d) for t in [-t_max, t_max]."""
    h = 2 * t_max / steps
    acc = 0.0
    for i in range(steps + 1):
        t = -t_max + i * h
        w = 1 if i in (0, steps) else (4 if i & 1 else 2)
        acc += w * rho(p[0] + t * d[0], p[1] + t * d[1], p[2] + t * d[2])
    return acc * h / 3.0


def random_ray(rng, b_lo=0.0, b_hi=8.0):
    """Unit direction d and closest point p (p . d = 0) with |p| in range."""
    while True:
        d = [rng.gauss(0, 1) for _ in range(3)]
        n = math.sqrt(sum(c * c for c in d))
        if n > 1e-6:
            d = [c / n for c in d]
            break
    while True:
        q = [rng.gauss(0, 1) for _ in range(3)]
        dot = sum(a * b for a, b in zip(q, d))
        q = [a - dot * b for a, b in zip(q, d)]
        n = math.sqrt(sum(c * c for c in q))
        if n > 1e-6:
            break
    b = rng.uniform(b_lo, b_hi)
    p = [c / n * b for c in q]
    return p, d, b


def col_from_tables(orbital, p, d, z=1.0):
    """Column density from tables for an orbital with effective charge z."""
    ps = [c * z for c in p]
    bs = math.sqrt(sum(c * c for c in ps))
    if orbital == "1s":
        return z * z * lookup(2, bs)
    if orbital == "2s":
        return z * z * lookup(3, bs)
    axis = {"2px": 0, "2py": 1, "2pz": 2}[orbital]
    return z * z * (ps[axis] ** 2 * lookup(0, bs) + d[axis] ** 2 * lookup(1, bs))


def scaled_rho(orbital, z):
    base = {"1s": gt.rho_1s, "2s": gt.rho_2s, "2px": gt.rho_2p(0),
            "2py": gt.rho_2p(1), "2pz": gt.rho_2p(2)}[orbital]
    return lambda x, y, w: z ** 3 * base(z * x, z * y, z * w)


class TestAnalyticValues(unittest.TestCase):
    def test_b_zero(self):
        a, b, s1, s2 = ROWS[0]
        self.assertAlmostEqual(a, 1 / (16 * math.pi), places=6)
        self.assertAlmostEqual(b, 1 / (8 * math.pi), places=6)
        self.assertAlmostEqual(s1, 1 / math.pi, places=6)
        self.assertAlmostEqual(s2, 4 / (32 * math.pi), places=6)

    def test_tail_is_negligible(self):
        for v in ROWS[-1]:
            self.assertLess(abs(v), 1e-9)

    def test_non_negative(self):
        for row in ROWS:
            for v in row:
                self.assertGreaterEqual(v, 0.0)


class TestNormalisation(unittest.TestCase):
    """Integrated over the image plane each orbital holds one electron."""

    def integrate(self, f):
        n = len(ROWS)
        h = B_MAX / (n - 1)
        acc = 0.0
        for i, row in enumerate(ROWS):
            w = 0.5 if i in (0, n - 1) else 1.0
            acc += w * f(i * h, row)
        return acc * h

    def test_s_orbitals(self):
        self.assertAlmostEqual(self.integrate(lambda b, r: 2 * math.pi * b * r[2]), 1, places=3)
        self.assertAlmostEqual(self.integrate(lambda b, r: 2 * math.pi * b * r[3]), 1, places=3)

    def test_p_along_axis(self):
        # ray parallel to the lobe: col = B(b)
        self.assertAlmostEqual(self.integrate(lambda b, r: 2 * math.pi * b * r[1]), 1, places=3)

    def test_p_across_axis(self):
        # ray perpendicular to the lobe: col = p_z^2 A(b), <p_z^2> = b^2 / 2
        self.assertAlmostEqual(self.integrate(lambda b, r: math.pi * b ** 3 * r[0]), 1, places=3)


class TestBruteForceRays(unittest.TestCase):
    """Table formula == direct numerical line integral of |psi|^2."""

    def test_random_rays(self):
        rng = random.Random(42)
        for orbital in ("1s", "2s", "2px", "2py", "2pz"):
            for z in (1.0, 2.7, 5.758):
                for _ in range(3):
                    p, d, _b = random_ray(rng, 0.0, 6.0 / z)
                    want = brute_line_integral(scaled_rho(orbital, z), p, d,
                                               t_max=60.0 / z)
                    got = col_from_tables(orbital, p, d, z)
                    tol = 2e-3 * max(want, 1e-3 * z * z) + 1e-7
                    self.assertAlmostEqual(
                        got, want, delta=tol,
                        msg=f"{orbital} z={z} p={p} d={d}")


class TestSymmetry(unittest.TestCase):
    def test_full_p_shell_is_isotropic(self):
        """px+py+pz summed depends only on b (Unsold's theorem)."""
        rng = random.Random(7)
        for b in (0.0, 0.5, 1.5, 3.0, 6.0):
            vals = []
            for _ in range(5):
                p, d, _ = random_ray(rng, b, b)
                vals.append(sum(col_from_tables(o, p, d) for o in ("2px", "2py", "2pz")))
            spread = max(vals) - min(vals)
            self.assertLess(spread, 1e-6 + 1e-5 * max(vals), msg=f"b={b} {vals}")
            # and it equals b^2 A + B
            self.assertAlmostEqual(vals[0], b * b * lookup(0, b) + lookup(1, b), places=6)


class TestMonteCarlo(unittest.TestCase):
    def test_histograms_match_tables(self):
        worst = gt.monte_carlo_verify(ROWS, n_samples=40000, seed=3, verbose=False)
        self.assertLess(worst, 0.01)


if __name__ == "__main__":
    unittest.main()
