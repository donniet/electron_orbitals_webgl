"""Validate data/elements.json (built by tools/build_elements.py).

Run: python3 -m unittest discover tests
"""
import json
import os
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOC = json.load(open(os.path.join(ROOT, "data", "elements.json")))
ELEMENTS = DOC["elements"]
CAP = {"s": 2, "p": 6, "d": 10, "f": 14}


class TestElements(unittest.TestCase):
    def test_sequential_and_complete(self):
        self.assertEqual([e["z"] for e in ELEMENTS], list(range(1, 87)))

    def test_occupancy(self):
        for e in ELEMENTS:
            self.assertEqual(sum(e["subshells"].values()), e["z"], e["symbol"])
            for k, c in e["subshells"].items():
                self.assertTrue(0 < c <= CAP[k[1]], (e["symbol"], k, c))

    def test_zeff_for_supported(self):
        for e in ELEMENTS:
            missing = [k for k in e["subshells"] if k not in e["zeff"]]
            self.assertEqual(e["supported"], not missing, e["symbol"])
            for v in e["zeff"].values():
                self.assertTrue(0 < v <= e["z"], (e["symbol"], v))

    def test_only_la_ce_gd_unsupported(self):
        bad = sorted(e["symbol"] for e in ELEMENTS if not e["supported"])
        self.assertEqual(bad, ["Ce", "Gd", "La"])
        for e in ELEMENTS:
            if not e["supported"]:
                self.assertIn("5d", e["unsupportedReason"])

    def test_groups_and_periods(self):
        by = {e["symbol"]: e for e in ELEMENTS}
        self.assertEqual((by["H"]["group"], by["He"]["group"], by["B"]["group"]), (1, 18, 13))
        self.assertEqual((by["Y"]["group"], by["Pd"]["group"], by["Xe"]["group"]), (3, 10, 18))
        self.assertEqual((by["Cs"]["group"], by["Hf"]["group"], by["Rn"]["group"]), (1, 4, 18))
        for e in ELEMENTS:
            if 57 <= e["z"] <= 71:
                self.assertIsNone(e["group"], e["symbol"])
        self.assertEqual([by[s]["period"] for s in ("H", "Li", "Na", "K", "Rb", "Cs")], [1, 2, 3, 4, 5, 6])

    def test_known_configurations(self):
        by = {e["symbol"]: e["subshells"] for e in ELEMENTS}
        self.assertEqual((by["Cr"]["3d"], by["Cr"]["4s"]), (5, 1))
        self.assertEqual((by["Cu"]["3d"], by["Cu"]["4s"]), (10, 1))
        self.assertNotIn("5s", by["Pd"])
        self.assertEqual((by["Pt"]["5d"], by["Pt"]["6s"]), (9, 1))
        self.assertEqual((by["Au"]["5d"], by["Au"]["6s"]), (10, 1))
        self.assertEqual((by["Gd"]["4f"], by["Gd"]["5d"]), (7, 1))

    def test_rahm_atomic_radius(self):
        by = {e["symbol"]: e for e in ELEMENTS}
        for e in ELEMENTS:
            self.assertTrue(100 <= e["atomicRadius"] <= 400, (e["symbol"], e["atomicRadius"]))
        self.assertEqual((by["H"]["atomicRadius"], by["He"]["atomicRadius"], by["Cs"]["atomicRadius"]),
                         (154, 134, 249))
        self.assertIn("Rahm", DOC["sources"]["atomicRadius"])

    def test_reviewed_values_kept_for_first_36(self):
        by = {e["symbol"]: e for e in ELEMENTS}
        self.assertEqual(by["Li"]["mass"], 6.94)
        self.assertEqual(by["Ar"]["mass"], 39.948)
        self.assertEqual(by["As"]["ionizationEnergy"], 9.789)


if __name__ == "__main__":
    unittest.main()
