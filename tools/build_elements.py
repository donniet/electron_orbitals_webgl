#!/usr/bin/env python3
"""Build data/elements.json (Z = 1..86) from cached, cited sources (stdlib only).

Sources (in data/sources/):

* pubchem_periodictable.json -- PubChem periodic table
  (https://pubchem.ncbi.nlm.nih.gov/rest/pug/periodictable/JSON):
  name, ground-state electron configuration, Pauling electronegativity,
  first ionization energy (eV), group block, atomic weight.
* clementi_zeff.csv -- effective nuclear charges Z_eff = n * zeta from
  Clementi & Raimondi, J. Chem. Phys. 38, 2686 (1963) and Clementi,
  Raimondi & Reinhardt, J. Chem. Phys. 47, 1300 (1967), as tabulated in the
  `mendeleev` package database (screening constants sigma = Z - Z_eff).
  Regenerate with --from-mendeleev-db PATH/elements.db.
* Elements whose ground-state configuration has a subshell without a
  Clementi value (La, Ce, Gd: 5d) are written with "supported": false and
  a reason; the app greys them out.
* overrides.json -- reviewed values that take precedence over PubChem
  (atomic weights, electronegativities and ionization energies for Z <= 36
  that the app used before this data pipeline; PubChem rounds several
  atomic weights more coarsely, e.g. Li 7.0, Ar 39.9).

Usage:
    python3 tools/build_elements.py            # writes data/elements.json
    python3 tools/build_elements.py --from-mendeleev-db elements.db
"""
import argparse
import csv
import json
import os
import re
import sqlite3
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "data", "sources")
PUBCHEM = os.path.join(SRC, "pubchem_periodictable.json")
ZEFF_CSV = os.path.join(SRC, "clementi_zeff.csv")
OVERRIDES = os.path.join(SRC, "overrides.json")
OUT = os.path.join(ROOT, "data", "elements.json")

MAX_Z = 86
CAPACITY = {"s": 2, "p": 6, "d": 10, "f": 14}
NOBLE = {"He": 2, "Ne": 10, "Ar": 18, "Kr": 36, "Xe": 54, "Rn": 86}
PERIOD_START = [1, 3, 11, 19, 37, 55, 87]

# PubChem GroupBlock -> category names used by the app (and its CSS).
CATEGORY = {
    "Nonmetal": "Reactive nonmetal",
    "Alkali metal": "Alkali metal",
    "Alkaline earth metal": "Alkaline earth metal",
    "Transition metal": "Transition metal",
    "Post-transition metal": "Post-transition metal",
    "Metalloid": "Metalloid",
    "Halogen": "Halogen",
    "Noble gas": "Noble gas",
    "Lanthanide": "Lanthanide",
}


def period_of(z):
    return sum(1 for s in PERIOD_START if z >= s)


def group_of(z):
    """IUPAC group 1-18; None for the f-block (La-Lu, shown in its own row)."""
    p = period_of(z)
    k = z - PERIOD_START[p - 1]
    if p == 1:
        return 1 if z == 1 else 18
    if p in (2, 3):
        return k + 1 if k < 2 else k + 11
    if p in (4, 5):
        return k + 1
    if p == 6:
        if k < 2:
            return k + 1
        if 57 <= z <= 71:
            return None
        return z - 72 + 4
    raise ValueError(z)


def parse_config(text, table):
    """'[Kr]5s2 4d10 5p1' -> {'1s': 2, ..., '5p': 1}; cores expand recursively."""
    out = {}
    m = re.match(r"\s*\[(\w+)\]", text)
    if m:
        out.update(parse_config(table[NOBLE[m.group(1)]]["ElectronConfiguration"], table))
        text = text[m.end():]
    for n, l, c in re.findall(r"(\d)([spdf])(\d+)", text):
        out[n + l] = out.get(n + l, 0) + int(c)
    return out


def load_pubchem():
    t = json.load(open(PUBCHEM))["Table"]
    cols = t["Columns"]["Column"]
    rows = [dict(zip(cols, r["Cell"])) for r in t["Row"]]
    return {int(r["AtomicNumber"]): r for r in rows}


def load_zeff():
    z = {}
    with open(ZEFF_CSV) as fh:
        for row in csv.DictReader(line for line in fh if not line.startswith("#")):
            z.setdefault(int(row["Z"]), {})[row["subshell"]] = float(row["zeff"])
    return z


def num(s):
    return float(s) if s not in ("", None) else None


def build():
    pub = load_pubchem()
    zeff = load_zeff()
    over = json.load(open(OVERRIDES))["elements"]
    out, notes = [], []
    for z in range(1, MAX_Z + 1):
        p = pub[z]
        sub = parse_config(p["ElectronConfiguration"], pub)
        assert sum(sub.values()) == z, (z, sub)
        for k, c in sub.items():
            assert 0 < c <= CAPACITY[k[1]], (z, k, c)
        # The model needs a Clementi Z_eff for every occupied subshell. Where
        # the 1967 tables computed a different configuration (La, Ce, Gd as
        # 4f^n 6s^2, so no 5d value), the element is marked unsupported and
        # shown greyed out rather than filled with an estimate.
        missing = [k for k in sub if k not in zeff.get(z, {})]
        rec = {
            "z": z,
            "symbol": p["Symbol"],
            "name": p["Name"],
            "mass": num(p["AtomicWeight"]),
            "group": group_of(z),
            "period": period_of(z),
            "category": CATEGORY[p["GroupBlock"]],
            "electronegativity": num(p["Electronegativity"]),
            "ionizationEnergy": num(p["IonizationEnergy"]),
            # order subshells by (n, l) for stable output
            "subshells": dict(sorted(sub.items(), key=lambda kv: (kv[0][0], "spdf".index(kv[0][1])))),
            "zeff": {k: round(zeff[z][k], 3)
                     for k in sorted(sub, key=lambda k: (k[0], "spdf".index(k[1])))
                     if k in zeff.get(z, {})},
            "supported": not missing,
        }
        if missing:
            rec["unsupportedReason"] = (
                f"No Clementi Z_eff for {', '.join(missing)}: the 1967 tables computed "
                f"{p['Symbol']} in a different configuration, so this model can't render "
                f"its ground state ({p['ElectronConfiguration']}).")
            notes.append(f"{rec['symbol']:>2} unsupported: no Clementi Z_eff for {', '.join(missing)}")
        o = over.get(str(z))
        if o:
            for key, val in o.items():
                if rec.get(key) != val:
                    notes.append(f"{rec['symbol']:>2} {key}: PubChem {rec.get(key)} -> override {val}")
                rec[key] = val
        out.append(rec)
    return out, notes


def write(elements):
    doc = {
        "generator": "tools/build_elements.py",
        "sources": {
            "pubchem": "PubChem periodic table, https://pubchem.ncbi.nlm.nih.gov/rest/pug/periodictable/JSON",
            "zeff": "Clementi & Raimondi, J. Chem. Phys. 38, 2686 (1963); Clementi, Raimondi & Reinhardt, "
                    "J. Chem. Phys. 47, 1300 (1967); via the mendeleev package database",
            "overrides": "data/sources/overrides.json (reviewed values for Z <= 36)",
        },
        "elements": elements,
    }
    with open(OUT, "w") as fh:
        json.dump(doc, fh, indent=1, ensure_ascii=False)
        fh.write("\n")


def extract_mendeleev(db_path):
    """Write clementi_zeff.csv from the mendeleev elements.db (read-only)."""
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    rows = con.execute(
        "select atomic_number, n, s, screening from screeningconstants "
        "where atomic_number <= ? order by atomic_number, n, s", (MAX_Z,)).fetchall()
    with open(ZEFF_CSV, "w", newline="") as fh:
        fh.write("# Z_eff = Z - sigma, sigma from the mendeleev package (screeningconstants table),\n")
        fh.write("# which tabulates Clementi & Raimondi (1963) and Clementi, Raimondi & Reinhardt (1967).\n")
        w = csv.writer(fh)
        w.writerow(["Z", "subshell", "zeff"])
        for z, n, s, sigma in rows:
            w.writerow([z, f"{n}{s}", f"{z - sigma:.4f}"])
    print(f"Wrote {ZEFF_CSV} ({len(rows)} rows)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--from-mendeleev-db", metavar="PATH",
                    help="regenerate data/sources/clementi_zeff.csv from mendeleev's elements.db")
    args = ap.parse_args(argv)
    if args.from_mendeleev_db:
        extract_mendeleev(args.from_mendeleev_db)
    elements, notes = build()
    write(elements)
    print(f"Wrote {OUT} ({len(elements)} elements)")
    if notes:
        print("Notes (overrides: PubChem value -> reviewed value):")
        for n in notes:
            print("  " + n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
