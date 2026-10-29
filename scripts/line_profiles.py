"""Which arms deviate in each line, and do they persist across reference panels?
Run from the repository root: python scripts/line_profiles.py"""
from pathlib import Path

import pandas as pd

from hpsc_integrity.annotation import load_reference
from hpsc_integrity.pipeline import score_arms

pd.set_option("display.width", 220)
ROOT = Path(__file__).resolve().parents[1]
counts = pd.read_csv(ROOT / "data/processed/gene_counts_all13.csv", index_col=0)
genes, arms = load_reference()
HESC = ["CSES9", "CSES7", "CSES22", "CSES21", "HUES14"]
HISEQ = [c for c in counts.columns if c not in HESC]
GAINED = {"CSES22": "12", "CSES21": "21"}


def arm_values(test, comps):
    t, _, _ = score_arms(counts[[test] + comps], genes, arms)
    t = t[(t["n_genes"] >= 50) & ~t.index.str.startswith(("X", "Y"))]
    t = t[t.index.str[:-1] != GAINED.get(test)]          # drop the arms of the known trisomy
    return t[test]


rows = []
for test in counts.columns:
    v = arm_values(test, [c for c in counts.columns if c != test])
    rows.append(dict(line=test, group="hESC" if test in HESC else "HiSeq", arms=len(v),
                     n_ge_0_12=int((v.abs() >= 0.12).sum()), n_ge_0_20=int((v.abs() >= 0.20).sum()),
                     n_ge_0_30=int((v.abs() >= 0.30).sum()), max_abs=round(v.abs().max(), 2),
                     robust_sd=round((v - v.median()).abs().median() * 1.4826, 3)))
print("Each line against all 12 others (arms of the known trisomy excluded):")
print(pd.DataFrame(rows).to_string(index=False))

print("\nhESC lines: arms with |log2| >= 0.15 under panel A (4 hESC) and panel D (8 HiSeq)  [A value, D value]")
key = lambda a: (int(a[:-1]), a[-1])
for test in HESC:
    a = arm_values(test, [c for c in HESC if c != test])
    d = arm_values(test, HISEQ)
    flagged = sorted(set(a[a.abs() >= 0.15].index) | set(d[d.abs() >= 0.15].index), key=key)
    print(f"{test}:", {k: (round(a.get(k, float("nan")), 2), round(d.get(k, float("nan")), 2)) for k in flagged})
