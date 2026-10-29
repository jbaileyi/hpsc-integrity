"""Robust per-arm summary: MAD-based noise, min-SNP filter, z-scores, mosaic detection limit.
    python scripts/allelic_summary_v2.py [--min-snps 50] [--z 3]"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--min-snps", type=int, default=50)
ap.add_argument("--z", type=float, default=3.0)
a = ap.parse_args()

ROOT = Path(__file__).resolve().parents[1]
LINES = ["CSES9", "CSES7", "CSES22", "CSES21", "HUES14"]
KNOWN = {"CSES22": ["12p", "12q"], "CSES21": ["21q"]}
WATCH = {"CSES7": ["17p", "17q"], "CSES9": ["20p", "20q"]}

def robust_sd(x):
    x = np.asarray(x); return 1.4826 * np.median(np.abs(x - np.median(x)))

rows = []
for line in LINES:
    t = pd.read_csv(ROOT / f"results/allelic/{line}_arms.csv", index_col=0)
    t = t[t["log2_vs_baseline"].notna() & (t["n_snps"] >= a.min_snps)]
    special = KNOWN.get(line, []) + WATCH.get(line, [])
    null = t.drop(index=[x for x in special if x in t.index])
    sd = robust_sd(null["log2_vs_baseline"])
    # drop arms >4 robust SDs from the null median, recompute (iterative, 3 rounds)
    for _ in range(3):
        med = null["log2_vs_baseline"].median()
        keep = null[(null["log2_vs_baseline"] - med).abs() <= 4 * sd]
        sd = robust_sd(keep["log2_vs_baseline"])
    # mosaic fraction f of cells with +1 copy: log2(1 + f/2) >= z*sd
    f_min = min(1.0, 2 * (2 ** (a.z * sd) - 1))
    outl = null[(null["log2_vs_baseline"] - med).abs() > 4 * sd]
    print(f"\n{line}: robust SD {sd:.3f} | min detectable gain fraction at z={a.z}: {f_min:.0%} | null outlier arms: "
          f"{ {k: round(v,3) for k,v in outl['log2_vs_baseline'].items()} }")
    for arm in special:
        if arm in t.index:
            v = t.loc[arm, "log2_vs_baseline"]
            print(f"  {arm}: log2 {v:+.3f}  z={(v-med)/sd:+.1f}  n_snps={int(t.loc[arm,'n_snps'])}")
            rows.append(dict(line=line, arm=arm, log2=v, z=(v-med)/sd, n_snps=t.loc[arm,'n_snps'], min_detectable_fraction=f_min))
pd.DataFrame(rows).to_csv(ROOT / "results/allelic/summary_robust.csv", index=False)
