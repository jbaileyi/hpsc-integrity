"""Numeric reading of eSNP-Karyotyping MajorMinorCalc output (rule pre-registered in docs/design.md, amendment 3).
python scripts/bench_esnp_numeric.py  -> results/benchmark/esnp_numeric.csv"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from hpsc_integrity.allelic import assign_snp_arms
from hpsc_integrity.annotation import load_reference

W = int(sys.argv[1]) if len(sys.argv) > 1 else 151
FACTOR, MIN_FRAC, MIN_SNPS = 1.25, 0.5, 30
_, arms = load_reference()
rows = []
for f in sorted(Path("results/benchmark/esnp_mm").glob("S*.csv")):
    d = pd.read_csv(f)
    d = d[d.chr <= 22].sort_values(["chr", "position"]).reset_index(drop=True)
    d["rm"] = d["MajorMinor"].rolling(W, center=True, min_periods=W).median()
    base = d["rm"].median()
    d["arm"] = assign_snp_arms(pd.DataFrame({"chrom": d.chr.astype(int).astype(str), "pos": d.position}), arms)
    d = d.dropna(subset=["arm"])
    for a, g in d.groupby("arm"):
        r = g["rm"].dropna()
        frac = float((r > FACTOR * base).mean()) if len(r) else np.nan
        rows.append(dict(id=f.stem, arm=a, n_snps=len(g), median_rm=round(float(r.median()), 3) if len(r) else np.nan,
                         frac_above=round(frac, 2), baseline=round(float(base), 3),
                         raised=bool(len(g) >= MIN_SNPS and frac >= MIN_FRAC)))
R = pd.DataFrame(rows)
R.to_csv("results/benchmark/esnp_numeric.csv" if W == 151 else f"results/benchmark/esnp_numeric_w{W}.csv", index=False)
for i, g in R.groupby("id"):
    print(i, "baseline", g.baseline.iloc[0], "| raised arms:", list(g.loc[g.raised, "arm"]) or "none")
