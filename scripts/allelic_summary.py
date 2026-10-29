"""Compare per-arm allelic results across lines. Run after allelic_profile.py has been run for each line:
    python scripts/allelic_summary.py
Each line's noise (SD) is taken from its own arms, leaving out the arms of its known trisomy and the arms
under test, so a real shift cannot inflate its own noise estimate."""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LINES = ["CSES9", "CSES7", "CSES22", "CSES21", "HUES14"]
KNOWN = {"CSES22": ["12p", "12q"], "CSES21": ["21q"]}         # G-banding trisomies (positive controls)
WATCH = {"CSES7": ["17p", "17q"], "CSES9": ["20p", "20q"]}     # arms with unexplained expression shifts

for line in LINES:
    f = ROOT / f"results/allelic/{line}_arms.csv"
    if not f.exists():
        print(f"\n{line}: no results yet")
        continue
    t = pd.read_csv(f, index_col=0)
    t = t[t["log2_vs_baseline"].notna()]
    special = KNOWN.get(line, []) + WATCH.get(line, [])
    null = t.drop(index=[a for a in special if a in t.index])
    sd = null["log2_vs_baseline"].std()
    print(f"\n{line}: {len(t)} arms assessed, null SD {sd:.3f}, largest |log2| in null {null['log2_vs_baseline'].abs().max():.3f}")
    for arm in special:
        if arm in t.index:
            r = t.loc[arm]
            print(f"  {arm}: log2 {r['log2_vs_baseline']:+.3f} ({r['log2_vs_baseline'] / sd:+.1f} null SDs), {int(r['n_snps'])} SNPs")
        else:
            print(f"  {arm}: not assessed (too few SNPs)")
    flagged = null[null["log2_vs_baseline"].abs() >= 0.10]
    if len(flagged):
        print("  other arms with |log2| >= 0.10:", {a: round(v, 3) for a, v in flagged["log2_vs_baseline"].items()})
