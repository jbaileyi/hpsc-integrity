"""Re-score arms after subtracting each arm's leave-one-out cross-line median (panel-of-normals style).
    python scripts/cohort_baseline.py"""
import numpy as np, pandas as pd
M = ["CSES9", "CSES7", "CSES22", "CSES21", "HUES14"]
SPECIAL = {"CSES22": ["12p","12q"], "CSES21": ["21q"], "CSES7": ["17p","17q"], "CSES9": ["20p","20q"]}
REFUSE = {("CSES9","9q"), ("HUES14","9q"), ("CSES22","16p")}
T = pd.DataFrame({l: pd.read_csv(f"results/allelic/{l}_arms.csv", index_col=0).query("n_snps >= 50")["log2_vs_baseline"] for l in M})
T = T.drop(index=["Xp", "Xq"], errors="ignore")
for l, a in REFUSE: T.loc[a, l] = np.nan
S = T.copy(); P = {}
for l in M:
    others = T.drop(columns=l)
    # leave-one-out panel median, excluding the other lines' own special arms
    for o in others.columns:
        for a in SPECIAL.get(o, []):
            if a in others.index: others.loc[a, o] = np.nan
    S[l] = T[l] - others.median(axis=1)
for l in M:
    null = S[l].drop(index=[a for a in SPECIAL.get(l, []) if a in S.index]).dropna()
    sd = 1.4826 * np.median(np.abs(null - null.median()))
    raw = T[l].drop(index=[a for a in SPECIAL.get(l, []) if a in T.index]).dropna()
    raw_sd = 1.4826 * np.median(np.abs(raw - raw.median()))
    print(f"{l}: robust SD raw {raw_sd:.3f} -> corrected {sd:.3f}")
    for a in SPECIAL.get(l, []):
        if a in S.index and pd.notna(S.loc[a, l]):
            print(f"   {a}: raw {T.loc[a, l]:+.3f} -> corrected {S.loc[a, l]:+.3f} (z {S.loc[a, l] / sd:+.1f})")
S.round(3).to_csv("results/allelic/cohort_corrected.csv")
