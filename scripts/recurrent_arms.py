"""Arms whose allelic log2 has the same sign across most lines (a sign of systematic artifact, not biology).
    python scripts/recurrent_arms.py"""
import numpy as np, pandas as pd
M = ["CSES9", "CSES7", "CSES22", "CSES21", "HUES14"]
EXCL = {("CSES22","12p"), ("CSES22","12q"), ("CSES21","21q"), ("CSES7","17p"), ("CSES7","17q"),
        ("CSES9","20p"), ("CSES9","20q"), ("CSES9","9q"), ("HUES14","9q")}
T = pd.DataFrame({l: pd.read_csv(f"results/allelic/{l}_arms.csv", index_col=0).query("n_snps >= 50")["log2_vs_baseline"] for l in M})
T = T.drop(index=["Xp", "Xq"], errors="ignore")
for l, a in EXCL:
    if a in T.index: T.loc[a, l] = np.nan
med = T.median(axis=1); n = T.notna().sum(axis=1)
same = T.apply(lambda c: np.sign(c) == np.sign(med)).sum(axis=1)
out = T.assign(median=med, n_lines=n, same_sign=same)
print(out[(med.abs() >= 0.05) & (n >= 4)].sort_values("median").round(3).to_string())
out.round(3).to_csv("results/allelic/recurrence.csv")
