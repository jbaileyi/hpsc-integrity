"""Per-arm het density relative to the cohort; flag arms the tool should refuse to call.
    python scripts/het_density.py"""
import pandas as pd
M = ["CSES9","CSES7","CSES22","CSES21","HUES14"]
A = pd.DataFrame({l: pd.read_csv(f"results/allelic/{l}_arms.csv", index_col=0)["n_snps"] for l in M}).drop(index=["Xp","Xq"], errors="ignore")
share = A / A.sum()
rows = []
for l in M:
    rel = share[l] / share.drop(columns=l).median(axis=1)
    rows += [dict(line=l, arm=a, n_snps=int(A.loc[a, l]), rel_het_density=round(v, 2), verdict="REFUSE" if v < 0.4 else "ok")
             for a, v in rel.items() if v < 0.6 and A.loc[a, l] >= 20]
df = pd.DataFrame(rows).sort_values("rel_het_density"); print(df.to_string(index=False))
df.to_csv("results/allelic/het_density_flags.csv", index=False)
