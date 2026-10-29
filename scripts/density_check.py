"""Single-sample vs leave-one-out arm-specific het density, from results/allelic/calls_<line>.csv. python scripts/density_check.py"""
import pandas as pd
M = ["CSES9", "CSES7", "CSES22", "CSES21", "HUES14"]
C = {l: pd.read_csv(f"results/allelic/calls_{l}.csv", index_col=0) for l in M}
R = pd.DataFrame({l: C[l]["rel_het_density"] for l in M})
for l in M:
    loo = R[l] / R.drop(columns=l).median(axis=1).where(lambda s: s > 0.3)
    t = pd.DataFrame({"single": R[l], "loo": loo, "n_snps": C[l].n_snps, "call": C[l].call})
    t = t[(t.n_snps >= 20) & ((t.single < 0.6) | (t.loo < 0.6))]
    print(f"\n{l}"); print(t.round(2).to_string())
print("\narms with single-sample density < 0.7 in >= 3 lines (systematic):")
print(R[(R < 0.7).sum(axis=1) >= 3].round(2).to_string())
