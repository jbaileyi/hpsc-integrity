"""Het-site depletion. (1) per-arm share of kept SNPs vs the other lines; (2) chr9 het/total sites per 5 Mb.
    python scripts/diag_het.py"""
import pandas as pd
M = {"CSES9":"SRR3501465","CSES7":"SRR3501466","CSES22":"SRR3501467","CSES21":"SRR3501468","HUES14":"SRR3501469"}
A = pd.DataFrame({l: pd.read_csv(f"results/allelic/{l}_arms.csv", index_col=0)["n_snps"] for l in M}).drop(index=["Xp","Xq"], errors="ignore")
share = A / A.sum()
print("arms with <40% of the typical (median of other lines) SNP share:")
for l in M:
    rel = share[l] / share.drop(columns=l).median(axis=1)
    print(f"  {l}: {rel[rel < 0.4].round(2).to_dict()}")
W, out = 5_000_000, {}
for l, s in M.items():
    d = pd.read_csv(f"data/bam/{s}.snps.tsv", sep="\t", header=None, names=["chrom","pos","ref","alt","c"], dtype={"chrom":str})
    d = d[d.chrom == "9"]
    ab = d.c.str.split(",", expand=True).astype(int)
    dep = ab.sum(axis=1); keep = dep >= 10
    d = d[keep].assign(het=(ab.min(axis=1)[keep] / dep[keep]) >= 0.2, win=lambda x: (x.pos // W) * 5)
    out[l] = d.groupby("win").het.agg(lambda h: f"{int(h.sum())}/{len(h)}")
print("\nchr9 het/total variant sites (depth>=10) per 5 Mb window (start Mb):")
print(pd.DataFrame(out).fillna("-").to_string())
