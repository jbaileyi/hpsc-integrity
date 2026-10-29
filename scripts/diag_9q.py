"""Windowed allelic ratio on chr9 per line + shared top SNPs. python scripts/diag_9q.py"""
import pandas as pd, numpy as np
M = {"CSES9":"SRR3501465","CSES7":"SRR3501466","CSES22":"SRR3501467","CSES21":"SRR3501468","HUES14":"SRR3501469"}
W = 5_000_000
res, tops = {}, {}
for line, srr in M.items():
    d = pd.read_csv(f"data/bam/{srr}.snps.tsv", sep="\t", header=None, names=["chrom","pos","ref","alt","c"], dtype={"chrom":str})
    d = d[d.chrom == "9"].copy()
    ab = d.c.str.split(",", expand=True).astype(int)
    d["depth"] = ab.sum(axis=1); d["maf"] = ab.min(axis=1) / d.depth
    d = d[(d.depth >= 20) & (d.maf >= 0.2)]
    d["ratio"] = ab.max(axis=1) / ab.min(axis=1).replace(0, np.nan)
    d["win"] = (d.pos // W) * W // 1_000_000
    res[line] = d.groupby("win").ratio.agg(["median", "size"]).apply(lambda r: f"{r['median']:.2f} (n={int(r['size'])})", axis=1)
    tops[line] = d.sort_values("ratio", ascending=False).head(15)
print(pd.DataFrame(res).fillna("-").to_string())
for l in ["CSES9", "HUES14"]:
    print(f"\n{l} top-ratio SNPs on chr9"); print(tops[l][["pos","ref","alt","c","ratio"]].to_string(index=False))
print("\nshared top positions:", sorted(set(tops["CSES9"].pos) & set(tops["HUES14"].pos)))
