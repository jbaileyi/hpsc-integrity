"""Per-arm allelic imbalance with locus-level bootstrap CIs and QC flags.
    python scripts/allelic_bootstrap.py data/bam/SRR3501465.snps.tsv CSES9 [--cluster-bp 100000] [--boot 2000]
SNPs within --cluster-bp of each other are collapsed to one locus (they share a gene, so they are not independent);
loci are then resampled within each arm. Writes results/allelic/<line>_boot.csv"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("snps"); ap.add_argument("line")
ap.add_argument("--cluster-bp", type=int, default=100_000)
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--min-snps", type=int, default=50)
a = ap.parse_args()
ROOT = Path(__file__).resolve().parents[1]
# approximate hg38 p/q boundaries; the CHECK line below tells you whether they match allelic_profile.py
CEN = {1:123_400_000, 2:93_900_000, 3:90_900_000, 4:50_000_000, 5:48_800_000, 6:59_800_000, 7:60_100_000, 8:45_200_000,
       9:43_300_000, 10:39_800_000, 11:53_400_000, 12:35_500_000, 13:17_700_000, 14:17_200_000, 15:19_000_000,
       16:36_800_000, 17:25_100_000, 18:18_500_000, 19:26_500_000, 20:28_100_000, 21:12_500_000, 22:15_000_000}
ACRO = {13, 14, 15, 21, 22}
SHOW = {"CSES22": ["12p","12q"], "CSES21": ["21q"], "CSES7": ["17p","17q"], "CSES9": ["20p","20q","9q"], "HUES14": ["9q"]}

d = pd.read_csv(a.snps, sep="\t", header=None, names=["chrom","pos","ref","alt","c"], dtype={"chrom": str})
d = d[d.chrom.isin([str(i) for i in range(1, 23)])].copy()
d["chrom"] = d.chrom.astype(int)
ab = d.c.str.split(",", expand=True).astype(int)
d["r"], d["a"] = ab[0].values, ab[1].values
dep = d.r + d.a
d = d[(dep >= 20) & (np.minimum(d.r, d.a) / dep >= 0.2)].copy()
d["ratio"] = np.maximum(d.r, d.a) / np.minimum(d.r, d.a)
d["refdom"] = (d.r > d.a).astype(float)
d["arm"] = d.chrom.astype(str) + np.where(d.pos < d.chrom.map(CEN), "p", "q")
d = d[~(d.arm.str.endswith("p") & d.chrom.isin(ACRO))]

n_snp = d.groupby("arm").size(); ok = n_snp[n_snp >= a.min_snps].index
snp_med = d.groupby("arm").ratio.median()
old = pd.read_csv(ROOT / f"results/allelic/{a.line}_arms.csv", index_col=0)["log2_vs_baseline"].dropna()
diff = (np.log2(snp_med[ok] / snp_med[ok].median()) - old).abs().dropna()
print(f"CHECK vs allelic_profile.py: {len(diff)} arms compared, max |diff| {diff.max():.3f} (should be ~0.000)")

d = d.sort_values(["chrom", "pos"])
d["locus"] = ((d.chrom.diff().fillna(1) != 0) | (d.pos.diff().fillna(0) > a.cluster_bp)).cumsum()
L = d.groupby(["arm", "locus"]).agg(ratio=("ratio", "median"), refdom=("refdom", "mean")).reset_index()
arms = [x for x in ok if (L.arm == x).sum() >= 3]
G = {x: L.loc[L.arm == x, "ratio"].to_numpy() for x in arms}
R = pd.Series({x: L.loc[L.arm == x, "refdom"].mean() for x in arms})

rng = np.random.default_rng(1)
pt = np.array([np.median(G[x]) for x in arms])
M = np.column_stack([np.median(G[x][rng.integers(0, len(G[x]), (a.boot, len(G[x])))], axis=1) for x in arms])
B = np.log2(M / np.median(M, axis=1, keepdims=True))
lo, hi = np.percentile(B, [2.5, 97.5], axis=0)

p0 = R.mean()
res = pd.DataFrame({"n_snps": n_snp[arms].values, "n_loci": [len(G[x]) for x in arms],
                    "log2": np.log2(pt / np.median(pt)), "ci_lo": lo, "ci_hi": hi}, index=arms)
res["refdom"] = R
res["ref_z"] = (res.refdom - p0) / np.sqrt(p0 * (1 - p0) / res.n_loci)
# crude model ratio ~ 1 + f/2: largest / smallest +1-copy clone fraction compatible with the CI
res["gain_frac_hi95"] = np.clip(2 * (2 ** res.ci_hi - 1), 0, 1)
res["gain_frac_lo95"] = np.clip(2 * (2 ** res.ci_lo - 1), 0, 1)
res["flags"] = [",".join(f for f, bad in [("FEW_LOCI", r.n_loci < 30), ("WIDE_CI", r.ci_hi - r.ci_lo > 0.30),
                                          ("REF_BIAS", abs(r.ref_z) > 3)] if bad) for _, r in res.iterrows()]
res.round(3).to_csv(ROOT / f"results/allelic/{a.line}_boot.csv")

print(f"{a.line}: {len(d)} SNPs -> {len(L)} loci (cluster {a.cluster_bp // 1000} kb); ref-dominant loci {p0:.2f}")
print(res.loc[[x for x in SHOW.get(a.line, []) if x in res.index]].round(3).to_string())
dev = res[((res.ci_lo > 0.10) | (res.ci_hi < -0.10)) & ~res.index.isin(SHOW.get(a.line, []))]
print("other arms whose CI lies fully beyond +/-0.10:", dev[["log2", "ci_lo", "ci_hi", "flags"]].round(2).to_dict("index"))
