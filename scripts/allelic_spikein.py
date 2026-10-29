"""Allelic spike-in: mimic a gain in fraction f of cells on one arm by moving reads between alleles at existing het SNPs
(depth unchanged, no reads invented; which allele sits on the gained haplotype is random per SNP), then run call_arms.
    python scripts/allelic_spikein.py [reps] [n_boot]     (try: 3 200 first)
Writes results/simulation/allelic_spikein_{runs,summary}.csv"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from hpsc_integrity.allelic import assign_snp_arms, filter_het_snps, read_allele_counts
from hpsc_integrity.allelic_calls import GAIN, call_arms
from hpsc_integrity.annotation import load_reference

reps = int(sys.argv[1]) if len(sys.argv) > 1 else 8
n_boot = int(sys.argv[2]) if len(sys.argv) > 2 else 400
LINES = {"CSES9": "SRR3501465", "HUES14": "SRR3501469", "CSES7": "SRR3501466"}
ARMS = ["1q", "12q", "20p"]
FRACS = [0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0]
genes, arms = load_reference()
rows = []
for li, (line, srr) in enumerate(LINES.items()):
    base = filter_het_snps(read_allele_counts(f"data/bam/{srr}.snps.tsv"))[["chrom", "pos", "ref_n", "alt_n"]]
    arm_of = assign_snp_arms(base, arms)
    for ai, arm in enumerate(ARMS):
        idx = np.flatnonzero((arm_of == arm).to_numpy())
        for fi, f in enumerate(FRACS):
            for rep in range(reps):
                rng = np.random.default_rng([rep, fi, ai, li])
                d = base.copy()
                ref, alt = d["ref_n"].to_numpy()[idx], d["alt_n"].to_numpy()[idx]
                alt_gained = rng.random(len(idx)) < 0.5
                moved = rng.binomial(np.where(alt_gained, ref, alt), f / (2 + f))   # f=1 -> allele fraction 2/3
                d.iloc[idx, d.columns.get_loc("ref_n")] = ref + np.where(alt_gained, -moved, moved)
                d.iloc[idx, d.columns.get_loc("alt_n")] = alt + np.where(alt_gained, moved, -moved)
                r = call_arms(filter_het_snps(d), genes, arms, n_boot=n_boot, seed=rep)
                rows.append(dict(line=line, arm=arm, f=f, rep=rep, call=r.at[arm, "call"], log2=r.at[arm, "log2"],
                                 false_gain_other=int((r.drop(index=arm)["call"] == GAIN).sum())))
    print(line, "done", flush=True)
runs = pd.DataFrame(rows)
summ = runs.groupby(["line", "arm", "f"]).agg(
    gain_rate=("call", lambda c: (c == GAIN).mean()),
    refused_or_few=("call", lambda c: c.isin(["REFUSED_HET_DEPLETED", "TOO_FEW_SNPS"]).mean()),
    mean_log2=("log2", "mean"), false_gain_other=("false_gain_other", "mean")).round(3).reset_index()
out = Path("results/simulation"); out.mkdir(parents=True, exist_ok=True)
runs.round(4).to_csv(out / "allelic_spikein_runs.csv", index=False)
summ.to_csv(out / "allelic_spikein_summary.csv", index=False)
pd.set_option("display.width", 200)
print("\nGAIN rate by f:"); print(summ.pivot_table(index=["line", "arm"], columns="f", values="gain_rate").to_string())
print("\nmean log2 by f:"); print(summ.pivot_table(index=["line", "arm"], columns="f", values="mean_log2").round(2).to_string())
print("\nfirst f with GAIN rate >= 0.8:", summ[summ.gain_rate >= 0.8].groupby(["line", "arm"]).f.min().to_dict())
print("refused/too-few rate at f=1:", summ[summ.f == 1.0].set_index(["line", "arm"]).refused_or_few.to_dict())
print("false GAINs on other arms, mean per run:", round(summ.false_gain_other.mean(), 3))
