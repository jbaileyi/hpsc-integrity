"""Validate clone-fraction CIs. python scripts/clonefrac_validate.py [reps]   (try 1 first; reps=3 takes roughly 15-20 min)
Part A: spike-in on CSES9, HUES14, CSES7. Part B: real trisomies. Part C: null arms in all five lines."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from hpsc_integrity.allelic import assign_snp_arms, filter_het_snps, read_allele_counts
from hpsc_integrity.allelic_calls import call_arms
from hpsc_integrity.annotation import load_reference
from hpsc_integrity.clone_fraction import apply_gain, estimate_clone_fractions

reps = int(sys.argv[1]) if len(sys.argv) > 1 else 1
SRR = {"CSES9": "SRR3501465", "CSES7": "SRR3501466", "CSES22": "SRR3501467", "CSES21": "SRR3501468", "HUES14": "SRR3501469"}
KNOWN = {"CSES22": ["12p", "12q"], "CSES21": ["21q"]}
ARMS, FRACS = ["1q", "12q", "20p"], [0.3, 0.5, 0.75, 1.0]
genes, arms = load_reference()
pd.set_option("display.width", 200)
load = lambda l: filter_het_snps(read_allele_counts(f"data/bam/{SRR[l]}.snps.tsv"))

rows = []
for li, line in enumerate(["CSES9", "HUES14", "CSES7"]):
    base = load(line)[["chrom", "pos", "ref_n", "alt_n"]]
    arm_of = assign_snp_arms(base, arms).to_numpy()
    for ai, arm in enumerate(ARMS):
        idx = np.flatnonzero(arm_of == arm)
        for fi, f in enumerate(FRACS):
            for rep in range(reps):
                mod = filter_het_snps(apply_gain(base, idx, f, np.random.default_rng([rep, fi, ai, li, 99])))
                res = call_arms(mod, genes, arms, n_boot=400, seed=rep)
                est, _ = estimate_clone_fractions(mod, genes, arms, res, exclude=[arm], seed=rep)
                if arm in est.index:
                    e = est.loc[arm]
                    rows.append(dict(line=line, arm=arm, f=f, rep=rep, call=e["call"], f_hat=e["f_hat"], f_lo=e["f_lo"], f_hi=e["f_hi"], note=e["note"]))
    print(line, "done", flush=True)
A = pd.DataFrame(rows)
A["resolved"] = A.f_hat.notna(); A["covered"] = (A.f_lo <= A.f) & (A.f <= A.f_hi); A["width"] = A.f_hi - A.f_lo; A["abs_err"] = (A.f_hat - A.f).abs()
out = Path("results/simulation"); out.mkdir(parents=True, exist_ok=True); A.to_csv(out / "clonefrac_spikein.csv", index=False)
print("\nPart A: by line and f (coverage, mean abs error, median CI width)")
R = A[A.resolved]
print(R.groupby(["line", "f"]).agg(coverage=("covered", "mean"), mae=("abs_err", "mean"), width=("width", "median"), n=("f", "size")).round(2).to_string())
print("fraction of runs resolved:", A.groupby("line").resolved.mean().round(2).to_dict())
print("\nCoverage among resolved runs:", round(R.covered.mean(), 3) if len(R) else "n/a", "| resolved:", len(R), "of", len(A))

print("\nPart B: real trisomies"); nulls = {}
for line in SRR:
    snps = load(line); res = call_arms(snps, genes, arms, n_boot=400)
    est, curve = estimate_clone_fractions(snps, genes, arms, res)
    if line in KNOWN: print(line, curve); print(est.loc[KNOWN[line]].to_string())
    nulls[line] = est[(est.call != "GAIN") & (est.f_lo > 0.1)]
print("\nPart C: null arms with f_lo > 0.1:", {l: dict(zip(d.index, d.f_lo)) for l, d in nulls.items() if len(d)}, "| total:", sum(len(d) for d in nulls.values()))
