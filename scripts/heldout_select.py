"""Pre-registered selection of HipSci lines (rule in docs/design.md, held-out protocol 2). python scripts/heldout_select.py"""
import numpy as np, pandas as pd

SEED, N, MIN_PAIRS = 20261030, 40, 20_000_000
m = pd.read_csv("data/heldout/hipsci_meta.tsv", sep="\t")
runs = pd.read_csv("data/heldout/hipsci_runs.tsv", sep="\t")[["run_accession", "fastq_ftp", "fastq_bytes"]]
m = m[(m.library_layout == "PAIRED") & (m.library_strategy == "RNA-Seq") & m.sample_alias.str.match(r"^HPSI\d{4}i-[a-z]{4}_\d+$")]
print(len(m), "paired RNA-seq runs;", m.sample_alias.nunique(), "iPSC lines")
best = m.sort_values("read_count", ascending=False).drop_duplicates("sample_alias")
best = best[best.read_count >= MIN_PAIRS].copy()
best["donor"] = best.sample_alias.str.split("_").str[0]
print(len(best), "lines with >=", MIN_PAIRS, "read pairs;", best.donor.nunique(), "donors")
rng = np.random.default_rng(SEED)
best = best.iloc[rng.permutation(len(best))]
pick = best.drop_duplicates("donor").iloc[:N].merge(runs, on="run_accession", how="left")
pick.to_csv("results/heldout/selection.csv", index=False)
gb = pick.fastq_bytes.fillna("0").str.split(";").apply(lambda x: sum(int(i) for i in x)).sum() / 1e9
print(len(pick), "selected | read pairs median", int(pick.read_count.median()), "| total FASTQ", round(gb, 1), "GB")
print("runs with no FASTQ link:", int(pick.fastq_ftp.isna().sum()))
print(pick[["sample_alias", "run_accession", "read_count"]].head(8).to_string(index=False))
