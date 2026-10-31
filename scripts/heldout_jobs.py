"""Build data/heldout/jobs.tsv (run, line, url1, url2, md5_1, md5_2) from the pre-registered selection. python scripts/heldout_jobs.py"""
import pandas as pd
sel = pd.read_csv("results/heldout/selection.csv")[["run_accession", "sample_alias"]]
runs = pd.read_csv("data/heldout/hipsci_runs.tsv", sep="\t")[["run_accession", "fastq_ftp", "fastq_md5"]]
rows = []
for r in sel.merge(runs, on="run_accession").itertuples():
    pairs = dict(zip(r.fastq_ftp.split(";"), r.fastq_md5.split(";")))
    u1 = [u for u in pairs if u.endswith("_1.fastq.gz")]
    u2 = [u for u in pairs if u.endswith("_2.fastq.gz")]
    assert len(u1) == 1 and len(u2) == 1, f"{r.run_accession}: unexpected files {list(pairs)}"
    rows.append((r.run_accession, r.sample_alias, u1[0], u2[0], pairs[u1[0]], pairs[u2[0]]))
pd.DataFrame(rows).to_csv("data/heldout/jobs.tsv", sep="\t", header=False, index=False)
print(len(rows), "jobs written")
