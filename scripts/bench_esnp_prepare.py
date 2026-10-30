"""Build blinded eSNP-Karyotyping input tables from our SNP calls. python scripts/bench_esnp_prepare.py
Writes results/benchmark/esnp_input/S1..S5.tsv; the S-id -> line key goes to ~/benchmark_key.csv (outside the repo). Unseeded shuffle."""
from pathlib import Path
import numpy as np, pandas as pd
from hpsc_integrity.allelic import read_allele_counts

SRR = {"CSES9": "SRR3501465", "CSES7": "SRR3501466", "CSES22": "SRR3501467", "CSES21": "SRR3501468", "HUES14": "SRR3501469"}
CODE = {**{str(i): i for i in range(1, 23)}, "X": 23, "Y": 24}
out = Path("results/benchmark/esnp_input"); out.mkdir(parents=True, exist_ok=True)
lines = list(SRR); np.random.default_rng().shuffle(lines)
key = []
for i, line in enumerate(lines, 1):
    d = read_allele_counts(f"data/bam/{SRR[line]}.snps.tsv")
    d = d[d.chrom.isin(CODE)]
    t = pd.DataFrame({"chr": d.chrom.map(CODE).astype(int), "position": d.pos, "AD1": d.ref_n, "AD2": d.alt_n,
                      "DP": d.ref_n + d.alt_n, "Karyotape": 0})   # Karyotape is unused by MajorMinorCalc
    t.to_csv(out / f"S{i}.tsv", sep="\t", index=False)
    key.append((f"S{i}", line))
pd.DataFrame(key, columns=["id", "line"]).to_csv(Path.home() / "benchmark_key.csv", index=False)
print("wrote", len(key), "blinded tables; key saved outside the repo")
