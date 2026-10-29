"""Per-arm allelic imbalance for one sample. Run from the repository root:
    python scripts/allelic_profile.py data/bam/SRR3501467.snps.tsv CSES22
Prints the arm table and writes results/allelic/<name>_arms.csv.
"""
import sys
from pathlib import Path

import pandas as pd

from hpsc_integrity.allelic import arm_allelic_table, filter_het_snps, read_allele_counts
from hpsc_integrity.annotation import load_reference

if len(sys.argv) != 3:
    sys.exit(__doc__)
path, name = sys.argv[1], sys.argv[2]
raw = read_allele_counts(path)
snps = filter_het_snps(raw)
_, arms = load_reference()
table, baseline = arm_allelic_table(snps, arms)

ROOT = Path(__file__).resolve().parents[1]
out = ROOT / "results/allelic"
out.mkdir(parents=True, exist_ok=True)
table.round(4).to_csv(out / f"{name}_arms.csv")

pd.set_option("display.width", 200)
print(f"{name}: {len(raw)} SNP rows read, {len(snps)} kept (depth >= 20, minor allele fraction >= 0.2)")
print(f"baseline (median of autosomal arm medians): {baseline:.3f}")
print(table.round(3).to_string())
