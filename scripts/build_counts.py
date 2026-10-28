"""Build one gene-by-sample count table from kallisto output for every run in an SRA run table.

Usage (from the repository root):
    python scripts/build_counts.py data/raw/SraRunTable.csv data/quant/primary \
        data/ref/cdna.primary.fa.gz data/processed/gene_counts_all13.csv
Sample names are the first word of the run table's source_name column (e.g. 'CSES9').
"""
import sys
from pathlib import Path

import pandas as pd

from hpsc_integrity.annotation import parse_cdna_headers
from hpsc_integrity.pipeline import gene_counts_from_kallisto


def main(run_table, quant_root, fasta, out):
    runs = pd.read_csv(run_table)
    names = runs["source_name"].str.split().str[0]
    if not names.is_unique:
        sys.exit(f"sample names are not unique: {names[names.duplicated()].tolist()}")
    dirs = {n: Path(quant_root) / r for n, r in zip(names, runs["Run"])}
    missing = [str(d) for d in dirs.values() if not (d / "abundance.tsv").exists()]
    if missing:
        sys.exit("missing kallisto output:\n" + "\n".join(missing))
    counts = gene_counts_from_kallisto(dirs, parse_cdna_headers(fasta))
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out)
    print(f"{counts.shape[0]} genes x {counts.shape[1]} samples written to {out}")
    print("samples:", ", ".join(counts.columns))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
