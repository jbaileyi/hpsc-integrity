"""Rebuild the bundled reference tables from the source files.

Usage (from the repository root):
    python scripts/build_reference.py data/ref/cdna.primary.fa.gz data/annotation/hg38_cytoBand.txt.gz
"""
import sys
from pathlib import Path

from hpsc_integrity.annotation import (
    arms_from_cytoband, build_gene_table, parse_cdna_headers, read_cytoband,
)

SOURCES = """Genes: protein-coding genes on chr1-22, X, Y from the headers of the Ensembl release 112
GRCh38 cDNA FASTA (Homo_sapiens.GRCh38.cdna.all.fa.gz). Chromosome, start, end, gene symbol.
Arms: first band start to last band end per arm from the UCSC hg38 cytoBand table.
Built by scripts/build_reference.py. Check each provider's terms before redistributing.
"""


def main(fasta, cytoband):
    out = Path(__file__).resolve().parents[1] / "src" / "hpsc_integrity" / "data"
    out.mkdir(parents=True, exist_ok=True)
    genes = build_gene_table(parse_cdna_headers(fasta))
    arms = arms_from_cytoband(read_cytoband(cytoband))
    genes.to_csv(out / "genes_ensembl112_grch38.csv.gz")
    arms.to_csv(out / "arms_hg38.csv", index=False)
    (out / "SOURCES.txt").write_text(SOURCES)
    print(f"{len(genes)} genes, {len(arms)} arms written to {out}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
