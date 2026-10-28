"""Reference annotation: Ensembl cDNA FASTA headers -> gene table; UCSC cytoBand -> chromosome arms."""
import gzip
import re
import warnings
from importlib.resources import files

import numpy as np
import pandas as pd

CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]

_HEADER = re.compile(
    r"^>(\S+) \S+ \S+:GRCh38:([^:]+):(\d+):(\d+):(-?1) "
    r"gene:(\S+) gene_biotype:(\S+) transcript_biotype:(\S+)"
    r"(?: gene_symbol:(\S+))?"
)
_TX_COLUMNS = ["transcript_id", "chrom", "start", "end", "strand",
               "gene_id", "gene_biotype", "transcript_biotype", "gene_symbol"]


def _open_text(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def parse_cdna_headers(path) -> pd.DataFrame:
    """Read the header lines of an Ensembl cDNA FASTA into a transcript table."""
    rows, unmatched = [], 0
    with _open_text(path) as f:
        for line in f:
            if line.startswith(">"):
                m = _HEADER.match(line)
                if m:
                    rows.append(m.groups())
                else:
                    unmatched += 1
    if not rows:
        raise ValueError(f"no Ensembl cDNA headers could be parsed from {path}")
    if unmatched:
        warnings.warn(f"{unmatched} FASTA headers did not match the expected Ensembl format")
    tx = pd.DataFrame(rows, columns=_TX_COLUMNS)
    tx["start"] = tx["start"].astype(int)
    tx["end"] = tx["end"].astype(int)
    return tx


def build_gene_table(tx: pd.DataFrame) -> pd.DataFrame:
    """One row per protein-coding gene on chromosomes 1-22, X, Y: chrom, start, end, symbol.

    Mitochondrial, alt-scaffold/patch, pseudogene and other non-coding genes are dropped
    (rationale in docs/design.md and docs/lessons.md).
    """
    keep = tx[tx["chrom"].isin(CHROMS) & (tx["gene_biotype"] == "protein_coding")]
    genes = keep.groupby("gene_id").agg(
        chrom=("chrom", "first"),
        start=("start", "min"),
        end=("end", "max"),
        symbol=("gene_symbol", "first"),
        n_chrom=("chrom", "nunique"),
    )
    if (genes["n_chrom"] != 1).any():
        bad = genes.index[genes["n_chrom"] != 1][:3].tolist()
        raise ValueError(f"genes found on more than one chromosome, e.g. {bad}")
    return genes.drop(columns="n_chrom")


def read_cytoband(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None,
                       names=["chrom", "start", "end", "band", "stain"])


def arms_from_cytoband(cb: pd.DataFrame) -> pd.DataFrame:
    """Chromosome arm boundaries: first band start to last band end for each p and q arm."""
    cb = cb.copy()
    cb["chrom"] = cb["chrom"].str.replace("^chr", "", regex=True)
    cb = cb[cb["chrom"].isin(CHROMS)].copy()
    cb["side"] = cb["band"].str[0]
    arms = cb.groupby(["chrom", "side"]).agg(start=("start", "min"), end=("end", "max")).reset_index()
    arms["arm"] = arms["chrom"] + arms["side"]
    arms["chrom_order"] = arms["chrom"].map({c: i for i, c in enumerate(CHROMS)})
    arms = arms.sort_values(["chrom_order", "side"])
    return arms[["chrom", "start", "end", "arm"]].reset_index(drop=True)


def assign_arms(genes: pd.DataFrame, arms: pd.DataFrame) -> pd.Series:
    """Arm name for each gene, using the gene midpoint. Genes outside every arm get NaN."""
    mid = (genes["start"] + genes["end"]) / 2
    out = pd.Series(np.nan, index=genes.index, dtype=object)
    for a in arms.itertuples():
        m = (genes["chrom"] == a.chrom) & (mid >= a.start) & (mid < a.end)
        out[m] = a.arm
    return out


def load_reference(genes_path=None, arms_path=None):
    """Load the gene and arm tables, from the bundled data unless paths are given."""
    data = files("hpsc_integrity") / "data"
    genes_path = genes_path or str(data / "genes_ensembl112_grch38.csv.gz")
    arms_path = arms_path or str(data / "arms_hg38.csv")
    genes = pd.read_csv(genes_path, index_col=0, dtype={"chrom": str})
    arms = pd.read_csv(arms_path, dtype={"chrom": str})
    return genes, arms
