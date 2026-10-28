"""End-to-end scoring: a gene-by-sample count table -> per-arm log2 ratios."""
from pathlib import Path

import pandas as pd

from .annotation import CHROMS, assign_arms
from .core import arm_medians, leave_one_out_ratio, log_cpm, remove_offset

MIN_COUNTS = 20


def filter_genes(counts: pd.DataFrame, min_counts: int = MIN_COUNTS) -> pd.DataFrame:
    """Keep genes with at least min_counts estimated reads in EVERY sample."""
    return counts[(counts >= min_counts).all(axis=1)]


def _strip_version(index: pd.Index) -> pd.Index:
    return index.astype(str).str.replace(r"\.\d+$", "", regex=True)


def match_genes(counts: pd.DataFrame, genes: pd.DataFrame):
    """Align a count table to the gene table by gene ID.

    Tries exact IDs first, then IDs without the Ensembl version suffix (ENSG...19 -> ENSG...).
    Raises a clear error if fewer than half of the count-table genes can be matched.
    """
    counts = counts.copy()
    counts.index = counts.index.astype(str)
    shared = counts.index.intersection(genes.index)
    if len(shared) / len(counts) >= 0.5:
        return counts.loc[shared], genes.loc[shared]

    c2, g2 = counts.copy(), genes.copy()
    c2.index, g2.index = _strip_version(c2.index), _strip_version(g2.index)
    if c2.index.is_unique and g2.index.is_unique:
        shared = c2.index.intersection(g2.index)
        if len(shared) / len(c2) >= 0.5:
            return c2.loc[shared], g2.loc[shared]

    raise ValueError(
        f"only {len(counts.index.intersection(genes.index))} of {len(counts)} gene IDs match the "
        "reference. The first column must hold Ensembl gene IDs (ENSG...), not symbols."
    )


def _sorted_arms(names):
    return sorted(names, key=lambda a: (CHROMS.index(a[:-1]), a[-1]))


def score_arms(counts: pd.DataFrame, genes: pd.DataFrame, arms: pd.DataFrame,
               min_counts: int = MIN_COUNTS):
    """Return (arm table, per-sample offsets, number of genes used)."""
    counts, genes = match_genes(counts, genes)
    counts = filter_genes(counts, min_counts)
    if len(counts) == 0:
        raise ValueError(f"no genes have >= {min_counts} estimated reads in every sample")
    genes = genes.loc[counts.index]

    ratio = leave_one_out_ratio(log_cpm(counts))
    ratio, offsets = remove_offset(ratio, genes["chrom"])
    table = arm_medians(ratio, assign_arms(genes, arms))
    return table.loc[_sorted_arms(table.index)], offsets, len(counts)


def gene_counts_from_kallisto(quant_dirs: dict, tx: pd.DataFrame) -> pd.DataFrame:
    """Sum kallisto est_counts over each gene's transcripts.

    quant_dirs: {sample name: folder containing abundance.tsv}; tx: transcript table from
    parse_cdna_headers. Only protein-coding genes on chr1-22, X, Y are kept (see docs/design.md).
    """
    keep = tx[tx["chrom"].isin(CHROMS) & (tx["gene_biotype"] == "protein_coding")]
    tx_to_gene = keep.set_index("transcript_id")["gene_id"]
    cols = {}
    for name, folder in quant_dirs.items():
        q = pd.read_csv(Path(folder) / "abundance.tsv", sep="\t", index_col=0)
        gene = tx_to_gene.reindex(q.index)                 # NaN for transcripts we do not keep
        cols[name] = q["est_counts"].groupby(gene).sum()   # NaN keys are dropped
    return pd.DataFrame(cols)
