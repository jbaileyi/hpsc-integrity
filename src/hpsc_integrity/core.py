"""Core scoring: gene counts -> per-arm log2 ratios -> tiered calls.

All functions take and return pandas objects. Gene tables are indexed by gene ID;
count tables have genes as rows and samples as columns.
"""
import numpy as np
import pandas as pd

# Provisional thresholds on |log2 ratio| (see docs/design.md). They come from the null
# distribution of five real lines and must be recalibrated on simulation.
HIGH = 0.20
POSSIBLE = 0.12
MIN_ARM_GENES = 50
MIN_SAMPLES = 3

AUTOSOMES = {str(i) for i in range(1, 23)}


def log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """log2(counts per million + 1): corrects for sequencing depth, then compresses range."""
    cpm = counts / counts.sum() * 1e6
    return np.log2(cpm + 1)


def leave_one_out_ratio(logcpm: pd.DataFrame) -> pd.DataFrame:
    """Each sample minus the per-gene median of all the OTHER samples.

    Leaving the sample out of its own baseline avoids the exact-zero artifact of a
    median that includes the sample itself.
    """
    if logcpm.shape[1] < MIN_SAMPLES:
        raise ValueError(
            f"need at least {MIN_SAMPLES} samples to build a leave-one-out baseline, "
            f"got {logcpm.shape[1]}"
        )
    out = pd.DataFrame(index=logcpm.index, columns=logcpm.columns, dtype=float)
    for s in logcpm.columns:
        out[s] = logcpm[s] - logcpm.drop(columns=s).median(axis=1)
    return out


def remove_offset(ratio: pd.DataFrame, chrom: pd.Series):
    """Subtract each sample's global offset: the median of its per-autosome median ratios.

    The median over chromosomes is robust to one or two abnormal chromosomes.
    Returns (corrected ratios, offsets per sample).
    """
    med = ratio.groupby(chrom).median()
    med = med.loc[[c for c in med.index if c in AUTOSOMES]]
    offset = med.median()
    return ratio - offset, offset


def arm_medians(ratio: pd.DataFrame, arm: pd.Series) -> pd.DataFrame:
    """Median ratio per chromosome arm, with the number of genes behind each value."""
    med = ratio.groupby(arm).median()
    med.insert(0, "n_genes", arm.value_counts().reindex(med.index))
    return med


def call_tiers(arm_table: pd.DataFrame, min_genes: int = MIN_ARM_GENES,
               high: float = HIGH, possible: float = POSSIBLE) -> pd.DataFrame:
    """Turn arm medians into one row per (sample, arm) with a direction and confidence tier.

    tier: high (|x| >= high), possible (|x| >= possible), none, or not_assessed
    (too few genes, or X/Y, which are not assessed in v0.x).
    """
    long = (
        arm_table.drop(columns="n_genes")
        .rename_axis("arm")
        .reset_index()
        .melt(id_vars="arm", var_name="sample", value_name="log2_ratio")
    )
    long["n_genes"] = long["arm"].map(arm_table["n_genes"]).astype(int)

    mag = long["log2_ratio"].abs()
    tier = np.where(mag >= high, "high", np.where(mag >= possible, "possible", "none"))
    direction = np.where(tier == "none", "none",
                         np.where(long["log2_ratio"] > 0, "gain", "loss"))
    long["tier"] = tier
    long["direction"] = direction

    skip = (long["n_genes"] < min_genes) | long["arm"].str.startswith(("X", "Y"))
    long.loc[skip, ["tier", "direction"]] = "not_assessed"
    return long[["sample", "arm", "n_genes", "log2_ratio", "direction", "tier"]]
