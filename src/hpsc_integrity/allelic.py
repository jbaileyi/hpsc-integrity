"""Allelic imbalance from RNA-seq: heterozygous SNPs -> per-arm major/minor allele ratio.

In a diploid region a heterozygous SNP is read ~50:50. In a region present in 3 copies (a gain) the
reads split ~2:1, so the major/minor ratio rises toward 2. Sampling noise keeps real diploid ratios
above 1 (about 1.3-1.4 in the source paper), so each arm is compared with the sample's own baseline.
"""
import numpy as np
import pandas as pd

from .annotation import CHROMS

MIN_DEPTH = 20      # reads at the SNP (as in the source paper)
MIN_MAF = 0.2       # minor allele fraction: keeps sites that look heterozygous (as in the source paper)
MIN_SNPS = 50       # SNPs needed before an arm is assessed


def read_allele_counts(path) -> pd.DataFrame:
    """Read `bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n'` output.

    Rows without exactly two allele depths (missing, or more than two alleles) are skipped.
    """
    df = pd.read_csv(path, sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "ad"],
                     dtype={"chrom": str, "ad": str})
    parts = df["ad"].str.split(",")
    df = df[parts.str.len() == 2].copy()
    parts = parts[df.index]
    df["ref_n"] = pd.to_numeric(parts.str[0], errors="coerce")
    df["alt_n"] = pd.to_numeric(parts.str[1], errors="coerce")
    df = df.dropna(subset=["ref_n", "alt_n"])
    df[["ref_n", "alt_n"]] = df[["ref_n", "alt_n"]].astype(int)
    return df[["chrom", "pos", "ref_n", "alt_n"]].reset_index(drop=True)


def filter_het_snps(snps: pd.DataFrame, min_depth: int = MIN_DEPTH, min_maf: float = MIN_MAF) -> pd.DataFrame:
    """Keep SNPs with enough reads that look heterozygous; add depth, maf and major/minor ratio."""
    depth = snps["ref_n"] + snps["alt_n"]
    minor = snps[["ref_n", "alt_n"]].min(axis=1)
    major = snps[["ref_n", "alt_n"]].max(axis=1)
    keep = (depth >= min_depth) & (minor / depth >= min_maf)
    out = snps[keep].copy()
    out["depth"] = depth[keep]
    out["maf"] = (minor / depth)[keep]
    out["ratio"] = (major / minor)[keep]
    return out.reset_index(drop=True)


def assign_snp_arms(snps: pd.DataFrame, arms: pd.DataFrame) -> pd.Series:
    """Arm name for each SNP from its position. SNPs outside every arm get NaN."""
    out = pd.Series(np.nan, index=snps.index, dtype=object)
    for a in arms.itertuples():
        m = (snps["chrom"] == a.chrom) & (snps["pos"] >= a.start) & (snps["pos"] < a.end)
        out[m] = a.arm
    return out


def arm_allelic_table(snps: pd.DataFrame, arms: pd.DataFrame, min_snps: int = MIN_SNPS):
    """Per-arm median major/minor ratio and its log2 relative to the sample's own baseline.

    The baseline is the median of the autosomal arm medians (robust to a few abnormal arms).
    X and Y are not assessed. Returns (table, baseline).
    """
    arm = assign_snp_arms(snps, arms)
    g = snps.groupby(arm)["ratio"]
    t = pd.DataFrame({"n_snps": g.size(), "median_ratio": g.median()})
    sex = t.index.str.startswith(("X", "Y"))
    usable = (t["n_snps"] >= min_snps) & ~sex
    baseline = t.loc[usable, "median_ratio"].median()
    t["log2_vs_baseline"] = np.log2(t["median_ratio"] / baseline)
    t.loc[~usable, "log2_vs_baseline"] = np.nan
    t = t.loc[sorted(t.index, key=lambda a: (CHROMS.index(a[:-1]), a[-1]))]
    return t, baseline
