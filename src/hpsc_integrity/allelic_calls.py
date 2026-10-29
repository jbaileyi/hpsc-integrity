"""Arm-level GAIN calls from allelic imbalance, with locus-level bootstrap CIs and explicit refusal states.

Allelic imbalance can support a gain only: a loss turns heterozygous sites homozygous, which the het filter
removes (it shows up as het depletion, not as a shifted ratio). Hence no LOSS state.

SNPs within `cluster_bp` of each other are collapsed to one locus (they sit in one gene, so they are not
independent), loci are resampled within each arm, and the baseline is recomputed in every replicate.
The gain-fraction columns use a crude linear model (ratio ~ 1 + f/2) fitted by assumption, not simulation.
"""
import numpy as np
import pandas as pd

from .allelic import MIN_SNPS, assign_snp_arms

GAIN, NORMAL, UNDERPOWERED = "GAIN", "NORMAL", "UNDERPOWERED"
REFUSED, TOO_FEW = "REFUSED_HET_DEPLETED", "TOO_FEW_SNPS"


def genes_per_arm(genes: pd.DataFrame, arms: pd.DataFrame) -> pd.Series:
    mid = (genes["start"] + genes["end"]) // 2
    return pd.Series({a.arm: int(((genes["chrom"] == a.chrom) & (mid >= a.start) & (mid < a.end)).sum())
                      for a in arms.itertuples()})


def call_arms(snps: pd.DataFrame, genes: pd.DataFrame, arms: pd.DataFrame, *, min_snps: int = MIN_SNPS,
              min_loci: int = 3, min_rel_density: float = 0.5, effect: float = 0.10,
              cluster_bp: int = 100_000, n_boot: int = 2000, seed: int = 1) -> pd.DataFrame:
    """`snps` is the output of filter_het_snps. Returns one row per autosomal arm."""
    d = snps.assign(arm=assign_snp_arms(snps, arms)).dropna(subset=["arm"])
    d = d[~d["arm"].str.startswith(("X", "Y"))].sort_values(["chrom", "pos"])
    new = (d["chrom"] != d["chrom"].shift()) | (d["arm"] != d["arm"].shift()) | (d["pos"].diff() > cluster_bp)
    d["locus"] = new.cumsum()
    loci = d.groupby("locus").agg(arm=("arm", "first"), ratio=("ratio", "median"))
    G = {a: g.to_numpy() for a, g in loci.groupby("arm")["ratio"]}

    ng = genes_per_arm(genes, arms)
    names = [a for a in ng.index if not a.startswith(("X", "Y"))]
    n_snps = d.groupby("arm").size().reindex(names, fill_value=0)
    dens = n_snps / ng[names].replace(0, np.nan)
    rel = dens / dens[ng[names] >= 30].median()

    state = {}
    for a in names:
        if n_snps[a] < min_snps or len(G.get(a, [])) < min_loci:
            state[a] = TOO_FEW
        elif rel[a] < min_rel_density:
            state[a] = REFUSED
    elig = [a for a in names if a not in state]
    if len(elig) < 5:
        raise ValueError(f"only {len(elig)} arms are assessable; need at least 5 to form a baseline")

    rng = np.random.default_rng(seed)
    M = np.column_stack([np.median(G[a][rng.integers(0, len(G[a]), (n_boot, len(G[a])))], axis=1) for a in elig])
    B = np.log2(M / np.median(M, axis=1, keepdims=True))
    lo, hi = np.percentile(B, [2.5, 97.5], axis=0)
    pt = np.array([np.median(G[a]) for a in elig])

    out = pd.DataFrame(index=names)
    out["n_snps"] = n_snps
    out["n_loci"] = [len(G.get(a, [])) for a in names]
    out["n_genes"] = ng[names]
    out["rel_het_density"] = rel
    out["log2"] = pd.Series(np.log2(pt / np.median(pt)), index=elig)
    out["ci_lo"] = pd.Series(lo, index=elig)
    out["ci_hi"] = pd.Series(hi, index=elig)
    out["gain_frac_hi95"] = np.clip(2 * (2 ** out["ci_hi"] - 1), 0, 1)
    call = pd.Series(state)
    for a in elig:
        call[a] = GAIN if out.at[a, "ci_lo"] > effect else NORMAL if out.at[a, "ci_hi"] < effect else UNDERPOWERED
    out["call"] = call.reindex(names)
    return out.round(3)
