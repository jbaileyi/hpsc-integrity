"""Spike-in simulation: how small a mosaic gain or loss can the method detect, at what depth?

A real sample is modified in silico so that a fraction f of its cells carry a gain (3 copies)
or loss (1 copy) of one chromosome arm. Reads are removed (binomial thinning), never invented,
so the real sampling noise of the sample is kept.
  gain, fraction f: mean copy number 2 + f  -> target/rest ratio rises by 1 + f/2 (others thinned)
  loss, fraction f: mean copy number 2 - f  -> target/rest ratio falls by 1 - f/2 (target thinned)
Assumes expression is proportional to copy number (the real trisomies reached ~91-94% of that).
"""
import numpy as np
import pandas as pd

from .annotation import assign_arms
from .core import call_tiers
from .pipeline import match_genes, score_arms


def thin(counts: pd.Series, p, rng) -> pd.Series:
    """Keep each read with probability p (scalar or per-gene). Counts are rounded to integers."""
    n = np.rint(counts.values).astype(np.int64)
    return pd.Series(rng.binomial(n, p), index=counts.index)


def spike_in(counts: pd.Series, gene_arm: pd.Series, arm: str, direction: str,
             f: float, depth: float, rng) -> pd.Series:
    in_arm = (gene_arm.reindex(counts.index) == arm).values
    if direction == "gain":
        keep = np.where(in_arm, 1.0, 1.0 / (1.0 + f / 2.0))
    elif direction == "loss":
        keep = np.where(in_arm, 1.0 - f / 2.0, 1.0)
    else:
        raise ValueError("direction must be 'gain' or 'loss'")
    return thin(counts, keep * depth, rng)      # thinning composes: spike-in x downsampling


def comparators_for(columns, test, chrom, known_abnormal):
    """Other samples, minus any known to be abnormal on the target chromosome."""
    known_abnormal = known_abnormal or {}
    return [c for c in columns if c != test and chrom not in known_abnormal.get(c, set())]


def run_grid(counts, genes, arms, test_samples, targets, fractions, depths, reps,
             seed=1, known_abnormal=None, min_counts=20) -> pd.DataFrame:
    """One row per simulated run. targets = [(arm, 'gain' | 'loss'), ...]."""
    counts, genes = match_genes(counts, genes)
    gene_arm = assign_arms(genes, arms)
    seeds = iter(np.random.SeedSequence(seed).spawn(10**6))
    rows = []
    for test in test_samples:
        for arm, direction in targets:
            comps = comparators_for(counts.columns, test, arm[:-1], known_abnormal)
            for f in fractions:
                for d in depths:
                    for rep in range(1 if (f == 0 and d == 1.0) else reps):   # no randomness: one run
                        rng = np.random.default_rng(next(seeds))
                        sim = counts[comps].copy()
                        sim[test] = spike_in(counts[test], gene_arm, arm, direction, f, d, rng)
                        table, _, n_genes = score_arms(sim, genes, arms, min_counts)
                        calls = call_tiers(table)
                        calls = calls[calls["sample"] == test]
                        hit = calls[calls["arm"] == arm].iloc[0]
                        rest = calls[(calls["arm"] != arm) & (calls["tier"] != "not_assessed")]
                        rows.append(dict(
                            test_sample=test, target=arm, direction=direction, f=f, depth=d, rep=rep,
                            n_comparators=len(comps), genes_used=n_genes, target_genes=int(hit["n_genes"]),
                            log2_ratio=hit["log2_ratio"], target_tier=hit["tier"],
                            target_direction=hit["direction"],
                            n_other_arms=len(rest),
                            other_high=int((rest["tier"] == "high").sum()),
                            other_possible=int((rest["tier"] == "possible").sum()),
                        ))
    return pd.DataFrame(rows)


def _wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return centre - half, centre + half


def summarize(runs: pd.DataFrame) -> pd.DataFrame:
    """Detection rate per (target, direction, f, depth): 'high' tier, and any call, in the right direction."""
    r = runs.copy()
    right = r["target_direction"] == r["direction"]
    r["det_high"] = (r["target_tier"] == "high") & right
    r["det_any"] = r["target_tier"].isin(["high", "possible"]) & right
    g = r.groupby(["target", "direction", "f", "depth"])
    out = g.agg(n_runs=("rep", "size"), k_high=("det_high", "sum"), k_any=("det_any", "sum"),
                mean_log2=("log2_ratio", "mean"), sd_log2=("log2_ratio", "std")).reset_index()
    out["rate_high"] = out["k_high"] / out["n_runs"]
    out["rate_any"] = out["k_any"] / out["n_runs"]
    ci = [_wilson(k, n) for k, n in zip(out["k_any"], out["n_runs"])]
    out["any_lo"], out["any_hi"] = [c[0] for c in ci], [c[1] for c in ci]
    return out


def detection_limits(summary: pd.DataFrame, level: float = 0.8) -> pd.DataFrame:
    """Smallest simulated mosaic fraction with detection rate >= level, per target and depth."""
    rows = []
    for (target, direction, depth), g in summary[summary["f"] > 0].groupby(["target", "direction", "depth"]):
        g = g.sort_values("f")
        ok_any = g.loc[g["rate_any"] >= level, "f"]
        ok_high = g.loc[g["rate_high"] >= level, "f"]
        rows.append(dict(target=target, direction=direction, depth=depth,
                         lod_any=ok_any.min() if len(ok_any) else np.nan,
                         lod_high=ok_high.min() if len(ok_high) else np.nan))
    return pd.DataFrame(rows)


def false_call_rates(runs: pd.DataFrame) -> pd.DataFrame:
    """False calls on arms OTHER than the spiked one, per run, by depth (uses all runs)."""
    g = runs.groupby("depth")
    return g.agg(runs=("rep", "size"), arms_per_run=("n_other_arms", "mean"),
                 false_high_per_run=("other_high", "mean"),
                 false_possible_per_run=("other_possible", "mean")).reset_index()
