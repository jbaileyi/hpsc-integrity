"""Clone-fraction estimates for allelic evidence: fraction of cells carrying a +1 copy of an arm, with a CI.

A per-sample response curve is built by simulating gains of known fraction f on the sample's own normal arms
(same read-shifting model as scripts/allelic_spikein.py: reads move between alleles so a full gain, f = 1,
gives allele fraction 2/3). The arm's log2 and its bootstrap CI (from call_arms) are inverted through the curve.
The CI covers arm sampling noise only, not uncertainty in the curve. A flat curve (noisy sample) gives 0-1.
"""
import numpy as np
import pandas as pd

from .allelic import assign_snp_arms, filter_het_snps
from .allelic_calls import NORMAL, UNDERPOWERED, call_arms

F_GRID = (0.25, 0.5, 0.75, 1.0)


def apply_gain(base: pd.DataFrame, idx, f: float, rng) -> pd.DataFrame:
    """Shift reads at the SNPs at positions `idx` as if a fraction f of cells carried one extra copy. Depth unchanged."""
    d = base.copy()
    ref, alt = d["ref_n"].to_numpy()[idx], d["alt_n"].to_numpy()[idx]
    alt_gained = rng.random(len(idx)) < 0.5
    moved = rng.binomial(np.where(alt_gained, ref, alt), f / (2 + f))
    d.iloc[idx, d.columns.get_loc("ref_n")] = ref + np.where(alt_gained, -moved, moved)
    d.iloc[idx, d.columns.get_loc("alt_n")] = alt + np.where(alt_gained, moved, -moved)
    return d


def response_curve(snps, genes, arms, res, *, f_grid=F_GRID, n_donors=6, n_boot=20, seed=1, exclude=()):
    """Median change in an arm's log2 when a gain of fraction f is simulated on it. Returns (f values, log2 shifts)."""
    base = snps[["chrom", "pos", "ref_n", "alt_n"]].reset_index(drop=True)
    arm_of = assign_snp_arms(base, arms).to_numpy()
    skip = set(exclude)
    donors = [a for a in res.index if res.at[a, "call"] in (NORMAL, UNDERPOWERED) and res.at[a, "ci_lo"] <= 0
              and res.at[a, "n_loci"] >= 10 and a not in skip]
    if len(donors) < 3:
        raise ValueError(f"only {len(donors)} normal arms available to calibrate on; need at least 3")
    rng = np.random.default_rng(seed)
    fs, ys = [0.0], [0.0]
    for f in f_grid:
        deltas = []
        for d in rng.choice(donors, size=min(n_donors, len(donors)), replace=False):
            idx = np.flatnonzero(arm_of == d)
            r = call_arms(filter_het_snps(apply_gain(base, idx, f, rng)), genes, arms, n_boot=n_boot, seed=1)
            if d in r.index and pd.notna(r.at[d, "log2"]):
                deltas.append(r.at[d, "log2"] - res.at[d, "log2"])
        if not deltas:
            raise ValueError(f"calibration failed at f={f}: simulated arms were no longer assessable")
        fs.append(f)
        ys.append(float(np.median(deltas)))
    return np.array(fs), np.maximum.accumulate(np.array(ys))


def _to_fraction(y, fs, ys):
    if y <= 0:
        return 0.0
    if y >= ys[-1]:
        return 1.0
    return float(np.interp(y, ys + np.arange(len(ys)) * 1e-9, fs))


def estimate_clone_fractions(snps, genes, arms, res=None, **kw):
    """Returns (table, curve): per assessed arm f_hat, f_lo, f_hi and a note; curve maps f to the simulated log2 shift."""
    res = call_arms(snps, genes, arms) if res is None else res
    fs, ys = response_curve(snps, genes, arms, res, **kw)
    flat = ys[-1] < 0.15
    rows = []
    for a, r in res.iterrows():
        if pd.isna(r["log2"]):
            continue
        notes = (["CURVE_FLAT"] if flat else []) + (["UPPER_UNBOUNDED"] if r["ci_hi"] >= ys[-1] else [])
        rows.append(dict(arm=a, call=r["call"], log2=r["log2"], f_hat=_to_fraction(r["log2"], fs, ys),
                         f_lo=_to_fraction(r["ci_lo"], fs, ys), f_hi=_to_fraction(r["ci_hi"], fs, ys),
                         note=",".join(notes)))
    return pd.DataFrame(rows).set_index("arm").round(3), {round(float(f), 2): round(float(y), 3) for f, y in zip(fs, ys)}
