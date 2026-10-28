import numpy as np
import pandas as pd
import pytest
from synthetic import make_synthetic

from hpsc_integrity.simulate import (
    comparators_for, detection_limits, false_call_rates, run_grid, spike_in, summarize, thin,
)


def test_thin_extremes_and_mean():
    rng = np.random.default_rng(0)
    c = pd.Series([1000, 2000, 3000], index=list("abc"))
    assert (thin(c, 1.0, rng) == c).all()
    assert (thin(c, 0.0, rng) == 0).all()
    big = pd.Series(np.full(2000, 1000))
    assert abs(thin(big, 0.5, rng).mean() - 500) < 5          # binomial mean


def test_spike_in_gain_and_loss_change_the_ratio_as_designed():
    rng = np.random.default_rng(0)
    counts = pd.Series(np.full(200, 100_000), index=[f"g{i}" for i in range(200)])
    arm = pd.Series(["1p"] * 100 + ["1q"] * 100, index=counts.index)
    ratio = lambda c: c[:100].sum() / c[100:].sum()
    gain = spike_in(counts, arm, "1p", "gain", f=1.0, depth=1.0, rng=rng)
    assert abs(ratio(gain) / ratio(counts) - 1.5) < 0.01       # f=1 gain: 3 copies vs 2
    loss = spike_in(counts, arm, "1p", "loss", f=0.5, depth=1.0, rng=rng)
    assert abs(ratio(loss) / ratio(counts) - 0.75) < 0.01      # half the cells lose one copy
    both = spike_in(counts, arm, "1p", "gain", f=0.0, depth=0.5, rng=rng)
    assert abs(both.sum() / counts.sum() - 0.5) < 0.01         # depth only


def test_comparators_exclude_known_abnormal():
    cols = ["A", "B", "C", "D"]
    assert comparators_for(cols, "A", "12", {"C": {"12"}}) == ["B", "D"]
    assert comparators_for(cols, "A", "5", {"C": {"12"}}) == ["B", "C", "D"]


def test_run_grid_detects_full_gain_not_small_or_null():
    counts, genes, arms = make_synthetic(trisomic_sample="S5", trisomic_chrom="4")
    runs = run_grid(counts, genes, arms, ["S1"], [("1p", "gain")], [0, 0.1, 1.0], [1.0],
                    reps=3, seed=1, known_abnormal={"S5": {"4"}})
    s = summarize(runs).set_index("f")
    assert s.loc[0.0, "rate_any"] == 0.0           # control: no false call on the target arm
    assert s.loc[0.1, "rate_high"] == 0.0          # expected log2 = 0.07, below the 0.12 tier
    assert s.loc[1.0, "rate_high"] == 1.0          # full trisomy, 0.58
    assert runs.loc[runs["f"] == 0, "rep"].nunique() == 1   # deterministic control runs once
    assert np.isclose(s.loc[1.0, "mean_log2"], np.log2(1.5), atol=0.05)


def test_summarize_limits_and_false_calls_hand_computed():
    base = dict(test_sample="S", target="1q", direction="gain", depth=1.0, n_comparators=4,
                genes_used=100, target_genes=100, n_other_arms=38)
    rows = []
    for f, tiers in [(0.1, ["none"] * 4), (0.3, ["possible", "possible", "none", "none"]),
                     (0.5, ["high", "high", "high", "possible"])]:
        for i, t in enumerate(tiers):
            rows.append({**base, "f": f, "rep": i, "log2_ratio": 0.0, "target_tier": t,
                         "target_direction": "gain" if t != "none" else "none",
                         "other_high": 0, "other_possible": 1 if i == 0 else 0})
    runs = pd.DataFrame(rows)
    s = summarize(runs).set_index("f")
    assert s.loc[0.3, "rate_any"] == 0.5 and s.loc[0.3, "rate_high"] == 0.0
    assert s.loc[0.5, "rate_any"] == 1.0 and s.loc[0.5, "rate_high"] == 0.75
    lod = detection_limits(summarize(runs), level=0.8).iloc[0]
    assert lod["lod_any"] == 0.5                     # first fraction with >= 80% any-call
    assert np.isnan(lod["lod_high"])                 # 0.75 < 0.8 at every fraction
    fc = false_call_rates(runs).iloc[0]
    assert np.isclose(fc["false_possible_per_run"], 3 / 12)   # one possible in 3 of 12 runs
