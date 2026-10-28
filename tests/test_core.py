import numpy as np
import pandas as pd
import pytest

from hpsc_integrity.core import (
    arm_medians, call_tiers, leave_one_out_ratio, log_cpm, remove_offset,
)


def test_log_cpm_normalizes_depth():
    # same proportions at two different depths must give identical values
    counts = pd.DataFrame({"a": [3, 1], "b": [6, 2]}, index=["g1", "g2"])
    out = log_cpm(counts)
    assert np.allclose(out["a"], out["b"])
    # 3 of 4 reads = 750,000 CPM
    assert np.isclose(out.loc["g1", "a"], np.log2(750000 + 1))


def test_leave_one_out_hand_computed():
    x = pd.DataFrame({"A": [1.0, 2.0], "B": [3.0, 4.0], "C": [5.0, 6.0]}, index=["g1", "g2"])
    out = leave_one_out_ratio(x)
    # A vs median(B, C) = 4 and 5; B vs median(A, C) = 3 and 4; C vs median(A, B) = 2 and 3
    assert np.allclose(out["A"], [-3, -3])
    assert np.allclose(out["B"], [0, 0])
    assert np.allclose(out["C"], [3, 3])


def test_leave_one_out_needs_three_samples():
    x = pd.DataFrame({"A": [1.0], "B": [2.0]}, index=["g1"])
    with pytest.raises(ValueError):
        leave_one_out_ratio(x)


def test_remove_offset_ignores_abnormal_chromosome_and_sex_chromosomes():
    # five autosomes, one of them abnormal (0.6), plus X at 5.0: offset must stay 0.1
    chrom = pd.Series(["1", "2", "3", "4", "5", "X"], index=list("abcdef"))
    ratio = pd.DataFrame({"S": [0.1, 0.1, 0.1, 0.1, 0.6, 5.0]}, index=list("abcdef"))
    corrected, offset = remove_offset(ratio, chrom)
    assert np.isclose(offset["S"], 0.1)
    assert np.isclose(corrected.loc["e", "S"], 0.5)
    assert np.isclose(corrected.loc["f", "S"], 4.9)


def test_arm_medians_and_counts():
    arm = pd.Series(["1p", "1p", "1p", "1q", "1q"], index=list("abcde"))
    ratio = pd.DataFrame({"S": [1.0, 3.0, 5.0, 2.0, 10.0]}, index=list("abcde"))
    out = arm_medians(ratio, arm)
    assert out.loc["1p", "S"] == 3.0 and out.loc["1p", "n_genes"] == 3
    assert out.loc["1q", "S"] == 6.0 and out.loc["1q", "n_genes"] == 2


def test_call_tiers_boundaries_and_exclusions():
    table = pd.DataFrame(
        {"n_genes": [100, 100, 100, 100, 10, 100, 100],
         "S": [0.25, 0.15, -0.30, 0.05, 0.50, 0.50, 0.20]},
        index=["1p", "1q", "2p", "2q", "3p", "Xp", "4p"],
    )
    out = call_tiers(table).set_index("arm")
    assert (out.loc["1p", "direction"], out.loc["1p", "tier"]) == ("gain", "high")
    assert (out.loc["1q", "direction"], out.loc["1q", "tier"]) == ("gain", "possible")
    assert (out.loc["2p", "direction"], out.loc["2p", "tier"]) == ("loss", "high")
    assert (out.loc["2q", "direction"], out.loc["2q", "tier"]) == ("none", "none")
    assert out.loc["3p", "tier"] == "not_assessed"   # too few genes
    assert out.loc["Xp", "tier"] == "not_assessed"   # X not assessed
    assert out.loc["4p", "tier"] == "high"           # exactly at the threshold counts


def test_end_to_end_synthetic_trisomy():
    # 4 chromosomes x 2 arms x 100 genes. Sample S5 has chromosome 1 at 1.5x (a trisomy).
    rng = np.random.default_rng(0)
    arms, chroms = [], []
    for c in ["1", "2", "3", "4"]:
        for a in ["p", "q"]:
            arms += [c + a] * 100
            chroms += [c] * 100
    genes = [f"g{i}" for i in range(len(arms))]
    arm = pd.Series(arms, index=genes)
    chrom = pd.Series(chroms, index=genes)
    means = rng.uniform(100, 1000, len(genes))

    counts = pd.DataFrame(index=genes, columns=[f"S{i}" for i in range(1, 6)], dtype=float)
    for s in counts.columns:
        scale = np.where((chrom == "1").values & (s == "S5"), 1.5, 1.0)
        counts[s] = rng.poisson(means * scale)

    ratio = leave_one_out_ratio(log_cpm(counts))
    ratio, _ = remove_offset(ratio, chrom)
    calls = call_tiers(arm_medians(ratio, arm))

    flagged = calls[calls["tier"] != "none"]
    assert set(zip(flagged["sample"], flagged["arm"])) == {("S5", "1p"), ("S5", "1q")}
    assert (flagged["direction"] == "gain").all()
    assert (flagged["tier"] == "high").all()
    # observed size should be close to log2(1.5) = 0.585
    assert np.allclose(flagged["log2_ratio"], np.log2(1.5), atol=0.05)
