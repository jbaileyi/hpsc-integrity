import numpy as np

from hpsc_integrity.allelic_calls import call_arms
from hpsc_integrity.clone_fraction import apply_gain, estimate_clone_fractions
from test_allelic_calls import _synthetic


def test_apply_gain_keeps_depth_and_raises_imbalance():
    snps, _, _ = _synthetic()
    base = snps[["chrom", "pos", "ref_n", "alt_n"]].reset_index(drop=True)
    idx = np.arange(300)
    mod = apply_gain(base, idx, 1.0, np.random.default_rng(0))
    assert np.array_equal((mod.ref_n + mod.alt_n).to_numpy(), (base.ref_n + base.alt_n).to_numpy())
    ratio = lambda d: (np.maximum(d.ref_n, d.alt_n) / np.minimum(d.ref_n, d.alt_n))[idx].mean()
    assert ratio(mod) > ratio(base)


def test_synthetic_full_gain_estimate():
    snps, genes, arms = _synthetic()
    res = call_arms(snps, genes, arms, n_boot=300)
    est, curve = estimate_clone_fractions(snps, genes, arms, res, f_grid=(0.5, 1.0), n_donors=3, n_boot=20)
    assert est.at["3q", "f_hi"] == 1.0 and est.at["3q", "f_lo"] > 0.3
    assert (est.drop(index="3q")["f_lo"] > 0.3).sum() == 0
    assert curve[1.0] > curve[0.5] > 0


def test_flat_curve_gives_no_estimate():
    snps, genes, arms = _synthetic()
    res = call_arms(snps, genes, arms, n_boot=300)
    est, _ = estimate_clone_fractions(snps, genes, arms, res, f_grid=(0.5, 1.0), n_donors=3, n_boot=20, min_curve=10.0)
    assert est["f_hat"].isna().all() and (est["note"] == "UNRESOLVABLE_NOISY_SAMPLE").all()
