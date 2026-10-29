import numpy as np
import pandas as pd

from hpsc_integrity.allelic import (
    arm_allelic_table, assign_snp_arms, filter_het_snps, read_allele_counts,
)


def test_read_allele_counts_skips_unusable_rows(tmp_path):
    f = tmp_path / "snps.tsv"
    f.write_text("1\t100\tA\tG\t12,8\n"       # fine
                 "1\t200\tC\tT\t5,3,1\n"      # three alleles: skipped
                 "2\t50\tG\tA\t.\n"           # missing depth: skipped
                 "X\t70\tT\tC\t0,30\n")
    out = read_allele_counts(f)
    assert out[["chrom", "pos", "ref_n", "alt_n"]].values.tolist() == [["1", 100, 12, 8], ["X", 70, 0, 30]]


def test_filter_het_snps_boundaries_and_ratio():
    snps = pd.DataFrame({"chrom": "1", "pos": range(5),
                         "ref_n": [16, 17, 10, 30, 15], "alt_n": [4, 3, 9, 30, 5]})
    out = filter_het_snps(snps)
    # (16,4): depth 20, maf 0.2 -> kept. (17,3): maf 0.15 -> dropped. (10,9): depth 19 -> dropped.
    assert out["pos"].tolist() == [0, 3, 4]
    assert np.allclose(out["ratio"], [4.0, 1.0, 3.0])      # major/minor
    assert np.allclose(out["maf"], [0.2, 0.5, 0.25])


def _arms():
    rows = []
    for c in ["1", "2", "3"]:
        rows += [(c, 0, 5000, c + "p"), (c, 5000, 10000, c + "q")]
    rows += [("X", 0, 5000, "Xp"), ("X", 5000, 10000, "Xq")]
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "arm"])


def test_assign_snp_arms_boundary():
    snps = pd.DataFrame({"chrom": ["1", "1", "1", "9"], "pos": [4999, 5000, 9999, 10]})
    out = assign_snp_arms(snps, _arms())
    assert out.iloc[0] == "1p" and out.iloc[1] == "1q" and out.iloc[2] == "1q"
    assert pd.isna(out.iloc[3])


def test_arm_allelic_table_hand_computed():
    rows = []
    for arm, ratio in [("1p", 1.0), ("1q", 2.0), ("2p", 1.0), ("2q", 1.0), ("3p", 1.0), ("3q", 1.0), ("Xq", 5.0)]:
        c, side = arm[:-1], arm[-1]
        for i in range(3):
            rows.append((c, (100 if side == "p" else 6000) + i, ratio))
    rows.append(("2", 9000, 1.0))                 # 2q has 4 SNPs; no arm below the minimum yet
    snps = pd.DataFrame(rows, columns=["chrom", "pos", "ratio"])
    t, baseline = arm_allelic_table(snps, _arms(), min_snps=3)
    assert baseline == 1.0                        # median of [1, 2, 1, 1, 1, 1]; X is excluded
    assert t.loc["1q", "median_ratio"] == 2.0 and np.isclose(t.loc["1q", "log2_vs_baseline"], 1.0)
    assert np.isclose(t.loc["1p", "log2_vs_baseline"], 0.0)
    assert t.loc["1q", "n_snps"] == 3 and t.loc["2q", "n_snps"] == 4
    assert np.isnan(t.loc["Xq", "log2_vs_baseline"])          # X not assessed


def test_arm_with_too_few_snps_is_not_assessed():
    rows = [("1", 100 + i, 1.0) for i in range(5)] + [("1", 6000 + i, 2.0) for i in range(2)]
    snps = pd.DataFrame(rows, columns=["chrom", "pos", "ratio"])
    t, _ = arm_allelic_table(snps, _arms(), min_snps=5)
    assert np.isnan(t.loc["1q", "log2_vs_baseline"])
    assert t.loc["1p", "log2_vs_baseline"] == 0.0
