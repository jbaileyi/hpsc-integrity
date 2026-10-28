import numpy as np
import pandas as pd
import pytest
from synthetic import make_synthetic

from hpsc_integrity.core import call_tiers
from hpsc_integrity.pipeline import filter_genes, gene_counts_from_kallisto, match_genes, score_arms


def test_filter_genes_boundary():
    counts = pd.DataFrame({"a": [20, 19, 500], "b": [20, 500, 500]}, index=["g1", "g2", "g3"])
    assert list(filter_genes(counts, 20).index) == ["g1", "g3"]   # 19 in one sample drops g2; 20 stays


def test_match_genes_exact_and_without_version():
    genes = pd.DataFrame({"chrom": ["1", "1"]}, index=["ENSG1.5", "ENSG2.7"])
    exact = pd.DataFrame({"s": [1, 2]}, index=["ENSG1.5", "ENSG2.7"])
    c, g = match_genes(exact, genes)
    assert list(c.index) == list(g.index) and len(c) == 2
    novers = pd.DataFrame({"s": [1, 2]}, index=["ENSG1", "ENSG2"])
    c, g = match_genes(novers, genes)
    assert sorted(c.index) == ["ENSG1", "ENSG2"] and sorted(g.index) == ["ENSG1", "ENSG2"]


def test_match_genes_rejects_symbols():
    genes = pd.DataFrame({"chrom": ["1", "1"]}, index=["ENSG1.5", "ENSG2.7"])
    symbols = pd.DataFrame({"s": [1, 2]}, index=["TP53", "MYC"])
    with pytest.raises(ValueError):
        match_genes(symbols, genes)


def test_score_arms_finds_synthetic_trisomy():
    counts, genes, arms = make_synthetic()
    table, offsets, n = score_arms(counts, genes, arms)
    assert n == 800
    assert list(table.index) == ["1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q"]
    calls = call_tiers(table)
    flagged = calls[calls["tier"] != "none"]
    assert set(zip(flagged["sample"], flagged["arm"])) == {("S5", "1p"), ("S5", "1q")}
    assert np.allclose(flagged["log2_ratio"], np.log2(1.5), atol=0.05)


def test_score_arms_needs_three_samples():
    counts, genes, arms = make_synthetic(n_samples=2)
    with pytest.raises(ValueError):
        score_arms(counts, genes, arms)


def test_gene_counts_from_kallisto_sums_transcripts_and_drops_unwanted(tmp_path):
    cols = ["transcript_id", "chrom", "start", "end", "strand", "gene_id",
            "gene_biotype", "transcript_biotype", "gene_symbol"]
    tx = pd.DataFrame([
        ("t1", "1", 1, 2, "1", "GA", "protein_coding", "protein_coding", "A"),
        ("t2", "1", 1, 2, "1", "GA", "protein_coding", "protein_coding", "A"),
        ("t3", "MT", 1, 2, "1", "GM", "protein_coding", "protein_coding", "M"),   # dropped
        ("t4", "2", 1, 2, "1", "GL", "lncRNA", "lncRNA", "L"),                   # dropped
        ("t5", "X", 1, 2, "1", "GX", "protein_coding", "protein_coding", "X"),
    ], columns=cols)
    for name, vals in {"s1": [1.5, 2.5, 100, 100, 7.0], "s2": [10, 20, 5, 5, 0.0]}.items():
        d = tmp_path / name
        d.mkdir()
        pd.DataFrame({"target_id": ["t1", "t2", "t3", "t4", "t5"], "length": 1, "eff_length": 1,
                      "est_counts": vals, "tpm": 0}).to_csv(d / "abundance.tsv", sep="\t", index=False)
    out = gene_counts_from_kallisto({"s1": tmp_path / "s1", "s2": tmp_path / "s2"}, tx)
    assert sorted(out.index) == ["GA", "GX"]
    assert out.loc["GA", "s1"] == 4.0 and out.loc["GA", "s2"] == 30.0     # 1.5 + 2.5; 10 + 20
    assert out.loc["GX", "s1"] == 7.0 and out.loc["GX", "s2"] == 0.0      # zero-count gene stays
