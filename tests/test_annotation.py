import warnings

import pandas as pd
import pytest

from hpsc_integrity.annotation import (
    arms_from_cytoband, assign_arms, build_gene_table, parse_cdna_headers,
)

FASTA = (
    ">ENST00000415118.1 cdna chromosome:GRCh38:14:22438547:22438554:1 gene:ENSG00000223997.1 "
    "gene_biotype:TR_D_gene transcript_biotype:TR_D_gene gene_symbol:TRDD1 "
    "description:T cell receptor delta diversity 1 [Source:HGNC Symbol;Acc:HGNC:12254]\nACGT\n"
    ">ENST00000999999.1 cdna chromosome:GRCh38:1:100:200:-1 gene:ENSG00000000001.1 "
    "gene_biotype:protein_coding transcript_biotype:protein_coding\nACGT\n"
)


def test_parse_headers(tmp_path):
    f = tmp_path / "x.fa"
    f.write_text(FASTA)
    tx = parse_cdna_headers(f)
    assert len(tx) == 2
    assert tx.loc[0, "chrom"] == "14" and tx.loc[0, "start"] == 22438547 and tx.loc[0, "end"] == 22438554
    assert tx.loc[0, "gene_symbol"] == "TRDD1"
    assert pd.isna(tx.loc[1, "gene_symbol"])           # gene without a symbol still parses
    assert tx.loc[1, "strand"] == "-1"


def test_parse_headers_rejects_other_files(tmp_path):
    f = tmp_path / "x.fa"
    f.write_text(">not an ensembl header\nACGT\n")
    with pytest.raises(ValueError):
        parse_cdna_headers(f)


def test_build_gene_table_filters_and_merges():
    cols = ["transcript_id", "chrom", "start", "end", "strand", "gene_id",
            "gene_biotype", "transcript_biotype", "gene_symbol"]
    tx = pd.DataFrame([
        ("t1", "1", 100, 200, "1", "A", "protein_coding", "protein_coding", "GA"),
        ("t2", "1", 150, 300, "1", "A", "protein_coding", "protein_coding", "GA"),   # same gene
        ("t3", "HSCHR6_MHC_COX_CTG1", 10, 20, "1", "B", "protein_coding", "protein_coding", "GB"),
        ("t4", "MT", 10, 20, "1", "C", "protein_coding", "protein_coding", "GC"),
        ("t5", "2", 10, 20, "1", "D", "lncRNA", "lncRNA", "GD"),
        ("t6", "X", 500, 900, "-1", "E", "protein_coding", "protein_coding", None),
    ], columns=cols)
    genes = build_gene_table(tx)
    assert list(genes.index) == ["A", "E"]             # alt scaffold, MT and lncRNA dropped
    assert genes.loc["A", "start"] == 100 and genes.loc["A", "end"] == 300   # spans both transcripts
    assert genes.loc["A", "chrom"] == "1" and genes.loc["E", "chrom"] == "X"


def test_build_gene_table_rejects_gene_on_two_chromosomes():
    cols = ["transcript_id", "chrom", "start", "end", "strand", "gene_id",
            "gene_biotype", "transcript_biotype", "gene_symbol"]
    tx = pd.DataFrame([
        ("t1", "1", 100, 200, "1", "A", "protein_coding", "protein_coding", "GA"),
        ("t2", "2", 100, 200, "1", "A", "protein_coding", "protein_coding", "GA"),
    ], columns=cols)
    with pytest.raises(ValueError):
        build_gene_table(tx)


def test_arms_from_cytoband():
    cb = pd.DataFrame([
        ("chr1", 0, 2_300_000, "p36.33", "gneg"),
        ("chr1", 2_300_000, 5_000_000, "p11.1", "acen"),
        ("chr1", 5_000_000, 6_000_000, "q11", "acen"),
        ("chr1", 6_000_000, 10_000_000, "q44", "gneg"),
        ("chrUn_KI270742v1", 0, 100, "q11", "gneg"),     # unplaced contig: ignored
    ], columns=["chrom", "start", "end", "band", "stain"])
    arms = arms_from_cytoband(cb)
    assert arms["arm"].tolist() == ["1p", "1q"]
    assert arms.loc[0, "start"] == 0 and arms.loc[0, "end"] == 5_000_000
    assert arms.loc[1, "start"] == 5_000_000 and arms.loc[1, "end"] == 10_000_000


def test_assign_arms_midpoints_and_boundary():
    arms = pd.DataFrame([("1", 0, 5_000_000, "1p"), ("1", 5_000_000, 10_000_000, "1q")],
                        columns=["chrom", "start", "end", "arm"])
    genes = pd.DataFrame({
        "chrom": ["1", "1", "1", "5"],
        "start": [1_000_000, 7_000_000, 4_999_000, 1],
        "end":   [2_000_000, 8_000_000, 5_001_000, 100],      # third gene midpoint = 5,000,000 exactly
    }, index=["a", "b", "c", "d"])
    out = assign_arms(genes, arms)
    assert out["a"] == "1p" and out["b"] == "1q"
    assert out["c"] == "1q"                    # boundary midpoint belongs to the q arm
    assert pd.isna(out["d"])                   # chromosome with no arm entry
