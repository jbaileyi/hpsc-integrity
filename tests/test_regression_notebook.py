"""The package must reproduce the notebook's arm table on the real GSE81402 data.

Skipped automatically when the local data files are absent (e.g. on CI)."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from hpsc_integrity.annotation import arms_from_cytoband, read_cytoband
from hpsc_integrity.pipeline import score_arms

ROOT = Path(__file__).resolve().parents[1]
COUNTS = ROOT / "data/processed/gene_counts.csv"
GENES = ROOT / "data/processed/genes.csv"
CYTO = ROOT / "data/annotation/hg38_cytoBand.txt.gz"
EXPECTED = ROOT / "results/arm_log2_ratios_leave_one_out.csv"


@pytest.mark.skipif(not all(p.exists() for p in [COUNTS, GENES, CYTO, EXPECTED]),
                    reason="local GSE81402 data not present")
def test_package_reproduces_notebook_arm_table():
    counts = pd.read_csv(COUNTS, index_col=0)
    genes = pd.read_csv(GENES, index_col=0, dtype={"chrom": str})
    arms = arms_from_cytoband(read_cytoband(CYTO))
    table, _, n_genes = score_arms(counts, genes, arms, min_counts=20)

    expected = pd.read_csv(EXPECTED, index_col=0)
    assert n_genes == 13174
    assert set(table.index) == set(expected.index)
    got = table.reindex(index=expected.index, columns=expected.columns)
    assert (got["n_genes"] == expected["n_genes"]).all()
    assert np.allclose(got.drop(columns="n_genes"), expected.drop(columns="n_genes"), atol=2e-4)
