from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from hpsc_integrity.allelic import filter_het_snps, read_allele_counts
from hpsc_integrity.allelic_calls import GAIN, REFUSED, call_arms
from hpsc_integrity.annotation import load_reference


def _synthetic(seed=0):
    rng = np.random.default_rng(seed)
    arms = pd.DataFrame([dict(arm=f"{c}{s}", chrom=str(c), start=st, end=en)
                         for c in range(1, 7) for s, st, en in [("p", 0, 50_000_000), ("q", 50_000_000, 150_000_000)]])
    genes = pd.DataFrame([dict(chrom=a.chrom, start=int(p), end=int(p) + 1000)
                          for a in arms.itertuples() for p in rng.integers(a.start, a.end - 1000, 800)])
    rows = []
    for a in arms.itertuples():
        n = 120 if a.arm == "5q" else 400            # 5q: het-depleted
        p = 2 / 3 if a.arm == "3q" else 0.5          # 3q: one extra copy
        depth = rng.integers(30, 120, n)
        alt = rng.binomial(depth, p)
        rows.append(pd.DataFrame(dict(chrom=a.chrom, pos=rng.integers(a.start, a.end, n), ref_n=depth - alt, alt_n=alt)))
    return filter_het_snps(pd.concat(rows, ignore_index=True)), genes, arms


def test_synthetic_gain_refusal_and_no_false_gains():
    snps, genes, arms = _synthetic()
    r = call_arms(snps, genes, arms, n_boot=500)
    assert r.at["3q", "call"] == GAIN
    assert r.at["5q", "call"] == REFUSED
    assert (r.drop(index="3q")["call"] == GAIN).sum() == 0


def _real(srr):
    p = Path("data/bam") / f"{srr}.snps.tsv"
    if not p.exists():
        pytest.skip("real SNP file not present")
    genes, arms = load_reference()
    return call_arms(filter_het_snps(read_allele_counts(p)), genes, arms)


def test_cses22_trisomy_12_is_gain():
    r = _real("SRR3501467")
    assert r.at["12p", "call"] == GAIN and r.at["12q", "call"] == GAIN


@pytest.mark.parametrize("srr", ["SRR3501465", "SRR3501469"])
def test_9q_het_depleted_is_refused(srr):
    assert _real(srr).at["9q", "call"] == REFUSED
