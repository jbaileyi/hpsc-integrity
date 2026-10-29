import pytest

from hpsc_integrity.allelic_report import parse_specs, run_allelic
from test_allelic_calls import _synthetic


def _write(snps, path):
    snps.assign(ref="A", alt="C", ad=snps.ref_n.astype(str) + "," + snps.alt_n.astype(str))[
        ["chrom", "pos", "ref", "alt", "ad"]].to_csv(path, sep="\t", header=False, index=False)


def test_run_allelic_synthetic(tmp_path):
    snps, genes, arms = _synthetic()
    f = tmp_path / "s.snps.tsv"; _write(snps, f)
    calls, qc = run_allelic([f"S1={f}"], genes, arms)
    c = calls.set_index("arm")["call"]
    assert c["3q"] == "GAIN" and c["5q"] == "REFUSED_HET_DEPLETED"
    assert qc.at[0, "sample"] == "S1" and qc.at[0, "flags"] == ""


def test_high_baseline_is_flagged(tmp_path):
    snps, genes, arms = _synthetic()
    f = tmp_path / "s.snps.tsv"; _write(snps, f)
    _, qc = run_allelic([f"S1={f}"], genes, arms, baseline_warn=1.0)
    assert qc.at[0, "flags"] == "ALLELIC_POWER_LOW"


def test_bad_spec_rejected():
    with pytest.raises(ValueError):
        parse_specs(["no_equals_sign"])
