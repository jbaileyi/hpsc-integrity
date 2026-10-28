import pandas as pd
import pytest
from synthetic import make_synthetic

from hpsc_integrity.cli import main


def _write(tmp_path, n_samples=5):
    counts, genes, arms = make_synthetic(n_samples=n_samples)
    counts.round().astype(int).to_csv(tmp_path / "counts.csv")
    genes.to_csv(tmp_path / "genes.csv")
    arms.to_csv(tmp_path / "arms.csv", index=False)
    return ["--genes", str(tmp_path / "genes.csv"), "--arms", str(tmp_path / "arms.csv")]


def test_cli_writes_outputs_and_flags_trisomy(tmp_path):
    ref = _write(tmp_path)
    out = tmp_path / "out"
    main([str(tmp_path / "counts.csv"), "-o", str(out)] + ref)
    calls = pd.read_csv(out / "arm_calls.csv")
    flagged = calls[calls["tier"] == "high"]
    assert set(zip(flagged["sample"], flagged["arm"])) == {("S5", "1p"), ("S5", "1q")}
    assert (out / "arm_log2_ratios.csv").exists() and (out / "qc.csv").exists()
    assert pd.read_csv(out / "qc.csv")["genes_used"].iloc[0] == 800


def test_cli_gives_clean_error_for_two_samples(tmp_path):
    ref = _write(tmp_path, n_samples=2)
    with pytest.raises(SystemExit) as e:
        main([str(tmp_path / "counts.csv"), "-o", str(tmp_path / "out")] + ref)
    assert e.value.code == 2


def test_cli_rejects_non_numeric_table(tmp_path):
    (tmp_path / "bad.csv").write_text("gene,a,b,c\nENSG1.1,x,y,z\nENSG2.1,x,y,z\n")
    with pytest.raises(SystemExit) as e:
        main([str(tmp_path / "bad.csv"), "-o", str(tmp_path / "out")])
    assert e.value.code == 2
