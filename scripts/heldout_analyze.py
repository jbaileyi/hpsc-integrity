"""Frozen-tool allelic calls on the held-out HipSci samples. python scripts/heldout_analyze.py
Writes results/heldout/{calls_long,sample_qc,failures}.csv"""
from pathlib import Path
import pandas as pd
from hpsc_integrity.allelic_report import run_allelic
from hpsc_integrity.annotation import load_reference

genes, arms = load_reference()
sel = pd.read_csv("results/heldout/selection.csv")
calls, qcs, failed = [], [], []
for r in sel.itertuples():
    p = Path(f"data/heldout/bam/{r.run_accession}.snps.tsv")
    if not p.exists():
        failed.append((r.sample_alias, "NO_SNP_FILE")); continue
    try:
        c, q = run_allelic([f"{r.sample_alias}={p}"], genes, arms)
    except ValueError as e:
        failed.append((r.sample_alias, str(e)[:100])); continue
    calls.append(c); qcs.append(q)
out = Path("results/heldout"); out.mkdir(parents=True, exist_ok=True)
pd.DataFrame(failed, columns=["sample", "reason"]).to_csv(out / "failures.csv", index=False)
if not calls:
    raise SystemExit(f"nothing analysed; failures: {failed[:5]}")
C, Q = pd.concat(calls, ignore_index=True), pd.concat(qcs, ignore_index=True)
C.to_csv(out / "calls_long.csv", index=False); Q.to_csv(out / "sample_qc.csv", index=False)
n = len(sel)
g = C[C.call == "GAIN"]
print(f"analysed {len(Q)} of {n} selected samples; not analysed: {len(failed)}")
print("samples with >= 1 GAIN:", g["sample"].nunique())
print(g[["sample", "arm", "log2", "ci_lo", "ci_hi", "n_snps"]].to_string(index=False) if len(g) else "no GAIN calls")
ref = C[C.call == "REFUSED_HET_DEPLETED"]
print("refused arms:", ref.arm.value_counts().to_dict(), "| samples with a refusal:", ref["sample"].nunique())
print("ALLELIC_POWER_LOW:", int((Q.flags != "").sum()), "of", len(Q))
print("kept het SNPs per sample: median", int(Q.het_snps_kept.median()), "min", int(Q.het_snps_kept.min()), "max", int(Q.het_snps_kept.max()))
print("arm-test mix:", C.call.value_counts(normalize=True).round(3).to_dict())
print("not analysed:", failed[:10])
