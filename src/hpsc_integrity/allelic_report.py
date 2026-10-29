"""Allelic arm calls for several samples, with a per-sample QC row (used by the CLI --snps option)."""
import pandas as pd

from .allelic import arm_allelic_table, filter_het_snps, read_allele_counts
from .allelic_calls import GAIN, NORMAL, REFUSED, UNDERPOWERED, call_arms

BASELINE_RATIO_WARN = 1.6   # provisional, set from one line (CSES7); a high baseline compresses the effect of a gain


def parse_specs(specs) -> dict:
    out = {}
    for s in specs:
        name, sep, path = s.partition("=")
        if not sep or not name or not path:
            raise ValueError(f"--snps expects SAMPLE=PATH, got {s!r}")
        if name in out:
            raise ValueError(f"sample {name!r} given twice in --snps")
        out[name] = path
    return out


def run_allelic(specs, genes, arms, baseline_warn: float = BASELINE_RATIO_WARN):
    """Returns (calls, qc): calls has one row per sample and arm; qc one row per sample."""
    calls, qc = [], []
    for name, path in parse_specs(specs).items():
        snps = filter_het_snps(read_allele_counts(path))
        res = call_arms(snps, genes, arms)
        baseline = float(arm_allelic_table(snps, arms)[1])
        assessed = res["call"].isin([GAIN, NORMAL, UNDERPOWERED])
        calls.append(res.rename_axis("arm").reset_index().assign(sample=name))
        qc.append(dict(sample=name, het_snps_kept=len(snps), baseline_ratio=round(baseline, 3),
                       arms_assessed=int(assessed.sum()),
                       frac_underpowered=round(float((res["call"] == UNDERPOWERED).sum() / max(assessed.sum(), 1)), 2),
                       flags="ALLELIC_POWER_LOW" if baseline > baseline_warn else ""))
    cols = ["sample", "arm"] + [c for c in calls[0].columns if c not in ("sample", "arm")]
    return pd.concat(calls, ignore_index=True)[cols], pd.DataFrame(qc)


def summarize(calls: pd.DataFrame, qc: pd.DataFrame) -> str:
    lines = ["Allelic evidence (arm-level; detects large clones only, see README):"]
    for r in qc.itertuples():
        c = calls[calls["sample"] == r.sample]
        gains = list(c.loc[c["call"] == GAIN, "arm"])
        refused = list(c.loc[c["call"] == REFUSED, "arm"])
        lines.append(f"  {r.sample}: GAIN {gains or 'none'}; refused (low het density) {refused or 'none'}; "
                     f"{r.frac_underpowered:.0%} of arms underpowered" + (f"  [{r.flags}]" if r.flags else ""))
    return "\n".join(lines)
