"""Does a larger or different reference panel help? Run from the repository root:
    python scripts/run_panel_experiment.py [reps]
Needs data/processed/gene_counts_all13.csv (scripts/build_counts.py).
Writes results/panel/*.csv.
"""
import sys
from pathlib import Path

import pandas as pd

from hpsc_integrity.annotation import load_reference
from hpsc_integrity.core import call_tiers
from hpsc_integrity.pipeline import score_arms
from hpsc_integrity.simulate import detection_limits, run_grid, summarize

ROOT = Path(__file__).resolve().parents[1]
reps = int(sys.argv[1]) if len(sys.argv) > 1 else 10
counts = pd.read_csv(ROOT / "data/processed/gene_counts_all13.csv", index_col=0)
genes, arms = load_reference()

HESC = ["CSES9", "CSES7", "CSES22", "CSES21", "HUES14"]       # NextSeq 500, G-banding truth
DIPLOID = ["CSES9", "CSES7", "HUES14"]                         # no whole-chromosome abnormality by G-banding
HISEQ = [c for c in counts.columns if c not in HESC]           # HiSeq 2000, parthenogenetic / SWAP lines
KNOWN_ABNORMAL = {"CSES22": {"12"}, "CSES21": {"21"}}
TRUTH_CHROM = {"CSES22": "12", "CSES21": "21"}                 # chromosome gained, from G-banding
PANELS = ["A_hESC_4", "B_all_12", "C_diploid", "D_HiSeq_8"]


def comparators(panel, test):
    if panel == "A_hESC_4":
        return [c for c in HESC if c != test]                  # the panel used so far
    if panel == "B_all_12":
        return [c for c in counts.columns if c != test]        # everything else, trisomic lines included
    if panel == "C_diploid":
        return [c for c in DIPLOID + HISEQ if c != test]       # diploid hESC lines + HiSeq lines
    return list(HISEQ)                                         # D: other platform only


# ---- 1. real data: each G-banding-labeled line against each panel ----
rows = []
for panel in PANELS:
    for test in HESC:
        comps = comparators(panel, test)
        table, _, n = score_arms(counts[[test] + comps], genes, arms)
        calls = call_tiers(table)
        calls = calls[(calls["sample"] == test) & (calls["tier"] != "not_assessed")]
        on_truth = calls["arm"].str[:-1] == TRUTH_CHROM.get(test)      # all arms of the gained chromosome
        truth = calls.loc[on_truth, "arm"].tolist()
        rest = calls[~on_truth]
        rows.append(dict(
            panel=panel, line=test, n_comparators=len(comps), genes_used=n,
            truth_arms=", ".join(f"{a} {calls.loc[calls['arm'] == a, 'log2_ratio'].iloc[0]:+.2f}" for a in truth),
            false_high=int((rest["tier"] == "high").sum()),
            false_weak=int((rest["tier"] == "possible").sum()),
            max_abs_other=round(rest["log2_ratio"].abs().max(), 3)))
real = pd.DataFrame(rows)

# ---- 2. spike-in simulation per panel (full depth) ----
TARGETS = [("1q", "gain"), ("12p", "gain"), ("17q", "gain"), ("10p", "loss"), ("18q", "loss")]
FRACTIONS = [0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0]
parts = []
for panel in PANELS:
    for test in DIPLOID:
        sub = counts[[test] + comparators(panel, test)]
        r = run_grid(sub, genes, arms, [test], TARGETS, FRACTIONS, [1.0], reps,
                     seed=1, known_abnormal=KNOWN_ABNORMAL)
        r.insert(0, "panel", panel)
        parts.append(r)
runs = pd.concat(parts, ignore_index=True)
lods = []
for panel, g in runs.groupby("panel"):
    l = detection_limits(summarize(g))
    l.insert(0, "panel", panel)
    lods.append(l)
lod = pd.concat(lods, ignore_index=True)
fc = runs.groupby(["panel", "test_sample"])[["other_high", "other_possible"]].mean().round(2)

out = ROOT / "results/panel"
out.mkdir(parents=True, exist_ok=True)
real.to_csv(out / "panel_real_calls.csv", index=False)
runs.round(4).to_csv(out / "panel_spikein_runs.csv", index=False)
lod.to_csv(out / "panel_lod.csv", index=False)

pd.set_option("display.width", 220)
print("REAL DATA: false calls exclude the arms of the gained chromosome (chr12 for CSES22, chr21 for CSES21)")
print(real.to_string(index=False))
print("\nSPIKE-IN: false calls per run on non-target arms")
print(fc.unstack("panel").to_string())
print("\nSPIKE-IN: smallest mosaic fraction detected in >= 80% of runs (NaN = never)")
print(lod.pivot_table(index=["target", "direction"], columns="panel",
                      values=["lod_any", "lod_high"]).round(2).to_string())
