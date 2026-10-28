"""Spike-in simulation on the real GSE81402 lines. Run from the repository root:
    python scripts/run_simulation.py [reps]
Writes results/simulation/*.csv and figures/spikein_sensitivity.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from hpsc_integrity.annotation import load_reference
from hpsc_integrity.simulate import detection_limits, false_call_rates, run_grid, summarize

ROOT = Path(__file__).resolve().parents[1]
reps = int(sys.argv[1]) if len(sys.argv) > 1 else 20

counts = pd.read_csv(ROOT / "data/processed/gene_counts.csv", index_col=0)
genes, arms = load_reference()

TEST = ["CSES9", "CSES7", "HUES14"]                      # G-banding diploid or near-diploid
TARGETS = [("1q", "gain"), ("12p", "gain"), ("17q", "gain"), ("10p", "loss"), ("18q", "loss")]
KNOWN_ABNORMAL = {"CSES22": {"12"}, "CSES21": {"21"}}    # from G-banding; dropped as comparators there
FRACTIONS = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0]
DEPTHS = [1.0, 0.5, 0.25, 0.1]

runs = run_grid(counts, genes, arms, TEST, TARGETS, FRACTIONS, DEPTHS, reps,
                seed=1, known_abnormal=KNOWN_ABNORMAL)
summary = summarize(runs)
lod = detection_limits(summary)
fc = false_call_rates(runs)

out = ROOT / "results/simulation"
out.mkdir(parents=True, exist_ok=True)
runs.round(4).to_csv(out / "spikein_runs.csv", index=False)
summary.round(4).to_csv(out / "spikein_summary.csv", index=False)
lod.to_csv(out / "spikein_lod.csv", index=False)
fc.round(4).to_csv(out / "false_calls_by_depth.csv", index=False)

fig, axes = plt.subplots(1, len(TARGETS), figsize=(15, 3.2), sharey=True)
for ax, (arm, direction) in zip(axes, TARGETS):
    s = summary[(summary["target"] == arm) & (summary["direction"] == direction)]
    for d, g in s.groupby("depth"):
        ax.plot(g["f"], g["rate_any"], marker="o", ms=3, label=f"{d:g} of reads")
    ax.axhline(0.8, color="grey", lw=0.5, ls="--")
    ax.set_title(f"{arm} {direction}")
    ax.set_xlabel("mosaic fraction")
axes[0].set_ylabel("detected (possible or high)")
axes[0].legend(fontsize=7)
plt.tight_layout()
(ROOT / "figures").mkdir(exist_ok=True)
plt.savefig(ROOT / "figures/spikein_sensitivity.png", dpi=150)

pd.set_option("display.width", 200)
print(f"{len(runs)} runs, {reps} replicates per cell, test lines {TEST}")
print("\nFalse calls on non-target arms, by depth:")
print(fc.round(3).to_string(index=False))
print("\nSmallest mosaic fraction detected in >= 80% of runs (NaN = never within the grid):")
print(lod.pivot_table(index=["target", "direction"], columns="depth",
                      values=["lod_any", "lod_high"]).round(2).to_string())
