"""Command-line interface: hpsc-integrity counts.csv -o outdir"""
import argparse
from pathlib import Path

import pandas as pd

from . import __version__
from .annotation import load_reference
from .core import call_tiers
from .allelic_report import run_allelic, summarize
from .pipeline import score_arms


def read_counts(path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python", index_col=0)
    if not df.index.is_unique:
        raise ValueError("gene IDs in the first column must be unique")
    try:
        return df.apply(pd.to_numeric, errors="raise")
    except (ValueError, TypeError):
        raise ValueError("the count table must be numeric apart from the gene ID column")


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="hpsc-integrity",
        description="Screen hPSC RNA-seq samples for chromosome-arm gains and losses.",
    )
    p.add_argument("counts", help="gene-by-sample count table (CSV/TSV); first column = Ensembl gene ID")
    p.add_argument("-o", "--out", default="hpsc_integrity_out", help="output directory")
    p.add_argument("--min-counts", type=int, default=20,
                   help="a gene needs this many estimated reads in every sample (default 20)")
    p.add_argument("--genes", help="gene table CSV (default: bundled Ensembl 112 / GRCh38)")
    p.add_argument("--arms", help="arm boundary CSV (default: bundled hg38)")
    p.add_argument("--snps", action="append", metavar="SAMPLE=PATH", default=[],
                   help="optional: per-sample heterozygous-SNP allele counts (bcftools query output); repeat per sample")
    p.add_argument("--clone-fraction", action="store_true",
                   help="with --snps: add clone-fraction bounds per arm (slow; simulates gains to calibrate each sample)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = p.parse_args(argv)

    try:
        counts = read_counts(args.counts)
        genes, arms = load_reference(args.genes, args.arms)
        table, offsets, n_genes = score_arms(counts, genes, arms, min_counts=args.min_counts)
    except (ValueError, FileNotFoundError) as e:
        p.error(str(e))

    calls = call_tiers(table)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    table.round(4).to_csv(out / "arm_log2_ratios.csv")
    calls.round({"log2_ratio": 4}).to_csv(out / "arm_calls.csv", index=False)
    qc = pd.DataFrame({"offset": offsets.round(4), "genes_used": n_genes})
    qc.index.name = "sample"
    qc.to_csv(out / "qc.csv")
    if args.snps:
        try:
            a_calls, a_qc = run_allelic(args.snps, genes, arms, clone_fraction=args.clone_fraction)
        except (ValueError, FileNotFoundError) as e:
            p.error(str(e))
        a_calls.to_csv(out / "allelic_calls.csv", index=False)
        a_qc.to_csv(out / "allelic_qc.csv", index=False)
        print(summarize(a_calls, a_qc))

    flagged = calls[calls["tier"].isin(["high", "possible"])].copy()
    flagged["abs"] = flagged["log2_ratio"].abs()
    flagged = flagged.sort_values(["tier", "abs"], ascending=[True, False]).drop(columns="abs")
    print(f"{counts.shape[1]} samples, {n_genes} genes used. Results in {out}/")
    if flagged.empty:
        print("No arm-level abnormalities above the provisional thresholds.")
    else:
        print(flagged.round(3).to_string(index=False))
