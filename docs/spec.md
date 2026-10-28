# hpsc-integrity: spec (draft v0)

## Problem
Cultured human pluripotent stem cells (hPSCs) acquire recurrent chromosomal
abnormalities (e.g. gains of 1q, 12p, 17q, 20q11.21; loss of 10p, 18q).
Standard checks (G-banding, SNP array, qPCR) are slow, costly, or targeted.
Many labs already have bulk RNA-seq from the same cultures.

## Goal
Detect recurrent hPSC abnormalities from bulk RNA-seq and report them
with a confidence level, and measure how low a mosaic fraction can be detected.

## Inputs
- v0.1: gene-level expression table (genes x samples) + gene coordinates
- v0.2: aligned reads (BAM) for allelic-imbalance evidence

## Outputs
- calls.tsv: per sample, per chromosome arm: gain / loss / neutral, score, confidence
- report.html: per-sample chromosome plot + call table + warnings

## Method
1. Normalize expression; filter low-expressed genes.
2. Compare each sample to a reference panel of diploid hPSC lines.
3. Smooth log-ratios along genomic position; summarize per arm and per hotspot.
4. (v0.2) Compute allelic imbalance at common heterozygous SNPs; combine
   with expression evidence.
5. Annotate overlap with curated hotspot table (data/hotspots.tsv).

## Validation
- Real data: public hPSC RNA-seq with known karyotype/array results
  (datasets listed in docs/datasets.md).
- Simulation: synthetic gains/losses at known arms, mosaic fractions
  (10/25/50/100%), several depths.
- Metrics: sensitivity, specificity, detection limit vs. mosaic fraction.
- Comparison: eSNP-Karyotyping, if it runs.

## Non-goals (v0.1)
Single-cell data, X-inactivation, focal CNVs, SNV/TP53 calling.

## Known limitations (expected)
Low-level mosaicism, low depth, regions with few expressed genes,
copy-neutral LOH, cell-type composition effects in differentiated samples.

## Open decisions (I decide, and record why in docs/design.md)
- Normalization method
- Reference panel: which diploid lines, how many
- Smoothing window size
- Call thresholds
- How to merge expression and allelic evidence
