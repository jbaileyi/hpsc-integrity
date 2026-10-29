- Version-stripped match to hg38 new-style uc IDs: 47,682 / 104,133 (46%). Too low,
  and same-ID-same-transcript is unverified after the UCSC annotation overhaul.
  Decision: do not use the GEO FPKM file for coordinates; re-quantify from raw reads
  (SRP075000) with a current annotation instead.

## Alt-haplotype contigs absorbed reads from main-chromosome genes (POU5F1)
- Marker check on undifferentiated lines: POU5F1 ~1-5 CPM in 4 of 5 lines (LIN28A, SOX2 normal);
  CSES7 ~1,100.
- Diagnostic: est_counts for POU5F1 transcripts summed by chrom/scaffold. In CSES9, CSES22,
  CSES21 and HUES14 the reads sat on HSCHR6_MHC_* alt contigs (e.g. HUES14 ~32.6k on each of two
  contigs), not on chr6. In CSES7 they were on chr6 (~19.4k).
- Cause: cdna.all FASTA contains near-identical MHC haplotype copies; kallisto assigns reads among
  them; the gene table (main chromosomes only) discarded those counts.
- Fix: index transcripts on chr1-22, X, Y, MT only (187,056) and re-quantify.
- Implication: any gene with alt/patch copies is affected, not just POU5F1.
- Result after primary-only re-index: POU5F1 760-1,898 CPM in all five lines (was 1.4-4.6 in
  four of five). Mapping rates unchanged (reads still pseudoaligned; they had been assigned to
  alt copies that the gene table dropped).

## Median baseline bias reproduced in synthetic data
- Synthetic test (5 samples, one with chr1 at 1.5x): trisomic sample offset -0.172; chr1 arms +0.576
  and +0.583 (expected 0.585); the four diploid samples read -0.03 to -0.08 on chr1.
- Same direction and size as real data (diploid lines -0.06 to -0.07 on chr12 when CSES22 trisomic).
- Cause: with 4 comparators, one abnormal sample shifts the median. Estimator property, not a data quirk.
- Mitigations to test in simulation: larger reference panel; two-pass baseline excluding flagged samples.

## Spike-in simulation: my back-of-envelope prediction was optimistic
- Predicted detection limits (full depth): ~20% mosaic for any call, ~33% for 'high'.
- Observed (80% detection): 0.30 and 0.50 on the grid (steps 0.2/0.3/0.5); 17q worse (0.75).
- Cause (post hoc): the arithmetic gave the fraction where the EXPECTED ratio reaches the threshold
  (~50% detection), not 80%. Adding arm noise (SD ~0.06; 0.84 SD for 80% power) gives ~0.28 and
  ~0.42, matching the grid. Lesson: state the power level when predicting a limit.

## Per-line heterogeneity in the simulation (depth 1.0 unless noted)
- False 'high' calls per run on non-target arms: CSES7 0.58, CSES9 0.00, HUES14 0.00 (0.12 at 10% depth).
  Pooled 0.2 per run is misleading. Probable source: CSES7 17p (0.199 in real run, 0.20 cutoff);
  arm names not recorded in spikein_runs.csv, so unverified.
- False 'possible' calls per run: CSES7 3.78, CSES9 2.05, HUES14 1.01 (about 10%, 5%, 3% of arms);
  30-100% more at 10% depth.
- 17q: CSES7 'detected' at f=0 (already a weak call at baseline); HUES14 needs 0.75 (baseline -0.17);
  CSES9 needs 0.5 (baseline -0.07). Pooled 17q limit is not meaningful. Line- and arm-specific
  baseline deviations shift the effective detection limit.
- Depth 1.0 -> 0.1: genes used 13,165 -> 10,223; target arm genes 312 -> 244; detection limits unchanged.
- Hypothesis recorded before checking: per-sample noise (robust SD of arm ratios) orders
  CSES7 > CSES9 > HUES14.

## Noise hypothesis result: robust SD does not explain the false calls
- Robust SD of arm ratios (real run): CSES9 0.045, CSES7 0.051, CSES22 0.044, CSES21 0.066, HUES14 0.042.
  Predicted order CSES7 > CSES9 > HUES14 held, but the spread among the three is ~20%, while false
  weak calls differ ~3.7x (3.78 / 2.05 / 1.01 per run).
- Simulated false weak calls equal the arms already flagged in the real run (CSES7 4, CSES9 2,
  HUES14 1): fixed, line-specific arm deviations, not random noise. This is why depth hardly matters.
  A robust SD ignores such outlier arms by design, so noise-scaled thresholds would not fix this.
- Consequence: the false-call rate measured here is a property of these three lines and this 4-line
  panel; it does not transfer to new samples. Independent data are needed.

## Prediction before the comparator-subset test
- CSES7 21q (-0.13) comes from CSES21's trisomy 21 in its comparator set: it should shrink toward 0
  when CSES21 is dropped. CSES7 chr17 (+0.16/+0.20) is line-specific: it should persist in every subset
  (or move only slightly when HUES14, whose 17q is -0.17, is dropped).

## Result of the comparator-subset test (real data, 3-comparator subsets)
- CSES7 21q: -0.13 (all 4) -> -0.08 when CSES21 (trisomy 21) is dropped: ~0.05 of the deficit comes from the
  trisomic comparator; the rest (-0.08 to -0.20 across subsets) persists. Prediction partly supported.
- CSES7 chr17: 17p 0.17-0.21, 17q 0.13-0.23 in every subset: line-specific. Prediction supported.
- Opposite chr17 deviations amplify each other: dropping CSES7 moves CSES9 17p/17q -0.04/-0.07 -> 0.00/-0.02
  and HUES14 17q -0.17 -> -0.13; dropping HUES14 lowers CSES7 17q 0.16 -> 0.13.
- HUES14 17q stays -0.13 to -0.21 in every subset; CSES9 19p/20q stay +0.13/+0.18 and +0.14/+0.19.
- A single line-arm value moves ~0.05-0.10 with comparator choice (panel-composition uncertainty).
- Unexplained: whether CSES7 (up) and HUES14 (down) chr17 shifts are biology (e.g., low-level mosaicism that
  G-banding and the allelic method could miss) or other line-to-line expression differences.

## Panel experiment: results against predictions (pooled over CSES9, CSES7, HUES14; my arithmetic)
- False weak calls per run: A 2.30, B 1.09, C 1.32, D 1.95. False strong calls per run: A 0.20,
  B 0.93, C 0.93, D 0.93 (about 4.6x higher than A). Trisomy values: A 0.50/0.55/0.53, B 0.40/0.58/0.45,
  C 0.40/0.57/0.45, D 0.41/0.56/0.40 (12p/12q/21q).
- P1 weak calls fall: true, but strong false calls rose (not predicted). P2 D worse than C: modestly.
  P3 limits stay ~0.3/0.5: true (17q uninterpretable).
- My decision rule mentioned weak calls only: an oversight, not a pass. Conclusion drawn without it.
- B, C, D nearly identical (8 of their comparators are the same HiSeq lines): other-platform comparators work
  about as well as same-platform ones.
- Per-line: CSES9 total false calls/run 2.06 (A) -> 3.18/4.02/4.87 (B/C/D); CSES7 4.44 -> 2.89/2.72/2.80;
  HUES14 1.01 -> 0.00/0.03/0.99. CSES21 false weak 2 (A) vs 11-12 (B/C/D): cause unknown (platform/library
  or line biology; not checked).
- Arm values move ~0.1 with panel choice (chr12: 0.50/0.55 under A vs 0.40/0.58 under B).
- Conclusion: a larger panel redistributes line-specific deviations; it does not remove them. Expression-only
  strong calls at 0.20 are not reliable alone. Independent evidence (allelic imbalance) is the planned fix (v0.2).

## 13-line profile: results and corrections
- Each line vs all 12 others (arms of known trisomies excluded; ~504 arm-line pairs): none >= 0.30;
  22 arms >= 0.12 (13 in CSES21 alone); 5 arms >= 0.20, all in NextSeq hESC lines (5/192 = 2.6%, ~1 per
  sample), none in the 8 HiSeq lines. HiSeq lines resemble each other, so their clean profiles are not a
  general null. Platform and line type are confounded.
- Robust SD: HiSeq 0.022-0.046; hESC CSES9 0.061, CSES7 0.039, CSES22 0.055, CSES21 0.100, HUES14 0.029.
  CSES21 is the clear noise outlier: a per-sample QC flag is useful (different use from the dropped
  threshold-scaling idea).
- CORRECTION: HUES14 17q is -0.17 vs 4 hESC lines but -0.04 vs the 8 HiSeq lines: an artifact of CSES7 being
  in its comparator set. HUES14 is the cleanest line.
- CSES7 chr17 persists and grows vs HiSeq (17p/17q +0.24/+0.26); CSES9 20q +0.17 -> +0.27. Both are arms
  commonly gained in hPSC culture. Unadjudicated, not "false": G-banding may be from another passage and
  cannot exclude clones < ~14%; the allelic method (Fig 3, image too coarse to read per arm) cannot
  exclude < ~1/3.
- PREDICTION (before computing): if these are real subclonal gains, allelic imbalance from our own reads
  will rise on CSES7 chr17 and CSES9 chr20 relative to each line's other chromosomes. If it does not,
  I will treat them as expression differences.

## Allelic pilot (CSES22) against predictions
- 12p +0.474 (ratio 1.950, 151 SNPs), 12q +0.555 (ratio 2.062, 490 SNPs); baseline 1.404; expected for a
  full trisomy log2(2.0/1.404) = 0.510. Prediction (+0.4 to +0.75): met.
- Other 37 autosomal arms: SD 0.050 (robust 0.047), max |x| 0.104 (21q). Prediction (within +/-0.15): met.
  Trisomy arms are 9.5 and 11.1 null SDs above the null mean.
- SNP counts: 64,718 rows, 9,971 kept. Prediction (low thousands; only larger arms reach 50): WRONG.
  I misapplied the paper's ~2,000 SNP figure; every assessed autosomal arm had >= 58 SNPs.
- X arms (male line) ratios 1.58 / 1.82: not assessed by design; probably artifacts of hemizygous X (not checked).
- Runtime: bcftools mpileup + call took 28.4 min on one core for 45.7M reads.
- Crude linear model (assumption from one trisomy): ratio(f) ~ 1.404 + 0.6 f. 0.104 (largest null) <-> f ~ 0.18;
  0.15 (3 SD) <-> f ~ 0.26. To be replaced by simulation.
