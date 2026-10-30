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

## Allelic results, all five lines, against pre-registered predictions
- CSES21 21q +0.418 (3.5 robust SDs; line SD 0.121): prediction met, marginal power.
- CSES9 20q -0.012, 20p -0.061: inside +/-0.10, so the rule gives "expression-only difference". Detection limit 38% (z=3).
- CSES7 17p -0.169, 17q +0.056: INCONCLUSIVE. Detection limit 58% exceeds the predicted 26-44% clone; line is underpowered.
- 9q: prediction (+/-0.15) WRONG. CSES9 +0.702 (65 SNPs), HUES14 +0.815 (85 SNPs), not seen in other 3 lines. Under investigation.
- Robust SD (MAD, outliers removed): CSES22 0.047, HUES14 0.065, CSES9 0.083, CSES21 0.121, CSES7 0.123.
- Earlier QC check of CSES7 alignment (50.9%) was invalid (sorted BAM head = chr1 only); redone with idxstats.

## 9q diagnosis and bootstrap results
- Hypothesis A (het depletion) MET for CSES9 and HUES14: 9q SNP share 0.25 / 0.27 of typical; chr9 ~70-110 Mb and ~135 Mb nearly
  homozygous in both lines (het/total 0-1 per 5 Mb window) while the other three lines are normal; 9q expression log2 +0.026 / -0.029.
- Interpretation: long homozygous / possible cnLOH on 9q, NOT a copy gain. Unconfirmed; allelic data alone cannot separate cnLOH from
  other causes of homozygosity.
- Lesson: median allelic ratio is invalid on het-depleted arms; the tool must flag het density and refuse to call such arms.
- Other het-depleted arms: CSES22 16p (0.29), CSES21 21p (0.31). Not assessed in bootstrap (<50 SNPs after filter).
- CSES7 unmapped reads: 49.0% of a 500k sample (matches BAM), diverse, not adapter/polyA. Likely non-human or low quality; not resolved.
  CSES9 unmapped reads are mostly Illumina adapter (trimming). MT is 9.2% of mapped CSES7 reads.
- Bootstrap (locus-level, 100 kb): CSES22 12q +0.529 [0.46, 0.60]; CSES21 21q +0.451 [0.26, 0.67];
  CSES9 20q +0.032 [-0.15, 0.27]; CSES7 17p -0.164 [-0.27, -0.01], 17q +0.064 [-0.08, 0.23].
- CSES9 20q: rule gives expression-only, but gain fraction up to ~40% not excluded. CSES7 chr17: no gain on both arms; clone <~35% on 17q not excluded.
- New: CSES21 19p -0.22 [-0.28, -0.11], unexplained.

## Corrections (after het-density run)
- CORRECTION: CSES22 16p has 58 SNPs and was assessed (log2 -0.090); the earlier line saying it was not assessed was wrong.
- CSES21 21p and HUES14 Yq have 1 SNP each: "too few SNPs", not het depletion. Refuse logic needs a separate TOO_FEW_SNPS state.
- GC content (0.495) of CSES7 reads is uninformative for human vs mouse. Unmapped reads may be rRNA, mouse, or low quality; unresolved.
- Observation: 19p negative in CSES7, CSES21, HUES14, CSES22; suspect a systematic artifact (see recurrence table).

## Recurrence across lines
- 19p negative in 5/5 lines (median -0.12), 21q 4/4, 8p 5/5, 7p 5/5, 15q positive 5/5. Modest (0.05-0.12) but same-sign across
  unrelated lines: likely systematic (reference bias / batch), n=5 only so cohort effects and batch are confounded.
- CSES21 19p (-0.22) is therefore not a line-specific loss. CSES21 21q (+0.418) sits on an arm with negative cohort median.
- Next: cohort-baseline rescoring (leave-one-out panel median).

## Cohort baseline test (negative result)
- Leave-one-out cross-line median subtraction, n=5. Robust SD raw -> corrected: CSES9 0.083->0.085, CSES7 0.123->0.135,
  CSES22 0.048->0.091 (worse), CSES21 0.121->0.101, HUES14 0.065->0.053. Mixed; hurts the cleanest line.
- CSES22 12q corrected +0.454 (raw +0.555, expected +0.51): correction worsened the best control.
- Decision: single-sample baseline is default; panel normalisation optional and off by default. n=5 is too small for a fair panel test.
- 19p allelic negative shift not mirrored in expression (CSES9 +0.156, CSES7 -0.102): likely technical, unproven.
- CSES7 unmapped-read BLAST: still pending.

## Allelic calls module: results against predictions
- Met: 12p/12q/21q GAIN; 9q (CSES9, HUES14) and CSES22 16p REFUSED; CSES7 17p NORMAL, 17q UNDERPOWERED; CSES9 20q UNDERPOWERED; 0 GAIN on ~190 null arm-tests.
- Partly failed: 19q refused in CSES7 (rel density 0.36) and CSES21 (0.46); 0.57-0.66 elsewhere. 19q is low in all five lines, so the single-sample
  density metric (SNPs per protein-coding gene) is biased for this arm. CSES21 19q refusal was not predicted.
- In-sample caveat: the CI rule was fixed after seeing the bootstrap CIs on these lines; 3 true positive arms only.
- 48-68% of arms are UNDERPOWERED: single-arm calls are limited to large clones. Module is arm-level; cannot see focal 20q11.21.
- gain_frac_hi95 uses the expression model (1 + f/2); pilot fit was 1.404 + 0.6 f. Unvalidated; to be replaced by simulation.
- Not duplicated by run_panel_experiment.py (expression); consistent with its conclusion that larger panels redistribute deviations.

## Density check and allelic spike-in (3-rep pass) against predictions
- Density: 19q single-sample rel density 0.36-0.66 in all five lines; arm-specific leave-one-out gives CSES9 1.15, HUES14 1.08, CSES21 0.79,
  CSES7 0.62. 9q x2 and CSES22 16p stay 0.23-0.30. Metric artifact for 19q confirmed in four lines; CSES7 stays 0.62 (not < 0.5). Fix needs a
  per-arm expected-density table (bundled), not five-line leave-one-out.
- Spike-in (reps=3, rough): f=0 and f<=0.3 never GAIN; false GAIN elsewhere 0. CSES9 and HUES14 reach GAIN at f=0.75 (f=0.5 gives 0-67%).
  PREDICTION FAILED: CSES7 never reaches GAIN even at f=1 (mean log2 +0.03 to +0.08 vs +0.38 to +0.50 elsewhere). Cause unchecked: noise or a bug.
- Caveat: allelic shift model (f/(2+f) read moves) validated against one real trisomy only.

## Allelic spike-in, 8 reps, and CSES7 diagnosis
- 8-rep results confirm the 3-rep pattern. CSES9 and HUES14: GAIN rate 1.0 at f>=0.75 (CSES9 20p 0.875), 0-0.5 at f=0.5, 0 at f<=0.3.
  CSES7: 0 at every f, mean log2 +0.01 to +0.06 at f=1 (others +0.40 to +0.46). 0 false GAINs elsewhere in any run, f=0 included.
- Limit is a grid point (0.75), true value is between 0.5 and 0.75; n=8 per cell, 2 usable lines, 3 arms.
- CSES7 is a real limit: depth and SNP counts match the other lines, but baseline median ratio is 1.7-1.9 (others 1.35-1.6) with wider IQR.
  A high baseline compresses the log2 effect of a gain and widens the bootstrap CI. Link to the 51% alignment / MT-rich library is a hypothesis, untested.
- Correction to my earlier note: I said a possible simulation bug could not be excluded; the identical code gives normal results in two lines
  and the CSES7 ratio data explain the result, so I now treat it as a limit.
- Caveat: shift model (moving f/(2+f) of reads) checked against one real trisomy only.

## Shift-model check against the real trisomy
- Real CSES22 12p/12q log2 +0.43/+0.53; simulated f=1 (CSES9, HUES14) +0.42/+0.43. Simulation is ~0.05 lower (conservative). Compared
  two arms of one real line with a mean over three arms in two lines; one trisomy only.
- CSES7 unmapped-read BLAST not run; CSES7 low alignment (50.9%) cause remains unresolved.

## CLI allelic step: QC flag against prediction
- Baseline median ratios: CSES9 1.434, CSES7 1.816, CSES22 1.404, CSES21 1.600, HUES14 1.402. Flag (> 1.6) hit only CSES7: prediction met,
  but CSES21 sits exactly on the cutoff (marginal). Threshold rests on one flagged line; CSES21 still calls its real 21q GAIN, so 1.6 may be conservative.
- CLI shows CSES7 and CSES21 19q as "refused (het-depleted)": known density-metric artifact (see density check), misleading to users until fixed.
- gain_frac_hi95 still in allelic_calls.csv; unvalidated.
- Allelic vs expression: expression flags CSES7 17p/17q, CSES9 20q as possible; allelic gives NORMAL/UNDERPOWERED. Complementary, not agreeing.

## Measured expression-tier limits (replaces back-of-envelope 20-33%)
- spikein_lod: lod_any 0.2-0.3, lod_high 0.5 for 1q/12p/10p/18q; 17q 0.5/0.75. Depth 0.1-1x changes little. Grid-resolution upper bounds.
- Null false calls per run: high 0.18-0.22, possible 2.3-3.3 (rises as depth falls). Prior estimate of 33% for 'high' was optimistic.
- datasets.md truth table was read by eye from Fig 3.

## CSES7 unmapped reads: closed without resolution
- Unmapped fraction 49.0% in a 500k sample (matches the BAM). Not adapter/polyA dominated; may be rRNA, mouse feeder, or low quality.
- 40-read FASTA for BLAST written (results/qc/SRR3501466_unmapped_sample.fa) but never classified. Decision: not pursued; documented as a limitation.

## Clone-fraction estimates, first validation (reps=1, n=3 per cell)
- P1 coverage 0.889 (n=36; too few runs to confirm calibration). P2 met (CSES9/HUES14 MAE 0.04-0.13 at f>=0.5). P3 met barely (width 0.40-0.49).
- P4 FAILED: CSES7 f=0.3 gave width 0.00, coverage 0.33 (confidently wrong, not uninformative).
- P5 met (12p 0.89, 12q 1.0, 21q 1.0). P6 met by count (6 <= 8) but 3 of 6 are CSES7 3q/14q/15q with f_lo = 1.0 on arms called NORMAL/UNDERPOWERED:
  likely a flat curve mapping any positive log2 to 1.0. Count alone hides this.
- CSES21 curve peaks at +0.27 at f=1 while real 21q is +0.45: simulated shift under-predicts real; CI ignores curve uncertainty, so f_lo=0.987 is overconfident.

## Clone-fraction fix, reps=3
- Q1 FAILED: CSES7 resolved in 4/36 spike-in runs, 2 of those missed the truth (one width 0.00). Q2 met (coverage 0.947, 72/76; CSES9+HUES14 70/72).
  Q3 mostly met (widths 0.27-0.56; f=1 slightly narrower than predicted). Q4 met (12p f_lo 0.63, 12q 0.90, 21q 0.75). Q5 met (1 null arm: CSES21 2q f_lo 0.19).
- Caveats: spike-in uses the same shift model as the calibration (partly self-testing); wide CIs cover easily; replicates share lines and arms; all real positives are full gains,
  so partial-fraction accuracy has no model-independent validation. CSES21 21q saturates at 1.0 (real +0.45 exceeds simulated f=1 +0.27).

## Clone-fraction gate diagnostic (10 seeds, 6 donors)
- Median-curve max: CSES9 0.36-0.49, HUES14 0.41-0.53, CSES22 0.41-0.50, CSES21 0.22-0.34 (all 10/10 >= 0.15); CSES7 0.001-0.158, 1/10 >= 0.15.
- Predictions met (CSES7 1/10). Spread overlaps threshold: the 4/36 CSES7 spike-in passes are gate noise from random donor choice.
- p20-curve: CSES7 max 0.049, CSES21 min 0.070 (separation too thin to rely on). Curve size tracks baseline ratio (1.82 -> 0.07, 1.60 -> 0.30, ~1.4 -> 0.45): 5 points, untested.

## Gate fix: 15 calibration donors (10 seeds)
- Median-curve max: CSES7 0.000-0.079 (0/10 >= 0.15; was 1/10, max 0.158), CSES21 0.248-0.355 (10/10), CSES9 0.412-0.486 (10/10). Predictions met.
- Only three lines re-run (HUES14, CSES22 not). 0.15 threshold remains in-sample (five lines).

## Clone-fraction validation, n_donors=15, reps=3
- R1 met (CSES7 resolved 0/36). R2 met (coverage 0.972, 70/72). R3 marginal fail (widths 0.32-0.67; f=0.5 cells grew to 0.67/0.58). R4 met (12p f_lo 0.65, 12q 0.95, 21q 0.78, f_hi 1).
- R5 half failed: 1 null arm but CSES21 2q f_lo = 0.36 (> 0.3). CSES21 curve flat at low f (0.03 at f=0.5): small log2 maps to a large fraction; arm-level offsets not in the CI. Untested hypothesis.
- Resolution: arms not called GAIN report only f_hi (what is excluded), not f_hat/f_lo.
- Intervals at f=0.5 are ~0.6 wide: coarse. Useful as lower bound for large gains and upper bound for exclusion.

## Clone fraction wired into CLI (--clone-fraction); display rule vs predictions
- Met: GAIN arms f_lo 0.67/0.95/0.77; 0 non-GAIN rows with f_lo; CSES7 UNRESOLVABLE on 37/37 arms. 45 tests pass.
- Partly failed: clean-line non-GAIN f_hi medians 0.50 (HUES14), 0.56 (CSES22), 0.68 (CSES9), not <= 0.5; CSES21 0.99.
- Same trisomies gave f_lo 0.65/0.95/0.78 in validation vs 0.67/0.95/0.77 in CLI: ~±0.02 calibration noise from random donor choice.
- Reading: arm-level allelic evidence excludes clones only above ~50-70% on clean lines, ~100% on CSES21; consistent with the ~75% spike-in limit.

## Head-to-head, eSNP-Karyotyping statistic with numeric rule (primary result, window 151)
- Key: S1 CSES22, S2 CSES21, S3 HUES14, S4 CSES9, S5 CSES7. Blinding was weak (baselines 1.40/1.60/1.41/1.44/1.86 identify the lines).
- B1 half met: chr12 raised in CSES22; 21q NOT raised in CSES21. B2 failed: 9q raised in HUES14 only, not CSES9. B4 failed: 0 other arms raised in any line.
- Cause of the misses: window 151 exceeds SNP count of the small arms (21q 141, CSES9 9q 65, HUES14 9q 85), so windows never sit inside the arm. A limit of our rule, not of eSNP-Karyotyping itself.
- Read: a ratio-only rule made no false calls elsewhere (matches our 0 false GAIN) and called HUES14 9q raised (het-depleted region read as a gain). Not shown: how eSNP's own LOH tools or CaSpER handle 9q (not run).
- This tests our rule on eSNP's statistic with identical allelic input; it does not test the published tool's own interpretation.

## Head-to-head exploratory window sweep (window 51; post-hoc, not primary)
- Predictions written after seeing the primary result and its explanation, so the first two were close to guaranteed. Window 151 rerun reproduces the primary result.
- Met: 21q raised in CSES21; 9q raised in both HUES14 and CSES9. Failed: no other arm raised in any line (0 vs predicted >= 1).
- Window 51 raised 12p, 12q (CSES22), 21q (CSES21), 9q (HUES14, CSES9), nothing in CSES7. Differences from hpsc-integrity: it refuses the two 9q arms; trisomy detection and specificity were equivalent on this input.
- Not shown: eSNP-Karyotyping's own interpretation or LOH tools, CaSpER.
