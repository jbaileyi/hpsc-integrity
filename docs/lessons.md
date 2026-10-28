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
