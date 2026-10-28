# Design decisions

## Reference (pinned)
- Ensembl release 112, GRCh38. FASTA: cdna.all; annotation: Homo_sapiens.GRCh38.112.gtf.gz.
- Reason: the GEO FPKM file used retired UCSC IDs (docs/lessons.md); re-quantifying from reads
  with a documented build makes coordinates unambiguous.
- Gene table will keep chromosomes 1-22, X, Y only (patches/alt scaffolds dropped).

## Single-end fragment length for kallisto
- Reads are single-end (~84 bp, TruSeq poly-A libraries per paper methods); fragment length
  is not measured, so kallisto needs an assumed mean and SD.
- Used -l 200 -s 20 (kallisto manual example). Assumption, not a measurement.
- TODO: sensitivity check (re-run one sample with a different -l, compare chromosome-arm results).

## Gene filter
- Keep protein-coding genes (Ensembl gene_biotype) on chr1-22, X, Y only.
- Dropped: MT (not nuclear copy number; few very abundant genes distort normalization),
  alt scaffolds/patches (duplicate genes), pseudogenes and other non-coding genes.
- Evidence: top transcript by TPM in CSES22 included a 177-nt processed pseudogene with
  effective length ~10 nt (inflated TPM; depends on the -l 200 assumption).
- Gene counts = sum of kallisto est_counts over a gene's transcripts.
- TODO (me): rewrite this reasoning in my own words.

## Index uses primary-assembly transcripts only
- See docs/lessons.md (POU5F1 / MHC alt haplotypes). MT transcripts stay in the index so
  mitochondrial reads are not forced onto nuclear lookalikes; they are dropped at the gene-table stage.

## First profile (expression), initial choices
- Gene filter: >= 20 estimated reads in every sample (MIN_COUNTS = 20). Starting value, not tuned.
- Normalization: counts per million, log2(CPM + 1).
- Baseline: per-gene median across the five hESC lines. Valid only because each chromosome is
  abnormal in at most one line; this is the batch-baseline dependence noted for e-Karyotyping.
- Smoothing: moving median, 200 genes, within chromosome (window size from the paper's
  expression method; not tuned).
- TODO (me): justify or change each choice after looking at results; rewrite in my own words.

## Baseline: leave-one-out instead of median-of-5
- Median-of-5 pins ~1/5 of values at exactly 0 (the line equal to the median), so flat regions show
  exact zeros and a moving median stays stuck until more than half a window deviates. It hides
  small shifts, which matters for a detection-limit claim.
- Leave-one-out (each line vs the median of the other four) removes this. Chosen for this structural
  reason, before any tuning; real-data truth is not used to pick settings.
- Result: trisomy 12 (CSES22) +0.55, trisomy 21 (CSES21) +0.53, vs 0.585 expected (94%, 91%).
- Known artifact: diploid lines dip slightly on a chromosome that is trisomic in another line.

## Per-sample offset and empirical null
- After leave-one-out, each line carries a global offset (median of 22 autosome medians):
  CSES9 -0.053, CSES7 +0.043, CSES22 +0.016, CSES21 0.000, HUES14 -0.040. Subtracted per line.
  Cause not established.
- Result: chr12 (CSES22) +0.53, chr21 (CSES21) +0.53 = 91% of log2(3/2) = 0.585.
- Empirical null (108 autosome medians excluding the two trisomies): SD 0.051, max |x| 0.177.
  Largest: CSES7 chr17 +0.18, HUES14 chr17 -0.16, CSES7 chr21 -0.13 (unexplained).
- Limitation: null comes from 5 lines of one study/platform; other labs' data will be noisier.
- Idea for simulation: two-pass baseline (call, then drop flagged samples from comparators).
  Evidence: all four diploid lines read -0.06 to -0.07 on chr12 because the trisomic line is
  in their comparator set.
- Back-of-envelope whole-chromosome detection limit from this null: ~25-30% mosaic. NOT a result.
