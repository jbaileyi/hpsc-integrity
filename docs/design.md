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

## Arm-level null and provisional call tiers
- Arms with >= 50 protein-coding genes (MIN_ARM_GENES = 50, provisional); X and Y excluded in v0.1.
- Null (39 arms x 5 lines minus 3 trisomic entries = 192): SD 0.060, 95th pct |x| 0.122, max 0.199.
- Provisional tiers: |log2| >= 0.20 high; 0.12-0.20 possible; < 0.12 no call. To be recalibrated
  on simulation, not on the five real lines.
- Back-of-envelope detectable mosaic fraction (assumes 91% attenuation holds for partial mosaics):
  ~20% at 0.12, ~33% at 0.20.
- 20q11.21 band (21 genes incl. BCL2L1) is not callable by expression; allelic method needed.
- TODO (me): confirm or change thresholds and write the reasoning in my own words.

## Package structure and decisions (v0.0.1)
- src layout; core.py holds pure functions (counts -> log2 ratios -> arm medians -> tiers).
- Call every chromosome arm, then annotate hotspots (trisomy 21 in the test data is not a hotspot).
- Tiers (provisional): high |log2| >= 0.20; possible >= 0.12; none. Arms need >= 50 genes. X/Y not assessed.
- Leave-one-out baseline needs >= 3 samples (code raises otherwise).
- OPEN: users with 1-2 samples need a bundled reference panel of diploid hPSC lines.
- TODO (me): re-read this in my own words; confirm or change each choice.

## Package v0.0.2
- annotation.py: Ensembl header parsing, gene table, UCSC cytoBand -> arms, gene-to-arm by midpoint.
- pipeline.py: score_arms (genes with >= 20 reads in every sample, log2 CPM, leave-one-out baseline,
  per-sample offset, arm medians). match_genes accepts Ensembl IDs with or without version suffix.
- cli.py: hpsc-integrity counts.csv -o outdir -> arm_log2_ratios.csv, arm_calls.csv, qc.csv.
- Reference bundled in the package (Ensembl 112 GRCh38 genes; hg38 arms), rebuilt by scripts/build_reference.py.
- Regression test: package reproduces results/arm_log2_ratios_leave_one_out.csv from the notebook.
- TODO before release: confirm Ensembl and UCSC redistribution terms for the bundled tables.

## Spike-in simulation: design and predictions (written BEFORE running)
- Real G-banding-diploid lines (CSES9, CSES7, HUES14) modified in silico: a gain in fraction f of
  cells raises the target/rest ratio by 1 + f/2 (non-target reads thinned); a loss lowers it by
  1 - f/2 (target reads thinned). Reads are only removed, never invented.
- Targets: 1q, 12p, 17q gains; 10p, 18q losses. f in {0, .05, .1, .2, .3, .5, .75, 1}; depth
  (fraction of reads kept for the test sample) in {1, .5, .25, .1}; 20 replicates.
- Comparators: the other four lines, minus any known abnormal on the target chromosome
  (CSES22 for chr12, CSES21 for chr21).
- Assumes expression is proportional to copy number (real trisomies reached 91-94% of that).
- Prediction from the real-data null (back-of-envelope): 'possible' (0.12) needs about 20% mosaic and
  'high' (0.20) about 33% at full depth. Depth should worsen both, because fewer genes pass the
  20-read filter; size not predicted.
- Limitations: thresholds derived from the same five lines (in-sample); panel is 4 comparators from
  one study; CSES7 and HUES14 carry unexplained chr17 shifts that affect their 17q results.

## Panel experiment: design and predictions (written BEFORE running)
- Panels: A = the 4 other hESC lines (as before); B = all 12 other lines; C = diploid hESC lines + the
  8 HiSeq lines (10-11 comparators); D = the 8 HiSeq lines only (pure platform mismatch: NextSeq 500,
  ~85 bp vs HiSeq 2000, 101 bp; parthenogenetic/SWAP lines).
- Test lines: real data = all five hESC lines; spike-in = CSES9, CSES7, HUES14 at full depth, 10 replicates.
- Predictions:
  1. C and B give fewer false weak calls than A, especially for CSES7 (line-specific deviations are
     diluted in a larger median); trisomy values stay near 0.5.
  2. D is worse than C if the platform shifts arm-level ratios. D matching C would mean platform mismatch
     is tolerable and a public reference panel is feasible.
  3. Detection limits stay near 0.3 (weak) and 0.5 (strong), set by the thresholds and the remaining
     noise; a drop to 0.2/0.3 would need baseline noise to fall.
- Decision rule (a judgment call): if C cuts false weak calls by roughly a third or more without hurting
  trisomy detection, a reference panel becomes the recommended design; if D fails clearly, the panel must
  be platform-matched.
- Caveat: the HiSeq lines have no G-banding truth here; if some carry abnormalities, a median over many
  lines is robust to a few, but I have not verified this.

## Allelic module v0.2: choices and pilot prediction (written BEFORE running)
- Alignment: HISAT2 2.2.3, prebuilt GRCh38 genome index (Ensembl chromosome names), single-end.
  CSES22: 86.3% overall alignment (81.0% unique), ~8 min wall time with 8 threads.
- SNPs: bcftools mpileup (MAPQ >= 20, base quality >= 20, max depth 2000) -> call -mv -> biallelic SNPs,
  QUAL >= 30 (arbitrary starting value). No known-SNP filter, so RNA editing and paralog artifacts can
  get through.
- Filters as in the source paper: depth >= 20, minor allele fraction >= 0.2. Statistic: per-arm median
  of the major/minor allele ratio; log2 relative to the sample's own baseline (median of autosomal arm
  medians); arms need >= 50 SNPs; X and Y not assessed.
- Expectation for a clonal gain: ratio ~ 1 + f (about 2 at f = 1) against a diploid baseline that sits
  above 1 because of sampling noise (~1.3-1.4 in the paper).
- PREDICTION, CSES22 pilot: 12p and 12q log2_vs_baseline between +0.4 and +0.75; most other arms
  within +/-0.15; kept SNPs in the low thousands (the paper reports ~2,000 SNPs at 15-20M reads), so
  only the larger arms reach 50 SNPs.
