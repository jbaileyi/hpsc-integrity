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

## Allelic predictions for the other four lines (written BEFORE running)
- CSES21 21q: log2 +0.4 to +0.75 (positive control). If 21q has fewer than 50 SNPs it is not assessed.
- CSES9, CSES7, HUES14: arms other than those under test within +/-0.15 of baseline.
- CSES7 chr17: expression shift +0.16 to +0.26 would imply a clone of roughly 26-44% if real
  (expression ~ copy number, ~91% attenuation); crude model then predicts allelic log2 +0.15 to +0.25
  on BOTH 17p and 17q.
- CSES9 20q: expression +0.17 to +0.27 implies roughly 28-46%; predicts allelic log2 +0.16 to +0.26.
- Decision rule: supports a real gain if allelic log2 exceeds +0.15 (~3 null SDs of that line) on both arms
  of chr17 (CSES7) or on 20q (CSES9); supports an expression-only difference if within +/-0.10; between
  is inconclusive. Each line's null SD is taken from its own other arms.
- Caveat: the linear model is fitted to one trisomy and is an assumption.

## 9q hypotheses (written BEFORE het-depletion and expression checks)
- A (het-depleted arm; leftover hets are artifacts): CSES9 and HUES14 9q hold <40% of the typical share of kept SNPs;
  het/total near 0 across much of 9q; no 9q expression shift.
- B (real imbalance): normal het density on 9q; 9q expression log2 > +0.1 in both lines.
- Note: ratios 2.33 / 2.47 exceed the full-trisomy value in CSES22 (2.06), which a single-copy gain should not reach.

## Allelic calls module (written BEFORE running)
- States: GAIN / NORMAL / UNDERPOWERED / REFUSED_HET_DEPLETED / TOO_FEW_SNPS. No LOSS (a loss depletes hets, it does not shift ratios).
- GAIN if bootstrap CI lower bound > +0.10; NORMAL if CI upper bound < +0.10; otherwise UNDERPOWERED. Refuse if SNPs-per-gene is < 0.5x the sample median.
- Predictions: CSES22 12p and 12q = GAIN; CSES21 21q = GAIN; CSES9 and HUES14 9q = REFUSED; CSES22 16p = REFUSED;
  CSES7 17p = NORMAL, 17q = UNDERPOWERED; CSES9 20q = UNDERPOWERED. Zero GAIN calls on any other arm in any line.
- Expected weakness: many arms in the noisy lines (CSES7, CSES21) will be UNDERPOWERED, not NORMAL. That is honest, not a bug.
- Risk: HUES14 9q relative density may sit near the 0.5 cutoff; CSES7 8p/19p/19q (0.53-0.58 vs other lines) may be wrongly refused.

## Density fix and allelic spike-in (written BEFORE running)
- Density test: compare single-sample relative density with leave-one-out arm-specific density (arm density / median of other 4 lines).
  Predict: 19q moves to ~1 in all lines under leave-one-out; 9q (x2) and CSES22 16p stay < 0.5; no new arm falls below 0.5.
  If 19q stays < 0.5 in CSES7/CSES21 after correction, treat as possibly real het depletion, not a metric artifact.
- Spike-in (reads moved between alleles at existing het SNPs, depth unchanged; targets 1q, 12q, 20p in CSES9, HUES14, CSES7):
  Predict f=0 gives <= 1% GAIN on the target arm and ~0 false GAINs elsewhere; f=1 gives GAIN in >= 90% of runs;
  80%-GAIN detection limit is >= 0.5 in all three lines (arm-level allelic is a large-clone detector), and worst in CSES7;
  f <= 0.3 is essentially never GAIN.

## CLI allelic step: QC flag (written BEFORE running)
- `hpsc-integrity counts.csv --snps NAME=snps.tsv` adds allelic calls (GAIN / refusals) and an allelic QC row per sample.
- Flag ALLELIC_POWER_LOW if the sample's baseline median ratio (median of autosomal arm medians) exceeds 1.6. Provisional: set from
  CSES7 only (arm medians 1.66-1.88 vs 1.35-1.6 in CSES9/HUES14; CSES22 baseline 1.404).
- Prediction: CSES7 flagged; CSES9, CSES22, CSES21, HUES14 not. If step 1 shows CSES21 above 1.6, the threshold is not specific to CSES7.
- Licence TODO resolved on 2026-10-29: terms recorded in SOURCES.txt.
- 2026-10-29: licence TODO resolved. Provider terms read and recorded in SOURCES.txt (Ensembl: no restrictions; UCSC: free, cytoBand not in the restricted list, cite UCSC).

## Clone-fraction estimates (allelic): written BEFORE running
- Method: per-sample response curve (shift f simulated on that sample's own normal arms, median log2 change per f); the arm's log2 and bootstrap CI are inverted through it.
- Predictions, spike-in on CSES9/HUES14 (arms 1q, 12q, 20p; f = 0.3, 0.5, 0.75, 1.0):
  P1 the 95% CI covers the true f in >= 85% of runs overall (nominal 95%; calibration noise is not in the CI).
  P2 at f = 1, estimate >= 0.8 and upper bound = 1.0; mean absolute error <= 0.2 for f >= 0.5.
  P3 expected weakness: median CI width >= 0.4 at f <= 0.5 (limited resolution at arm level).
  P4 CSES7: median width >= 0.8 (curve nearly flat), reported as uninformative, not wrong.
- Predictions, real data:
  P5 CSES22 12p/12q and CSES21 21q (G-banding trisomies): estimate >= 0.7, upper bound = 1.0.
  P6 null arms (about 190, five lines, known gains excluded): at most 8 with lower bound > 0.1 (95% CIs allow ~2.5% one-sided exceedance by chance).
- Failure rule: overall coverage < 75% means the CI is wrong, and it will not be added to the CLI.

## Clone-fraction fix (written BEFORE rerunning; reps=3)
- Changes: (1) curve max < 0.15 -> f_hat/f_lo/f_hi = NaN, note UNRESOLVABLE_NOISY_SAMPLE; (2) CI includes curve spread: f_lo from the 80th-percentile curve, f_hi from the 20th.
- Predictions:
  Q1 CSES7: every spike-in run unresolved (its curve max was 0.007). Zero numeric estimates, so no wrong ones.
  Q2 coverage among resolved runs >= 85% (CSES9, HUES14; n about 72).
  Q3 median CI width at f >= 0.5 in CSES9/HUES14 grows (previously 0.2-0.4); expect 0.3-0.6.
  Q4 real trisomies: CSES22 12p/12q f_lo >= 0.5 and f_hi = 1; CSES21 21q resolved with f_hi = 1.
  Q5 null arms: at most 4 with f_lo > 0.1 across the five lines, none above 0.3.
- Failure rule: coverage < 75% -> the interval is wrong; do not add it to the CLI.

## Clone-fraction gate diagnostic (written BEFORE running)
- Question: is the response-curve max (f=1) stable across calibration seeds?
- Predictions over 10 seeds per line: CSES9, HUES14, CSES22: curve max >= 0.15 in 10/10; CSES21 (curve 0.27): 10/10;
  CSES7: >= 0.15 in at most 2/10 (but not 0/10, since 4/36 spike-in runs passed the gate).
- If CSES7 seed spread overlaps 0.15, the gate needs more donors and a stricter rule (e.g. also require the 20th-percentile curve to exceed a floor).
  Any threshold fitted on five lines is in-sample and must be labelled so.

## Gate fix attempt: more calibration donors (written BEFORE running)
- Change only n_donors 6 -> 15 (or all eligible arms if fewer). Threshold stays 0.15 (provisional, in-sample).
- Predictions over 10 seeds: CSES7 median-curve max stays < 0.15 in 10/10 (its median is 0.074; spread should shrink); CSES9, HUES14, CSES22, CSES21 stay >= 0.15 in 10/10.
- If CSES7 still reaches 0.15 in any seed, the gate needs a second condition, labelled in-sample.

## Clone-fraction validation with n_donors=15 (written BEFORE running; reps=3)
- Change: response_curve default n_donors 6 -> 15. Threshold 0.15 unchanged.
- Predictions:
  R1 CSES7 resolved in <= 1 of 36 spike-in runs (previously 4/36), and any resolved run covers the truth.
  R2 coverage among resolved runs >= 85% (previously 0.947).
  R3 CSES9/HUES14 widths at f >= 0.5 not larger than before (0.27-0.56), because curves are less noisy.
  R4 real trisomies: 12p/12q f_lo >= 0.5 with f_hi = 1; 21q resolved.
  R5 null arms: at most 4 with f_lo > 0.1, none above 0.3.
- Failure rule: coverage < 75% means do not wire the estimates into the CLI.

## Clone-fraction display rule and CLI flag (written BEFORE running)
- f_hat and f_lo shown only for arms called GAIN; all assessed arms show f_hi. CLI flag --clone-fraction (off by default; slow).
- Predictions on the five real lines: GAIN arms show f_lo >= 0.5 (12p, 12q, 21q); no non-GAIN arm shows an f_lo; CSES7 shows UNRESOLVABLE_NOISY_SAMPLE on all arms.
  Expected non-GAIN f_hi: wide (>= 0.8) on most arms in CSES7/CSES21; at most 0.5 on many arms in CSES9 and HUES14.

## Head-to-head: eSNP-Karyotyping and CaSpER vs hpsc-integrity (written BEFORE running)
- Data: the five GSE81402 lines. Truth (G-banding, Weissbein Fig 3, read by eye): CSES22 +12, CSES21 +21; CSES9, CSES7, HUES14 normal.
  Unadjudicated, not scored as true or false: CSES7 chr17, CSES9 chr20 (expression shifts), 9q (het-depleted; possible cnLOH).
- Each tool run with default settings, no tuning on these lines. Bias to state: hpsc-integrity was developed on these lines, the others were not.
- CaSpER controls: for each test line, the three G-banding-diploid lines minus itself (CSES9, CSES7, HUES14; the same line is never its own control). Decided before seeing output.
- Measured per tool: (a) the three known trisomic arms called (12p, 12q, 21q); (b) arm-sized (>= 1 arm) gain/loss events on other chromosomes, counted per line; (c) what is reported at 9q in CSES9 and HUES14; (d) whether any clone-fraction or cell-fraction output exists.
- Predictions (low confidence on B2 and B4):
  B1 both tools call chr12 in CSES22 and chr21 in CSES21.
  B2 eSNP-Karyotyping flags 9q as aberrant in at least one of CSES9, HUES14; CaSpER reports it as LOH or no CNV, not a gain.
  B3 neither tool reports a clone fraction with an interval. If either does, our clone-fraction feature is not distinctive and the README says so.
  B4 CaSpER makes >= 3 arm-sized events on other chromosomes across the five lines; hpsc-integrity makes 0 GAIN calls there.
- Decision rules: (1) a tool that handles 9q without a gain call means refusal logic is not a differentiator; say so. (2) a tool that misses a trisomy we call is reported with its settings. (3) results go in the README either way, in-sample caveat attached.

## Head-to-head amendment 1 (written BEFORE any tool is run)
- eSNP-Karyotyping, route B: use only its analysis functions on our own SNP calls (HISAT2 GRCh38 reads, bcftools calls), not its Tophat2/GATK pipeline.
  MajorMinorCalc minDP=20, maxDP=very large, minAF=0.2 (same filter as ours); PlotGenome Window=151, PValue=TRUE.
- eSNP gives plots, not calls. Measure (b) is replaced: the five PlotGenome images are renamed with random IDs, shuffled, and read blind by the project owner,
  who lists chromosomes that look gained. The key is opened after the readings are recorded. Ours-vs-theirs is then compared on chromosomes read as gained.
- eSNP's LOH tool (DeletionTable, Plot_Zygosity_*) needs dbSNP 142 files (build likely hg19, unverified) and is NOT run. Decision rule (1) is therefore restricted:
  9q can only be judged on the allelic-ratio plot, and the README will say eSNP's LOH tools exist and were not run.
- Bias to state: the allelic input is nearly identical (same depth and minor-allele thresholds), so differences are in summarisation and calling, not data.
- CaSpER: BAFExtract on our BAMs (hg38), controls and settings as pre-registered. Calls read from its segment output at arm level.

## Head-to-head amendment 2 (written BEFORE any plot is made)
- eSNP input built directly from our bcftools SNP calls (columns chr, position, AD1, AD2, DP = AD1+AD2; X=23, Y=24; other contigs dropped), bypassing EditVCF.
  MajorMinorCalc(minDP=20, maxDP=1e9, minAF=0.2) is then the same filter as ours, so allelic input is effectively identical.
- PlotGenome uses hg19 chromosome sizes/centromeres hardcoded; our positions are GRCh38. X-axis offsets are slightly off. Recorded deviation.
- Blind reading: lines are shuffled with an unseeded RNG into S1..S5; the key is written outside the repo (~/benchmark_key.csv).
  The reader lists, per image, chromosomes with a visibly raised moving median across most of the chromosome or one arm ("none" allowed), with no threshold.
  Readings are committed to results/benchmark/esnp_readings.txt BEFORE the key is opened. Limitation: the reader knows two of the five lines are trisomic, so blinding to line identity only partly removes bias.
- Predictions: the reader marks chr12 in one image and chr21 in another (B1); also marks at least one other chromosome as raised somewhere (hand-reading is generous).
- CaSpER: not started. Needs hg38 annotation/cytoband and BAFExtract pileups; set up only after the code is read.

## Head-to-head amendment 3: numeric reading replaces visual reading (written BEFORE any numeric output)
- eSNP-Karyotyping's MajorMinorCalc (minDP=20, maxDP=1e9, minAF=0.2) is run in R and its output saved; the visual reading plan is dropped.
- Numeric rule, applied by script to autosomes only: moving median of the major/minor ratio over genome-ordered SNPs, window 151 (as PlotGenome);
  baseline = median of all window values; each window is assigned to the arm of its centre SNP (GRCh38 arms from the package, not PlotGenome's hg19 table);
  an arm is "raised" if it has >= 30 SNPs and >= 50% of its windows exceed 1.25 x baseline. The thresholds are fixed now and are not tuned afterwards.
- Lines stay shuffled (S1..S5; key ~/benchmark_key.csv unopened) until the numeric results are committed.
- Predictions (B2 and B4 are low confidence):
  B1 arms 12p and 12q raised in one image, 21q in another, and no other trisomy arms.
  B2 9q raised in two images (CSES9 and HUES14): this ratio-only reading has no het-density check, so het-depleted arms look like gains.
  B4 at least 2 other arms raised across the five images (the rule has no confidence interval, so noisy lines should produce some).
- Limitations: this tests a calling rule written by us on eSNP's statistic, not the published tool's own interpretation. hpsc-integrity was developed on these lines; the rule here was not.
- CaSpER is still not run.

## Head-to-head amendment 4: exploratory window sweep (written BEFORE running; NOT the primary result)
- Primary result (window 151) is final and logged. This run changes only the window (51) and is labelled post-hoc exploratory, since it was chosen after seeing the primary result.
- Predictions: 21q raised in CSES21 (S2); 9q raised in both HUES14 and CSES9; at least 1 other arm raised across the five images (a smaller window is noisier).
- Failure of these predictions is reported the same way.

## Held-out protocol: iPSC passage time course (written BEFORE data access; tool frozen at v0.0.4)
- Data: Zhang et al. 2026, GSA-Human HRA012852/HRA012814 (access pending). iPSC#1-3, P10-P50, PE150 RNA-seq. Use only passages present in the deposit; report the full set.
- Truth (karyotype, 50 metaphases, as published): iPSC#3 +12: P10 0, P20 14/50, P30 46/50, P40 100%, P50 100%; iPSC#3 dup 20q11.2-q13.3: P30 0, P40 39/52, P50 60/60;
  iPSC#2 dup chr17 3/50 at P40 (OGM: 17q 40.56-74.17 Mb); iPSC#1 dup chr5 3/51 at P50. The line carrying the 3.4 Mb chr20 gain is inconsistent in the paper: not scored.
- Processing deviations (preprocessing only, package unchanged): paired-end reads aligned with hisat2 -1/-2 on the same GRCh38 index, same bcftools calls as before. No threshold or code change.
- Predictions (T3 and T7 low confidence):
  T1 iPSC#3 12p and 12q are GAIN at P30, P40, P50.
  T2 iPSC#3 chr12 at P20 (28%) is not GAIN, and its f_hi is >= 0.28 (the tool does not exclude the true fraction).
  T3 iPSC#3 chr12 at P10 is not GAIN.
  T4 chr12 clone-fraction interval [f_lo, f_hi] contains the karyotype fraction at >= 2 of P30, P40, P50.
  T5 No GAIN on any arm other than 12p/12q (iPSC#3 P30-P50), 20q (iPSC#3 P40-P50) and 17q (iPSC#2), across all samples; zero false GAINs on documented-normal arms.
  T6 iPSC#3 20q (dup covers about half the arm, so the arm-level ratio is diluted) is not GAIN at P40; the P50 outcome is not predicted.
  T7 No sample is flagged ALLELIC_POWER_LOW (baseline ratio <= 1.6).
- Failure rules: T1 failing means high-clone sensitivity does not replicate on real data. Any false GAIN under T5 means the specificity claim fails. All results go in the README either way, with exact binomial intervals on the karyotype fractions.
- Limits: 3 lines, a handful of passages; RNA-seq and karyotype cultures may not be the same split; karyotype counts have their own sampling error.

## Held-out protocol 2: HipSci open RNA-seq, specificity-first (written BEFORE selection or any alignment; tool frozen at v0.0.4)
- Data: ENA PRJEB7388 (HipSci iPSC RNA-seq, paired-end 75 bp). Lines matching ^HPSI\d{4}i-[a-z]{4}_\d+$; one run per line (the one with the most read pairs);
  lines with < 20,000,000 read pairs excluded (power criterion, set before seeing any result); one line per donor, chosen at random.
  N = 40 donors, numpy default_rng seed 20261030. Selection script: scripts/heldout_select.py; output results/heldout/selection.csv.
- Processing deviation (preprocessing only, package and thresholds unchanged): paired-end alignment (hisat2 -1/-2, same GRCh38 index), same bcftools calls as before.
- Truth: HipSci portal copy-number check against the donor's primary tissue. It is consulted ONLY for lines where the tool calls GAIN. A GAIN call counts as supported if the portal
  lists a copy-number gain overlapping >= 50% of the called arm's length (length used as a proxy if coordinates are not shown); otherwise it is a false positive.
- Limits: this measures specificity (and refusal and power behaviour); sensitivity on held-out data is not measured unless a line with an arm-sized gain happens to be drawn.
  The array check has its own mosaic limit, and gives regions, not clone fractions.
- Predictions:
  U1 at most 2 of 40 samples carry an unsupported GAIN call (false-positive sample rate <= 5%).
  U2 at least 1 sample has an arm refused as het-depleted (long runs of homozygosity in some donors); low confidence.
  U3 ALLELIC_POWER_LOW (baseline ratio > 1.6) is raised in <= 20% of samples.
  U4 median kept het SNPs per sample >= 5,000 (replaces an earlier shallower-reads guess, dropped).
- Failure rule: >= 5 of 40 samples with an unsupported GAIN (> 10%) means the specificity claim does not replicate; reported as such.
