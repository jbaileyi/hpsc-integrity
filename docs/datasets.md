# Validation datasets

## GSE81402 (Weissbein et al. 2016, Nat Commun 7:12144)
- BioProject PRJNA321469 / SRA study SRP075000
- RNA-seq, single-end ~84 bp. hESC lines on NextSeq 500; parthenogenetic/SWAP on HiSeq 2000
  (platform is perfectly confounded with sample type: do not pool in one baseline)
- GEO FPKM file NOT usable (retired UCSC IDs, see docs/lessons.md). Re-quantify from reads.

## Copy-number truth: five hESC lines (karyotype from SRA metadata; VERIFY vs paper Fig 3)
| Line | GSM | Run | Karyotype (metadata) | Bases |
|---|---|---|---|---|
| CSES9 | GSM2151889 | SRR3501465 | Diploid | 2.90e9 |
| CSES7 | GSM2151890 | SRR3501466 | Diploid | 3.25e9 |
| CSES22 | GSM2151891 | SRR3501467 | Trisomy 12 | 3.81e9 |
| CSES21 | GSM2151892 | SRR3501468 | Trisomy 21 | 5.49e9 |
| HUES14 | GSM2151893 | SRR3501469 | "d20q21.11" (paper text says 20q11.21: CHECK Fig 3) | 3.82e9 |

## Not copy-number truth (HiSeq 2000): GSM2151894-901 = SRR3501470-477
Parthenogenetic / SWAP samples. Labels inconsistent (swaPS2 labeled Parthenogenetic): verify before use.

## Mixing / depth experiments from the paper (for simulation benchmark)
- SRR1561108: diploid (same study as SRR1561105)
- SRR1561105: trisomy 12
- SRR1693240: trisomies 12 and 17, >50M mapped reads

## Comparison tool
- eSNPKaryotyping (R): https://github.com/BenvenLab/eSNPKaryotyping
- Last commit / installs? TODO

## Truth table verified against Figure 3 (G-banding panels, read by eye)
| Line | Run | G-banding (Fig 3) | eSNP signal (Fig 3, approx.) |
|---|---|---|---|
| CSES9 | SRR3501465 | no abnormality marked | none significant |
| CSES7 | SRR3501466 | no abnormality marked | none significant (small spike at start of chr1: note) |
| CSES22 | SRR3501467 | trisomy 12 (red box) | allelic ratio ~2.0-2.2 on chr12 |
| CSES21 | SRR3501468 | trisomy 21 (red box) | peak ~1.9 at chr21 |
| HUES14 | SRR3501469 | abnormal chr20 homolog (red arrowhead) | small signal at 20q11.21/BCL-XL; needed 51-SNP windows + common SNPs |

Notes
- SRA metadata labels HUES14 "d20q21.11"; figure inset and paper text say 20q11.21. Use 20q11.21; gain vs loss not stated in metadata.
- All five panels appear to be XY: no X-status test possible here.
- Expectation written BEFORE running: arm-level expression calls should find chr12 (CSES22) and
  probably chr21 (CSES21); HUES14 focal 20q11.21 is expected to be a miss.

## kallisto quantification (Ensembl 112, -l 200 -s 20)
| Line | Run | Reads processed | Pseudoaligned | % |
|---|---|---|---|---|
| CSES9 | SRR3501465 | 34,372,963 | 26,439,916 | 76.9 |
| CSES7 | SRR3501466 | 38,449,904 | 19,481,278 | 50.7 |
| CSES22 | SRR3501467 | 45,666,767 | 37,615,635 | 82.4 |
| CSES21 | SRR3501468 | 65,859,863 | 45,762,736 | 69.5 |
| HUES14 | SRR3501469 | 45,162,849 | 35,928,017 | 79.6 |
- CSES7 is the low outlier (50.7%). Cause not investigated; first suspect if its arm-level profile is noisy.

## First expression profile (median-of-5 baseline, MIN_COUNTS=20, window=200 genes; no tuning)
- CSES22 chr12 median +0.48 (trisomy 12 detected); CSES21 chr21 +0.43 (trisomy 21 detected).
- Expected for one extra copy: log2(3/2) = 0.585; observed ~82% and ~74% of that.
- HUES14 chr20: 0.00 (focal 20q11.21 not detected, as predicted beforehand).
- Candidates to check, not findings: CSES7 chr17 +0.11, HUES14 chr17 -0.09, CSES21 chrX -0.12.
- Caveat: with 5 lines the median baseline pins ~20% of values at exactly 0 (flat lines, sharp steps).

## Leave-one-out profile (MIN_COUNTS=20, window=200)
- CSES22 chr12 +0.55; CSES21 chr21 +0.53. Next-largest autosome deviations: CSES7 chr17 +0.22,
  HUES14 chr17 -0.20 (both also present, smaller, under median-of-5). CSES21 chrX -0.23 (X not
  assessed in v0.1).
- G-banding calls CSES7 and HUES14 normal. Chr17 deviations unexplained: candidate false positives
  until checked against Fig 3 allelic panels.

## Arm-level results (leave-one-out, per-sample offset removed; 12p/12q/21q from cytoBand hg38)
- Detected: CSES22 12p +0.50, 12q +0.55; CSES21 21q +0.53.
- Acrocentric short arms (13p, 14p, 15p, 21p, 22p) have no protein-coding genes: not assessable.
- Y arms (5 and 9 genes) give nonsense values (Yq +1.41 CSES7, -0.75 CSES22): minimum-genes rule needed.
- 20q11.21 band (hg38 chr20:30.4-33.5 Mb, 21 genes): CSES9 +0.38, CSES7 -0.29, CSES21 -0.25,
  CSES22 +0.17, HUES14 +0.13 => not callable by expression (HUES14 inside the noise).
- Candidate false positives: CSES7 chr17 (both arms +0.16/+0.20), HUES14 17q -0.17, CSES21 18q -0.17.
  If expression ~ copy number, these imply ~25-33% subclones, which a 20-metaphase karyotype would
  almost certainly detect (number of metaphases scored here unknown).

## Figure 3 check of chr17 (read by eye from the figure image)
- CSES7 and HUES14 allelic traces are flat across chr17 (~1.3); no significant marks there.
  Only red mark in CSES7 panel: start of chr1.
- Caveat: the paper's own mixing test saw a 1/3 mixture as visible but not significant, so a 25-33%
  subclone would be near the allelic method's limit too. Chr17 expression shifts remain unexplained.
