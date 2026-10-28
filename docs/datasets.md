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
