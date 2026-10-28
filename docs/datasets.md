# Validation datasets

## GSE81402 (Weissbein et al. 2016, Nat Commun 7:12144)
- BioProject: PRJNA321469 / SRA: SRP075000
- Platforms: GPL11154 (HiSeq 2000), GPL18573 (NextSeq 500) -- TODO: which samples are on which
- Processed file: GSE81402_Samples_FPKM.txt.gz (FPKM, 13 samples)
- Samples with G-banding truth: GSM2151889-893 = CSES9, CSES7, CSES22, CSES21, HUES14
- Other samples (parthenogenetic / SWAP, not copy-number truth): GSM2151894-901
  = swaPS2, pES2, pES7, pES10, swaPS1, swaPS6, swaPS9, swaPS10

| Line | G-banding result | eSNP-Karyotyping result | Source |
|---|---|---|---|
| CSES9 | TODO | TODO | Fig 3 |
| CSES7 | TODO | TODO | Fig 3 |
| CSES22 | TODO | TODO | Fig 3 |
| CSES21 | TODO | TODO | Fig 3 |
| HUES14 | TODO (small chr20 region) | TODO | Fig 3 |

## Mixing / depth experiments from the paper (for simulation benchmark)
- SRR1561108: diploid (same study as SRR1561105)
- SRR1561105: trisomy 12
- SRR1693240: trisomies 12 and 17, >50M mapped reads

## Comparison tool
- eSNPKaryotyping (R): https://github.com/BenvenLab/eSNPKaryotyping
- Last commit / installs? TODO
