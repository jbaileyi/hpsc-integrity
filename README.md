# hpsc-integrity

Screen human pluripotent stem cell (hPSC) RNA-seq for chromosome-arm gains and losses, with optional allelic-imbalance evidence from heterozygous SNPs. Research prototype (v0.0.x); not a replacement for karyotyping or SNP arrays.

## What it does

- **Expression tiers** (`hpsc-integrity counts.csv -o out/`): per-arm log2 ratios against a leave-one-out baseline, called as `high` / `possible` / none. Gains and losses. Needs at least 3 samples.
- **Allelic evidence** (add `--snps SAMPLE=snps.tsv`, repeatable): per-arm calls from major/minor allele ratios at heterozygous SNPs, with locus-level bootstrap CIs. States: `GAIN`, `NORMAL`, `UNDERPOWERED`, `REFUSED_HET_DEPLETED`, `TOO_FEW_SNPS`. Gains only; a loss depletes heterozygous sites instead of shifting ratios, so there is no `LOSS` state.
- **Clone-fraction bounds** (add `--clone-fraction` with `--snps`; slow, because each sample simulates gains on its own normal arms to calibrate): for `GAIN` arms, lower and upper bounds on the fraction of cells carrying the extra copy; for other assessed arms, an upper bound on the clone size that is not excluded. Noisy samples (e.g. CSES7) get `UNRESOLVABLE_NOISY_SAMPLE` and no numbers.

## Validation (five hESC lines, GSE81402; three gains per the G-banding in Weissbein et al. 2016 Fig 3)

| Check | Result |
|---|---|
| Positive controls (G-banding trisomies) | CSES22 12p, 12q and CSES21 21q called `GAIN` (3 of 3) |
| Null arms | 0 `GAIN` calls in about 190 null arm-tests; 0 false GAINs elsewhere in every spike-in run |
| Het-depleted regions | 9q in CSES9 and HUES14, and 16p in CSES22, refused |
| Allelic spike-in detection limit | `GAIN` reliable at about 75% clone fraction in CSES9 and HUES14; 50% gives 0-50%; nothing at 30% or below |
| Expression-tier mosaic limit (spike-in, 3 lines, grid-resolution) | `possible` call from about 20-30% clone fraction, `high` from about 50% (1q, 12p, 10p, 18q); 17q needs about 50% / 75%. Depth 0.1-1x changes this little |
| Expression-tier false calls (spike-in null) | about 0.2 `high` and 2.3-3.3 `possible` arms per sample (more at low depth) |

## Known limitations and failures

- **All validation is in-sample:** five lines from one study, three true positives, and the CI rule was fixed after seeing the CIs on these lines.
- **Allelic calls are a large-clone detector.** 55-81% of arms are `UNDERPOWERED` at this depth. Do not read `NORMAL` or `UNDERPOWERED` as "no abnormality".
- **CSES7 never reaches `GAIN`, even at a simulated full trisomy.** Its baseline allelic ratio (1.82 vs 1.40-1.43) compresses the effect. The `ALLELIC_POWER_LOW` flag (baseline ratio > 1.6) is set from this one line; CSES21 sits exactly at 1.60.
- **19q refusals are an artifact** of the single-sample density metric (19q is low in all five lines). A per-arm expected-density table is not yet bundled.
- **9q is likely long homozygosity (possible copy-neutral LOH), not a gain.** Unconfirmed.
- **Arm-level only:** focal events such as 20q11.21 are not detectable. Losses are not called by the allelic module.
- **Expression `high` calls are not reliable alone** (about 1 false `high` per 5 samples in the null simulation); treat them as leads to confirm.
- **Clone-fraction bounds are coarse and only simulation-validated.** Spike-in coverage 0.97 (70 of 72 resolved runs; two clean lines; one read-shifting model) with interval widths of about 0.3-0.7. All real positives are full trisomies, so partial fractions have no model-independent test. Estimates saturate near 100% when the real shift exceeds the simulated one (CSES21 21q). The 0.15 threshold for "resolvable" was set from five lines. On a clean line, a typical arm not called a gain still cannot exclude a clone of about 50-70%.
- **Expression and allelic evidence are not combined** into one call; they can disagree (e.g. CSES7 17p/17q).
- The spike-in shift model matches one real trisomy (simulated f=1 about +0.43 vs real +0.43/+0.53).
- CSES7's 51% alignment rate is unexplained. Unmapped reads were not adapter, polyA or all-N dominated (a sample was written to `results/qc/` for BLAST but not classified). Cause not investigated further.

## Development record

Predictions were written to `docs/design.md` before each experiment, and results, including failures, are in `docs/lessons.md`.

## License and data

Code: MIT (see LICENSE). Bundled reference tables are derived coordinate tables: gene coordinates from Ensembl release 112 (GRCh38) and chromosome arms from the UCSC hg38 cytoBand table. Provider terms and links are recorded in `src/hpsc_integrity/data/SOURCES.txt`. Please cite Ensembl and the UCSC Genome Browser if you use this in your work.
