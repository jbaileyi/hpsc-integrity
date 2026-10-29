#!/usr/bin/env bash
# Call SNPs from an aligned BAM and write a tab-separated table of allele depths.
# Usage (from the repository root): bash scripts/call_snps.sh SRR3501467
# Needs data/bam/<run>.bam and data/hisat2/genome.fa (with .fai).
set -euo pipefail
r=$1
bcftools mpileup -f data/hisat2/genome.fa -q 20 -Q 20 -d 2000 -a FORMAT/AD -Ou "data/bam/$r.bam" \
  | bcftools call -mv -Ou \
  | bcftools view -v snps -m2 -M2 -i 'QUAL>=30' -Ob -o "data/bam/$r.snps.bcf"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "data/bam/$r.snps.bcf" > "data/bam/$r.snps.tsv"
echo "$r: $(wc -l < data/bam/$r.snps.tsv) SNP rows written to data/bam/$r.snps.tsv"
