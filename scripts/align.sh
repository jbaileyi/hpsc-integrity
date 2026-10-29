#!/usr/bin/env bash
# Align single-end FASTQs with HISAT2 and write sorted, indexed BAMs.
# Usage (from the repository root): bash scripts/align.sh SRR3501465 SRR3501466 ...
set -euo pipefail
IDX=data/hisat2/grch38/genome
mkdir -p data/bam
for r in "$@"; do
  hisat2 -p 8 -x "$IDX" -U "data/fastq/$r.fastq" --no-unal --summary-file "data/bam/$r.hisat2.txt" \
    | samtools sort -@ 4 -o "data/bam/$r.bam" -
  samtools index "data/bam/$r.bam"
  echo "$r: $(grep 'overall alignment rate' data/bam/$r.hisat2.txt)"
done
