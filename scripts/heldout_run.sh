#!/usr/bin/env bash
# One held-out run: download (md5-checked, resumable), align paired-end to GRCh38, call SNPs as in scripts/call_snps.sh.
# Usage (from the repository root): bash scripts/heldout_run.sh ERR1555117
set -euo pipefail
r=$1
D=data/heldout; IDX=data/hisat2/grch38/genome
mkdir -p $D/fastq $D/bam
if [ -s "$D/bam/$r.snps.tsv" ]; then echo "$r: already done"; exit 0; fi
IFS=$'\t' read -r run name u1 u2 m1 m2 < <(awk -F'\t' -v r="$r" '$1==r' $D/jobs.tsv)
[ -n "${run:-}" ] || { echo "$r: not in jobs.tsv" >&2; exit 1; }
for i in 1 2; do
  eval u=\$u$i; eval m=\$m$i
  f="$D/fastq/${run}_${i}.fastq.gz"
  if [ ! -s "$f" ] || [ "$(md5sum "$f" | cut -d' ' -f1)" != "$m" ]; then
    curl -sS -L -C - --retry 5 -o "$f" "https://$u" || curl -sS -L -C - --retry 5 -o "$f" "ftp://$u"
  fi
  [ "$(md5sum "$f" | cut -d' ' -f1)" = "$m" ] || { echo "$run: md5 mismatch on mate $i" >&2; exit 1; }
done
hisat2 -p 2 -x "$IDX" -1 "$D/fastq/${run}_1.fastq.gz" -2 "$D/fastq/${run}_2.fastq.gz" --no-unal \
  --summary-file "$D/bam/$run.hisat2.txt" | samtools sort -@ 1 -o "$D/bam/$run.bam" -
samtools index "$D/bam/$run.bam"
bcftools mpileup -f data/hisat2/genome.fa -q 20 -Q 20 -d 2000 -a FORMAT/AD -Ou "$D/bam/$run.bam" \
  | bcftools call -mv -Ou \
  | bcftools view -v snps -m2 -M2 -i 'QUAL>=30' -Ob -o "$D/bam/$run.snps.bcf"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "$D/bam/$run.snps.bcf" > "$D/bam/$run.snps.tsv.tmp"
mv "$D/bam/$run.snps.tsv.tmp" "$D/bam/$run.snps.tsv"
rm -f "$D/fastq/${run}_1.fastq.gz" "$D/fastq/${run}_2.fastq.gz"
echo "$run ($name): $(grep 'overall alignment rate' $D/bam/$run.hisat2.txt) | $(wc -l < $D/bam/$run.snps.tsv) SNP rows"
