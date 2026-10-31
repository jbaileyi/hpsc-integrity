#!/usr/bin/env bash
# One held-out run: resumable download (retries until md5 matches), paired-end alignment to GRCh38, SNP calls as in scripts/call_snps.sh.
# Usage (from the repository root): bash scripts/heldout_run.sh ERR1555117
set -euo pipefail
r=$1
D=data/heldout; IDX=data/hisat2/grch38/genome
mkdir -p $D/fastq $D/bam
if [ -s "$D/bam/$r.snps.tsv" ]; then echo "$r: already done"; exit 0; fi
IFS=$'\t' read -r run name u1 u2 m1 m2 b1 b2 < <(awk -F'\t' -v r="$r" '$1==r' $D/jobs.tsv)
[ -n "${run:-}" ] || { echo "$r: not in jobs.tsv" >&2; exit 1; }
CURL="curl -sS -L -C - --connect-timeout 30 --speed-limit 100000 --speed-time 60"
sum() { md5sum "$1" 2>/dev/null | cut -d' ' -f1; }
fetch() {   # url md5 bytes outfile; delete only a complete-size file with a bad md5, or an oversized file
  local url=$1 md5=$2 want=$3 f=$4 n have
  for n in $(seq 1 400); do
    have=$(stat -c %s "$f" 2>/dev/null || echo 0)
    if [ "$have" -eq "$want" ]; then
      if [ "$(sum "$f")" = "$md5" ]; then return 0; fi
      echo "$run: full size but md5 wrong; restarting $f" >&2; rm -f "$f"
    elif [ "$have" -gt "$want" ]; then
      echo "$run: file larger than expected; restarting $f" >&2; rm -f "$f"
    fi
    $CURL -o "$f" "https://$url" || $CURL -o "$f" "ftp://$url" || sleep 5
  done
  echo "$run: gave up after 400 attempts: $f" >&2; return 1
}
fetch "$u1" "$m1" "$b1" "$D/fastq/${run}_1.fastq.gz"
fetch "$u2" "$m2" "$b2" "$D/fastq/${run}_2.fastq.gz"
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
