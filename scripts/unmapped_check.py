"""Characterise unmapped reads among the first N FASTQ records. python scripts/unmapped_check.py SRR3501466 [N]"""
import subprocess, sys, collections
r = sys.argv[1]; N = int(sys.argv[2]) if len(sys.argv) > 2 else 500_000
reads = {}
with open(f"data/fastq/{r}.fastq") as fh:
    for _ in range(N):
        h = fh.readline()
        if not h: break
        s = fh.readline().strip(); fh.readline(); fh.readline()
        reads[h[1:].split()[0].split("/")[0]] = s
p = subprocess.Popen(["samtools", "view", "-F", "4", f"data/bam/{r}.bam"], stdout=subprocess.PIPE, text=True)
mapped = 0
for line in p.stdout:
    q = line.split("\t", 1)[0]
    if reads.get(q) is not None:
        reads[q] = None; mapped += 1
if mapped == 0: sys.exit("no read names matched between FASTQ and BAM; check name formats")
un = [s for s in reads.values() if s is not None]
print(f"{r}: {len(un)}/{len(reads)} sampled reads unmapped ({len(un)/len(reads):.1%}); mean length {sum(map(len, un))/len(un):.0f}")
print("top 32-mer prefixes:", collections.Counter(s[:32] for s in un).most_common(4))
for name, motif in [("adapter", "AGATCGGAAGAGC"), ("polyA", "A"*20), ("polyT", "T"*20)]:
    print(f"  {name}: {sum(motif in s for s in un)/len(un):.1%} of unmapped")
print("most common full sequences (paste 2-3 into NCBI BLAST):")
for s, c in collections.Counter(un).most_common(3): print(c, s)
