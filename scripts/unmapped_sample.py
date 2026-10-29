"""Write 40 random unmapped reads (no Ns, len>=60) from the first N FASTQ records to FASTA for NCBI BLAST.
    python scripts/unmapped_sample.py SRR3501466"""
import random, subprocess, sys
r = sys.argv[1]; N = 500_000; random.seed(1)
reads = {}
with open(f"data/fastq/{r}.fastq") as fh:
    for _ in range(N):
        h = fh.readline()
        if not h: break
        s = fh.readline().strip(); fh.readline(); fh.readline()
        reads[h[1:].split()[0].split("/")[0]] = s
p = subprocess.Popen(["samtools", "view", "-F", "4", f"data/bam/{r}.bam"], stdout=subprocess.PIPE, text=True)
for line in p.stdout:
    reads.pop(line.split("\t", 1)[0], None)
un = [(k, s) for k, s in reads.items() if "N" not in s and len(s) >= 60]
pick = random.sample(un, 40)
with open(f"results/qc/{r}_unmapped_sample.fa", "w") as o:
    for k, s in pick: o.write(f">{k}\n{s}\n")
print(len(un), "eligible unmapped reads; wrote 40 to results/qc/")
