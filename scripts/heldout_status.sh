#!/usr/bin/env bash
# Progress of the held-out queue. bash scripts/heldout_status.sh
D=data/heldout
date
echo "done (SNP file written): $(ls $D/bam/*.snps.tsv 2>/dev/null | wc -l) of $(wc -l < $D/queue.txt)"
echo "disk used by held-out data: $(du -sh $D 2>/dev/null | cut -f1) | free: $(df -h . | tail -n 1 | awk '{print $4}')"
echo "logs with real problems (curl retry noise ignored):"
for f in $D/logs/ERR*.log; do grep -v '^curl:' "$f" 2>/dev/null | grep -qi 'error\|fail\|gave up\|md5 wrong\|larger than' && echo "  $f"; done
echo "latest completed:"; grep -h "SNP rows" $D/logs/ERR*.log 2>/dev/null | tail -n 3
echo "running:"; pgrep -fa "heldout_run.sh ERR" | grep -v "sh -c" | cut -c1-60
