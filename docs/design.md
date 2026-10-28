# Design decisions

## Reference (pinned)
- Ensembl release 112, GRCh38. FASTA: cdna.all; annotation: Homo_sapiens.GRCh38.112.gtf.gz.
- Reason: the GEO FPKM file used retired UCSC IDs (docs/lessons.md); re-quantifying from reads
  with a documented build makes coordinates unambiguous.
- Gene table will keep chromosomes 1-22, X, Y only (patches/alt scaffolds dropped).
