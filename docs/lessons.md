- Version-stripped match to hg38 new-style uc IDs: 47,682 / 104,133 (46%). Too low,
  and same-ID-same-transcript is unverified after the UCSC annotation overhaul.
  Decision: do not use the GEO FPKM file for coordinates; re-quantify from raw reads
  (SRP075000) with a current annotation instead.
