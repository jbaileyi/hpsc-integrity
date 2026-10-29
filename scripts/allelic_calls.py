"""Allelic arm calls for all five lines -> results/allelic/calls_<line>.csv. python scripts/allelic_calls.py"""
import pandas as pd
from hpsc_integrity.allelic import filter_het_snps, read_allele_counts
from hpsc_integrity.allelic_calls import call_arms
from hpsc_integrity.annotation import load_reference

M = {"CSES9": "SRR3501465", "CSES7": "SRR3501466", "CSES22": "SRR3501467", "CSES21": "SRR3501468", "HUES14": "SRR3501469"}
WATCH = {"CSES22": ["12p", "12q"], "CSES21": ["21q"], "CSES7": ["17p", "17q"], "CSES9": ["20p", "20q", "9q"], "HUES14": ["9q"]}
genes, arms = load_reference()
pd.set_option("display.width", 200)
for line, srr in M.items():
    r = call_arms(filter_het_snps(read_allele_counts(f"data/bam/{srr}.snps.tsv")), genes, arms)
    r.to_csv(f"results/allelic/calls_{line}.csv")
    print(f"\n{line}: {r.call.value_counts().to_dict()}")
    print(r.loc[[a for a in WATCH[line] if a in r.index], ["n_snps", "n_loci", "rel_het_density", "log2", "ci_lo", "ci_hi", "call"]].to_string())
    print("  lowest relative het density:", r.rel_het_density.nsmallest(6).round(2).to_dict())
    print("  GAIN calls outside the watch arms:", [a for a in r.index[r.call == "GAIN"] if a not in WATCH[line]])
