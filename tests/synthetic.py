"""Small synthetic dataset: 4 chromosomes x 2 arms x 100 genes, 5 samples, S5 has chr1 at 1.5x."""
import numpy as np
import pandas as pd


def make_synthetic(seed=0, trisomic_sample="S5", trisomic_chrom="1", n_samples=5):
    rng = np.random.default_rng(seed)
    rows = []
    for c in ["1", "2", "3", "4"]:
        for j in range(100):                       # p arm: positions 1000-1995
            rows.append((c, 1000 + 10 * j, 1005 + 10 * j))
        for j in range(100):                       # q arm: positions 10000-10995
            rows.append((c, 10000 + 10 * j, 10005 + 10 * j))
    genes = pd.DataFrame(rows, columns=["chrom", "start", "end"],
                         index=[f"ENSG{i:05d}.1" for i in range(len(rows))])
    genes["symbol"] = [f"G{i}" for i in range(len(rows))]
    genes.index.name = "gene_id"
    arms = pd.DataFrame(
        [(c, 0, 5000, c + "p") for c in "1234"] + [(c, 5000, 20000, c + "q") for c in "1234"],
        columns=["chrom", "start", "end", "arm"],
    )
    means = rng.uniform(100, 1000, len(genes))
    cols = [f"S{i}" for i in range(1, n_samples + 1)]
    counts = pd.DataFrame(index=genes.index, columns=cols, dtype=float)
    for s in cols:
        scale = np.where((genes["chrom"] == trisomic_chrom).values & (s == trisomic_sample), 1.5, 1.0)
        counts[s] = rng.poisson(means * scale)
    return counts, genes, arms
