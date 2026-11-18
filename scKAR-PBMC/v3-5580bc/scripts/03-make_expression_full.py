#!/usr/bin/env python3
"""05-make_expression_full.py — gene x cell matrix for the full 5,580-cell run.

scKAR clusters cells from a gene expression matrix (MODE=GENE_EXPRESSION); the
k-mer matrix is never used for clustering. This reads the CellRanger filtered h5,
subsets to the QC-passing barcodes in meta/cells_full.tsv, and writes
data/expression_matrix/expression.csv in the genes-as-rows, cells-as-columns
layout the clustering step expects (column names are bare 16 nt barcodes).
"""
import os, sys
import numpy as np, pandas as pd, scanpy as sc

F   = sys.argv[1] if len(sys.argv) > 1 else "/u/scratch/m/mhuang/claude/compare3/sckar_full"
H5  = ("/u/scratch/m/mhuang/std_scripts/counts/"
       "5k_Human_Donor4_PBMC_3p_gem-x_5k_Human_Donor4_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.h5")
OUT = os.path.join(F, "data", "expression_matrix")
os.makedirs(OUT, exist_ok=True)

cells = pd.read_csv(os.path.join(F, "meta", "cells_full.tsv"), sep="\t")
ad = sc.read_10x_h5(H5)
ad.var_names_make_unique()
print("h5:", ad.shape, flush=True)

present = [b for b in cells.barcode if b in set(ad.obs_names)]
missing = len(cells) - len(present)
print(f"requested={len(cells)} found_in_h5={len(present)} missing={missing}", flush=True)
if missing:
    raise SystemExit(f"{missing} barcodes absent from the h5 — aborting rather than silently dropping cells")

ad = ad[present].copy()
X = ad.X.toarray().T.astype(np.float32)      # genes x cells
df = pd.DataFrame(X, index=ad.var_names, columns=[b.split("-")[0] for b in present])
print("matrix:", df.shape, flush=True)
df.to_csv(os.path.join(OUT, "expression.csv"))
print("wrote", os.path.getsize(os.path.join(OUT, "expression.csv")), "bytes", flush=True)
