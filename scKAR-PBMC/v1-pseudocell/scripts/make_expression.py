
import sys, os
import numpy as np, pandas as pd, scanpy as sc
H5  = "/u/scratch/m/mhuang/std_scripts/counts/5k_Human_Donor4_PBMC_3p_gem-x_5k_Human_Donor4_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.h5"
MAP = "/u/scratch/m/mhuang/claude/compare3/pseudocell_map.csv"
OUT = "/u/scratch/m/mhuang/claude/compare3/sckar_pbmc/data/expression_matrix"
os.makedirs(OUT, exist_ok=True)

ad = sc.read_10x_h5(H5)
ad.var_names_make_unique()
print("h5:", ad.shape, flush=True)

m = pd.read_csv(MAP)
m["bc"] = m.Barcode
common = [b for b in m.bc if b in set(ad.obs_names)]
print("map barcodes:", len(m), " found in h5:", len(common), flush=True)
m = m[m.bc.isin(set(ad.obs_names))].copy()

ad = ad[m.bc.values].copy()
X = ad.X.tocsc() if hasattr(ad.X, "tocsc") else ad.X
pcs = sorted(m.pseudocell.unique())
idx = {pc: np.where(m.pseudocell.values == pc)[0] for pc in pcs}

import scipy.sparse as sp
Xc = ad.X.tocsr()
mat = np.zeros((len(pcs), ad.n_vars), dtype=np.float64)
for i, pc in enumerate(pcs):
    mat[i] = np.asarray(Xc[idx[pc]].sum(axis=0)).ravel()
counts = pd.DataFrame(mat, index=pcs, columns=ad.var_names)
print("pseudocell x gene:", counts.shape, flush=True)

tot = counts.sum(axis=1)
tpm = counts.div(tot, axis=0) * 1e6
# scKAR's graph-clustering.py does read_csv(index_col=0) then .T -> needs GENE-MAJOR
gene_major = tpm.T
gene_major.index.name = "gene"
gene_major.to_csv(os.path.join(OUT, "expression.csv"))
counts.T.to_csv(os.path.join(OUT, "pseudocell_counts.csv"))

summ = pd.DataFrame({
    "pseudocell": pcs,
    "n_barcodes": [len(idx[pc]) for pc in pcs],
    "total_umi": [int(tot[pc]) for pc in pcs],
    "genes_detected": [int((counts.loc[pc] > 0).sum()) for pc in pcs],
    "celltype": [m.loc[m.pseudocell == pc, "GeneralCellType"].iloc[0] for pc in pcs],
})
summ.to_csv("/u/scratch/m/mhuang/claude/compare3/sckar_pbmc/pseudocell_summary.csv", index=False)
print(summ.to_string(), flush=True)
print("WROTE", os.path.join(OUT, "expression.csv"), gene_major.shape, flush=True)
