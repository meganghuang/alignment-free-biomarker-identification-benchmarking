
import os, pandas as pd, numpy as np, scanpy as sc, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
SK="/u/scratch/m/mhuang/claude/compare3/sckar_pbmc"
OUT=SK+"/summary"; os.makedirs(OUT, exist_ok=True)
sc.settings.figdir = OUT

df = pd.read_csv(SK+"/data/expression_matrix/expression.csv", index_col=0)   # gene-major
ad = sc.AnnData(df.T)                                                        # cells x genes
ad.obs["pseudocell"] = ad.obs_names
pref = {"B":"B cells","DC":"Dendritic cells","Mono":"Monocytes","NK":"NK","Plt":"Platelets","T":"T cells"}
ad.obs["celltype"] = [pref[n.rsplit("_",1)[0]] for n in ad.obs_names]
cl = pd.read_csv(SK+"/data/clustering_results/cluster.csv").set_index("cell")["cluster"]
ad.obs["leiden_scKAR"] = [str(cl[n]) for n in ad.obs_names]
bp = pd.read_csv(SK+"/data/clustering_results/bipartitions/bipartition_0.csv", sep="\t", header=None, names=["cell","cond"]).set_index("cell")["cond"]
ad.obs["condition"] = [bp[n] for n in ad.obs_names]

sc.pp.normalize_total(ad); sc.pp.log1p(ad)
sc.pp.pca(ad, n_comps=20)
sc.pp.neighbors(ad, n_neighbors=10, n_pcs=20)
sc.tl.umap(ad)
coords = pd.DataFrame(ad.obsm["X_umap"], columns=["UMAP_1","UMAP_2"], index=ad.obs_names)
coords = coords.join(ad.obs[["celltype","leiden_scKAR","condition"]])
coords.index.name="pseudocell"
coords.to_csv(OUT+"/umap_coords_sckar.csv")

for key, fname in [("celltype","umap_sckar_celltype"),("leiden_scKAR","umap_sckar_leiden"),("condition","umap_sckar_condition")]:
    fig, ax = plt.subplots(figsize=(6,5))
    for lev in sorted(coords[key].unique()):
        m = coords[key]==lev
        ax.scatter(coords.loc[m,"UMAP_1"], coords.loc[m,"UMAP_2"], s=70, label=str(lev), edgecolor="k", linewidth=0.3)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_title(f"scKAR pseudo-cells (n=52) — {key}")
    ax.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.01,0.5))
    fig.tight_layout()
    fig.savefig(f"{OUT}/{fname}.png", dpi=150); fig.savefig(f"{OUT}/{fname}.pdf"); plt.close(fig)

ct = pd.crosstab(coords["celltype"], coords["leiden_scKAR"])
ct.to_csv(OUT+"/celltype_by_leiden.csv")
pd.crosstab(coords["celltype"], coords["condition"]).to_csv(OUT+"/celltype_by_condition.csv")
print(ct.to_string())
print(pd.crosstab(coords["celltype"], coords["condition"]).to_string())
print("WROTE umap", coords.shape)
