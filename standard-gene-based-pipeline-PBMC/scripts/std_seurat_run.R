
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(ggplot2); library(patchwork)})
set.seed(7)
DATA <- "/u/scratch/m/mhuang/std_scripts/counts"
OUT  <- "/u/scratch/m/mhuang/claude/compare3/std"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
sink(file.path(OUT,"run_log.txt"), split=TRUE)
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("hdf5r:", tryCatch(as.character(packageVersion("hdf5r")), error=function(e) "ABSENT"), "\n")

h5 <- file.path(DATA,"5k_Human_Donor4_PBMC_3p_gem-x_5k_Human_Donor4_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.h5")
m <- Read10X_h5(h5); if (is.list(m)) m <- m[[1]]
cat("raw matrix:", nrow(m), "genes x", ncol(m), "barcodes\n")

sample <- CreateSeuratObject(counts=m, project="StdCounts", min.cells=1, min.features=100)
cat("after CreateSeuratObject:", ncol(sample), "cells,", nrow(sample), "genes\n")

cm <- GetAssayData(sample, layer="counts")
mito.features <- grep("^MT-", rownames(sample), value=TRUE)
cat("n mito features:", length(mito.features), "\n")
sample[["percent.mito"]] <- Matrix::colSums(cm[mito.features,]) / Matrix::colSums(cm)
rm(cm); gc()

p <- VlnPlot(sample, features=c("nFeature_RNA","nCount_RNA","percent.mito"), group.by="orig.ident", ncol=3, pt.size=0.4)
ggsave(file.path(OUT,"qc_violin.pdf"), p, width=15, height=6)
ggsave(file.path(OUT,"qc_violin.png"), p, width=15, height=6, dpi=150)

sample <- subset(sample, subset = nCount_RNA <= 5e4 & percent.mito <= 0.2)
cat("after QC filter:", ncol(sample), "cells\n")
sample <- NormalizeData(sample, normalization.method="LogNormalize", scale.factor=7000)
sample <- FindVariableFeatures(sample, selection.method="vst", nfeatures=2000)
sample <- ScaleData(sample)
sample <- RunPCA(sample, verbose=FALSE)
p <- ElbowPlot(sample, ndims=50); ggsave(file.path(OUT,"elbow.pdf"), p, width=6, height=5)

set.seed(7)
dims_use <- 1:30
sample <- RunUMAP(sample, reduction="pca", dims=dims_use, verbose=FALSE)
sample <- FindNeighbors(sample, reduction="pca", dims=dims_use, verbose=FALSE)
sample <- FindClusters(sample, resolution=0.4, verbose=FALSE)
cat("n clusters:", nlevels(sample$seurat_clusters), "\n")
print(table(sample$seurat_clusters))

# attach the authoritative per-barcode labels from the user's own standard run
annot <- read.csv(file.path(DATA,"barcode_celltype_STD.csv"), stringsAsFactors=FALSE)
rownames(annot) <- annot$Barcode
sample$CellType        <- annot[colnames(sample), "CellType"]
sample$GeneralCellType <- annot[colnames(sample), "GeneralCellType"]
cat("labelled cells:", sum(!is.na(sample$CellType)), "of", ncol(sample), "\n")
print(table(sample$GeneralCellType, useNA="ifany"))

mk <- function(p, f, w=7, h=6){ ggsave(file.path(OUT,paste0(f,".pdf")), p, width=w, height=h)
                                 ggsave(file.path(OUT,paste0(f,".png")), p, width=w, height=h, dpi=150) }
mk(DimPlot(sample, label=TRUE, reduction="umap", pt.size=0.1, shuffle=TRUE, raster=FALSE), "umap_cluster")
mk(DimPlot(sample, reduction="umap", group.by="CellType", label=TRUE, pt.size=0.1, shuffle=TRUE, raster=FALSE), "umap_CellType", 8, 6)
mk(DimPlot(sample, reduction="umap", group.by="GeneralCellType", label=TRUE, pt.size=0.1, shuffle=TRUE, raster=FALSE), "umap_GeneralCellType", 8, 6)
p <- VlnPlot(sample, features=c("nFeature_RNA","nCount_RNA","percent.mito"), ncol=3, pt.size=0)
ggsave(file.path(OUT,"qc_by_cluster.pdf"), p, width=18, height=6, limitsize=FALSE)
ggsave(file.path(OUT,"qc_by_cluster.png"), p, width=18, height=6, dpi=150, limitsize=FALSE)

# export embedding + metadata for cross-pipeline figures
emb <- as.data.frame(Embeddings(sample, "umap"))
colnames(emb) <- c("UMAP_1","UMAP_2")
emb$Barcode <- rownames(emb)
emb$cluster <- as.character(sample$seurat_clusters)
emb$CellType <- sample$CellType
emb$GeneralCellType <- sample$GeneralCellType
emb$nCount_RNA <- sample$nCount_RNA
emb$nFeature_RNA <- sample$nFeature_RNA
write.csv(emb[,c("Barcode","UMAP_1","UMAP_2","cluster","CellType","GeneralCellType","nCount_RNA","nFeature_RNA")],
          file.path(OUT,"umap_coords_std.csv"), row.names=FALSE)
write.csv(as.data.frame(table(sample$GeneralCellType)), file.path(OUT,"cells_per_generaltype.csv"), row.names=FALSE)
write.csv(as.data.frame(table(sample$CellType)),        file.path(OUT,"cells_per_celltype.csv"),    row.names=FALSE)

cat("=== FindAllMarkers ===\n")
cluster_markers <- FindAllMarkers(sample, only.pos=TRUE, logfc.threshold=0.25)
write.csv(cluster_markers, file.path(OUT,"cluster_markers_std.csv"), row.names=FALSE)
cat("marker rows:", nrow(cluster_markers), " unique genes:", length(unique(cluster_markers$gene)), "\n")
sig <- cluster_markers[cluster_markers$p_val_adj < 0.05, ]
write.csv(sig, file.path(OUT,"cluster_markers_std_sig.csv"), row.names=FALSE)
cat("significant rows:", nrow(sig), " unique genes:", length(unique(sig$gene)), "\n")

# markers per CELL TYPE (the unit the three pipelines are compared on)
Idents(sample) <- "GeneralCellType"
ct_markers <- FindAllMarkers(sample, only.pos=TRUE, logfc.threshold=0.25)
write.csv(ct_markers, file.path(OUT,"celltype_markers_std.csv"), row.names=FALSE)
cat("celltype marker rows:", nrow(ct_markers), "unique genes:", length(unique(ct_markers$gene)), "\n")
print(table(ct_markers$cluster))
Idents(sample) <- "seurat_clusters"

d <- sig[order(sig$cluster, -sig$avg_log2FC), ]
top <- do.call(rbind, lapply(split(d, d$cluster), function(x) head(x,3)))
genes <- unique(as.character(top$gene))
p <- DotPlot(sample, features=rev(genes), cols=c("grey90","red")) + coord_flip() + labs(x=NULL,y=NULL)
ggsave(file.path(OUT,"dotplot_cluster_markers.pdf"), p, width=12, height=12)
ggsave(file.path(OUT,"dotplot_cluster_markers.png"), p, width=12, height=12, dpi=150)

saveRDS(sample, file.path(OUT,"std_seurat5.rds"))
cat("=== sessionInfo ===\n"); print(sessionInfo())
sink()
