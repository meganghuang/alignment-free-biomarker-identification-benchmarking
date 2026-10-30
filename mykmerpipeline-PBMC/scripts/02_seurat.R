# ============================================================
# Script: 2_kmer_seuratv2.R
# Author: Megan Huang
#
# Description: [CURRENTLY UNDER REVISION]
# Performs single-cell analysis using Seurat on a k-mer count
# matrix. Includes quality control, normalization, dimensionality
# reduction (PCA, UMAP), clustering, marker detection, and
# visualization. Generates plots and cluster-specific marker files.
#
# Inputs:
#   From DATA_DIR:
#     - kmer_matrix.mtx        : Sparse k-mer count matrix
#     - kmer_ids.txt           : Row names (k-mers)
#     - cell_barcodes.txt      : Column names (cells)
#
#   Optional external inputs:
#     - annotation.csv                 : Cluster → cell type mapping
#     - barcode_celltype_STD.csv       : Barcode → cell type annotations
#     - (optional) cellranger_matrix.h5: For standard (non-kmer) pipeline
#
# Outputs:
#   In DATA_DIR/plots/:
#     - quality_by_origin.pdf
#     - quality_by_cluster.pdf
#     - umap_cluster.pdf
#     - umap_origin.pdf
#     - umap_celltype.pdf
#     - umap_celltype_STD.pdf
#     - umap_celltype_from_csv.pdf
#     - dotplot_cluster_markers.pdf
#     - Azimuth_FeaturePlots_AllCellTypes.pdf
#     - sample_seurat4.rds          : Saved Seurat object
#     - <cluster>_cluster_markers.csv : Marker genes per cluster
#
#   Also outputs:
#     - sample_seurat4.rds (working directory)
#
# TODO:
#   - parameter testing
#     - features (20K)
#     - logfc threshold (TODO: try lower values)
#   - speed up clustering (issue is that it checks all features for 
#     clustering, not just those found by nfeatures)
#     
# ============================================================

library(Seurat)
if (requireNamespace("harmony", quietly = TRUE)) library(harmony)  # absent in R/4.2.2-BIO; RunHarmony is commented out below
library(plyr)
library(dplyr)
library(Matrix)
require("ggrepel")

rm(list=ls())

# RUN COPY (not the repo version): paths come from the environment so a run can
# read the real inputs without writing over previously saved results.
data_dir  <- Sys.getenv("KMER_DATA_DIR",  "/u/scratch/m/mhuang/scripts/bkc_counts")
plots_dir <- Sys.getenv("KMER_PLOTS_DIR", file.path(data_dir, "plots"))
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}
cat("data_dir  =", data_dir, "\nplots_dir =", plots_dir, "\n")

#______________________________________________________________________________

# For STANDARD PIPELINE:
# h5_file <- file.path(data_dir, "cellranger_matrix.h5")
# 
# kmer_mat <- Read10X_h5(h5_file)
# if (is.list(kmer_mat)) {
#   kmer_mat <- kmer_mat[[1]]
# }
# if (!dir.exists("plots_STD")) {
#   dir.create("plots_STD")
# }

#For KMER PIPELINE:
# Read sparse matrix
kmer_mat <- readMM(file.path(data_dir, "kmer_matrix.mtx"))
kmer_ids <- readLines(file.path(data_dir, "kmer_ids.txt"))
cell_barcodes <- readLines(file.path(data_dir, "cell_barcodes.txt"))
rownames(kmer_mat) <- kmer_ids
colnames(kmer_mat) <- cell_barcodes

#______________________________________________________________________________

# Create Seurat object (don't overwrite the function name!)
sample <- CreateSeuratObject(counts = kmer_mat, project = "KmerCounts", min.cells = 1, min.features = 100)

# Add in metadata
# sample$donor = as.character(sapply(as.character(sample$orig.ident), function(x){
#   unlist(strsplit(x,'_'))[1]
# }))
# sample$cycle = as.character(sapply(as.character(sample$orig.ident), function(x){
#   unlist(strsplit(x,'_'))[2]
# }))

# Count Mito genes and make plots
# BUGFIX: on the k-mer branch rownames are k-mer sequences, not gene symbols, so
# "^MT-" never matches and percent.mito came out identically 0 for every cell --
# a real-looking QC column that measured nothing. Compute it only when MT- features
# actually exist; otherwise record NA and drop it from the QC panels.
mito.features = grep(pattern="^MT-", x=rownames(x=sample), value=TRUE)
if (length(mito.features) > 0) {
  percent.mito = Matrix::colSums(x=GetAssayData(object=sample, slot="counts")[mito.features, , drop=FALSE]) /
    Matrix::colSums(x=GetAssayData(object=sample, slot="counts"))
  sample[["percent.mito"]] = percent.mito
} else {
  warning("No ^MT- features found -- expected on the k-mer branch (rownames are k-mer ",
          "sequences). percent.mito set to NA and omitted from QC plots.")
  sample[["percent.mito"]] = NA_real_
}
qc_features = c("nFeature_RNA", "nCount_RNA")
if (!all(is.na(sample$percent.mito))) qc_features = c(qc_features, "percent.mito")

VlnPlot(object=sample, features=qc_features,
        group.by="orig.ident", ncol=length(qc_features), pt.size=0.4)
ggsave(file.path(plots_dir, "quality_by_origin.pdf"), width=30, height=9)
if ("percent.mito" %in% qc_features) {
  FeatureScatter(object=sample, feature1="nCount_RNA", feature2="percent.mito")
}
FeatureScatter(object=sample, feature1="nCount_RNA", feature2="nFeature_RNA")

# STD: Filter cells and Normalize data
# sample = subset(x=sample, subset=nCount_RNA <= 5e4 & percent.mito <= 0.2)
# sample = NormalizeData(object = sample, normalization.method = "LogNormalize", scale.factor = 7000)

# KMER: Filter cells and Normalize data
sample = subset(x=sample, subset=nCount_RNA <= 1.7e6)
sample = NormalizeData(object = sample, normalization.method = "LogNormalize",
                       scale.factor = 150000)

gc()

# Find variable Genes and scale data by number of UMIs and Mito gene percentage
sample = FindVariableFeatures(object=sample, selection.method="vst", nfeatures = 20000)
VariableFeatures(sample)
sample = ScaleData(object=sample)

gc()

# Run PCA
sample = RunPCA(object=sample, verbose=FALSE)
DimPlot(sample, reduction="pca", group.by="ident")
ElbowPlot(object=sample, ndims=50)

# Run Harmony
# pca_dims = 1:30
# sample = RunHarmony(sample, group.by.vars="orig.ident", dims.use=pca_dims)
# ElbowPlot(object=sample, ndims=50, reduction="harmony")

# Run UMAP and clustering
set.seed(7)
dims_use = 1:30
sample = RunUMAP(object=sample, reduction="pca", dims=dims_use, verbose=FALSE)
sample = FindNeighbors(object=sample, reduction="pca", dims=dims_use, verbose=FALSE)
sample = FindClusters(object=sample, resolution=0.2, verbose=FALSE)

# Make DM plots
set.seed(7)
DimPlot(sample, label=TRUE, reduction="umap", pt.size=0.1, shuffle=TRUE, raster = F)
ggsave(file.path(plots_dir, "umap_cluster.pdf"), width=7, height=6)
DimPlot(sample, label=FALSE, reduction="umap", group.by="orig.ident", pt.size=0.1, shuffle=TRUE, raster = F) + labs(title=NULL)
ggsave(file.path(plots_dir, "umap_origin.pdf"), width=7, height=6)

# Save the Seurat object
saveRDS(file="sample_seurat4.rds", sample)
saveRDS(file = file.path(plots_dir, "sample_seurat4.rds"), sample)
# sample = readRDS("sample_seurat4.rds")

# Remove object from memory
# BUGFIX: this was rm(list=ls()), which also deleted data_dir/plots_dir; the reload
# below then pointed at a hardcoded *different* dataset (output_files_28), so the
# object read back was not the one just saved. Keep the paths across the cleanup.
rm(list = setdiff(ls(), c("data_dir", "plots_dir")))
gc()  # optional: force garbage collection

# Reload Seurat object from plots folder (same data_dir the save block wrote to)
sample <- readRDS(file.path(plots_dir, "sample_seurat4.rds"))

# Quality by cluster
qc_features = c("nFeature_RNA", "nCount_RNA")
if ("percent.mito" %in% colnames(sample[[]]) && !all(is.na(sample$percent.mito))) {
  qc_features = c(qc_features, "percent.mito")
}
VlnPlot(object=sample, features=qc_features,
        ncol=length(qc_features), pt.size=0)
ggsave(file.path(plots_dir, "quality_by_cluster.pdf"), width=30, height=9, limitsize=FALSE)

# Find cluster markers
library(Seurat)
library(future)
plan(strategy="multicore", workers=16)
options(future.globals.maxSize=6.4e10)

#cluster_markers = FindAllMarkers(object=sample, only.pos=TRUE, logfc.threshold=0.25)
#write.csv(cluster_markers, "cluster_markers.csv")

clusters <- levels(sample$seurat_clusters)
marker_list <- list()
for (cl in clusters) {
  
  message("Finding markers for cluster ", cl)
  
  markers <- FindMarkers(
    object = sample,
    ident.1 = cl,
    only.pos = TRUE,
    logfc.threshold = 0.25
  )
  
  # Add gene column explicitly (since rownames contain gene names)
  markers$gene <- rownames(markers)
  # BUGFIX: FindMarkers() output carries no cluster label, so the combined table
  # had nothing for group_by(cluster) downstream. Stamp it here.
  markers$cluster <- cl
  marker_list[[cl]] <- markers
  
  # Save file
  filename <- paste0(cl, "_cluster_markers.csv")
  write.csv(markers, file.path(plots_dir, filename), row.names = FALSE)
}

# BUGFIX: the combined cluster_markers.csv was only ever written by the commented-out
# FindAllMarkers() call above, and was read from the working directory rather than
# plots_dir -- so this block failed on a missing file. Build it from the loop instead.
cluster_markers <- do.call(rbind, marker_list)
write.csv(cluster_markers, file.path(plots_dir, "cluster_markers.csv"), row.names = FALSE)
cluster_markers = read.csv(file.path(plots_dir, "cluster_markers.csv"))
top = cluster_markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
genes = unique(as.character(top$gene))
DotPlot(sample, features=rev(genes), cols=c("grey90","red")) + 
  coord_flip() + labs(x=NULL, y=NULL)
ggsave(file.path(plots_dir, "dotplot_cluster_markers.pdf"), width=12, height=12)

# Feature plots and Violin plots
#FeaturePlot(sample, features=c("MLLT3-L","MLLT3-S","MLLT3-C"), cols=c("grey90","red"))

# change label
annotation = read.csv('annotation.csv')
sample$celltype = plyr::mapvalues(sample$seurat_clusters, from = annotation$cluster, to=annotation$annotation)
set.seed(1)
DimPlot(sample, label=TRUE, reduction="umap",  group.by="celltype",pt.size=0.1, shuffle=TRUE, raster = F)
ggsave(file.path(plots_dir, "umap_celltype.pdf"), width=7, height=6)

# Save the Seurat object
saveRDS(file="sample_seurat4.rds", sample)


#______________________________________


# ACCESSORY FUNCTIONS

#=======================================
# Azimuth references
#=======================================

azimuth_markers <- list(
  "B intermediate" = c("MS4A1","TNFRSF13B","IGHM","IGHD","AIM2","CD79A","LINC01857","RALGPS2","BANK1","CD79B"),
  "B memory" = c("MS4A1","COCH","AIM2","BANK1","SSPN","CD79A","TEX9","RALGPS2","TNFRSF13C","LINC01781"),
  "B naive" = c("IGHM","IGHD","CD79A","IL4R","MS4A1","CXCR4","BTG1","TCL1A","CD79B","YBX3"),
  "Plasmablast" = c("IGHA2","MZB1","TNFRSF17","DERL3","TXNDC5","TNFRSF13B","POU2AF1","CPNE5","HRASLS2","NT5DC2"),
  "CD4 CTL" = c("GZMH","CD4","FGFBP2","ITGB1","GZMA","CST7","GNLY","B2M","IL32","NKG7"),
  "CD4 Naive" = c("TCF7","CD4","CCR7","IL7R","FHIT","LEF1","MAL","NOSIP","LDHB","PIK3IP1"),
  "CD4 Proliferating" = c("MKI67","TOP2A","PCLAF","CENPF","TYMS","NUSAP1","ASPM","PTTG1","TPX2","RRM2"),
  "CD4 TCM" = c("IL7R","TMSB10","CD4","ITGB1","LTB","TRAC","AQP3","LDHB","IL32","MAL"),
  "CD4 TEM" = c("IL7R","CCL5","FYB1","GZMK","IL32","GZMA","KLRB1","TRAC","LTB","AQP3"),
  "Treg" = c("RTKN2","FOXP3","AC133644.2","CD4","IL2RA","TIGIT","CTLA4","FCRL3","LAIR2","IKZF2"),
  "CD8 Naive" = c("CD8B","S100B","CCR7","RGS10","NOSIP","LINC02446","LEF1","CRTAM","CD8A","OXNAD1"),
  "CD8 Proliferating" = c("MKI67","CD8B","TYMS","TRAC","PCLAF","CD3D","CLSPN","CD3G","TK1","RRM2"),
  "CD8 TCM" = c("CD8B","ANXA1","CD8A","KRT1","LINC02446","YBX3","IL7R","TRAC","NELL2","LDHB"),
  "CD8 TEM" = c("CCL5","GZMH","CD8A","TRAC","KLRD1","NKG7","GZMK","CST7","CD8B","TRGC2"),
  "ASDC" = c("PPP1R14A","LILRA4","AXL","IL3RA","SCT","SCN9A","LGMN","DNASE1L3","CLEC4C","GAS6"),
  "cDC1" = c("CLEC9A","DNASE1L3","C1orf54","IDO1","CLNK","CADM1","FLT3","ENPP1","XCR1","NDRG2"),
  "cDC2" = c("FCER1A","HLA-DQA1","CLEC10A","CD1C","ENHO","PLD4","GSN","SLC38A1","NDRG2","AFF3"),
  "pDC" = c("ITM2C","PLD4","SERPINF1","LILRA4","IL3RA","TPM2","MZB1","SPIB","IRF4","SMPD3"),
  "CD14 Mono" = c("S100A9","CTSS","S100A8","LYZ","VCAN","S100A12","IL1B","CD14","G0S2","FCN1"),
  "CD16 Mono" = c("CDKN1C","FCGR3A","PTPRC","LST1","IER5","MS4A7","RHOC","IFITM3","AIF1","HES4"),
  "NK" = c("GNLY","TYROBP","NKG7","FCER1G","GZMB","TRDC","PRF1","FGFBP2","SPON2","KLRF1"),
  "NK Proliferating" = c("MKI67","KLRF1","TYMS","TRDC","TOP2A","FCER1G","PCLAF","CD247","CLSPN","ASPM"),
  "NK_CD56bright" = c("XCL2","FCER1G","SPINK2","TRDC","KLRC1","XCL1","SPTSSB","PPP1R9A","NCAM1","TNFRSF11A"),
  "Eryth" = c("HBD","HBM","AHSP","ALAS2","CA1","SLC4A1","IFIT1B","TRIM58","SELENBP1","TMCC2"),
  "HSPC" = c("SPINK2","PRSS57","CYTL1","EGFL7","GATA2","CD34","SMIM24","AVP","MYB","LAPTM4B"),
  "ILC" = c("KIT","TRDC","TTLL10","LINC01229","SOX4","KLRB1","TNFRSF18","TNFRSF4","IL1R1","HPGDS"),
  "Platelet" = c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","HIST1H2AC","RGS18","GP9"),
  "dnT" = c("PTPN3","MIR4422HG","NUCB2","CAV1","DTHD1","GZMA","MYB","FXYD2","GZMK","AC004585.1"),
  "gdT" = c("TRDC","TRGC1","TRGC2","KLRC1","NKG7","TRDV2","CD7","TRGV9","KLRD1","KLRG1"),
  "MAIT" = c("KLRB1","NKG7","GZMK","IL7R","SLC4A10","GZMA","CXCR6","PRSS35","RBM24","NCR3")
)

library(ggplot2)
library(patchwork)

pdf(file.path(plots_dir, "Azimuth_FeaturePlots_AllCellTypes.pdf"), width = 20, height = 10)

for (ct in names(azimuth_markers)) {
  
  genes <- azimuth_markers[[ct]]
  
  message("Plotting: ", ct)
  
  p <- FeaturePlot(sample, features = genes)
  
  p_final <- p +
    plot_annotation(
      title = ct,
      theme = theme(
        plot.title = element_text(size = 26, face = "bold", hjust = 0.5)
      )
    )
  
  print(p_final)
}

dev.off()

#---------------------------------------------
# UMAP plot colored by cell type
#---------------------------------------------

set.seed(7)
DimPlot(sample, label=TRUE, reduction="umap", group.by="celltype",
        pt.size=0.1, shuffle=TRUE, raster = FALSE)
ggsave(file.path(plots_dir, "umap_celltype_STD.pdf"), width=7, height=6)

#KMER
library(ggplot2)

annot <- read.csv("barcode_celltype_STD.csv", stringsAsFactors = FALSE)
rownames(annot) <- annot$Barcode

sample$celltype_from_csv <- annot[colnames(sample), "CellType"]

p <- DimPlot(sample, reduction = "umap", group.by = "celltype_from_csv",
             label = TRUE, pt.size = 0.1, shuffle = TRUE)
p

ggsave(file.path(plots_dir, "umap_celltype_from_csv.pdf"), p, width = 7, height = 6)