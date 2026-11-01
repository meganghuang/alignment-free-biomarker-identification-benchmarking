library(Seurat)
library(harmony)
library(plyr)
library(dplyr)
library(Matrix)
require("ggrepel")

rm(list=ls())

data_dir <- "/u/scratch/m/mhuang/std_scripts/counts"

# create plots folder inside data_dir
plots_dir <- file.path(data_dir, "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

#______________________________________________________________________________

#For STANDARD PIPELINE:
h5_file <- file.path(data_dir, "5k_Human_Donor4_PBMC_3p_gem-x_5k_Human_Donor4_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.h5")

kmer_mat <- Read10X_h5(h5_file)
if (is.list(kmer_mat)) {
  kmer_mat <- kmer_mat[[1]]
}
if (!dir.exists("plots")) {
  dir.create("plots")
}

#For KMER PIPELINE:
# Read sparse matrix
#kmer_mat <- readMM(file.path(data_dir, "kmer_matrix.mtx"))
#kmer_ids <- readLines(file.path(data_dir, "kmer_ids.txt"))
#cell_barcodes <- readLines(file.path(data_dir, "cell_barcodes.txt"))
#rownames(kmer_mat) <- kmer_ids
#colnames(kmer_mat) <- cell_barcodes

#______________________________________________________________________________

# Create Seurat object (don't overwrite the function name!)
sample <- CreateSeuratObject(counts = kmer_mat, project = "StdCounts", min.cells = 1, min.features = 100)

# Add in metadata
# sample$donor = as.character(sapply(as.character(sample$orig.ident), function(x){
#   unlist(strsplit(x,'_'))[1]
# }))
# sample$cycle = as.character(sapply(as.character(sample$orig.ident), function(x){
#   unlist(strsplit(x,'_'))[2]
# }))

# Count Mito genes and make plots
mito.features = grep(pattern="^MT-", x=rownames(x=sample), value=TRUE)
percent.mito = Matrix::colSums(x=GetAssayData(object=sample, slot="counts")[mito.features,]) /
  Matrix::colSums(x=GetAssayData(object=sample, slot="counts"))
sample[["percent.mito"]] = percent.mito
VlnPlot(object=sample, features=c("nFeature_RNA","nCount_RNA","percent.mito"),
        group.by="orig.ident", ncol=3, pt.size=0.4)
ggsave(file.path(plots_dir, "quality_by_origin.pdf"), width=30, height=9)
FeatureScatter(object=sample, feature1="nCount_RNA", feature2="percent.mito")
FeatureScatter(object=sample, feature1="nCount_RNA", feature2="nFeature_RNA")

# STD: Filter cells and Normalize data
sample = subset(x=sample, subset=nCount_RNA <= 5e4 & percent.mito <= 0.2)
sample = NormalizeData(object = sample, normalization.method = "LogNormalize", scale.factor = 7000)

# KMER: Filter cells and Normalize data
#sample = subset(x=sample, subset=nCount_RNA <= 1.7e6)
#sample = NormalizeData(object = sample, normalization.method = "LogNormalize",
                       scale.factor = 1.5e5)

gc()

# Find variable Genes and scale data by number of UMIs and Mito gene percentage
sample = FindVariableFeatures(object=sample, selection.method="vst", nfeatures = 2000)
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
sample = FindClusters(object=sample, resolution=0.4, verbose=FALSE)

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



# Quality by cluster
VlnPlot(object=sample, features=c("nFeature_RNA","nCount_RNA","percent.mito"), 
        ncol=3, pt.size=0)
ggsave(file.path(plots_dir, "quality_by_cluster.pdf"), width=30, height=9, limitsize=FALSE)

# Remove object from memory
rm(list=ls())
gc()  # optional: force garbage collection

# Reload Seurat object from plots folder
data_dir <- "/u/scratch/m/mhuang/std_scripts/counts"
plots_dir <- file.path(data_dir, "plots")
sample <- readRDS(file.path(plots_dir, "sample_seurat4.rds"))

# Find cluster markers
library(Seurat)
library(future)
plan(strategy="multicore", workers=16)
options(future.globals.maxSize=6.4e10)

#cluster_markers = FindAllMarkers(object=sample, only.pos=TRUE, logfc.threshold=0.25)
#write.csv(cluster_markers, "cluster_markers.csv")

# GET CLUSTERS
clusters <- levels(sample$seurat_clusters)
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
  
  # Save file
  filename <- paste0(cl, "_cluster_markers.csv")
  write.csv(markers, file.path(plots_dir, filename), row.names = FALSE)
}

#CLUSTER DOT PLOT
# Read all cluster marker files
files <- list.files(plots_dir, pattern="*_cluster_markers.csv", full.names=TRUE)

cluster_markers <- do.call(rbind, lapply(files, function(f) {
  df <- read.csv(f)
  df$cluster <- gsub("_cluster_markers.csv", "", basename(f))
  return(df)
}))

# Get top 3 genes per cluster
top <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3)

genes <- unique(as.character(top$gene))

# DotPlot
DotPlot(sample, features=rev(genes), cols=c("grey90","red")) +
  coord_flip() +
  labs(x=NULL, y=NULL)

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

#___________________________________________________________________________

# Reload Seurat object from plots folder
data_dir <- "/u/scratch/m/mhuang/std_scripts/counts"
plots_dir <- file.path(data_dir, "plots")
sample <- readRDS(file.path(plots_dir, "sample_seurat4.rds"))

# Read cluster-to-celltype annotation
annotation <- read.csv(
  file.path(data_dir, "annotation_STD.csv"),
  stringsAsFactors = FALSE
)

# Make sure columns are named correctly
colnames(annotation) <- c("Cluster", "CellType", "GeneralCellType")

# Build barcode → cluster table
df <- data.frame(
  Barcode = colnames(sample),
  Cluster = as.numeric(as.character(sample$seurat_clusters)),
  stringsAsFactors = FALSE
)

# Merge to assign annotations
df_annot <- merge(df, annotation, by = "Cluster", all.x = TRUE)

# Keep original cell order
df_annot <- df_annot[match(colnames(sample), df_annot$Barcode), ]

# Save barcode + detailed + general cell type
write.csv(
  df_annot[, c("Barcode", "CellType", "GeneralCellType")],
  file.path(data_dir, "barcode_celltype_STD.csv"),
  row.names = FALSE
)

# Save the Seurat object
saveRDS(sample, file = "sample_seurat4.rds")
#_____________________________________________________________________________

# Azimuth references

azimuth_markers <- list(
  "B intermediate" = c("MS4A1","TNFRSF13B","IGHM","IGHD","AIM2","CD79A","LINC01857","RALGPS2","BANK1","CD79B"),
  "B memory" = c("MS4A1","COCH","AIM2","BANK1","SSPN","CD79A","TEX9","RALGPS2","TNFRSF13C","LINC01781"),
  "B naive" = c("IGHM","IGHD","CD79A","IL4R","MS4A1","CXCR4","BTG1","TCL1A","CD79B","YBX3"),
  "Plasmablast" = c("IGHA2","MZB1","TNFRSF17","DERL3","TXNDC5","TNFRSF13B","POU2AF1","CPNE5","NT5DC2"),
  "CD4 CTL" = c("GZMH","CD4","FGFBP2","ITGB1","GZMA","CST7","GNLY","B2M","IL32","NKG7"),
  "CD4 Naive" = c("TCF7","CD4","CCR7","IL7R","FHIT","LEF1","MAL","NOSIP","LDHB","PIK3IP1"),
  "CD4 Proliferating" = c("MKI67","TOP2A","PCLAF","CENPF","TYMS","NUSAP1","ASPM","PTTG1","TPX2","RRM2"),
  "CD4 TCM" = c("IL7R","TMSB10","CD4","ITGB1","LTB","TRAC","AQP3","LDHB","IL32","MAL"),
  "CD4 TEM" = c("IL7R","CCL5","FYB1","GZMK","IL32","GZMA","KLRB1","TRAC","LTB","AQP3"),
  "Treg" = c("RTKN2","FOXP3","CD4","IL2RA","TIGIT","CTLA4","FCRL3","LAIR2","IKZF2"),
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
  "Platelet" = c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","RGS18","GP9"),
  "dnT" = c("PTPN3","MIR4422HG","NUCB2","CAV1","DTHD1","GZMA","MYB","FXYD2","GZMK"),
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

library(ggplot2)

# Read annotation file
annot <- read.csv(
  file.path(data_dir, "barcode_celltype_STD.csv"),
  stringsAsFactors = FALSE
)

# Set barcode rownames
rownames(annot) <- annot$Barcode

# Add annotations to Seurat object
sample$CellType <- annot[colnames(sample), "CellType"]
sample$GeneralCellType <- annot[colnames(sample), "GeneralCellType"]

# Optional: clean encoding issues
sample$CellType <- iconv(sample$CellType, from = "", to = "UTF-8", sub = "")
sample$GeneralCellType <- iconv(sample$GeneralCellType, from = "", to = "UTF-8", sub = "")

# ----- CellType UMAP -----
p1 <- DimPlot(
  sample,
  reduction = "umap",
  group.by = "CellType",
  label = TRUE,
  pt.size = 0.1,
  shuffle = TRUE
)

p1

ggsave(
  file.path(plots_dir, "umap_CellType.pdf"),
  p1,
  width = 7,
  height = 6
)

# ----- GeneralCellType UMAP -----
p2 <- DimPlot(
  sample,
  reduction = "umap",
  group.by = "GeneralCellType",
  label = TRUE,
  pt.size = 0.1,
  shuffle = TRUE
)

p2

ggsave(
  file.path(plots_dir, "umap_GeneralCellType.pdf"),
  p2,
  width = 7,
  height = 6
)