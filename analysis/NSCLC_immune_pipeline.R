suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(tidyr)
  library(readr)
  library(stringr)
  library(Matrix)
})

set.seed(1)

# 0) Paths
setwd("~/Desktop/10Xdata")
stopifnot(dir.exists(getwd()))
cat("Working directory:", getwd(), "\n")

# 0b) Make a unique run folder so outputs NEVER mix
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTDIR <- file.path(getwd(), paste0("RUN_", run_id))
dir.create(OUTDIR, showWarnings = FALSE)
cat("This run output folder:", OUTDIR, "\n")

# helper to save plots cleanly (optional)
save_plot <- function(p, name, w=10, h=6) {
  ggsave(filename = file.path(OUTDIR, name), plot = p, width = w, height = h)
  cat("Saved plot:", file.path(OUTDIR, name), "\n")
}

# 0c) Quick check input files
h5_files <- list.files(getwd(), pattern = "donor_[1-7]_count_sample_feature_bc_matrix\\.h5$", full.names = FALSE)
cat("Found donor h5 files:", length(h5_files), "\n")
print(h5_files)
stopifnot(length(h5_files) == 7)
# ============================
# STEP 1: Load donors + merge
# ============================

donor_list <- list()

for (i in 1:7) {
  cat("\n--- Loading Donor", i, "---\n")
  
  filename <- paste0(
    "20k_NSCLC_DTC_3p_nextgem_intron_donor_",
    i,
    "_count_sample_feature_bc_matrix.h5"
  )
  stopifnot(file.exists(filename))
  
  data <- Read10X_h5(filename)
  
  # RNA required
  stopifnot("Gene Expression" %in% names(data))
  rna_counts <- data$`Gene Expression`
  
  obj <- CreateSeuratObject(
    counts = rna_counts,
    project = "NSCLC",
    min.cells = 3,
    min.features = 200
  )
  
  # ADT optional (but expected in this dataset)
  if ("Antibody Capture" %in% names(data)) {
    adt_counts <- data$`Antibody Capture`
    common <- intersect(colnames(obj), colnames(adt_counts))
    cat("  RNA cells:", ncol(obj), " | ADT cells:", ncol(adt_counts),
        " | overlap:", length(common), "\n")
    
    if (length(common) > 0) {
      adt_subset <- adt_counts[, common, drop = FALSE]
      adt_subset <- adt_subset[, colnames(obj)[colnames(obj) %in% common], drop = FALSE]
      obj[["ADT"]] <- CreateAssayObject(counts = adt_subset)
    } else {
      warning("No overlapping barcodes between RNA and ADT for Donor ", i)
    }
  } else {
    warning("No Antibody Capture modality found for Donor ", i)
  }
  
  obj$donor <- paste0("Donor", i)
  donor_list[[paste0("Donor", i)]] <- obj
  
  cat("  Final cells kept:", ncol(obj), "\n")
}

cat("\nAll donors loaded. Donors in list:", length(donor_list), "\n")
print(sapply(donor_list, ncol))

# Merge donors
nsclc <- merge(
  x = donor_list[[1]],
  y = donor_list[2:7],
  add.cell.ids = names(donor_list),
  project = "NSCLC"
)

cat("\nMerged object cells:", ncol(nsclc), "\n")
cat("Donor distribution:\n")
print(table(nsclc$donor))

# Save checkpoint (inside this RUN folder)
saveRDS(nsclc, file = file.path(OUTDIR, "STEP1_nsclc_merged_raw.rds"))
cat("Saved checkpoint:", file.path(OUTDIR, "STEP1_nsclc_merged_raw.rds"), "\n")
# ============================
# STEP 2: QC + filtering
# ============================

# Load checkpoint (safety)
nsclc <- readRDS(file.path(OUTDIR, "STEP1_nsclc_merged_raw.rds"))

DefaultAssay(nsclc) <- "RNA"

# Join layers ONCE after merge (Seurat v5 best practice)
nsclc <- JoinLayers(nsclc)

# QC metrics
nsclc[["percent.mt"]]   <- PercentageFeatureSet(nsclc, pattern = "^MT-")
nsclc[["percent.ribo"]] <- PercentageFeatureSet(nsclc, pattern = "^RP[SL]")

cat("QC summary (all cells):\n")
print(summary(nsclc$nFeature_RNA))
print(summary(nsclc$nCount_RNA))
print(summary(nsclc$percent.mt))
print(summary(nsclc$percent.ribo))

# QC violin plot by donor
p_qc <- VlnPlot(
  nsclc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
  ncol = 4,
  group.by = "donor",
  pt.size = 0
) + plot_annotation(title = "QC metrics by donor")

print(p_qc)
save_plot(p_qc, "STEP2_QC_violin_by_donor.png", w=14, h=4.5)

# Apply your filter thresholds (same as before)
nsclc_filt <- subset(
  nsclc,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 7500 &
    percent.mt < 20
)

cat("\nCells before filtering:", ncol(nsclc), "\n")
cat("Cells after filtering :", ncol(nsclc_filt), "\n")
cat("Cells removed         :", ncol(nsclc) - ncol(nsclc_filt), "\n")

# Re-join layers after filtering (keeps v5 layers consistent)
nsclc_filt <- JoinLayers(nsclc_filt)

# Save checkpoint
saveRDS(nsclc_filt, file = file.path(OUTDIR, "STEP2_nsclc_filtered.rds"))
cat("Saved checkpoint:", file.path(OUTDIR, "STEP2_nsclc_filtered.rds"), "\n")
# ============================
# STEP 3: RNA preprocessing + PCA
# ============================

nsclc <- readRDS(file.path(OUTDIR, "STEP2_nsclc_filtered.rds"))
DefaultAssay(nsclc) <- "RNA"

# Normalize (creates "data" layer in Seurat v5)
nsclc <- NormalizeData(nsclc, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features
nsclc <- FindVariableFeatures(nsclc, selection.method = "vst", nfeatures = 2000)

# Scale (only HVGs)
nsclc <- ScaleData(nsclc, features = VariableFeatures(nsclc))

# PCA
nsclc <- RunPCA(nsclc, features = VariableFeatures(nsclc), npcs = 50)

# Plots
p_elbow <- ElbowPlot(nsclc, ndims = 50) + ggtitle("RNA PCA elbow (filtered)")
print(p_elbow)
save_plot(p_elbow, "STEP3_RNA_PCA_elbow.png", w=7, h=5)

# Save checkpoint
saveRDS(nsclc, file = file.path(OUTDIR, "STEP3_nsclc_RNA_PCA.rds"))
cat("Saved checkpoint:", file.path(OUTDIR, "STEP3_nsclc_RNA_PCA.rds"), "\n")

# Quick sanity
cat("Top HVGs:", paste(head(VariableFeatures(nsclc), 10), collapse = ", "), "\n")
cat("PCA dims:", ncol(Embeddings(nsclc, 'pca')), "\n")
# ============================
# STEP 4: ADT preprocessing + APCA
# ============================

nsclc <- readRDS(file.path(OUTDIR, "STEP3_nsclc_RNA_PCA.rds"))
stopifnot("ADT" %in% names(nsclc@assays))

DefaultAssay(nsclc) <- "ADT"

# Use all proteins as "features"
VariableFeatures(nsclc) <- rownames(nsclc[["ADT"]])

# CLR normalize across cells
nsclc <- NormalizeData(nsclc, normalization.method = "CLR", margin = 2)

# Scale ADT
nsclc <- ScaleData(nsclc, features = VariableFeatures(nsclc))

# APCA (keep small because ADT has few proteins)
adt_npcs <- min(5, nrow(nsclc[["ADT"]]) - 1)
nsclc <- RunPCA(
  nsclc,
  features = VariableFeatures(nsclc),
  npcs = adt_npcs,
  reduction.name = "apca",
  reduction.key = "APCA_"
)

# Quick check + save
print(nsclc[["apca"]], dims = 1:3, nfeatures = 5)

saveRDS(nsclc, file = file.path(OUTDIR, "STEP4_nsclc_ADT_APCA.rds"))
cat("Saved checkpoint:", file.path(OUTDIR, "STEP4_nsclc_ADT_APCA.rds"), "\n")
cat("ADT npcs used:", adt_npcs, "\n")
cat("ADT features:", nrow(nsclc[["ADT"]]), "\n")
# ============================
# STEP 5: WNN (RNA + ADT) → UMAP → clusters
# ============================

nsclc <- readRDS(file.path(OUTDIR, "STEP4_nsclc_ADT_APCA.rds"))

# CRITICAL: switch back to RNA before WNN
DefaultAssay(nsclc) <- "RNA"

# Use fixed dimensions (stable + standard)
rna_dims <- 1:30
adt_dims <- 1:5

# Build WNN graph
nsclc <- FindMultiModalNeighbors(
  nsclc,
  reduction.list = list("pca", "apca"),
  dims.list = list(rna_dims, adt_dims)
)

# UMAP on weighted nearest neighbors
nsclc <- RunUMAP(
  nsclc,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_",
  seed.use = 1
)

# Clustering on wsnn graph
nsclc <- FindClusters(
  nsclc,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = 0.5,
  verbose = FALSE
)

# Sanity + plots
cat("Cells:", ncol(nsclc), "\n")
cat("Clusters:", length(unique(nsclc$seurat_clusters)), "\n")
print(table(nsclc$seurat_clusters))

p1 <- DimPlot(nsclc, reduction = "wnn.umap", label = TRUE, repel = TRUE) +
  ggtitle("WNN clusters (RNA + ADT)")
p2 <- DimPlot(nsclc, reduction = "wnn.umap", group.by = "donor") +
  ggtitle("Cells by donor")
print(p1 + p2)

# Save checkpoint
saveRDS(nsclc, file = file.path(OUTDIR, "STEP5_nsclc_WNN.rds"))
cat("Saved checkpoint:", file.path(OUTDIR, "STEP5_nsclc_WNN.rds"), "\n")
# ============================
# STEP 6: Immune vs non-immune marker scan + immune subset
# ============================

nsclc <- readRDS(file.path(OUTDIR, "STEP5_nsclc_WNN.rds"))
DefaultAssay(nsclc) <- "RNA"

# Quick marker panel (broad, not over-specific)
immune_markers   <- c("PTPRC","CD3D","CD3E","TRAC","CD8A","CD4","NKG7","GNLY","MS4A1","CD79A","LYZ","S100A8","FCGR3A")
epithelial_markers <- c("EPCAM","KRT8","KRT18","KRT19","ALDH1A1","MSLN")
endothelial_markers <- c("PECAM1","VWF","KDR")
fibroblast_markers  <- c("COL1A1","COL1A2","DCN","LUM","COL3A1")

p_dot <- DotPlot(
  nsclc,
  features = c(immune_markers, epithelial_markers, endothelial_markers, fibroblast_markers),
  group.by = "seurat_clusters"
) + RotatedAxis() + ggtitle("Step 6: Marker scan (immune / epi / endo / fibro)")

print(p_dot)
save_plot(p_dot, "STEP6_marker_dotplot.png", w=16, h=6)

# OPTIONAL: also a quick FeaturePlot sanity for PTPRC and EPCAM
p_fp <- FeaturePlot(nsclc, features = c("PTPRC","EPCAM"), reduction = "wnn.umap", ncol = 2)
print(p_fp)
save_plot(p_fp, "feature_PTPRC_EPCAM.png", w=10, h=4)

nsclc <- readRDS(file.path(OUTDIR, "STEP5_nsclc_WNN.rds"))
DefaultAssay(nsclc) <- "RNA"

# Immune clusters from Step 6 dotplot
immune_clusters <- c(0, 1, 2, 3, 4, 8, 11, 12, 15, 19, 21)

# Safety checks
immune_clusters <- intersect(immune_clusters, levels(nsclc$seurat_clusters))
stopifnot(length(immune_clusters) > 0)

immune <- subset(nsclc, seurat_clusters %in% immune_clusters)

cat("Immune clusters:", paste(immune_clusters, collapse = ", "), "\n")
cat("Immune cells extracted:", ncol(immune), "out of", ncol(nsclc), "\n")

saveRDS(immune, file = file.path(OUTDIR, "STEP6_immune_rawsubset.rds"))
cat("Saved checkpoint:", file.path(OUTDIR, "STEP6_immune_rawsubset.rds"), "\n")
# ============================
# STEP 7) Immune reprocessing (pre-doublet) with WNN
# ============================

immune <- readRDS(file.path(OUTDIR, "STEP6_immune_rawsubset.rds"))

# keep Seurat v5 layers consistent
immune <- JoinLayers(immune)

# ---- RNA side ----
DefaultAssay(immune) <- "RNA"
immune <- NormalizeData(immune)
immune <- FindVariableFeatures(immune, nfeatures = 2000)
immune <- ScaleData(immune, features = VariableFeatures(immune))
immune <- RunPCA(immune, npcs = 30)

# ---- ADT side ----
stopifnot("ADT" %in% names(immune@assays))
DefaultAssay(immune) <- "ADT"
VariableFeatures(immune) <- rownames(immune[["ADT"]])
immune <- NormalizeData(immune, normalization.method = "CLR", margin = 2)
immune <- ScaleData(immune, features = VariableFeatures(immune))

adt_npcs <- min(5, nrow(immune[["ADT"]]) - 1)
immune <- RunPCA(
  immune,
  features = VariableFeatures(immune),
  npcs = adt_npcs,
  reduction.name = "apca",
  reduction.key  = "APCA_"
)

# ---- WNN ----
DefaultAssay(immune) <- "RNA"   # IMPORTANT
immune <- FindMultiModalNeighbors(
  immune,
  reduction.list = list("pca", "apca"),
  dims.list      = list(1:30, 1:adt_npcs)
)

immune <- RunUMAP(
  immune,
  nn.name        = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key  = "wnnUMAP_",
  seed.use       = 1
)

immune <- FindClusters(
  immune,
  graph.name = "wsnn",
  algorithm  = 3,
  resolution = 0.8,
  verbose    = FALSE
)

# ---- plot + checkpoint ----
p_immune <- DimPlot(immune, reduction = "wnn.umap", label = TRUE, repel = TRUE) +
  ggtitle("Immune subclusters (pre-doublet)")
print(p_immune)
save_plot(p_immune, "STEP7_immune_preDoublet_umap.png", w=10, h=6)

saveRDS(immune, file = file.path(OUTDIR, "STEP7_immune_preDoublet.rds"))
cat("Saved checkpoint:", file.path(OUTDIR, "STEP7_immune_preDoublet.rds"), "\n")

cat("Immune cells:", ncol(immune), "\n")
cat("Immune clusters:", length(unique(immune$seurat_clusters)), "\n")
print(table(immune$seurat_clusters))
# ============================
# STEP 8) Doublet detection (immune-only) with scDblFinder
# ============================

immune <- readRDS(file.path(OUTDIR, "STEP7_immune_preDoublet.rds"))

# Install/load scDblFinder if needed
if (!requireNamespace("scDblFinder", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("scDblFinder")
}
suppressPackageStartupMessages({
  library(scDblFinder)
  library(SingleCellExperiment)
})

# Convert Seurat -> SingleCellExperiment (use RNA assay)
sce <- as.SingleCellExperiment(immune, assay = "RNA")

set.seed(1)
sce <- scDblFinder(sce)

# Bring results back into Seurat metadata
immune$doublet_score <- sce$scDblFinder.score
immune$doublet_class <- sce$scDblFinder.class

# Summary
cat("\nDoublet class counts:\n")
print(table(immune$doublet_class))

cat("\nDoublet % per immune cluster:\n")
dbl_tab <- table(immune$seurat_clusters, immune$doublet_class)
dbl_pct <- round(prop.table(dbl_tab, margin = 1) * 100, 2)
print(dbl_pct)

# Subset singlets
immune_singlets <- subset(immune, subset = doublet_class == "singlet")

cat("\nImmune cells before doublet removal:", ncol(immune), "\n")
cat("Immune cells after  doublet removal:", ncol(immune_singlets), "\n")
cat("Doublets removed:", ncol(immune) - ncol(immune_singlets), "\n")

# Save checkpoints
saveRDS(immune, file = file.path(OUTDIR, "STEP8_immune_withDoublets.rds"))
saveRDS(immune_singlets, file = file.path(OUTDIR, "STEP8_immune_singlets_raw.rds"))
cat("Saved:\n  -", file.path(OUTDIR, "STEP8_immune_withDoublets.rds"),
    "\n  -", file.path(OUTDIR, "STEP8_immune_singlets_raw.rds"), "\n")

# Quick plot (optional)
p_dbl <- DimPlot(immune, reduction = "wnn.umap", group.by = "doublet_class") +
  ggtitle("Immune doublets (scDblFinder)")
print(p_dbl)
save_plot(p_dbl, "STEP8_immune_doublet_class_umap.png", w=10, h=6)
# ============================
# STEP 9) Reprocess immune singlets (RNA + ADT WNN)
# ============================

immune_singlets <- readRDS(file.path(OUTDIR, "STEP8_immune_singlets_raw.rds"))

# keep Seurat v5 layers consistent
immune_singlets <- JoinLayers(immune_singlets)

# ---- RNA preprocessing ----
DefaultAssay(immune_singlets) <- "RNA"
immune_singlets <- NormalizeData(immune_singlets)
immune_singlets <- FindVariableFeatures(immune_singlets, nfeatures = 2000)
immune_singlets <- ScaleData(immune_singlets, features = VariableFeatures(immune_singlets))
immune_singlets <- RunPCA(immune_singlets, npcs = 30)

# ---- ADT preprocessing ----
stopifnot("ADT" %in% names(immune_singlets@assays))
DefaultAssay(immune_singlets) <- "ADT"
VariableFeatures(immune_singlets) <- rownames(immune_singlets[["ADT"]])
immune_singlets <- NormalizeData(immune_singlets, normalization.method = "CLR", margin = 2)
immune_singlets <- ScaleData(immune_singlets, features = VariableFeatures(immune_singlets))

adt_npcs <- min(5, nrow(immune_singlets[["ADT"]]) - 1)
immune_singlets <- RunPCA(
  immune_singlets,
  features = VariableFeatures(immune_singlets),
  npcs = adt_npcs,
  reduction.name = "apca",
  reduction.key  = "APCA_"
)

# ---- WNN + UMAP + clustering ----
DefaultAssay(immune_singlets) <- "RNA"  # IMPORTANT

immune_singlets <- FindMultiModalNeighbors(
  immune_singlets,
  reduction.list = list("pca", "apca"),
  dims.list      = list(1:30, 1:adt_npcs)
)

immune_singlets <- RunUMAP(
  immune_singlets,
  nn.name        = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key  = "wnnUMAP_",
  seed.use       = 1
)

immune_singlets <- FindClusters(
  immune_singlets,
  graph.name = "wsnn",
  algorithm  = 3,
  resolution = 0.8,
  verbose    = FALSE
)

# ---- outputs ----
cat("Immune singlets cells:", ncol(immune_singlets), "\n")
cat("Immune singlets clusters:", length(unique(immune_singlets$seurat_clusters)), "\n")
print(table(immune_singlets$seurat_clusters))

p_clean <- DimPlot(immune_singlets, reduction = "wnn.umap", label = TRUE, repel = TRUE) +
  ggtitle("Immune singlets (clean WNN)")
print(p_clean)
save_plot(p_clean, "STEP9_immune_singlets_clean_umap.png", w=10, h=6)

saveRDS(immune_singlets, file = file.path(OUTDIR, "STEP9_immune_singlets_clean.rds"))
cat("Saved checkpoint:", file.path(OUTDIR, "STEP9_immune_singlets_clean.rds"), "\n")
# ============================
# STEP 10) Immune annotation markers (RNA + ADT validation)
# ============================

immune_singlets <- readRDS(file.path(OUTDIR, "STEP9_immune_singlets_clean.rds"))
immune_singlets <- JoinLayers(immune_singlets)
Idents(immune_singlets) <- "seurat_clusters"

# ---- Marker panels ----
# T cells
tcell <- c("CD3D","CD3E","TRAC","IL7R","CCR7","LTB","MAL","TCF7")
cd8   <- c("CD8A","CD8B","NKG7","GZMK","GZMB","PRF1","GNLY")
treg  <- c("FOXP3","IL2RA","CTLA4","IKZF2","TIGIT")

# NK
nk    <- c("NKG7","GNLY","PRF1","GZMB","KLRD1","FCGR3A","TRDC")

# B / Plasma
bcell <- c("MS4A1","CD79A","CD74","HLA-DRA","CD37","CD19")
plasma<- c("MZB1","JCHAIN","XBP1","SDC1","IGHG1","IGKC")

# Myeloid / Mono / Mac
mono  <- c("LYZ","S100A8","S100A9","FCN1","CTSS","LGALS3")
fcgr3a_mono <- c("FCGR3A","MS4A7","LST1","IFITM3","CTSD")
macro <- c("APOE","C1QA","C1QB","C1QC","TREM2","MSR1","MRC1","VSIG4")

# Dendritic
dc1   <- c("FCER1A","CST3","CLEC10A")
dc2   <- c("ITGAX","FCER1A","CD1C","CLEC10A")
pdc   <- c("GZMB","IRF7","IL3RA","TCF4")

# ---- RNA DotPlot ----
DefaultAssay(immune_singlets) <- "RNA"
rna_markers <- list(
  Tcell = tcell,
  CD8_NK = cd8,
  Treg = treg,
  NK = nk,
  Bcell = bcell,
  Plasma = plasma,
  Mono = mono,
  FCGR3A_mono = fcgr3a_mono,
  Macrophage = macro,
  DC = c(dc1, dc2),
  pDC = pdc
)

p_rna_dot <- DotPlot(
  immune_singlets,
  features = unique(unlist(rna_markers)),
  group.by = "seurat_clusters"
) + RotatedAxis() + ggtitle("Immune RNA marker validation (DotPlot)")

print(p_rna_dot)
save_plot(p_rna_dot, "STEP10_immune_RNA_marker_dotplot.png", w=18, h=7)

# ---- ADT DotPlot ----
DefaultAssay(immune_singlets) <- "ADT"
adt_markers <- c("CD45","CD3","CD4.1","CD8","CD14.1","CD11c","CD16","CD56","CD19.1")

p_adt_dot <- DotPlot(
  immune_singlets,
  features = adt_markers,
  group.by = "seurat_clusters"
) + RotatedAxis() + ggtitle("Immune ADT marker validation (DotPlot)")

print(p_adt_dot)
save_plot(p_adt_dot, "STEP10_immune_ADT_marker_dotplot.png", w=14, h=5)

# ---- Quick UMAP FeaturePlots (RNA) ----
DefaultAssay(immune_singlets) <- "RNA"
p_fp <- FeaturePlot(
  immune_singlets,
  features = c("CD3D","CD8A","NKG7","MS4A1","LYZ","FCGR3A","APOE","FCER1A"),
  reduction = "wnn.umap",
  ncol = 4
)
print(p_fp)
save_plot(p_fp, "STEP10_immune_featureplots_RNA.png", w=14, h=8)

# Save checkpoint
saveRDS(immune_singlets, file = file.path(OUTDIR, "STEP10_immune_markers_checked.rds"))
cat("Saved checkpoint:", file.path(OUTDIR, "STEP10_immune_markers_checked.rds"), "\n")
# ----------------------------
# STEP 11) Immune cluster markers + labeling (EVIDENCE-BASED)
# ----------------------------
immune <- readRDS(file.path(OUTDIR, "STEP10_immune_markers_checked.rds"))
immune <- JoinLayers(immune)
DefaultAssay(immune) <- "RNA"
Idents(immune) <- "seurat_clusters"

# 11A) Find markers per immune cluster
markers_imm <- FindAllMarkers(
  immune,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Save markers table
write.csv(markers_imm, file.path(OUTDIR, "STEP11_immune_FindAllMarkers.csv"), row.names = FALSE)
cat("Saved:", file.path(OUTDIR, "STEP11_immune_FindAllMarkers.csv"), "\n")

# 11B) Print top markers per cluster (quick evidence view)
top10 <- markers_imm %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

print(top10 %>% select(cluster, gene, avg_log2FC, pct.1, pct.2))

# 11C) Create a draft label map (YOU will finalize after seeing top10)
cluster_to_label <- c(
  "0"  = "Inflammatory/Classical Monocytes (S100A8/S100A9/FCN1/VCAN+)",
  "1"  = "B cells (FCRL/BANK1/MS4A1/EBF1+; memory-like)",
  "2"  = "Cytotoxic CD8 T (CD8A/CD8B/NKG7/CCL5; LAG3+)",
  "3"  = "FCGR3A+ Mono/Mac (CD14/FCGR3A/CCL18/MARCO+)",
  "4"  = "Activated CD4 / Treg-like (TNFRSF4/TNFRSF18/CTLA4/ICOS/IL2RA+)",
  "5"  = "CXCL13+ CD4 T (CXCL13/TOX2+; TLS-like)",
  "6"  = "FOLR2+ Macrophages (FOLR2/SELENOP/APOE+)",
  "7"  = "IL7R+ CD4 T (IL7R/KLRB1/RORA; memory-like)",
  "8"  = "Activated/Class-switched B (IGHA1/TNFRSF13B/C+)",
  "9"  = "Cytotoxic CD8/NK-like T (CD8/IFNG/NKG7/CCL5+)",
  "10" = "Activated T (TNFRSF4/BATF/ICOS/LEF1+)",
  "11" = "NON-IMMUNE contamination (endothelial/stromal: HSPG2/COL4A1/CAV1+)",
  "12" = "MARCO+ Macrophages (MARCO/NR1H3/APOC2/MSR1+)",
  "13" = "B cells (BANK1/BLK/VPREB3+ subset)",
  "14" = "NON-IMMUNE contamination (epithelial/secretory: SLPI/WFDC2/AKR1C+)",
  "15" = "NON-IMMUNE contamination (neuroendocrine: PHOX2B/OTX2/NEUROD1+)",
  "16" = "NK cells (KIR/KLRC/NCAM1/GNLY+)",
  "17" = "pDC (SCT/CLEC4C/LILRA4+)",
  "18" = "FCGR3A+ MARCO+ myeloid (SIGLEC10/LILRA5+)",
  "19" = "γδ T cells (TRDC/TRGC2+)"
)
missing <- setdiff(levels(immune$seurat_clusters), names(cluster_to_label))
cat("Clusters missing labels:", paste(missing, collapse=", "), "\n")

# assign a placeholder label instead of crashing
if (length(missing) > 0) {
  cluster_to_label[missing] <- paste0("Unknown cluster ", missing, " (inspect markers)")
}

immune$cell_type <- unname(cluster_to_label[as.character(immune$seurat_clusters)])

# Plot labeled UMAP
p_lab <- DimPlot(immune, reduction="wnn.umap", group.by="cell_type", label=TRUE, repel=TRUE) +
  ggtitle("Immune cell types (labeled from markers)")
print(p_lab)
save_plot(p_lab, "STEP11_immune_labeled_umap.png", w=12, h=7)

# Save checkpoint
saveRDS(immune, file.path(OUTDIR, "STEP11_immune_labeled.rds"))
cat("Saved:", file.path(OUTDIR, "STEP11_immune_labeled.rds"), "\n")



# ----------------------------
# STEP 12) Remove non-immune contamination (10/14/15) + save clean immune object
# ----------------------------
immune <- readRDS(file.path(OUTDIR, "STEP11_immune_labeled.rds"))
immune <- JoinLayers(immune)
Idents(immune) <- "seurat_clusters"

# Remove these clusters (you labeled them as non-immune)
non_immune_clusters <- c("11","14","15")

cat("Cells BEFORE:", ncol(immune), "\n")
cat("Cluster sizes BEFORE:\n")
print(table(Idents(immune)))

immune_clean <- subset(immune, idents = setdiff(levels(Idents(immune)), non_immune_clusters))

cat("\nCells AFTER :", ncol(immune_clean), "\n")
cat("Cluster sizes AFTER:\n")
print(table(immune_clean$seurat_clusters))

cat("\nCell-type counts AFTER:\n")
print(sort(table(immune_clean$cell_type), decreasing = TRUE))

# Save checkpoint
saveRDS(immune_clean, file = file.path(OUTDIR, "STEP12_immune_clean_labeled.rds"))
cat("Saved:", file.path(OUTDIR, "STEP12_immune_clean_labeled.rds"), "\n")

# Plot
p_clean <- DimPlot(
  immune_clean,
  reduction = "wnn.umap",
  group.by  = "cell_type",
  label = TRUE, repel = TRUE
) + ggtitle("Immune microenvironment (clean, labeled; clusters 11/14/15 removed)")

print(p_clean)
save_plot(p_clean, "STEP12_immune_clean_labeled_umap.png", w=12, h=7)
# ----------------------------
# STEP 13) Immune composition by donor (counts + proportions + plot)
# ----------------------------
immune_clean <- readRDS(file.path(OUTDIR, "STEP12_immune_clean_labeled.rds"))

# Sanity checks
stopifnot("donor" %in% colnames(immune_clean@meta.data))
stopifnot("cell_type" %in% colnames(immune_clean@meta.data))

# Count table: donor x cell_type
comp_counts <- table(immune_clean$donor, immune_clean$cell_type)
comp_props  <- prop.table(comp_counts, margin = 1) * 100

# Save tables
write.csv(as.data.frame(comp_counts),
          file.path(OUTDIR, "STEP13_immune_composition_counts.csv"),
          row.names = FALSE)
write.csv(as.data.frame(round(comp_props, 3)),
          file.path(OUTDIR, "STEP13_immune_composition_percent.csv"),
          row.names = FALSE)

cat("Saved:\n",
    " -", file.path(OUTDIR, "STEP13_immune_composition_counts.csv"), "\n",
    " -", file.path(OUTDIR, "STEP13_immune_composition_percent.csv"), "\n")

# Long format for ggplot
comp_df <- as.data.frame(comp_props)
colnames(comp_df) <- c("Donor", "CellType", "Percent")

# Plot
p_comp <- ggplot(comp_df, aes(x = Donor, y = Percent, fill = CellType)) +
  geom_col(width = 0.85) +
  labs(title = "Immune cell composition by donor",
       x = "Donor", y = "Proportion (%)", fill = "celltype") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_comp)
save_plot(p_comp, "STEP13_immune_composition_by_donor.png", w=12, h=7)

# Quick summary (optional, but useful)
cat("\nCells per donor (immune_clean):\n")
print(table(immune_clean$donor))

# ============================
# STEP14_A) Option A: Donor immune state scores (RNA) + plots
# ============================


# Load your clean immune object
immune_clean <- readRDS(file.path(OUTDIR, "STEP12_immune_clean_labeled.rds"))
immune_clean <- JoinLayers(immune_clean)
DefaultAssay(immune_clean) <- "RNA"

# Safety checks (NO assumptions)
stopifnot("donor" %in% colnames(immune_clean@meta.data))
stopifnot("cell_type" %in% colnames(immune_clean@meta.data))
stopifnot("wnn.umap" %in% Reductions(immune_clean))

cat("Cells:", ncol(immune_clean), "\n")
cat("Donors:", paste(unique(immune_clean$donor), collapse = ", "), "\n")
cat("Cell types:", length(unique(immune_clean$cell_type)), "\n")

# ----------------------------
# 14A) Define immune gene programs (RNA)
# ----------------------------
gene_sets <- list(
  CYTOTOXICITY = c("NKG7","GNLY","PRF1","GZMB","GZMH","CTSW","FCGR3A","CCL5"),
  EXHAUSTION   = c("PDCD1","LAG3","HAVCR2","TIGIT","CTLA4","TOX"),
  IFN_ISG      = c("ISG15","IFIT1","IFIT2","IFIT3","MX1","OAS1","STAT1","IRF7"),
  INFLAM_MYELOID = c("S100A8","S100A9","FCN1","IL1B","LGALS3","LYZ","CTSS"),
  TAM_MACRO    = c("APOE","C1QA","C1QB","C1QC","TREM2","MSR1","MRC1","MARCO")
)

present_missing <- lapply(names(gene_sets), function(nm){
  genes <- gene_sets[[nm]]
  tibble::tibble(
    program = nm,
    gene = genes,
    present = genes %in% rownames(immune_clean)
  )
}) |> dplyr::bind_rows()

print(present_missing)

cat("\nMissing genes per program:\n")
print(present_missing |> dplyr::filter(!present) |> dplyr::count(program, name="n_missing"))

cat("\nPresent genes per program:\n")
print(present_missing |> dplyr::filter(present) |> dplyr::count(program, name="n_present"))

expr_check <- lapply(names(gene_sets), function(nm){
  g <- gene_sets[[nm]]
  avg <- Matrix::rowMeans(GetAssayData(immune_clean, layer="data")[g, , drop=FALSE])
  data.frame(program=nm, gene=g, avg_expr=as.numeric(avg))
}) |> dplyr::bind_rows()

# show lowest-expressed genes per program
expr_check %>%
  dplyr::group_by(program) %>%
  dplyr::arrange(avg_expr) %>%
  dplyr::slice_head(n=3) %>%
  print(n=200)

pct_check <- lapply(names(gene_sets), function(nm){
  g <- gene_sets[[nm]]
  mat <- GetAssayData(immune_clean, layer="data")[g, , drop=FALSE]
  pct <- Matrix::rowMeans(mat > 0) * 100
  data.frame(program=nm, gene=g, pct_cells=as.numeric(pct))
}) |> dplyr::bind_rows()

pct_check %>%
  dplyr::group_by(program) %>%
  dplyr::arrange(pct_cells) %>%
  dplyr::slice_head(n=3) %>%
  print(n=200)

# Keep only genes that exist in dataset (important and honest)
gene_sets <- lapply(gene_sets, function(g) intersect(g, rownames(immune_clean)))
print(sapply(gene_sets, length))

# Require at least 4 genes per set to score (prevents junk)
keep_sets <- names(gene_sets)[sapply(gene_sets, length) >= 4]
stopifnot(length(keep_sets) >= 3)
gene_sets <- gene_sets[keep_sets]

cat("Scoring sets:", paste(names(gene_sets), collapse=", "), "\n")

# Remove old columns if re-running
old <- c(paste0(names(gene_sets), "_score"))
old <- intersect(old, colnames(immune_clean@meta.data))
if (length(old) > 0) immune_clean@meta.data[, old] <- NULL

# ----------------------------
# 14B) Add module scores
# ----------------------------
# AddModuleScore expects a LIST of feature lists; name creates *_1 columns
immune_clean <- AddModuleScore(immune_clean, features = gene_sets, name = "STATE_")

# Map STATE_1, STATE_2 ... into readable columns
for (i in seq_along(gene_sets)) {
  nm <- names(gene_sets)[i]
  immune_clean[[paste0(nm, "_score")]] <- immune_clean[[paste0("STATE_", i)]]
}

# Save checkpoint
saveRDS(immune_clean, file.path(OUTDIR, "STEP14A_immune_clean_withStateScores.rds"))
cat("Saved:", file.path(OUTDIR, "STEP14A_immune_clean_withStateScores.rds"), "\n")

# ----------------------------
# 14C) UMAP feature plots (nice for capstone)
# ----------------------------
features_to_plot <- paste0(names(gene_sets), "_score")
p_umap <- FeaturePlot(immune_clean, features = features_to_plot, reduction="wnn.umap", ncol=3) &
  theme(plot.title = element_text(size=10))
print(p_umap)
ggsave(file.path(OUTDIR, "STEP14A_stateScores_umap.png"), p_umap, width=14, height=8)

# ----------------------------
# 14D) Donor-level plots (GLOBAL)
# ----------------------------
df <- immune_clean@meta.data %>%
  select(donor, cell_type, all_of(features_to_plot)) %>%
  tidyr::pivot_longer(cols = all_of(features_to_plot),
                      names_to = "program",
                      values_to = "score")

p_donor_global <- ggplot(df, aes(x=donor, y=score)) +
  geom_violin(trim=TRUE, scale="width") +
  geom_boxplot(width=0.15, outlier.size=0.2) +
  facet_wrap(~program, scales="free_y", ncol=3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Immune state scores by donor (all immune cells)",
       x="Donor", y="Module score")
print(p_donor_global)
ggsave(file.path(OUTDIR, "STEP14A_stateScores_byDonor_GLOBAL.png"),
       p_donor_global, width=14, height=8)

# ----------------------------
# 14E) Donor-level plots WITHIN major compartments (most important)
# ----------------------------
# Define major compartments from your cell_type labels (NO assumptions: use pattern matching)
df$compartment <- dplyr::case_when(
  df$cell_type %in% c(
    "Naive CD4 T (LEF1/TCF7/CCR7-like)",
    "Activated CD4 T (TNFRSF4/TNFRSF18/CTLA4/BATF+)",
    "CXCL13+ CD4 T (CXCL13/TOX2+)",
    "Cytotoxic T (CD8/NKG7/CCL5+)",
    "Activated/Exhausted CD8 T (CD8A/B, GZMH/LAG3-like)",
    "γδ T cells (TRDC/TRGC2+)"
  ) ~ "T_cells",
  
  df$cell_type %in% c(
    "NK cells (KIR/KLRC/GNLY/NCAM1+)"
  ) ~ "NK_cells",
  
  df$cell_type %in% c(
    "Naive B (MS4A1/CD79A/PAX5+)",
    "Activated / class-switched B (IGHA1/TNFRSF13B/C+)",
    "B cell subset (BLK/BANK1/FCRL1+)"
  ) ~ "B_cells",
  
  df$cell_type %in% c(
    "Classical Monocytes (S100A8/S100A9/FCN1/VCAN+)",
    "Intermediate Mono/Mac (CD14/FCGR3A/CCL18+)",
    "FOLR2+ Macrophages (FOLR2/SELENOP/APOE+)",
    "MARCO+ Macrophages (MARCO/NR1H3/APOC2+)"
  ) ~ "Myeloid",
  
  df$cell_type %in% c(
    "pDC (CLEC4C/LILRA4/LAMP5+)"
  ) ~ "DCs",
  
  TRUE ~ "Other"
)
# Keep main ones
df2 <- df %>% filter(compartment %in% c("T_cells","NK_cells","Myeloid","B_cells","DCs"))

p_by_comp <- ggplot(df2, aes(x=donor, y=score)) +
  geom_boxplot(outlier.size=0.2) +
  facet_grid(program ~ compartment, scales="free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Immune state scores by donor within compartments",
       x="Donor", y="Module score")
print(p_by_comp)
ggsave(file.path(OUTDIR, "STEP14A_stateScores_byDonor_WITHINcompartments.png"),
       p_by_comp, width=16, height=10)

# ----------------------------
# 14F) Stats: donor differences (Kruskal-Wallis per program)
# ----------------------------
stats_kw <- df2 %>%
  group_by(program, compartment) %>%
  summarise(
    n_cells = n(),
    p_kw = tryCatch(kruskal.test(score ~ donor)$p.value, error=function(e) NA_real_),
    .groups="drop"
  ) %>%
  arrange(p_kw)

write.csv(stats_kw, file.path(OUTDIR, "STEP14A_Kruskal_byDonor_withinCompartments.csv"), row.names=FALSE)
print(stats_kw)
cat("Saved stats:", file.path(OUTDIR, "STEP14A_Kruskal_byDonor_withinCompartments.csv"), "\n")

# 1) Read Step14F stats
stats_file <- file.path(OUTDIR, "STEP14A_Kruskal_byDonor_withinCompartments.csv")
stopifnot(file.exists(stats_file))
stats_kw <- readr::read_csv(stats_file, show_col_types = FALSE)

# 2) Clean names + compute -log10(p)
stats_kw2 <- stats_kw %>%
  mutate(
    program = str_replace(program, "_score$", ""),   # remove suffix
    neglog10p = -log10(p_kw),
    neglog10p_cap = pmin(neglog10p, 50)              # cap so extreme p-values don't blow up the color scale
  )

# 3) Order programs by strongest signal overall (nice heatmap ordering)
program_order <- stats_kw2 %>%
  group_by(program) %>%
  summarise(max_sig = max(neglog10p, na.rm=TRUE), .groups="drop") %>%
  arrange(desc(max_sig)) %>%
  pull(program)

comp_order <- c("T_cells","NK_cells","Myeloid","B_cells","DCs")

stats_kw2$program <- factor(stats_kw2$program, levels = program_order)
stats_kw2$compartment <- factor(stats_kw2$compartment, levels = comp_order)

# 4) Plot heatmap
p_heat <- ggplot(stats_kw2, aes(x = compartment, y = program, fill = neglog10p_cap)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(is.na(p_kw), "", sprintf("%.1f", neglog10p))),
            size = 3) +
  theme_classic() +
  labs(
    title = "Donor differences in immune state programs (Kruskal–Wallis)",
    subtitle = "Heatmap shows -log10(p) per program × compartment (higher = stronger donor effect)",
    x = "Compartment", y = "Program", fill = "-log10(p) (capped)"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face="bold")
  )

print(p_heat)
ggsave(file.path(OUTDIR, "STEP14F_KW_heatmap_neglog10p.png"), p_heat, width=9, height=5)
cat("Saved heatmap:", file.path(OUTDIR, "STEP14F_KW_heatmap_neglog10p.png"), "\n")


obj_file <- file.path(OUTDIR, "STEP14A_immune_clean_withStateScores.rds")
stopifnot(file.exists(obj_file))
immune_clean <- readRDS(obj_file)

# Programs we scored (columns end with _score)
score_cols <- grep("_score$", colnames(immune_clean@meta.data), value = TRUE)

# Build long table with compartment assignment (same logic you used)
df <- immune_clean@meta.data %>%
  dplyr::select(donor, cell_type, all_of(score_cols)) %>%
  tidyr::pivot_longer(cols = all_of(score_cols),
                      names_to = "program",
                      values_to = "score") %>%
  mutate(
    program = str_replace(program, "_score$", ""),
    compartment = case_when(
      grepl("NK", cell_type, ignore.case=TRUE) ~ "NK_cells",
      grepl("T",  cell_type, ignore.case=TRUE) ~ "T_cells",
      grepl("Mono|Mac", cell_type, ignore.case=TRUE) ~ "Myeloid",
      grepl("B", cell_type, ignore.case=TRUE) ~ "B_cells",
      grepl("pDC|DC", cell_type, ignore.case=TRUE) ~ "DCs",
      TRUE ~ "Other"
    )
  ) %>%
  filter(compartment %in% c("T_cells","NK_cells","Myeloid","B_cells","DCs"))

# Median score per donor within each program×compartment
med_by_donor <- df %>%
  group_by(program, compartment, donor) %>%
  summarise(med = median(score, na.rm=TRUE), .groups="drop")

# Effect direction summary: which donor highest/lowest + effect size (max median - min median)
effect_summary <- med_by_donor %>%
  group_by(program, compartment) %>%
  summarise(
    donor_high = donor[which.max(med)],
    donor_low  = donor[which.min(med)],
    med_high = max(med),
    med_low  = min(med),
    effect_range = med_high - med_low,
    .groups="drop"
  )

# Join with KW p-values from Step14F
stats_file <- file.path(OUTDIR, "STEP14A_Kruskal_byDonor_withinCompartments.csv")
stats_kw <- readr::read_csv(stats_file, show_col_types = FALSE) %>%
  mutate(program = str_replace(program, "_score$", ""))

top_hits <- stats_kw %>%
  left_join(effect_summary, by = c("program","compartment")) %>%
  arrange(p_kw) %>%
  mutate(neglog10p = -log10(p_kw))

# Save + display top 15
write.csv(top_hits, file.path(OUTDIR, "STEP14F_TOP_hits_with_direction.csv"), row.names = FALSE)
cat("Saved:", file.path(OUTDIR, "STEP14F_TOP_hits_with_direction.csv"), "\n")

print(top_hits %>%
        select(program, compartment, p_kw, neglog10p, effect_range, donor_high, donor_low, med_high, med_low) %>%
        head(15))


# ----------------------------
# STEP 15) Single-gene validation of donor differences (NO assumptions)
# Goal: confirm module-score donor effects with actual marker genes
# ----------------------------


# Load state-scored object (same one used for Step14)
obj_file <- file.path(OUTDIR, "STEP14A_immune_clean_withStateScores.rds")
stopifnot(file.exists(obj_file))
immune_clean <- readRDS(obj_file)
immune_clean <- JoinLayers(immune_clean)
DefaultAssay(immune_clean) <- "RNA"

# Sanity
stopifnot(all(c("donor","cell_type") %in% colnames(immune_clean@meta.data)))
cat("Cells:", ncol(immune_clean), "\n")
cat("Donors:", paste(sort(unique(immune_clean$donor)), collapse=", "), "\n")
cat("Cell types:", length(unique(immune_clean$cell_type)), "\n")

# Helper: safe plotting for genes that exist
genes_present <- function(g) intersect(g, rownames(immune_clean))

# ----------------------------
# 15A) Define key genes to validate (RNA) - matches Step14 programs
# ----------------------------
gene_panels <- list(
  CYTOTOXICITY = c("NKG7","GNLY","PRF1","GZMB","GZMH","CTSW","CCL5"),
  EXHAUSTION   = c("PDCD1","LAG3","HAVCR2","TIGIT","CTLA4","TOX"),
  IFN_ISG      = c("ISG15","IFIT1","IFIT2","IFIT3","MX1","OAS1","STAT1","IRF7"),
  INFLAM_MYELOID = c("S100A8","S100A9","FCN1","IL1B","LGALS3","LYZ","CTSS"),
  TAM_MACRO    = c("APOE","C1QA","C1QB","C1QC","TREM2","MSR1","MRC1","MARCO")
)

# Keep only genes actually present
gene_panels <- lapply(gene_panels, genes_present)
cat("\nGenes present per panel:\n")
print(sapply(gene_panels, length))
stopifnot(all(sapply(gene_panels, length) >= 3))

# ----------------------------
# 15B) Define compartments STRICTLY from YOUR cell_type labels
# (No regex guessing; uses exact labels you created in Step11)
# ----------------------------
meta <- immune_clean@meta.data

meta$compartment <- dplyr::case_when(
  meta$cell_type %in% c(
    "Naive CD4 T (LEF1/TCF7/CCR7-like)",
    "Activated CD4 T (TNFRSF4/TNFRSF18/CTLA4/BATF+)",
    "CXCL13+ CD4 T (CXCL13/TOX2+)",
    "Cytotoxic T (CD8/NKG7/CCL5+)",
    "Activated/Exhausted CD8 T (CD8A/B, GZMH/LAG3-like)",
    "γδ T cells (TRDC/TRGC2+)"
  ) ~ "T_cells",
  
  meta$cell_type %in% c("NK cells (KIR/KLRC/GNLY/NCAM1+)") ~ "NK_cells",
  
  meta$cell_type %in% c(
    "Naive B (MS4A1/CD79A/PAX5+)",
    "Activated / class-switched B (IGHA1/TNFRSF13B/C+)",
    "B cell subset (BLK/BANK1/FCRL1+)"
  ) ~ "B_cells",
  
  meta$cell_type %in% c(
    "Classical Monocytes (S100A8/S100A9/FCN1/VCAN+)",
    "Intermediate Mono/Mac (CD14/FCGR3A/CCL18+)",
    "FOLR2+ Macrophages (FOLR2/SELENOP/APOE+)",
    "MARCO+ Macrophages (MARCO/NR1H3/APOC2+)"
  ) ~ "Myeloid",
  
  meta$cell_type %in% c("pDC (CLEC4C/LILRA4/LAMP5+)") ~ "DCs",
  
  TRUE ~ "Other"
)

immune_clean$compartment <- meta$compartment
cat("\nCells per compartment:\n")
print(sort(table(immune_clean$compartment), decreasing = TRUE))

# ----------------------------
# 15C) Quick reality-check:
# Some genes should be compartment-enriched (not everywhere).
# We'll print % expressing for a few representative genes.
# ----------------------------
DefaultAssay(immune_clean) <- "RNA"
expr_mat <- GetAssayData(immune_clean, layer = "data")

rep_genes <- unique(c("PDCD1","LAG3","TOX","NKG7","GNLY","FCN1","S100A8","APOE","C1QC"))
rep_genes <- intersect(rep_genes, rownames(immune_clean))
cat("\nRepresentative genes used for sanity:", paste(rep_genes, collapse=", "), "\n")

pct_by_comp <- lapply(rep_genes, function(g){
  v <- expr_mat[g, ]
  data.frame(
    gene = g,
    compartment = immune_clean$compartment,
    expr = as.numeric(v > 0)
  ) %>%
    group_by(gene, compartment) %>%
    summarise(pct_expr = mean(expr)*100, .groups="drop")
}) %>% bind_rows()

cat("\n% cells expressing (expr>0) by compartment (selected genes):\n")
print(pct_by_comp %>% arrange(gene, desc(pct_expr)))

# ----------------------------
# 15D) Donor comparison plots WITHIN each compartment for key genes
# Output: 1 plot per (panel × compartment) that has enough cells
# ----------------------------
dir.create(OUTDIR, showWarnings = FALSE)

plot_gene_violin <- function(obj, genes, compartment_name, filename_prefix){
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) return(NULL)
  
  sub <- subset(obj, subset = compartment == compartment_name)
  if (ncol(sub) < 200) {
    cat("Skipping", filename_prefix, "-> too few cells:", ncol(sub), "\n")
    return(NULL)
  }
  
  # Long format expression from RNA "data" layer
  mat <- GetAssayData(sub, layer="data")[genes, , drop=FALSE]
  df_long <- as.data.frame(t(as.matrix(mat)))
  df_long$donor <- sub$donor
  
  df_long <- df_long %>%
    pivot_longer(cols = all_of(genes), names_to="gene", values_to="expr")
  
  p <- ggplot(df_long, aes(x = donor, y = expr)) +
    geom_violin(trim=TRUE, scale="width") +
    geom_boxplot(width=0.15, outlier.size=0.2) +
    facet_wrap(~gene, scales="free_y", ncol=4) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(
      title = paste0("Donor comparison (", compartment_name, "): single-gene expression"),
      subtitle = paste0("Genes: ", paste(genes, collapse=", ")),
      x="Donor", y="Log-normalized expression"
    )
  
  out_png <- file.path(OUTDIR, paste0(filename_prefix, "_", compartment_name, ".png"))
  ggsave(out_png, p, width=14, height=8)
  cat("Saved:", out_png, "\n")
  return(p)
}

# Run panels within the compartments where they biologically make sense.
# (Still NO assumptions: we just run and let the data show patterns.)
compartments_to_run <- c("T_cells","NK_cells","Myeloid","B_cells","DCs")

# For each compartment, run a focused set of panels
for (comp in compartments_to_run) {
  if (comp %in% c("T_cells","NK_cells")) {
    plot_gene_violin(immune_clean, gene_panels$CYTOTOXICITY, comp, "STEP15_genePanel_CYTOTOXICITY")
    plot_gene_violin(immune_clean, gene_panels$EXHAUSTION,   comp, "STEP15_genePanel_EXHAUSTION")
    plot_gene_violin(immune_clean, gene_panels$IFN_ISG,      comp, "STEP15_genePanel_IFNISG")
  }
  if (comp == "Myeloid") {
    plot_gene_violin(immune_clean, gene_panels$INFLAM_MYELOID, comp, "STEP15_genePanel_INFLAM_MYELOID")
    plot_gene_violin(immune_clean, gene_panels$TAM_MACRO,      comp, "STEP15_genePanel_TAM_MACRO")
    plot_gene_violin(immune_clean, gene_panels$IFN_ISG,        comp, "STEP15_genePanel_IFNISG")
  }
  if (comp == "B_cells") {
    plot_gene_violin(immune_clean, gene_panels$IFN_ISG, comp, "STEP15_genePanel_IFNISG")
    plot_gene_violin(immune_clean, gene_panels$EXHAUSTION, comp, "STEP15_genePanel_EXHAUSTION")
  }
  if (comp == "DCs") {
    plot_gene_violin(immune_clean, gene_panels$IFN_ISG, comp, "STEP15_genePanel_IFNISG")
  }
}

# Save checkpoint (with compartment column)
saveRDS(immune_clean, file.path(OUTDIR, "STEP15_immune_clean_with_compartment.rds"))
cat("Saved checkpoint:", file.path(OUTDIR, "STEP15_immune_clean_with_compartment.rds"), "\n")

# ----------------------------
# STEP 16) Check myeloid contamination inside T_cells (NO assumptions)
# ----------------------------

immune_clean <- readRDS(file.path(OUTDIR, "STEP15_immune_clean_with_compartment.rds"))
immune_clean <- JoinLayers(immune_clean)
DefaultAssay(immune_clean) <- "RNA"

# Focus only on T cells
tobj <- subset(immune_clean, subset = compartment == "T_cells")
cat("T_cells:", ncol(tobj), "\n")

# A strict myeloid-marker panel (if these are high in T cells, it's contamination/doublets/ambient)
myeloid_leak <- c("LYZ","LST1","TYROBP","FCER1G","AIF1","S100A8","S100A9","FCN1","C1QA","C1QB","C1QC","APOE")
myeloid_leak <- intersect(myeloid_leak, rownames(tobj))
cat("Myeloid-leak genes present:", length(myeloid_leak), "\n")

# Score myeloid-leak program inside T cells
tobj <- AddModuleScore(tobj, features = list(myeloid_leak), name = "MYELOID_LEAK_")
tobj$MYELOID_LEAK <- tobj$MYELOID_LEAK_1

# Summaries by donor
df <- tobj@meta.data %>%
  select(donor, MYELOID_LEAK)

cat("\nMYELOID_LEAK score summary (T cells):\n")
print(df %>% group_by(donor) %>% summarise(
  n = n(),
  med = median(MYELOID_LEAK),
  mean = mean(MYELOID_LEAK),
  sd = sd(MYELOID_LEAK),
  .groups="drop"
))

# Plot distribution by donor (T cells only)
p <- ggplot(df, aes(x=donor, y=MYELOID_LEAK)) +
  geom_violin(trim=TRUE, scale="width") +
  geom_boxplot(width=0.15, outlier.size=0.2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Myeloid-leak score within T cells (quality check)",
       x="Donor", y="Module score")
print(p)
ggsave(file.path(OUTDIR, "STEP16_Tcells_MYELOID_LEAK_byDonor.png"), p, width=10, height=6)

# How many T cells look suspiciously myeloid-like?
thr <- median(tobj$MYELOID_LEAK) + 2*sd(tobj$MYELOID_LEAK)
cat("\nThreshold used (median + 2*sd):", thr, "\n")
cat("T cells above threshold:", sum(tobj$MYELOID_LEAK > thr), "out of", ncol(tobj), "\n")

# ----------------------------
# STEP 17) Remove high-myeloid-leak T cells + rerun donor comparisons
# ----------------------------
immune_clean <- readRDS(file.path(OUTDIR, "STEP15_immune_clean_with_compartment.rds"))
immune_clean <- JoinLayers(immune_clean)
DefaultAssay(immune_clean) <- "RNA"

# Recompute leak on T cells (so it's reproducible)
tobj <- subset(immune_clean, subset = compartment == "T_cells")

myeloid_leak <- c("LYZ","LST1","TYROBP","FCER1G","AIF1",
                  "S100A8","S100A9","FCN1","C1QA","C1QB","C1QC","APOE")
myeloid_leak <- intersect(myeloid_leak, rownames(tobj))

tobj <- AddModuleScore(tobj, features = list(myeloid_leak), name = "MYELOID_LEAK_")
tobj$MYELOID_LEAK <- tobj$MYELOID_LEAK_1

thr <- median(tobj$MYELOID_LEAK) + 2*sd(tobj$MYELOID_LEAK)
cat("Leak threshold:", thr, "\n")
cat("T cells removed:", sum(tobj$MYELOID_LEAK > thr), "out of", ncol(tobj), "\n")

# Keep clean T cells
t_keep <- colnames(tobj)[tobj$MYELOID_LEAK <= thr]

# Keep ALL non-T cells as they are
nonT_keep <- colnames(immune_clean)[immune_clean$compartment != "T_cells"]

# Combined kept cells
keep_cells <- c(nonT_keep, t_keep)

immune_clean2 <- immune_clean[, keep_cells]
cat("Cells before:", ncol(immune_clean), "\n")
cat("Cells after :", ncol(immune_clean2), "\n")
cat("T cells before:", sum(immune_clean$compartment=="T_cells"), "\n")
cat("T cells after :", sum(immune_clean2$compartment=="T_cells"), "\n")

saveRDS(immune_clean2, file.path(OUTDIR, "STEP17_immune_clean_noLeakTcells.rds"))

# ---- Rerun Step 14F stats quickly on the filtered object ----
score_cols <- grep("_score$", colnames(immune_clean2@meta.data), value = TRUE)
df <- immune_clean2@meta.data %>%
  dplyr::select(donor, cell_type, compartment, all_of(score_cols)) %>%
  tidyr::pivot_longer(cols = all_of(score_cols), names_to="program", values_to="score") %>%
  mutate(program = stringr::str_replace(program, "_score$", "")) %>%
  filter(compartment %in% c("T_cells","NK_cells","Myeloid","B_cells","DCs"))

stats_kw2 <- df %>%
  group_by(program, compartment) %>%
  summarise(
    n_cells = n(),
    p_kw = tryCatch(kruskal.test(score ~ donor)$p.value, error=function(e) NA_real_),
    .groups="drop"
  ) %>%
  arrange(p_kw)

write.csv(stats_kw2, file.path(OUTDIR, "STEP17_Kruskal_byDonor_withinCompartments_noLeakTcells.csv"),
          row.names = FALSE)

print(stats_kw2)

# ----------------------------
# STEP 18) FINAL PLOTS (post Step17 cleanup)
# ----------------------------

# 0) Load cleaned object (post leak removal)
obj_file <- file.path(OUTDIR, "STEP17_immune_clean_noLeakTcells.rds")
stopifnot(file.exists(obj_file))
immune_final <- readRDS(obj_file)
immune_final <- JoinLayers(immune_final)
DefaultAssay(immune_final) <- "RNA"

# 0b) Make an output folder for FINAL plots
FINALDIR <- file.path(OUTDIR, "FINAL_PLOTS_STEP18")
dir.create(FINALDIR, showWarnings = FALSE)
cat("Saving final plots to:", FINALDIR, "\n")

save_plot2 <- function(p, name, w=12, h=7) {
  out <- file.path(FINALDIR, name)
  ggsave(out, p, width=w, height=h)
  cat("Saved:", out, "\n")
}

# Sanity
stopifnot(all(c("donor","cell_type","compartment") %in% colnames(immune_final@meta.data)))
stopifnot("wnn.umap" %in% Reductions(immune_final))
cat("Cells:", ncol(immune_final), " | Donors:", length(unique(immune_final$donor)),
    " | Cell types:", length(unique(immune_final$cell_type)), "\n")

# ----------------------------
# 18A) UMAP labeled (cell types)
# ----------------------------
p_umap_label <- DimPlot(
  immune_final,
  reduction = "wnn.umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) + ggtitle("Immune microenvironment (final; post leak-removal)")

save_plot2(p_umap_label, "A1_UMAP_celltype_FINAL.png", w=13, h=7)

# ----------------------------
# 18B) Composition by donor (percent)
# ----------------------------
comp_counts <- table(immune_final$donor, immune_final$cell_type)
comp_props  <- prop.table(comp_counts, margin = 1) * 100
comp_df <- as.data.frame(comp_props)
colnames(comp_df) <- c("Donor","CellType","Percent")

p_comp <- ggplot(comp_df, aes(x=Donor, y=Percent, fill=CellType)) +
  geom_col(width=0.85) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Immune cell-type composition by donor (final)",
       x="Donor", y="Proportion (%)", fill="Cell type")

save_plot2(p_comp, "B1_Composition_byDonor_FINAL.png", w=13, h=7)

# Also save tables (final)
write.csv(as.data.frame(comp_counts),
          file.path(FINALDIR, "B1_Composition_counts_FINAL.csv"),
          row.names = FALSE)
write.csv(as.data.frame(round(comp_props, 3)),
          file.path(FINALDIR, "B1_Composition_percent_FINAL.csv"),
          row.names = FALSE)

# ----------------------------
# 18C) State score UMAPs (if scores exist)
# ----------------------------
score_cols <- setdiff(
  grep("_score$", colnames(immune_final@meta.data), value = TRUE),
  "doublet_score"
)
cat("Score columns found:", paste(score_cols, collapse=", "), "\n")
stopifnot(length(score_cols) >= 3)

p_state_umap <- FeaturePlot(
  immune_final,
  features = score_cols,
  reduction = "wnn.umap",
  ncol = 3
) & theme(plot.title = element_text(size=10))

ggsave(file.path(FINALDIR, "C1_StateScores_UMAP_FINAL.png"),
       p_state_umap, width=14, height=8)
cat("Saved:", file.path(FINALDIR, "C1_StateScores_UMAP_FINAL.png"), "\n")

# ----------------------------
# 18D) Donor comparison within compartments (boxplots) - KEY FIGURE
# ----------------------------
df_long <- immune_final@meta.data %>%
  select(donor, compartment, all_of(score_cols)) %>%
  pivot_longer(cols = all_of(score_cols), names_to="program", values_to="score") %>%
  mutate(program = str_replace(program, "_score$", "")) %>%
  filter(compartment %in% c("T_cells","NK_cells","Myeloid","B_cells","DCs"))

p_by_comp <- ggplot(df_long, aes(x=donor, y=score)) +
  geom_boxplot(outlier.size=0.2) +
  facet_grid(program ~ compartment, scales="free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Immune state scores by donor within compartments (final)",
       x="Donor", y="Module score")

save_plot2(p_by_comp, "D1_StateScores_byDonor_WITHINcompartments_FINAL.png", w=16, h=10)

# ----------------------------
# 18E) Kruskal–Wallis heatmap (final)
# ----------------------------
stats_kw <- df_long %>%
  group_by(program, compartment) %>%
  summarise(
    n_cells = n(),
    p_kw = tryCatch(kruskal.test(score ~ donor)$p.value, error=function(e) NA_real_),
    .groups="drop"
  ) %>%
  mutate(
    neglog10p = -log10(p_kw),
    neglog10p_cap = pmin(neglog10p, 50)
  )

write.csv(stats_kw, file.path(FINALDIR, "E1_Kruskal_byDonor_withinCompartments_FINAL.csv"),
          row.names = FALSE)

program_order <- stats_kw %>%
  group_by(program) %>%
  summarise(max_sig = max(neglog10p, na.rm=TRUE), .groups="drop") %>%
  arrange(desc(max_sig)) %>%
  pull(program)

comp_order <- c("T_cells","NK_cells","Myeloid","B_cells","DCs")

stats_kw$program <- factor(stats_kw$program, levels = program_order)
stats_kw$compartment <- factor(stats_kw$compartment, levels = comp_order)

p_heat <- ggplot(stats_kw, aes(x=compartment, y=program, fill=neglog10p_cap)) +
  geom_tile(color="white", linewidth=0.4) +
  geom_text(aes(label = ifelse(is.na(p_kw), "", sprintf("%.1f", neglog10p))),
            size=3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=30, hjust=1),
        plot.title = element_text(face="bold")) +
  labs(title="Donor differences in immune state programs (Kruskal–Wallis, final)",
       subtitle="Heatmap shows -log10(p) per program × compartment (higher = stronger donor effect)",
       x="Compartment", y="Program", fill="-log10(p) (capped)")

save_plot2(p_heat, "E2_KW_heatmap_neglog10p_FINAL.png", w=9, h=5)

# ----------------------------
# 18F) Direction table: which donor high/low (final)
# ----------------------------
med_by_donor <- df_long %>%
  group_by(program, compartment, donor) %>%
  summarise(med = median(score, na.rm=TRUE), .groups="drop")

effect_summary <- med_by_donor %>%
  group_by(program, compartment) %>%
  summarise(
    donor_high = donor[which.max(med)],
    donor_low  = donor[which.min(med)],
    med_high = max(med),
    med_low  = min(med),
    effect_range = med_high - med_low,
    .groups="drop"
  )

top_hits <- stats_kw %>%
  left_join(effect_summary, by=c("program","compartment")) %>%
  arrange(p_kw)

write.csv(top_hits, file.path(FINALDIR, "F1_TOP_hits_with_direction_FINAL.csv"),
          row.names = FALSE)

cat("\nDONE. Final plots + tables saved in:\n", FINALDIR, "\n")











