# Detailed Methods

Extended technical documentation for each pipeline step. For the high-level overview, see the main [README](../README.md).

---

## Step 1 – Data Loading & Merging

Seven donor HDF5 files are loaded with `Seurat::Read10X_h5()`. Both `Gene Expression` (RNA) and `Antibody Capture` (ADT) modalities are extracted. Per-donor Seurat objects are created (`min.cells = 3`, `min.features = 200`), tagged with donor identity, and merged.

## Step 2 – Quality Control

Four QC metrics are computed per cell: `nFeature_RNA`, `nCount_RNA`, `percent.mt` (MT- genes), and `percent.ribo` (RPS/RPL genes). Thresholds applied:

| Filter | Threshold | Rationale |
|--------|-----------|-----------|
| nFeature_RNA | > 200 | Remove empty/dying cells |
| nFeature_RNA | < 7,500 | Remove likely doublets |
| percent.mt | < 20% | Remove damaged cells |

## Step 3 – RNA Normalization & PCA

LogNormalize (scale factor 10,000) → 2,000 HVGs (VST) → ScaleData → PCA (50 PCs). Elbow plot guides selection of 30 PCs for downstream use.

## Step 4 – ADT Normalization & APCA

All ADT features set as variable → CLR normalization (margin = 2, across cells) → ScaleData → PCA with 5 components (`apca` reduction).

## Step 5 – WNN Integration

`FindMultiModalNeighbors` with RNA PCs 1:30 and ADT PCs 1:5 produces a per-cell weighted neighborhood. UMAP and Louvain clustering (SLM algorithm, resolution 0.5) are run on the WNN graph.

## Step 6 – Immune Subset Extraction

A curated marker panel (PTPRC, CD3D, CD3E, TRAC, CD8A, MS4A1, CD79A, LYZ, S100A8, EPCAM, KRT18, PECAM1, COL1A1, etc.) is evaluated via DotPlot. PTPRC-high / EPCAM-low clusters are classified as immune and extracted.

## Steps 7–9 – Immune Reprocessing & Doublet Removal

Three-pass approach: (7) full WNN reprocessing at resolution 0.8, (8) scDblFinder doublet detection on SingleCellExperiment conversion, (9) singlet-only WNN reprocessing. This ensures doublets don't distort the immune embedding.

## Steps 10–12 – Annotation & Contamination Removal

Annotation uses three convergent evidence streams:
- **FindAllMarkers**: Unbiased DE per cluster (30,673 total marker genes)
- **RNA DotPlot**: 60+ curated markers across all immune lineages
- **ADT DotPlot**: CD45, CD3, CD4, CD8, CD14, CD11c, CD16, CD56, CD19

Three contamination clusters removed: endothelial (HSPG2/COL4A1/CAV1), epithelial (SLPI/WFDC2/AKR1C), neuroendocrine (PHOX2B/OTX2/NEUROD1).

## Step 14 – Immune State Scoring

`AddModuleScore` computes background-subtracted enrichment for 5 gene programs (CYTOTOXICITY, EXHAUSTION, IFN_ISG, INFLAM_MYELOID, TAM_MACRO). Genes are binned by expression and compared to matched control genes.

## Steps 16–17 – Myeloid-Leak Correction

A myeloid-leak module (12 genes: LYZ, LST1, TYROBP, FCER1G, AIF1, S100A8/A9, FCN1, C1QA/B/C, APOE) is scored in T cells. Cells above median + 2×SD are removed. Statistical results shift after correction, validating this step.


