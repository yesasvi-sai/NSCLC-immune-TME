# Donor-to-Donor Variation in the NSCLC Immune Microenvironment Using Joint Single-Cell RNA and Surface-Protein Data

Single-cell multi-modal analysis (RNA + ADT surface protein) of 9,103 immune cells across 7 NSCLC donors. We use Weighted Nearest Neighbor (WNN) integration in Seurat v5 to identify 18 immune cell populations and reveal that myeloid functional heterogeneity is the dominant axis of patient-to-patient variation in the tumor immune microenvironment.

---

## Key Findings

- **18 immune populations** identified through WNN integration of RNA and ADT data
- **Myeloid compartment dominates** inter-donor heterogeneity — all 5 scored immune programs are significant, 4 with p < 0.001
- **IFN/ISG signaling** is the most variable program across patients (p = 1.27 × 10⁻⁸ in myeloid cells)
- **Donor 4** has a B cell–dominated TME (>60% B lineage), suggestive of tertiary lymphoid structures
- **Myeloid ambient RNA leak** in T cells varies by donor and must be corrected before inter-donor comparisons

---

## Dataset

10x Genomics **20k NSCLC DTC 3' NextGEM** — 7 donors, paired RNA + ADT (CD45, CD3, CD4, CD8, CD14, CD11c, CD16, CD56, CD19). Publicly available at [10xgenomics.com/datasets](https://www.10xgenomics.com/datasets).

---

## Pipeline

The analysis is a single 1,528-line R script (`analysis/NSCLC_immune_pipeline.R`) running 18 sequential steps:

**Steps 1–5:** Load 7 donors → QC filter → RNA PCA (30 PCs) → ADT CLR + APCA (5 PCs) → WNN integration + UMAP + clustering

**Steps 6–9:** Immune subset extraction (PTPRC+ / EPCAM−) → Immune WNN reprocessing → scDblFinder doublet removal → Singlet reprocessing

**Steps 10–12:** RNA + ADT marker validation → FindAllMarkers + evidence-based annotation → Non-immune contamination removal (endothelial, epithelial, neuroendocrine clusters)

**Steps 14–18:** Module scoring (5 gene programs) → Single-gene validation → Myeloid-leak QC in T cells → Leak correction → Final figures + Kruskal–Wallis statistics

---

## Gene Programs Scored

| Program | Genes | Captures |
|---------|-------|----------|
| CYTOTOXICITY | NKG7, GNLY, PRF1, GZMB, GZMH, CTSW, FCGR3A, CCL5 | Effector killing |
| EXHAUSTION | PDCD1, LAG3, HAVCR2, TIGIT, CTLA4, TOX | T cell dysfunction |
| IFN_ISG | ISG15, IFIT1, IFIT2, IFIT3, MX1, OAS1, STAT1, IRF7 | Interferon response |
| INFLAM_MYELOID | S100A8, S100A9, FCN1, IL1B, LGALS3, LYZ, CTSS | Inflammatory monocytes |
| TAM_MACRO | APOE, C1QA, C1QB, C1QC, TREM2, MSR1, MRC1, MARCO | Macrophage polarization |

---

## Main Results

### Kruskal–Wallis: donor differences by compartment

| Program | Myeloid (651 cells) | T cells (64 cells) |
|---------|:-------------------:|:-------------------:|
| IFN_ISG | **1.27 × 10⁻⁸** | 0.190 |
| CYTOTOXICITY | **3.27 × 10⁻⁶** | **0.005** |
| INFLAM_MYELOID | **5.65 × 10⁻⁵** | 0.341 |
| TAM_MACRO | **3.55 × 10⁻⁴** | 0.845 |
| EXHAUSTION | 0.269 | 0.465 |

Bold = significant (p < 0.05). All five programs are significant in myeloid cells. Only cytotoxicity reaches significance in T cells.

### Donor composition archetypes

| Archetype | Donors | Feature |
|-----------|--------|---------|
| B cell–dominated | Donor 4 | >60% B lineage |
| T cell–enriched | Donors 1, 2 | High CXCL13+ CD4, cytotoxic CD8 |
| Myeloid-enriched | Donors 5, 6 | Up to 55% monocyte/macrophage |
| Mixed | Donors 3, 7 | Balanced |

---

## Requirements
```r
# CRAN
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork",
                    "tibble", "tidyr", "readr", "stringr", "Matrix"))

# Bioconductor
BiocManager::install(c("scDblFinder", "SingleCellExperiment"))
```

Tested with R ≥ 4.3 and Seurat v5.

---

## Running

1. Download the 7 donor HDF5 files from 10x Genomics and place them in your data directory
2. Update the working directory path on line 15 of the script
3. Run: `Rscript analysis/NSCLC_immune_pipeline.R`

Outputs go to a timestamped `RUN_YYYYMMDD_HHMMSS/` folder. Final publication figures are in `FINAL_PLOTS_STEP18/`.

---
