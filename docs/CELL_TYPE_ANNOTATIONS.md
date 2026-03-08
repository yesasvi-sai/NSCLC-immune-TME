# Cell Type Annotation 

Per-cluster annotation rationale for the 18 immune populations and 3 removed contamination clusters. All labels are based on convergent evidence from FindAllMarkers, RNA DotPlot, ADT DotPlot, and FeaturePlot localization.

---

## Myeloid Compartment

| Cluster | Label | Top RNA Markers | ADT Validation |
|:-------:|-------|-----------------|----------------|
| 0 | Inflammatory/Classical Monocytes | S100A8, S100A9, FCN1, VCAN, PLXDC2, AIF1 | CD14↑ CD45↑ CD3↓ |
| 3 | FCGR3A+ Mono/Mac | FCGR3A, CD14, CCL18, MARCO, MS4A7, LST1 | CD16↑ CD14↑ |
| 6 | FOLR2+ Macrophages | FOLR2, SELENOP, APOE, C1QA, C1QB, C1QC | CD45↑ CD14mod |
| 12 | MARCO+ Macrophages | MARCO, NR1H3, APOC2, MSR1, MRC1 | CD45↑ CD14mod |
| 18 | FCGR3A+ MARCO+ Myeloid | SIGLEC10, LILRA5, FCGR3A, MARCO | CD16↑ |

## T Cell Compartment

| Cluster | Label | Top RNA Markers | ADT Validation |
|:-------:|-------|-----------------|----------------|
| 2 | Cytotoxic CD8 T (LAG3+) | CD8A, CD8B, NKG7, CCL5, GZMK, LAG3 | CD8↑ CD3↑ |
| 9 | Cytotoxic CD8/NK-like T | CD8A, IFNG, NKG7, GZMB, PRF1, GNLY | CD8↑ CD56mod |
| 7 | IL7R+ CD4 T (memory) | IL7R, KLRB1, RORA, LTB, MAL, TCF7 | CD4↑ CD3↑ |
| 5 | CXCL13+ CD4 T (TLS-like) | CXCL13, TOX2, PDCD1, ICOS | CD4↑ CD3↑ |
| 4 | Activated CD4 / Treg-like | TNFRSF4, TNFRSF18, CTLA4, ICOS, IL2RA, FOXP3 | CD4↑ CD3↑ |
| 10 | Activated T | TNFRSF4, BATF, ICOS, LEF1 | CD3↑ |
| 19 | γδ T cells | TRDC, TRGC1, TRGC2, NKG7 | CD3mod CD4↓ CD8↓ |

## B Cell Compartment

| Cluster | Label | Top RNA Markers | ADT Validation |
|:-------:|-------|-----------------|----------------|
| 1 | B cells (memory-like) | FCRL1, BANK1, MS4A1, EBF1, CD79A, PAX5 | CD19↑ |
| 13 | B cells (BANK1/BLK subset) | BANK1, BLK, VPREB3, MS4A1 | CD19↑ |
| 8 | Activated/Class-switched B | IGHA1, TNFRSF13B, TNFRSF13C, JCHAIN | CD19 variable |

## Other Populations

| Cluster | Label | Top RNA Markers | ADT Validation |
|:-------:|-------|-----------------|----------------|
| 16 | NK cells | KIR2DL4, KLRC1, KLRC2, NCAM1, GNLY, PRF1 | CD56↑ CD16↑ CD3↓ |
| 17 | pDC | SCT, CLEC4C, LILRA4, LAMP5, IRF7, TCF4 | Low for lineage markers |

## Removed: Non-Immune Contamination

| Cluster | Identity | Markers | Reason |
|:-------:|----------|---------|--------|
| 11 | Endothelial/Stromal | HSPG2, COL4A1, CAV1, VWF | Non-immune lineage |
| 14 | Epithelial/Secretory | SLPI, WFDC2, AKR1C1/2/3 | Non-immune lineage |
| 15 | Neuroendocrine | PHOX2B, OTX2, NEUROD1, INSM1 | Non-immune lineage |
