# Mac review request: GSE205506 annotation strategy

## What is complete

- Table S1 response mapping passed exactly for all 40 GEO samples and 19 subjects.
- The R-only Seurat workflow produced 238,934 retained cells and 28 working clusters at resolution 0.6.
- Top10/Top50 markers, reviewed broad annotations, cluster QC, author marker overlap, figures, and reproducibility files are included in this upload.
- The author marker workbook supplied by the user was parsed exactly: 45 clusters from five compartment-specific analyses.

## Interpretation of the author marker workbook

The workbook is useful as an external marker reference, but its 45 labels should not be transferred directly to the current 28 whole-dataset clusters. The author labels come from separate T/I/NK, B, Myeloid, Endothelial, and Fibroblast reclustering analyses. The current object was clustered once across all cells, and the workbook contains no epithelial compartment.

Same-compartment overlap supports several broad calls, including Bn, pB IgG, FOLR2+Mac, CXCL12+Fibro, Myofibroblast, and LYVE1+LEC. However, cluster 3 is canonically CD4 T and its strongest author overlap is CD8+MAIT, so overlap rank alone is not sufficient for subtype transfer.

## QC concern

Nine working clusters contain at least five mitochondrial genes among their Top10 markers: 2, 6, 12, 15, 19, 21, 23, 24, and 27. Together they contain 53,708 cells (22.478%). Several have extreme median mitochondrial fractions, including 92.06% (cluster 2), 96.76% (cluster 12), and 99.82% (cluster 27).

The cells are retained and flagged in `cluster_quality_reviewed`; no hidden deletion was performed.

## Requested Mac decision

1. Approve excluding or separately reprocessing the nine mitochondrial-marker-enriched clusters before subtype-level reclustering.
2. After that QC decision, approve splitting the clean cells into the same five author compartments and reclustering each compartment before assigning the 45 author subtype labels.

Recommended next step: perform the QC cleanup first, then compartment-specific reclustering and marker validation. Do not adopt the 45 subtype labels directly from the current 28-cluster object.

