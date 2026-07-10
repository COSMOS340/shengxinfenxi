# PC Note: GSE189926 UMAP Repair Policy

Created: 2026-07-10

## Decision

Do not restart the whole Pan-GI ICI workflow from scratch.

The current issue is local to the `GSE189926` object construction step:

- the existing `gse189926_umap_overview.png` is a failed QC embedding
- the embedding is strongly structured by sample and patient
- the labels are too broad for downstream communication analysis
- the marker dotplot is too compressed for label review

Use the PC-local object from:

`uploads/20260709_response_object_construction/object_inventory.tsv`

Then run the current request:

- `20260710_gse189926_annotation_qc_refinement_manifest.tsv`
- `20260710_gse189926_annotation_qc_refinement_request.md`

## Software Route

Scanpy is acceptable. The problem is not the software name; the problem is that the previous route did not sufficiently audit QC, sample structure, HVG/PCA drivers, and refined immune labels.

If using Scanpy, record exact choices for:

- QC-pass thresholds
- HVG rule
- technical gene exclusion rule for PCA-driving features
- PCA dimensions
- neighbor parameters
- UMAP parameters
- Leiden resolution
- sample-aware correction method

If PC prefers Seurat, that is also acceptable, but it should not become a full project restart. Use the same biological and QC requirements and record the Seurat route in `gse189926_umap_repair_strategy_audit.tsv`.

## Required UMAP Comparison

Upload three UMAPs:

1. `gse189926_umap_raw_qc.png`
2. `gse189926_umap_qc_pass.png`
3. `gse189926_umap_sample_aware_corrected.png`

Also upload:

- `gse189926_umap_repair_strategy_audit.tsv`
- `gse189926_umap_neighbor_mixing_metrics.tsv`
- `gse189926_sample_patient_mixing_audit.tsv`

The aim is to determine whether the embedding is dominated by technical or sample-level structure. Do not erase all patient biology blindly.

## Stop Condition

Do not run CellChat, LIANA, NicheNet, differential communication, final statistics, or manuscript figures until these repaired UMAPs and refined labels are reviewed.
