# GSE189926 resolution and UMAP topology sweep

Status: `COMPLETED_AWAITING_MAC_REVIEW`

- Exact input path, size, SHA-256, cell count, assays, layers, reductions, commands, and metadata fields were verified.
- The saved object had no graph; `harmony_snn_sweep` was rebuilt from Harmony dimensions 1:30 with `k.param=15`.
- Seven resolutions were calculated on the same graph and displayed on one unchanged baseline UMAP.
- Seven UMAP routes were evaluated separately on the fixed Harmony reduction.
- All 89,587 accepted QC-pass cells were retained; local linearity was used only as a topology diagnostic.
- Evidence-ranked resolution for Mac review: `0.5`, which supports retaining the existing formal resolution; formal object unchanged.
- Evidence-ranked UMAP route for Mac review: `neighbor30`; formal object unchanged.
- Baseline q99 filament audit contains 896 cells across 17 resolution-0.5 clusters.
- Every PNG/PDF comparison sheet uses a white background, fixed seed 340 plotting order, alpha 0.90, and the approved bright palette.
- The sweep object is local-only and does not overwrite the prior object.
- Seurat 5.4.0 called uwot through the default `future::nbrOfWorkers()=1` plan; the actual route table records one worker.
- The post-UMAP summary resumed from a verified checkpoint after a grouped-median type-consistency error; the corrected source casts all grouped medians to numeric, and the successful resume log is included.
