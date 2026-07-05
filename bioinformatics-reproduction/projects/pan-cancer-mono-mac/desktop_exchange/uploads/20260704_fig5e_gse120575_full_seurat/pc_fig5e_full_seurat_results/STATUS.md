# Fig5E GSE120575 full Seurat PC run status

Run completed on 2026-07-05 Asia/Shanghai.

## Scope

Request: `20260704_fig5e_gse120575_full_seurat`

Workflow executed:

1. Downloaded GSE120575 TPM and patient metadata from GEO.
2. Converted dense TPM text to sparse MatrixMarket locally.
3. Built a Seurat object from TPM values without additional normalization.
4. Clustered all cells at resolution 1.0.
5. Subset paper-stated myeloid clusters `8`, `12`, and `15`.
6. Reclustered the subset at resolution 0.4.
7. Annotated myeloid subclusters by correlation to the supplied Fig. 2 monocyte/macrophage average-expression reference.
8. Computed sample-level proportions and Wilcoxon summaries for `C1QC+ TAMs` and `THBS1+ MDSCs + SPP1+ TAMs`.

## Results summary

- Myeloid subset cells: 1,637
- Reclustering resolution: 0.4
- Myeloid subclusters: 4 (`g0`, `g1`, `g2`, `g3`)
- Subcluster sizes: `g0` 566, `g1` 517, `g2` 426, `g3` 128
- Best Fig. 2 reference matches:
  - `g0` -> Fig. 2 cluster 8, `AMs`, transfer score 0.181738173666529
  - `g1` -> Fig. 2 cluster 8, `AMs`, transfer score 0.0437478707402778
  - `g2` -> Fig. 2 cluster 31, `S100A12+ MDSCs`, transfer score -0.0429996693256586
  - `g3` -> Fig. 2 cluster 0, `AMs`, transfer score -0.0185957288075855

## Fig5E outcome

The requested target lineages were not assigned to any reclustered myeloid subcluster under this exact full-Seurat/correlation workflow.

- `C1QC+ TAMs`: total assigned cells = 0
- `THBS1+ MDSCs + SPP1+ TAMs`: total assigned cells = 0
- All sample-level proportions are 0.
- Wilcoxon p-values are `NA` because both response groups have constant zero proportions.

This was not forced or manually relabeled. The marker audit and subcluster annotation table are included for review.

## Visual QA

The generated PNG was opened locally and visually inspected. Text is not clipped, panels do not overlap, and axis labels are visible. The plotted values are all zero with `NA` labels, matching the generated tables.

## Local execution notes

The PC runner used:

- Rscript: `H:\R\bin\Rscript.exe`
- Python: `H:\python312\python.exe`
- R library: `H:\R\library`
- Temporary directory: `I:\codex-run-tmp\fig5e_gse120575_full_seurat`

Local script adjustments made for this PC run:

- Resolved `request_dir` from the script path instead of `sys.frame(1)$ofile`.
- Used `FIG5E_PYTHON` to call the H-drive Python executable.
- Avoided duplicate `cell_id` metadata names by storing Seurat row names as `seurat_cell_id`.
- Set myeloid subcluster labels to `g0/g1/g2/g3` before average-expression annotation so Seurat group names match metadata group names.

## GitHub handoff contents

This GitHub handoff includes text outputs only:

- `tables/fig5e_full_cluster_marker_audit.tsv`
- `tables/fig5e_full_myeloid_subcluster_annotation.tsv`
- `tables/fig5e_full_sample_proportions.tsv`
- `tables/fig5e_full_wilcoxon.tsv`
- `logs/run_fig5e_gse120575_full_seurat.log`
- `files/gse120575_myeloid_reclustered_seurat.rds.not_uploaded.txt`
- `upload_manifest.tsv`

Generated locally but not uploaded in this handoff:

- `outputs/figures/fig5e_full_seurat_proportions.png`
- `outputs/figures/fig5e_full_seurat_proportions.pdf`
- `outputs/files/gse120575_myeloid_reclustered_seurat.rds`
- Raw GEO downloads under `outputs/raw/`
- Sparse MatrixMarket intermediates under `outputs/sparse/`

Reason: the Codex GitHub connector can commit text files, but it does not stream the generated 39.8 MB local RDS. Local PATs available in this session returned 404 for this private repository through the GitHub REST contents API.

See `upload_manifest.tsv` for uploaded byte sizes and SHA256 hashes. The RDS local path, byte size, and SHA256 are recorded in `files/gse120575_myeloid_reclustered_seurat.rds.not_uploaded.txt`.