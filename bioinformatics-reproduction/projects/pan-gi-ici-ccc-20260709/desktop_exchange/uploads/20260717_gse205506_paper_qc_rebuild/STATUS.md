# GSE205506 source-paper QC rebuild

Status: `COMPLETE_GSE205506_PAPER_QC_REBUILD`

## Formal route

- Re-imported all 40 author matrices: 324,020 raw barcodes from 40 samples and 19 subjects.
- Article basic-QC pass: 200,663 cells; scDblFinder was sensitivity-only and excluded zero cells from the primary route.
- Compartment mitochondrial-QC pass and final atlas: 166,624 cells.
- Article reported total: 155,397; observed difference: +11227 cells. Equality was not forced.
- First-round broad clusters: 34; final broad clusters: 35; exactly six broad compartments represented in both.
- No whole broad or subtype cluster was deleted.

## Compartment mitochondrial QC

- Lymphoid: 1322/57464 removed (2.30%); article comparison 9.20%.
- Myeloid: 468/9972 removed (4.69%); article comparison 12.84%.
- Fibroblast: 112/8783 removed (1.28%); article comparison 8.11%.
- Endothelial: 190/15420 removed (1.23%); article comparison 8.50%.
- Epithelial: 31947/109024 removed (29.30%); article comparison 29.75%.

## Compartment review

- T/I/NK: 15 clusters; 10 author-overlap labels supported; 5 retained as `Unresolved_review`.
- B: 13 clusters; 7 author-overlap labels supported; 6 retained as `Unresolved_review`.
- Myeloid: 16 clusters; 10 author-overlap labels supported; 6 retained as `Unresolved_review`.
- Endothelial: 17 clusters; 11 author-overlap labels supported; 6 retained as `Unresolved_review`.
- Fibroblast: 18 clusters; 10 author-overlap labels supported; 8 retained as `Unresolved_review`.
- Epithelial was saved as a separate sixth object and received no mmc3-derived subtype label.

## Verification

- All required PNG/PDF figures were generated with the approved bright palette.
- Visual inspection passed for 16 PNG figures; details are in `gse205506_paper_qc_figure_visual_qc.tsv`.
- Exact local paths, byte sizes, and SHA-256 hashes are in `gse205506_paper_qc_local_object_inventory.tsv`.
- The prior QC-version-1 RDS files were preserved and remain listed as rejected in `gse205506_rejected_qc_v1_object_inventory.tsv`.

