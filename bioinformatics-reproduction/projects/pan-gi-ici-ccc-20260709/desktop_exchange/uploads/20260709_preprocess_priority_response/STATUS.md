# Pan-GI ICI priority response preprocessing

Completed: 2026-07-09

Local matrix directory: `I:\shengxinfenxi\bioinformatics-reproduction\projects\pan-gi-ici-ccc-20260709\03_matrices`
Upload directory: `bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260709_preprocess_priority_response`

Summary:
- Checksum verification: 32/32 verified.
- GSE189926: 22 text matrices shape-checked and joined to GEO outcome/treatment/MMR metadata; no dense import or filtering.
- GSE235863 CD45/CD8: h5ad obs metadata summarized, response joined by patient, and P18 ambiguity preserved in join audit.
- ICB_Zenodo_Liver_Ma: Seurat metadata summarized; HCC and iCCA kept separate.
- OMIX001073: zip matrix dimensions and annotation counts audited without extracting full matrices.
- GSE236581: metadata summarized; no response-related column detected, so it was kept out of response comparisons.

No CellChat, LIANA, NicheNet, differential communication, final statistical testing, integration, downsampling, or full expression export was run.
No large processed object was created or uploaded.
