# Mac review request: GSE189926 resolution and UMAP topology sweep

## Separate decisions

- Resolution evidence rank 1: `0.5` (clusters=19, mean adjacent ARI=0.9494, broad-marker conflict flags=2).
- UMAP evidence rank 1: `neighbor30` (q90 local linearity=2.507, mean same-sample neighbor fraction=0.1405, centroid/dispersion ratio=8.093).
- These ranks use separate, fully exported evidence tables. No choice was made from visual compactness alone.

## Filament review

- Baseline all-cell q99 local-linearity rule flagged 896 cells; no cells were removed.
- `gse189926_filament_region_review.tsv` records sample/patient dominance, QC associations, broad-label mixing, and affected cluster Top50 markers.
- The local-linearity metric is a geometry diagnostic and is not used as a biological quality or exclusion score.

## Requested decision

1. Review whether resolution `0.5` should be retained as the formal resolution after inspecting marker support and transitions.
2. Review whether UMAP route `neighbor30` should replace the current display embedding after inspecting biology, sample mixing, QC overlays, and filament evidence.
3. Do not update the formal GSE189926 object until both decisions are approved.
