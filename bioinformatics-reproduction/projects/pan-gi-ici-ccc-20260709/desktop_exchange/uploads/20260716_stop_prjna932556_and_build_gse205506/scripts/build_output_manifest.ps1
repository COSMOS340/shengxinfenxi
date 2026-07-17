param(
    [Parameter(Mandatory = $true)]
    [string]$UploadDirectory
)

$ErrorActionPreference = 'Stop'
$uploadDirectory = (Resolve-Path -LiteralPath $UploadDirectory).Path

$githubFiles = @(
    'STATUS.md',
    'MAC_REVIEW_REQUEST.md',
    'prjna932556_deletion_manifest_before.tsv',
    'prjna932556_deletion_audit_after.tsv',
    'prjna932556_deletion_summary.tsv',
    'prjna932556_storage_audit_before.tsv',
    'prjna932556_storage_audit_after.tsv',
    'prjna932556_storage_audit_final.tsv',
    'gse205506_input_inventory.tsv',
    'gse205506_tar_inventory.tsv',
    'gse205506_geo_sample_metadata.tsv',
    'gse205506_response_mapping.tsv',
    'gse205506_response_join_audit.tsv',
    'gse205506_response_join_summary.tsv',
    'gse205506_matrix_structure_audit.tsv',
    'gse205506_qc_summary.tsv',
    'gse205506_qc_sensitivity_audit.tsv',
    'gse205506_qc_filter_summary.tsv',
    'gse205506_qc_filter_sensitivity_audit.tsv',
    'gse205506_doublet_audit.tsv',
    'gse205506_environment_rna_audit.tsv',
    'gse205506_batch_structure_knn_audit.tsv',
    'gse205506_batch_structure_pc_eta_squared.tsv',
    'gse205506_clustering_resolution_summary.tsv',
    'gse205506_cluster_top10_markers.tsv',
    'gse205506_cluster_top50_markers.tsv',
    'gse205506_annotation_evidence.tsv',
    'gse205506_cluster_qc_and_reviewed_annotation.tsv',
    'gse205506_author_marker_cluster_summary.tsv',
    'gse205506_author_marker_reference_top10.tsv',
    'gse205506_author_marker_reference_top50.tsv',
    'gse205506_author_marker_overlap_top50_non_nuisance_top3.tsv',
    'gse205506_author_marker_overlap_expected_compartment_top3.tsv',
    'gse205506_patient_sample_cluster_cell_counts.tsv',
    'gse205506_analysis_parameters.tsv',
    'gse205506_qc_singlet_object_summary.tsv',
    'gse205506_integrated_object_summary.tsv',
    'gse205506_final_object_summary.tsv',
    'gse205506_reviewed_object_summary.tsv',
    'sessionInfo_inventory_basic_qc.txt',
    'sessionInfo_response_and_author_markers.txt',
    'sessionInfo_qc_singlet_object.txt',
    'sessionInfo_normalize_integrate_cluster.txt',
    'sessionInfo_markers_annotation_figures.txt',
    'sessionInfo_review_annotation_and_author_markers.txt',
    '02_parse_table_s1_and_author_markers.log',
    '03a_build_qc_singlet_seurat.log',
    '03b_normalize_integrate_cluster.stderr.log',
    '03b_normalize_integrate_cluster.stdout.log',
    '03c_markers_annotation_figures.stderr.log',
    '03d_review_annotation_and_author_markers.stderr.log',
    'scripts/01_gse205506_inventory_basic_qc.R',
    'scripts/02_parse_table_s1_and_author_markers.R',
    'scripts/03a_build_qc_singlet_seurat.R',
    'scripts/03b_normalize_integrate_cluster.R',
    'scripts/03c_markers_annotation_figures.R',
    'scripts/03d_review_annotation_and_author_markers.R',
    'scripts/manage_prjna932556_deletion.ps1',
    'scripts/build_output_manifest.ps1',
    'scripts/upload_github_handoff.ps1',
    'figures/gse205506_qc_sample_retention.png',
    'figures/gse205506_qc_sample_retention.pdf',
    'figures/gse205506_umap_working_clusters_reviewed.png',
    'figures/gse205506_umap_working_clusters_reviewed.pdf',
    'figures/gse205506_umap_broad_cell_type_reviewed.png',
    'figures/gse205506_umap_broad_cell_type_reviewed.pdf',
    'figures/gse205506_umap_cluster_quality_reviewed.png',
    'figures/gse205506_umap_cluster_quality_reviewed.pdf',
    'figures/gse205506_umap_timepoint_reviewed.png',
    'figures/gse205506_umap_timepoint_reviewed.pdf',
    'figures/gse205506_umap_response_reviewed.png',
    'figures/gse205506_umap_response_reviewed.pdf',
    'figures/gse205506_marker_dotplot_working_clusters_reviewed.png',
    'figures/gse205506_marker_dotplot_working_clusters_reviewed.pdf'
)

$rows = foreach ($relativePath in $githubFiles) {
    $absolutePath = Join-Path $uploadDirectory ($relativePath -replace '/', '\')
    if (-not (Test-Path -LiteralPath $absolutePath -PathType Leaf)) {
        throw "Manifest input is missing: $absolutePath"
    }
    $item = Get-Item -LiteralPath $absolutePath
    [pscustomobject]@{
        relative_path = $relativePath
        absolute_path = $item.FullName
        size_bytes = $item.Length
        sha256 = (Get-FileHash -LiteralPath $absolutePath -Algorithm SHA256).Hash
        delivery = 'github_upload'
        notes = 'lightweight reproducibility output'
    }
}

$localOnly = @(
    @{ Path = 'F:\pan-gi-ici-ccc-20260709\GSE205506\raw\GSE205506_RAW.tar'; Note = 'raw GEO archive; local only' },
    @{ Path = 'F:\pan-gi-ici-ccc-20260709\GSE205506\supplement\mmc2.xlsx'; Note = 'user-provided Table S1 workbook; local only' },
    @{ Path = 'F:\pan-gi-ici-ccc-20260709\GSE205506\supplement\mmc3.xlsx'; Note = 'user-provided author marker workbook; local only' },
    @{ Path = 'F:\pan-gi-ici-ccc-20260709\GSE205506\r_analysis\gse205506_seurat_final_reviewed.rds'; Note = 'final reviewed Seurat object; local only' }
)

foreach ($entry in $localOnly) {
    $item = Get-Item -LiteralPath $entry.Path
    $rows += [pscustomobject]@{
        relative_path = ''
        absolute_path = $item.FullName
        size_bytes = $item.Length
        sha256 = (Get-FileHash -LiteralPath $item.FullName -Algorithm SHA256).Hash
        delivery = 'local_only'
        notes = $entry.Note
    }
}

$manifestPath = Join-Path $uploadDirectory 'output_manifest.tsv'
$rows | Export-Csv -LiteralPath $manifestPath -Delimiter "`t" -NoTypeInformation -Encoding utf8

$manifest = Import-Csv -LiteralPath $manifestPath -Delimiter "`t"
if (($manifest | Where-Object delivery -eq 'github_upload').Count -ne $githubFiles.Count) {
    throw 'Manifest GitHub file count verification failed'
}

[pscustomobject]@{
    manifest_path = $manifestPath
    github_files = $githubFiles.Count
    local_only_files = $localOnly.Count
    manifest_rows = $manifest.Count
}
