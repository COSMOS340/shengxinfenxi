# Reproduction Project Scaffold

Use this reference when starting a new bioinformatics paper reproduction project, building a public dataset analysis project, or turning analysis outputs into a manuscript package.

The goal is to create the project scaffold before analysis starts, so data sources, sample labels, validation roles, scripts, figures, and limitations remain traceable.

## Template Assets

Copy templates from:

`assets/reproduction-project-scaffold/`

Required templates:

1. `analysis_plan.md`
2. `data_manifest.tsv`
3. `sample_metadata.tsv`
4. `methods_log.tsv`
5. `figure_manifest.tsv`
6. `validation_log.tsv`
7. `limitations.md`
8. `handoff_checklist.md`

## Startup Rule

Create the scaffold before downloading, processing, or modeling data:

1. Select one route.
2. Create the project directory tree.
3. Fill `data_manifest.tsv`.
4. Fill `sample_metadata.tsv`.
5. Write `analysis_plan.md`.
6. Lock discovery, training, test, external validation, and mechanism support roles.
7. Initialize `methods_log.tsv`, `figure_manifest.tsv`, and `validation_log.tsv`.
8. Start data download and analysis only after the previous records exist.

## Standard Project Tree

```text
project/
  00_metadata/
    data_manifest.tsv
    sample_metadata.tsv
    gene_set_manifest.tsv
    genetic_resource_manifest.tsv
    target_source_manifest.tsv
    structure_manifest.tsv
  01_plan/
    analysis_plan.md
    figure_plan.md
    limitations.md
  02_scripts/
    00_setup/
    01_data_download/
    02_preprocessing/
    03_discovery/
    04_modeling/
    05_validation/
    06_mechanism/
    07_figures/
  03_data_raw/
    geo/
    tcga/
    single_cell/
    spatial/
    gwas/
    qtl/
    compound_targets/
    disease_targets/
    structures/
    gene_sets/
  04_data_processed/
    expression/
    clinical/
    single_cell/
    spatial/
    target_lists/
    model_inputs/
  05_results/
    batch_correction/
    differential_expression/
    wgcna/
    enrichment/
    models/
    immune/
    subtype/
    mr/
    smr/
    colocalization/
    network_toxicology/
    docking/
    md/
    single_cell/
    spatial/
  06_figures/
    main/
    supplementary/
    review/
  07_tables/
    main/
    supplementary/
  99_logs/
    methods_log.tsv
    figure_manifest.tsv
    validation_log.tsv
    package_session_info.txt
    handoff_checklist.md
```

Keep the top-level folders. Remove unused subdirectories after selecting the route.

## Required Files

| File | Create before | Purpose |
|---|---|---|
| `00_metadata/data_manifest.tsv` | Data download | Record every data source, role, version, path, and group definition |
| `00_metadata/sample_metadata.tsv` | Expression analysis | Lock sample groups, batches, platforms, inclusion, and labels |
| `01_plan/analysis_plan.md` | Feature selection | State the question, route, success criteria, figure plan, and risks |
| `01_plan/limitations.md` | Manuscript drafting | Record data boundaries, method boundaries, and unsupported claims |
| `99_logs/methods_log.tsv` | Each run | Record scripts, inputs, outputs, thresholds, seeds, packages, and status |
| `99_logs/figure_manifest.tsv` | Figure generation | Bind each panel to a claim, source table, script, file, and QA note |
| `99_logs/validation_log.tsv` | Validation reporting | Record external validation, ROC, MR sensitivity, HEIDI, immune tests, and other checks |
| `99_logs/handoff_checklist.md` | Delivery | State what is complete, blocked, and ready for review |

## Manifest Schemas

### `data_manifest.tsv`

Use one row per source dataset.

Columns:

`dataset_id`, `role`, `source_database`, `accession_or_id`, `url_or_doi`, `organism`, `tissue_or_cell_type`, `platform`, `data_type`, `raw_path`, `processed_path`, `download_date`, `version_or_release`, `sample_count`, `group_definition`, `notes`

Allowed role values:

`discovery`, `training`, `test`, `external_validation`, `mechanism`, `genetic_support`, `structure_support`

### `sample_metadata.tsv`

Use one row per analyzed sample. Match `sample_id` exactly to the expression matrix column name.

Columns:

`sample_id`, `dataset_id`, `role`, `group`, `batch`, `platform`, `tissue`, `organism`, `include`, `exclude_reason`, `source_label`, `verified_label`, `notes`

Do not infer case and control labels only from file name suffixes. Read the metadata source and record the source label.

### `methods_log.tsv`

Columns:

`step_id`, `step_name`, `script_path`, `input_files`, `output_files`, `package_versions`, `parameters`, `thresholds`, `random_seed`, `start_time`, `end_time`, `status`, `notes`

Allowed status values:

`planned`, `running`, `done`, `blocked`, `review`

### `figure_manifest.tsv`

Columns:

`figure_id`, `panel_id`, `claim`, `source_table`, `source_script`, `output_file`, `route_module`, `data_role`, `qa_status`, `qa_notes`, `final_width`, `final_height`, `format`

Allowed QA status values:

`planned`, `review`, `fixed`, `pass`

No panel enters the main manuscript without a source table and a source script.

### `validation_log.tsv`

Columns:

`validation_id`, `validation_type`, `dataset_id`, `input_feature_set`, `model_or_test`, `metric`, `effect_size`, `p_value`, `adjusted_p_value`, `sample_count`, `pass_fail`, `notes`

Use this file to record every validation check that supports a Results claim.

## Route Selector

| Route | Use when | Required references |
|---|---|---|
| R01 MR/GWAS plus GEO plus machine learning | The paper needs genetic support and expression-based validation | `method-routes.md`, `bulk-geo-tcga-enrichment-immune.md`, `wgcna-machine-learning.md` |
| R02 Single-cell plus bulk diagnostic | The paper needs cell type localization for bulk biomarkers | `single-cell-spatial-trajectory.md`, `bulk-geo-tcga-enrichment-immune.md`, `wgcna-machine-learning.md` |
| R03 Comorbidity WGCNA plus theme genes | The paper starts from disease overlap, aging, pathway, or theme gene sets | `bulk-geo-tcga-enrichment-immune.md`, `wgcna-machine-learning.md` |
| R04 Treatment response WGCNA or ANN | The paper compares responders, nonresponders, baseline, or post-treatment states | `bulk-geo-tcga-enrichment-immune.md`, `wgcna-machine-learning.md`, `single-cell-spatial-trajectory.md` |
| R05 TCGA prognosis plus single-cell and spatial | The paper needs tumor prognosis and localization support | `bulk-geo-tcga-enrichment-immune.md`, `wgcna-machine-learning.md`, `single-cell-spatial-trajectory.md` |
| R06 Network toxicology plus docking, MD, and AOP | The paper focuses on a compound, exposure, target network, and structural support | `network-toxicology-docking-md.md`, `bulk-geo-tcga-enrichment-immune.md` |
| R07 Network toxicology plus multi-GEO WGCNA and SHAP | The paper combines compound targets, expression datasets, feature ranking, and docking | `network-toxicology-docking-md.md`, `bulk-geo-tcga-enrichment-immune.md`, `wgcna-machine-learning.md` |
| R08 Bulk WGCNA plus ML plus immune plus SMR | The paper needs bulk expression discovery, model validation, immune analysis, and genetic gene support | `bulk-geo-tcga-enrichment-immune.md`, `wgcna-machine-learning.md`, `method-routes.md` |

## Route Scaffolds

### R01: MR/GWAS plus GEO plus Machine Learning

Required directories:

```text
03_data_raw/gwas/
03_data_raw/geo/
03_data_raw/qtl/
04_data_processed/expression/
05_results/mr/
05_results/colocalization/
05_results/differential_expression/
05_results/models/
```

Lock before analysis:

1. Exposure GWAS and outcome GWAS.
2. Ancestry and sample overlap notes.
3. GEO discovery and validation roles.
4. Feature selection boundary.

### R02: Single-Cell plus Bulk Diagnostic

Required directories:

```text
03_data_raw/geo/
03_data_raw/single_cell/
04_data_processed/expression/
04_data_processed/single_cell/
05_results/differential_expression/
05_results/models/
05_results/single_cell/
05_results/immune/
```

Lock before analysis:

1. Single-cell sample identity and disease state.
2. Cell type annotation review rule.
3. Bulk training and validation data roles.
4. Whether single-cell evidence is used for localization, marker intersection, or trajectory.

### R03: Comorbidity WGCNA plus Theme Genes

Required directories:

```text
03_data_raw/geo/
03_data_raw/gene_sets/
04_data_processed/expression/
05_results/differential_expression/
05_results/wgcna/
05_results/enrichment/
05_results/models/
```

Lock before analysis:

1. Disease dataset and comorbidity dataset roles.
2. Theme gene set source and version.
3. DEG threshold and WGCNA module selection rule.
4. Model validation boundary.

### R04: Treatment Response WGCNA or ANN

Required directories:

```text
03_data_raw/geo/
03_data_raw/clinical/
03_data_raw/single_cell/
04_data_processed/expression/
05_results/wgcna/
05_results/models/
05_results/immune/
05_results/single_cell/
05_results/colocalization/
05_results/docking/
```

Lock before analysis:

1. Clinical response definition.
2. Baseline and post-treatment time points.
3. Prediction task and differential task as separate analyses.
4. ANN input features and scaling rule.

### R05: TCGA Prognosis plus Single-Cell and Spatial

Required directories:

```text
03_data_raw/tcga/
03_data_raw/single_cell/
03_data_raw/spatial/
04_data_processed/expression/
04_data_processed/clinical/
05_results/prognosis/
05_results/single_cell/
05_results/spatial/
```

Lock before analysis:

1. TCGA project, workflow type, and expression scale.
2. Clinical endpoint, follow-up time, and status columns.
3. Tumor or normal sample type code.
4. Single-cell and spatial roles as localization support.

### R06: Network Toxicology plus Docking, MD, and AOP

Required directories:

```text
03_data_raw/compound_targets/
03_data_raw/disease_targets/
03_data_raw/geo/
03_data_raw/structures/
04_data_processed/target_lists/
05_results/network_toxicology/
05_results/enrichment/
05_results/docking/
05_results/md/
05_results/aop/
```

Lock before analysis:

1. Compound target database list and query date.
2. Disease gene database list and score or filter rule.
3. GEO validation dataset role.
4. Docking protein and ligand source.
5. MD force field and simulation settings.
6. AOP event to evidence table.

### R07: Network Toxicology plus Multi-GEO WGCNA and SHAP

Required directories:

```text
03_data_raw/compound_targets/
03_data_raw/geo/
04_data_processed/expression/
05_results/batch_correction/
05_results/differential_expression/
05_results/wgcna/
05_results/models/
05_results/shap/
05_results/docking/
```

Lock before analysis:

1. Discovery and validation GEO datasets.
2. Batch correction inputs and outputs.
3. Feature selection data boundary.
4. SHAP model object and feature order.

### R08: Bulk WGCNA plus ML plus Immune plus SMR

Required directories:

```text
03_data_raw/geo/
03_data_raw/qtl/
03_data_raw/gwas/
04_data_processed/expression/
05_results/differential_expression/
05_results/wgcna/
05_results/models/
05_results/immune/
05_results/smr/
05_results/colocalization/
```

Lock before analysis:

1. Discovery, test, and external validation roles.
2. WGCNA trait metadata.
3. Immune scoring method and reference.
4. eQTL or pQTL tissue, GWAS ancestry, and SMR binary path.

## Script Numbering

Use one numbered script per analysis step:

| Script | Purpose |
|---|---|
| `00_setup/00_check_environment.R` | Package, path, seed, and session setup |
| `01_data_download/01_download_geo.R` | GEO or source download |
| `01_data_download/02_download_tcga.R` | TCGA download |
| `02_preprocessing/01_build_metadata.R` | Sample metadata |
| `02_preprocessing/02_process_expression.R` | Matrix and annotation |
| `02_preprocessing/03_batch_correction.R` | Multi-dataset merge and PCA |
| `03_discovery/01_differential_expression.R` | Differential expression |
| `03_discovery/02_wgcna.R` | WGCNA |
| `04_modeling/01_feature_selection.R` | LASSO, random forest, SVM-RFE, or model panel |
| `05_validation/01_external_validation.R` | Locked validation |
| `06_mechanism/01_enrichment.R` | GO, KEGG, or GSEA |
| `06_mechanism/02_immune_scores.R` | ssGSEA, GSVA, or CIBERSORT |
| `07_figures/01_main_figures.R` | Main figures |
| `07_figures/02_supplementary_figures.R` | Supplementary figures |

Each script header must state:

1. Input files.
2. Output files.
3. Dataset role.
4. Thresholds.
5. Random seed.
6. Package versions or session file.

## Copy-On-Start Checklist

1. Copy templates from `assets/reproduction-project-scaffold/`.
2. Create the standard project tree.
3. Fill `data_manifest.tsv`.
4. Fill `sample_metadata.tsv`.
5. Fill `analysis_plan.md`.
6. Initialize `methods_log.tsv`, `figure_manifest.tsv`, and `validation_log.tsv`.
7. Read the selected route in `figure-narrative-templates.md`.
8. Read the required method references.

## Delivery Checklist

Before reporting completion:

1. `data_manifest.tsv` has one row per source dataset.
2. `sample_metadata.tsv` matches expression matrix columns.
3. Raw data remain untouched.
4. Processed matrices have reproducible scripts.
5. Train, test, and external validation roles are locked.
6. Feature selection did not use external validation.
7. Each figure panel has a row in `figure_manifest.tsv`.
8. Each statistical table includes sample count and adjusted P value when relevant.
9. Final figures have visual QA notes.
10. `limitations.md` lists data and method boundaries.
11. `handoff_checklist.md` is updated for the next agent or manual review.
