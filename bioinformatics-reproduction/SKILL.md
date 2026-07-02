---
name: bioinformatics-reproduction
description: Plan, reproduce, and audit publication-oriented bioinformatics analyses. Use when Codex needs to design or execute workflows for GEO, TCGA, single-cell RNA-seq, WGCNA, machine learning diagnostic models, Mendelian randomization, eQTL, SMR, colocalization, immune infiltration, enrichment analysis, network toxicology, molecular docking, molecular dynamics, spatial transcriptomics, cell-cell communication, trajectory analysis, scTenifoldKnk virtual perturbation, or paper reproduction from a manuscript-style bioinformatics study.
---

# Bioinformatics Reproduction

## Role

Use this skill to turn a bioinformatics paper idea or manuscript workflow into a reproducible analysis plan, implementation scaffold, and scientific QC checklist.

The skill is designed for publication-facing work: it prioritizes traceable data sources, explicit sample definitions, validation separation, figure quality, and method choices that can be defended in a paper or response letter.

## First Reads

Open only the reference files needed for the task:

| Task | Read |
|---|---|
| Build a full paper reproduction plan | `references/article-reproduction-playbook.md` |
| Start a new reproduction project scaffold with manifests, logs, templates, and route-specific folders | `references/reproduction-project-scaffold.md` |
| Download and cache public PMC supplementary files, including CloudPMC proof-of-work assets | `references/literature-data-downloads.md` |
| Run or audit a dense `csv.gz` single-cell expression atlas first pass | `references/dense-expression-single-cell-atlas.md` |
| Design manuscript main figures, supplementary figures, or a Figure 1-7 narrative | `references/figure-narrative-templates.md` |
| Choose a method route | `references/method-routes.md` |
| Design or audit bulk GEO/TCGA expression, enrichment, immune scoring, or subtyping workflows | `references/bulk-geo-tcga-enrichment-immune.md` |
| Design or audit WGCNA plus machine learning diagnostic workflows | `references/wgcna-machine-learning.md` |
| Design or audit network toxicology, molecular docking, or molecular dynamics workflows | `references/network-toxicology-docking-md.md` |
| Design or audit single-cell, spatial transcriptomics, trajectory, communication, or virtual perturbation workflows | `references/single-cell-spatial-trajectory.md` |
| Audit analysis validity | `references/scientific-qc.md` |
| Prepare final figures and tables | `references/figure-table-qc.md` |

## Operating Rules

1. Treat every course script, prior notebook, or online tutorial as implementation evidence, not as scientific authority.
2. Verify exact dataset accessions, phenotype labels, group definitions, genome builds, gene identifiers, and platform annotations before analysis.
3. Keep raw inputs immutable. Write cleaned data, derived matrices, figures, logs, and reports into a project-specific output folder.
4. Record package versions, random seeds, filtering thresholds, model formulas, feature lists, and output file paths.
5. Separate discovery, training, validation, and external validation datasets before feature selection unless the study design explicitly justifies a different order.
6. Do not report diagnostic or prognostic performance until leakage checks, class-direction checks, and external validation checks are complete.
7. Visually inspect every generated figure before reporting completion.

## Workflow

1. Define the biological question.
   - Disease or phenotype.
   - Exposure, treatment, comorbidity, compound, pathway, or cell type focus.
   - Primary claim the paper will support.

2. Lock the data design.
   - List datasets with accession, organism, platform, tissue, sample count, and case/control or outcome definition.
   - Decide which dataset is discovery, training, internal validation, external validation, single-cell support, spatial support, or genetic-association support.

3. Choose the route.
   - Read `references/method-routes.md`.
   - Pick the smallest route that supports the claim.
   - Add method modules only when they answer a defined biological or validation question.

4. Build the reproducible scaffold.
   - Read `references/reproduction-project-scaffold.md`.
   - Run `scripts/create_reproduction_project.py` when creating a project on disk.
   - Copy the required templates from `assets/reproduction-project-scaffold/` when manual setup is requested.
   - Use directories such as `00_metadata/`, `01_plan/`, `02_scripts/`, `03_data_raw/`, `04_data_processed/`, `05_results/`, `06_figures/`, `07_tables/`, and `99_logs/`.
   - Keep one numbered script per major step.
   - Make script inputs and outputs explicit at the top of each file.

5. Run analysis with checkpoints.
   - Validate sample metadata after every merge.
   - Save intermediate objects and tables.
   - Write a short log for thresholds, removed samples, removed genes, and model settings.

6. Audit scientific validity.
   - Read `references/scientific-qc.md`.
   - Fix data leakage, weak validation, label mismatch, batch artifacts, and unsupported claims before preparing figures.

7. Prepare publication outputs.
   - Read `references/figure-table-qc.md`.
   - Export vector figures when practical.
   - Inspect labels, legends, panel sizes, colors, and statistical annotations.

## Report Back

For completed work, report:

1. Data sources used and exact accessions.
2. Route selected and method modules executed.
3. Scripts, outputs, figures, and tables created.
4. Validation design and leakage checks.
5. Scientific risks that remain.
6. Files that should be opened for manual review before manuscript use.
