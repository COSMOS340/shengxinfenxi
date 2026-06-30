# Article Reproduction Playbook

Use this reference when the task is to reproduce a manuscript-style bioinformatics paper or design a new paper route from public datasets.

## Paper Logic First

Write a one-paragraph analysis thesis before coding:

1. What biological or clinical question is being answered.
2. Which public datasets provide discovery evidence.
3. Which independent datasets validate the result.
4. Which orthogonal evidence layer supports mechanism, causality, localization, or therapeutic relevance.
5. Which final claim the figures will defend.

## Common Storyline Patterns

### GEO plus Machine Learning Diagnosis

Use for disease classification or biomarker discovery in bulk expression data.

Core chain:

1. GEO discovery data preprocessing.
2. Metadata-verified case/control labels.
3. Differential expression and optional WGCNA.
4. Feature selection using LASSO, random forest, SVM-RFE, or multi-model ranking.
5. Diagnostic ROC in training data.
6. Independent validation dataset ROC.
7. Nomogram or calibration only when a clinical prediction framing is justified.

### Single-Cell plus Bulk Validation

Use when the paper needs cell-type localization for bulk-derived genes.

Core chain:

1. Bulk expression identifies disease-associated genes.
2. Single-cell data identifies cell types and marker genes.
3. Intersect bulk genes with cell-type markers or test expression in annotated cell types.
4. Validate expression and ROC in independent bulk datasets.
5. Add pathway, immune, communication, or trajectory analysis only when tied to the cell-type mechanism.

### WGCNA plus Theme Gene Set

Use when the paper starts from a pathway set, aging genes, comorbidity genes, compound targets, or disease-gene database.

Core chain:

1. Define the theme gene set from a named source.
2. Run differential expression on phenotype groups.
3. Run WGCNA on expression data with metadata-verified traits.
4. Intersect differential genes, module genes, and the theme gene set.
5. Perform enrichment and diagnostic validation.
6. Keep WGCNA module selection and intersection logic fixed before model evaluation.

### Mendelian Randomization and Colocalization

Use when the paper claims genetic support for exposure, gene, protein, metabolite, or microbiome effects.

Core chain:

1. Select exposure GWAS, eQTL, pQTL, mQTL, or microbiome QTL sources.
2. Select outcome GWAS with matched ancestry and clear phenotype.
3. Clump instruments, harmonize alleles, and filter weak instruments.
4. Run MR, heterogeneity, pleiotropy, leave-one-out, and sensitivity analyses.
5. Add colocalization or SMR/HEIDI for gene-level support.
6. Report ancestry, sample overlap, instrument strength, and multiple-testing control.

### Network Toxicology plus Docking or Molecular Dynamics

Use for compound-disease mechanism papers.

Core chain:

1. Collect compound targets from named databases.
2. Collect disease genes from named disease resources and expression datasets.
3. Intersect compound targets, disease genes, differential genes, and WGCNA module genes when available.
4. Build PPI and hub-gene ranking.
5. Run enrichment and pathway narrative.
6. Dock compound-protein pairs and add molecular dynamics only for a small final set.
7. Treat docking and simulation as structural support, not proof of clinical efficacy.

## Figure Narrative

For detailed Figure 1-7 templates and source-table contracts, read `figure-narrative-templates.md`.

Use figures to advance one claim at a time:

1. Dataset and workflow overview.
2. Discovery differential analysis or module discovery.
3. Feature-selection evidence.
4. Validation evidence.
5. Mechanism evidence.
6. External or orthogonal support.
7. Final model or integrated schematic.

## Minimum Handoff Package

For a complete directory scaffold and copy-ready templates, read `reproduction-project-scaffold.md` and use `assets/reproduction-project-scaffold/`.

For a complete reproduction package, create:

1. `data_manifest.tsv` with accession, source, platform, tissue, groups, and role.
2. `analysis_plan.md` with the chosen route and success criteria.
3. Numbered scripts with fixed inputs and outputs.
4. `methods_log.tsv` with thresholds, seeds, packages, and model settings.
5. `figure_manifest.tsv` with source tables and scripts for each panel.
6. `limitations.md` with validation gaps and scientific boundaries.
