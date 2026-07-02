# Shengxinfenxi Bioinformatics Skill

Open-source Codex skill for planning, reproducing, and auditing publication-oriented bioinformatics analyses.

The skill lives in:

`bioinformatics-reproduction/`

## What It Covers

- GEO and TCGA expression analysis
- Bulk expression preprocessing, identifier mapping, batch correction, differential expression, enrichment, immune scoring, and subtyping guardrails
- Single-cell RNA-seq support for bulk findings
- WGCNA
- Machine learning diagnostic models
- WGCNA plus machine learning diagnostic workflow guardrails
- MR, eQTL, SMR, and colocalization workflows
- Immune infiltration and enrichment analysis
- Network toxicology
- Molecular docking and molecular dynamics
- Network toxicology plus docking/MD workflow guardrails
- Spatial transcriptomics, cell-cell communication, trajectory analysis, and scTenifoldKnk virtual perturbation
- Dense `csv.gz` single-cell expression atlas first-pass audits with stratified sampling
- Figure 1-7 manuscript narrative templates and source-table contracts
- Reproduction project scaffolds with a generator script, manifest, methods log, figure manifest, validation log, limitations, and handoff templates
- Public PMC supplementary file caching, including CloudPMC proof-of-work handling for allowed public assets
- Scientific QC for manuscript-facing outputs

## Install

Copy or symlink the skill folder into your Codex skills directory:

```bash
mkdir -p ~/.codex/skills
cp -R bioinformatics-reproduction ~/.codex/skills/
```

Restart Codex after installation so the skill metadata is loaded.

## Use

Ask Codex for tasks such as:

```text
Use the bioinformatics-reproduction skill to design a GEO plus single-cell diagnostic model workflow.
```

```text
Use the bioinformatics-reproduction skill to audit a GEO/TCGA expression, enrichment, immune infiltration, and ConsensusClusterPlus subtyping workflow.
```

```text
Use the bioinformatics-reproduction skill to design the Figure 1-7 narrative for a GEO plus single-cell diagnostic biomarker manuscript.
```

```text
Use the bioinformatics-reproduction skill to create a project scaffold for a GEO plus WGCNA plus machine learning reproduction.
```

Direct script usage:

```bash
python3 bioinformatics-reproduction/scripts/create_reproduction_project.py --route-id R03 --project-name geo_wgcna_project --output-root ./projects
```

```bash
python3 bioinformatics-reproduction/scripts/download_pmc_supplementary_pow.py --url-file pmc_supplementary_urls.tsv --output-dir ./project/03_data_raw/paper/supplementary
```

```bash
Rscript bioinformatics-reproduction/scripts/stratified_dense_expression_atlas.R ./project/03_data_raw/single_cell/GSE154763 ./project/04_data_processed/GSE154763_expression_first_pass GSE154763_
```

```text
Use the bioinformatics-reproduction skill to audit a WGCNA plus machine learning paper reproduction plan.
```

```text
Use the bioinformatics-reproduction skill to prepare QC checks for MR, SMR, and colocalization results.
```

```text
Use the bioinformatics-reproduction skill to audit a network toxicology plus molecular docking and GROMACS workflow.
```

```text
Use the bioinformatics-reproduction skill to audit a single-cell plus spatial transcriptomics workflow with Monocle3 trajectory and scTenifoldKnk simulation.
```

## Scope

This repository contains generalized workflow and QC guidance. It does not include private datasets, proprietary course materials, manuscript PDFs, or local machine paths.

## License

MIT License. See `LICENSE`.
