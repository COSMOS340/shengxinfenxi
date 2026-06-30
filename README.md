# Shengxinfenxi Bioinformatics Skill

Open-source Codex skill for planning, reproducing, and auditing publication-oriented bioinformatics analyses.

The skill lives in:

`bioinformatics-reproduction/`

## What It Covers

- GEO and TCGA expression analysis
- Single-cell RNA-seq support for bulk findings
- WGCNA
- Machine learning diagnostic models
- MR, eQTL, SMR, and colocalization workflows
- Immune infiltration and enrichment analysis
- Network toxicology
- Molecular docking and molecular dynamics
- Spatial transcriptomics, cell-cell communication, and trajectory analysis
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
Use the bioinformatics-reproduction skill to audit a WGCNA plus machine learning paper reproduction plan.
```

```text
Use the bioinformatics-reproduction skill to prepare QC checks for MR, SMR, and colocalization results.
```

## Scope

This repository contains generalized workflow and QC guidance. It does not include private datasets, proprietary course materials, manuscript PDFs, or local machine paths.

## License

MIT License. See `LICENSE`.
