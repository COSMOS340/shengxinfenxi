# Network Toxicology, Docking, and Molecular Dynamics

Use this reference when a reproduction project studies a compound, drug, metabolite, environmental exposure, food additive, contaminant, or small molecule and connects it to a disease mechanism through target databases, transcriptomics, docking, molecular dynamics, or an AOP narrative.

## Route Selection

| Study aim | Route |
|---|---|
| Compound-disease mechanism | compound identity -> compound targets -> disease targets -> target integration -> enrichment -> validation |
| Structural validation of hub genes | receptor/ligand preparation -> pocket definition -> AutoDock Vina -> interaction figure |
| Dynamic stability after docking | docked complex -> GROMACS setup -> equilibration -> production simulation -> trajectory analysis |
| Toxicology story | toxicity/ADMET -> targets -> pathways -> immune or cell context -> adverse outcome pathway |
| Stronger biological support | add GEO/TCGA, single-cell, MR, or experimental validation after target integration |

## Minimum Project Scaffold

```text
00_manifest/
01_compound_identity/
02_toxicity_admet/
03_compound_targets/
04_disease_targets/
05_target_integration/
06_ppi_hub/
07_enrichment/
08_transcriptome_validation/
09_structural_docking/
10_molecular_dynamics/
11_aop_narrative/
12_figures_tables/
99_logs/
```

## Required Manifests

| File | Required fields |
|---|---|
| `compound_manifest.tsv` | compound name, CAS, PubChem CID, SMILES or InChIKey, structure source, retrieval date |
| `compound_target_manifest.tsv` | database, version or retrieval date, species, query field, score field, threshold, exported file |
| `disease_target_manifest.tsv` | disease term, synonyms, database, score or relevance rule, retrieval date, exported file |
| `ppi_manifest.tsv` | network database, organism, confidence cutoff, interaction evidence types, software version |
| `docking_manifest.tsv` | receptor ID, chain, resolution, ligand file, receptor preparation, grid center, box size, Vina version, exhaustiveness |
| `md_manifest.tsv` | GROMACS version, force field, ligand parameter method, water model, ion settings, MDP files, production length |
| `figure_manifest.tsv` | source table, plotting script, output file, visual inspection status |

## Compound and Disease Target Integration

Minimum route:

1. Lock compound identity using CAS, PubChem CID and canonical structure string.
2. Export compound-related targets from named databases such as ChEMBL, STITCH, Swiss Target Prediction, PharmMapper, CTD, DrugBank or BindingDB.
3. Export disease-related genes from named databases such as GeneCards, OMIM, DisGeNET, CTD or curated disease resources.
4. Normalize gene symbols using an authoritative annotation source.
5. Remove blank identifiers, synonyms that map to the same gene, and duplicated records.
6. Save each source-specific target list before merging.
7. Use union within compound-target sources when the aim is broad mechanism discovery.
8. Use intersection between compound targets and disease genes only after the source lists are frozen.
9. Store every intermediate list and the final integrated target list.

QC focus:

1. Database names alone are not enough. Record query date, organism, thresholds and exported columns.
2. Disease target resources use different scoring rules. Do not mix score columns without documenting the rule.
3. Gene-symbol conversion can merge synonyms and remove older identifiers. Save the mapping table.
4. Directory-wide file loading is fragile. Use an explicit input manifest instead of reading every `.txt` or `.csv` in a folder.
5. A database intersection is not mechanism proof. Treat it as a prioritization step.

## PPI and Hub Ranking

Minimum route:

1. Submit the frozen integrated target list to STRING, GeneMANIA or another declared interaction resource.
2. Record organism, confidence cutoff and evidence channels.
3. Export node and edge tables.
4. If using Cytoscape, record Cytoscape version and app versions.
5. If using cytoHubba, record each algorithm and the final intersection or ranking rule.
6. If using MCODE, record all MCODE parameters.
7. Save the hub-gene table before enrichment, docking or machine learning.

QC focus:

1. Different hub algorithms answer different network questions.
2. A top-degree hub is not automatically the best docking receptor or biological driver.
3. PPI confidence thresholds change the hub list. Report them.

## Enrichment Analysis

Minimum route:

1. Use the frozen target or hub list as input.
2. Convert symbols to ENTREZ or Ensembl IDs with a saved mapping table.
3. Run GO and KEGG with declared packages and database dates.
4. Apply multiple-testing correction.
5. Save full enrichment results before selecting plotted terms.
6. Record the rule used to choose terms for figures.

QC focus:

1. Enrichment requires a defined tested gene universe.
2. Small target lists create unstable top terms.
3. Pathway enrichment is pathway annotation evidence, not direct pathway activation.
4. Do not add new pathway claims in the text that are absent from the enrichment table.

## Transcriptomic and Genetic Validation Layers

Use validation layers when the paper needs more than a database-derived mechanism.

| Layer | Use | Guardrail |
|---|---|---|
| GEO/TCGA differential analysis | Show disease expression context | verify group labels from metadata |
| WGCNA | Link genes to disease traits or modules | choose modules before downstream model evaluation |
| Machine learning | Build diagnostic or prioritization models | split training and validation before feature selection |
| Immune infiltration | Interpret microenvironment signatures | correct across many immune cell scores |
| Single-cell | Localize genes to cell types | verify annotation with markers |
| MR/pQTL-MR | Add genetic direction evidence | check instruments, ancestry, harmonization and pleiotropy |

Read the dedicated WGCNA/machine-learning and MR/SMR references when using those layers.

## Molecular Docking

Minimum route:

1. Lock the ligand identity and structure source.
2. Generate or download 3D ligand structure.
3. Prepare ligand file formats such as `mol2` and `pdbqt`.
4. Lock receptor source: PDB ID, chain, resolution, mutations, cofactors and missing residues.
5. Prepare receptor structure and record water, ion, cofactor and native-ligand handling.
6. Define binding pocket and grid center from a known active site, co-crystal ligand, prediction tool or documented residue set.
7. Run docking with recorded Vina version and parameters.
8. Save config, log, input structures and output structures.
9. Extract docking score from the docking log or output file.
10. Draw interaction diagrams from the selected pose.

AutoDock Vina fields to record:

```text
receptor =
ligand =
center_x =
center_y =
center_z =
size_x =
size_y =
size_z =
exhaustiveness =
energy_range =
num_modes =
```

QC focus:

1. Record protonation, tautomer, charge and rotatable-bond handling.
2. Search boxes that are too large reduce interpretability and trigger Vina warnings.
3. Low exhaustiveness reduces search confidence.
4. A docking score is not an experimental binding affinity.
5. Do not compare scores across receptors without controlling preparation and grid settings.
6. Preserve raw docking logs. Never type scores manually into the final table without traceable extraction.

## Molecular Dynamics

Minimum GROMACS route:

1. Build receptor-ligand complex.
2. Build topology for protein and ligand.
3. Place the complex in a simulation box.
4. Solvate the system.
5. Add ions and neutralize.
6. Run energy minimization.
7. Run NVT equilibration.
8. Run NPT equilibration.
9. Run production simulation.
10. Correct periodic boundary conditions.
11. Analyze RMSD, RMSF, radius of gyration, SASA and hydrogen bonds.
12. Build a free-energy landscape only after defining the variables used.

Command families to preserve:

```bash
gmx editconf
gmx solvate
gmx grompp
gmx genion
gmx mdrun
gmx genrestr
gmx make_ndx
gmx trjconv
gmx rms
gmx rmsf
gmx gyrate
gmx sasa
gmx hbond
gmx sham
```

MDP fields to record:

```text
integrator
nsteps
dt
constraints
coulombtype
rcoulomb
rvdw
tcoupl
pcoupl
gen_vel
```

Report production time from:

```text
production_time_ps = nsteps * dt
```

QC focus:

1. Force-field choice and ligand parameterization dominate MD validity.
2. Save all warnings from `grompp`; do not use `-maxwarn` without explaining each warning.
3. Energy minimization, NVT and NPT stability must be checked before production.
4. Production length must match the numerical MDP settings, not comments in the file.
5. RMSD/RMSF/gyration curves are stability diagnostics, not disease mechanism proof.
6. A free-energy landscape needs the variables, input file and conversion script reported.

## AOP Narrative

AOP can organize toxicology logic after the analytical evidence is complete.

Minimum AOP fields:

1. Compound or exposure.
2. Molecular initiating event.
3. Key genes or proteins.
4. Pathways.
5. Immune or cellular events.
6. Adverse outcome or disease phenotype.
7. Evidence source for each arrow.

QC focus:

1. AOP is a narrative evidence map, not a statistical test.
2. Each arrow needs a literature source or a project result.
3. Do not introduce genes, pathways or immune cells that are absent from upstream tables.
4. Keep speculative wording out of final claims.

## Common Failure Modes

1. Database exports have no retrieval date.
2. Compound identity is ambiguous.
3. Disease gene thresholds are not recorded.
4. Gene symbol conversion is not saved.
5. Target intersections are treated as causal evidence.
6. PPI and hub rankings are described without parameters.
7. Enrichment results are cherry-picked without the full table.
8. Docking scores are typed into spreadsheets by hand.
9. Vina warnings are ignored.
10. Protein preparation steps are not recorded.
11. Ligand charge and protonation are not recorded.
12. GROMACS `-maxwarn` is used without warning logs.
13. MD duration is copied from comments rather than calculated from `nsteps * dt`.
14. AOP diagrams add unsupported claims.
15. Figures are exported without visual inspection.

## Completion Checklist

Before reporting a network toxicology, docking or MD reproduction as complete:

1. Data manifests exist for compound, target databases, disease targets, PPI, docking and MD.
2. Every source table has a retrieval date or version.
3. Integrated target lists are reproducible from saved inputs.
4. Enrichment has full tables and correction columns.
5. Docking config, logs and structures are preserved.
6. Vina warnings have been resolved or documented.
7. MD command history and MDP files are preserved.
8. Production simulation length is calculated from numeric parameters.
9. AOP arrows point to upstream tables or literature sources.
10. Final figures pass visual inspection.

