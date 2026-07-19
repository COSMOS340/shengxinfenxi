# Mac review request: GSE205506 source-paper QC rebuild

Please review these items before using this object for CRC Figure 1:

1. Compare observed mitochondrial-removal fractions with the article values in `gse205506_paper_qc_removal_fraction_comparison.tsv`; parameters were not tuned to the article total.
2. Review first-round and final broad evidence, especially conflicts listed in `gse205506_paper_qc_hybrid_review.tsv` and the two broad-cluster evidence tables.
3. Review every `Unresolved_review` compartment cluster in the six annotation-evidence tables; no author subtype was forced.
4. Review the bright-palette UMAPs and marker dotplots in `figures/`.
5. Confirm the final cell total and sample/subject traceability from `gse205506_paper_qc_final_object_summary.tsv` and `gse205506_paper_qc_final_counts.tsv`.

The RDS objects are local-only on F: and are fully identified by path, size, and SHA-256 in `gse205506_paper_qc_local_object_inventory.tsv`.

