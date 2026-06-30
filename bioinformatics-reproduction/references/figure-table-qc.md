# Figure and Table QC

Use this reference before delivering manuscript-facing bioinformatics figures or tables.

## Figure Export

1. Export vector PDF or SVG for line art, heatmaps, enrichment plots, forest plots, and network diagrams when practical.
2. Export high-resolution PNG or TIFF for raster-heavy plots.
3. Use consistent fonts, point sizes, and panel labels.
4. Keep color palettes interpretable in grayscale and color-vision-deficiency views when the figure carries categorical information.
5. Save each panel with a reproducible script path and source table path.

## Visual Inspection

Open each figure and check:

1. Axis labels are not clipped.
2. Legends stay inside the canvas or have enough margin.
3. Text is readable at final manuscript size.
4. Long gene names and pathway names do not overlap.
5. Heatmap row names and column annotations remain aligned.
6. ROC curves show class direction, AUC, and confidence interval when used.
7. Boxplots show statistical test, sample count, and group labels.
8. Network diagrams have readable labels or a separate node table.

## Table Export

For each table, include:

1. Source script.
2. Source dataset.
3. Filtering rule.
4. Column definitions.
5. Identifier namespace.
6. Units where relevant.
7. Full unfiltered table when a filtered table is shown in the manuscript.

## Manuscript Figure Set

A bioinformatics reproduction article usually needs:

1. Workflow overview.
2. Dataset QC and grouping evidence.
3. Discovery result.
4. Feature or module selection result.
5. Validation result.
6. Mechanistic interpretation.
7. Orthogonal evidence when used.
8. Final schematic only after the data panels support it.

## Completion Rule

Do not report a figure package as complete until:

1. Every figure opens.
2. Every panel maps to a source table.
3. Every statistical annotation is traceable.
4. Every key label is readable.
5. The final figure directory contains no placeholder plots.
