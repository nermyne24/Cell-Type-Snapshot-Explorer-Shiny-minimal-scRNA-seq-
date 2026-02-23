Overview

This Shiny web application allows users to explore cell-type specific marker genes from a simulated single-cell RNA-seq dataset.

The app:

Computes per-gene specificity statistics

Automatically selects the most specific marker gene

Visualizes marker expression on a UMAP embedding

Displays gene-level statistics dynamically

Gene Specificity Score

For a selected cell type, the specificity score is defined as:

diff = mean_in − mean_out

Where:

mean_in = average expression inside selected cell type

mean_out = average expression outside selected cell type

Marker gene selection rule:

Maximum diff

Tie-breaker → highest detection rate (det_in)

This ensures deterministic and reproducible results.
