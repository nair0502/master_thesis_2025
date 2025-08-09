# master_thesis_2025
TRR356 Project title: Data Mining for Pathogenic and Symbiotic Plant-Microbe Interaction signatures Using Co-Expression Networks
## Aim: Pathway Enrichment Analysis in *Solanum lycopersicum* Under Different Stress and Interaction Conditions

## Overview

This project explores the gene expression and pathway enrichment patterns in *Solanum lycopersicum* (tomato) under various biological conditions including:

- **Biotrophic interaction**
- **Symbiotic interaction**
- **Necrotrophic stress** (control)
- **Drought stress** (control)

The goal is to identify both **unique** and **shared** molecular pathways activated during biotrophic and symbiotic interactions, providing insight into plant responses and potential cross-talk between beneficial and pathogenic relationships.

## Research Approach

To uncover biologically relevant gene modules and their functional significance, two co-expression inference strategies were employed:

- **Weighted Gene Co-expression Network Analysis (WGCNA)** – identifies modules of co-expressed genes based on correlation and hierarchical clustering.
- **GENIE3** – a machine learning-based method using random forests to infer gene regulatory relationships. Since GENIE3 does not include clustering, **Louvain community detection** was applied to the inferred network to group genes into modules.

To assess the similarity between modules inferred by WGCNA and GENIE3, **cosine similarity** was computed on their respective adjacency matrices, allowing for direct comparison of co-expression structures.

Functional annotation and pathway enrichment were performed using tools such as **Mercator** and **clusterProfiler**, enabling biological interpretation of gene modules.

## Tools and Workflow

The analysis pipeline included the following tools:

- **Quality Control & Preprocessing**
  - [`fastp`](https://github.com/OpenGene/fastp): read trimming and quality filtering
  - [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): read quality checks

- **Transcript Quantification**
  - [`Kallisto`](https://pachterlab.github.io/kallisto/): pseudo-alignment and transcript quantification

- **Normalization & Co-expression Analysis**
  - Variance Stabilizing Transformation (**VST**)
  - [`WGCNA`](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/): module detection and correlation
  - [`GENIE3`](https://bioconductor.org/packages/release/bioc/html/GENIE3.html): random forest-based co-expression inference
  - [`igraph`](https://igraph.org/r/): Louvain clustering for GENIE3 network

- **Functional Annotation**
  - [`clusterProfiler`](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html): pathway enrichment analysis
  - [`Mercator`](https://www.plabipd.de/portal/mercator-sequence-annotation): for assigning genes to functional categories

- **Comparison & Evaluation**
  - Cosine similarity analysis between WGCNA and GENIE3 adjacency matrices
