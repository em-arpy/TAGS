# TAGS

This repository provides R code to calculate Transcriptomic Aging Gene Scores (TAGS) as RNA versions of 5 epigenetic clocks, calculated from gene expression data.

See LICENSE.txt.


Please cite the source when using these scores:

Arpawong, T.E., Cole, S., Badhesha, H. et al. How epigenetic clocks tick: unpacking the black box by deciphering biological pathways and transcriptomic signatures of accelerated aging. npj Aging (2026). https://doi.org/10.1038/s41514-026-00446-x 


**************************************************************

Notes: 
• The code uses a gene expression matrix to calculate TAGS, where values are numeric counts
• Code assumes DEG lists and expression matrix are in the same directory

Code steps:
1. Read in data matrix of gene expression counts (gene x sample ID), then log2-transform, cpm-normalize, transform to where sample IDs are in rows and genes (by Ensemble Gene ID) are in columns
2. Extract from matrix the genes needed for each clock by ENSG ID
   -Note: if your gene matrix does not include all of the genes for each clock, we suggest imputing the base_mean value, provided in the gene lists, for each individual
3. Calculate each clock & output
