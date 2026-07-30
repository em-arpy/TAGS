# TAGS

This repository provides R code and gene names to calculate Transcriptomic Aging Gene Scores (TAGS) as RNA versions of 5 epigenetic clocks, calculated from gene expression data.

For license information, see LICENSE.


Please cite the source when using these scores:

Arpawong, Thalida Em, Steve Cole, Harshanna Badhesha, Jung Ki Kim, Christopher R. Beam, Eric T. Klopack, Kimberly Siegmund, Bharat Thyagarajan, and Eileen M. Crimmins. "How epigenetic clocks tick: Unpacking the black box by deciphering biological pathways and transcriptomic signatures of accelerated aging." npj Aging (2026). https://doi.org/10.1038/s41514-026-00446-x 


**************************************************************

Notes: <br>
• The code uses a gene expression matrix to calculate TAGS, where values are numeric counts <br>
• Code assumes the folder containing 5 DEG lists and expression matrix are in the same directory<br>
<br>
Code steps:
1. Read in data matrix of gene expression counts (gene x sample ID), log2-transform, cpm-normalize, then transform the matrix to where sample IDs are in rows and genes (by Ensemble Gene ID) are in columns (sample ID x gene)
2. Extract from matrix the genes needed for each clock by ENSG ID
   -Note: if your gene matrix does not include all of the genes for each clock, we suggest imputing the base_mean value, provided in the gene lists, for each individual
3. Calculate each clock & output scores to a csv file
