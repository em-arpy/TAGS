#Calculate Transcriptomic Aging Gene Scores (TAGS) as gene expression proxies
#of Age Acceleration Epigenetic Clocks

#Source: Arpawong, T.E., Cole, S., Badhesha, H. et al. How epigenetic clocks tick: unpacking the black box by deciphering biological pathways and transcriptomic signatures of accelerated aging. npj Aging (2026). https://doi.org/10.1038/s41514-026-00446-x 

###############################################################

#Important Notes:
#Gene expression is in matrix format, values are numeric counts
#Code assumes "gene_lists" folder containing DEG lists and expression matrix is in the same directory
#Change 2 things below, indicated by ***: 
#1. the path to your working directory
#2. name of your gene expression file


#Code steps:
#1. Read in data matrix of gene expression counts (gene x sample ID), then log2-transform, cpm-normalize, transform to where sample IDs are in rows and genes (by Ensemble Gene ID) are in columns
#2. Extract from matrix the genes needed for each clock by ENSG ID
#   -Note: if your gene matrix does not include all of the genes for each clock, we suggest imputing 
#   the base_mean value, provided in the gene lists, for each individual
#3. Calculate each clock & output


###############################################################

#Install libraries if necessary

library(dplyr)
library(edgeR)

###############################################################
#### Read in matrix data of gene expression

#***Change directory
setwd("<your directory>")

#load matrix of expression counts with rownames as gene IDs
#***Change name of file
rsem<-read.table('<name of gene expression matrix file>', header=T)


# normalize
dge <- DGEList(counts = rsem)
dge <- calcNormFactors(dge)
# compute log2CPM
log2cpm <- cpm(dge, log = TRUE, prior.count = 2)
log2cpm<-as.data.frame(log2cpm)
rownames(log2cpm)<-rownames(rsem)
tlog2cpm<-t(log2cpm)
rm(rsem,dge,log2cpm)


##############################################################
# Read DEG files
##############################################################
clocks <- c(
  GrimAge = "GrimAge",
  Horvath = "Horvath",
  Hannum = "Hannum",
  PhenoAge = "PhenoAge",
  DunedinPACE = "DunedinPACE"
)

deg_list <- lapply(clocks, function(x) {
  read.csv(
    file.path("./gene_lists/",
              paste0("Significant_DEG_", x, "_fdr01.csv")),
    header = TRUE
  )
})

##############################################################
# Function to calculate weighted TAGS
##############################################################
calc_score <- function(deg, expr) {
  
  idx <- match(deg$geneID, colnames(expr))
  keep <- !is.na(idx)
  
  weighted <- sweep(
    expr[, idx[keep], drop = FALSE],
    2,
    deg$coeff_log2[keep],
    `*`
  )
  
  rowMeans(weighted, na.rm = TRUE)
}

##############################################################
# Calculate scores & output
##############################################################
scores <- lapply(names(deg_list), function(clock) {
  calc_score(deg_list[[clock]], tlog2cpm)
})

names(scores) <- c(
  "TAGSGrimAge",
  "TAGSHorvath",
  "TAGSHannum",
  "TAGSPhenoAge",
  "TAGSDunedinPACE"
)

tscores <- as.data.frame(scores)
write.csv(tscores, file='TAGS_scores.csv', row.names = T)



