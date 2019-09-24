if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DOSE")
library(DOSE)

BiocManager::install("clusterProfiler")
library(clusterProfiler)


data(DO2EG)

enrichDGO <- function(gene, Gont = "MF", Dont = "DO", 
                      pvalueCutoff = 0.05, pAdjustMethod = "BH", universe,
                      qvalueCutoff = 0.2, minGSSize = 10, maxSize = 500,
                      readable = FALSE, pool = FALSE){
  
  # perform enrichDO and enrichGO on input of genes.
  # returns a vector of both analyses.
  
  DOanalysis <- enrichDO(gene, Dont, pvalueCutoff, pAdjustMethod,
           universe, minGSSize, maxGSSize, qvalueCutoff, readable)
  
  GOanalysis <- enrichGO(gene, OrgDb, Gont, pvalueCutoff, pAdjustMethod, universe,
                         qvalueCutoff, minGSSize, maxGSSize,
                         readable, pool)
  
  DGOanalysis <- c(DOanalysis, GOanalysis)
  names(DGOanalysis) <- c("DO", "GO")
  
  return(DGOanalysis)
  
}

DGObarplot < function(enrichResult) {
  # plot GO and DO analysis 
  
}

DGOnetplot <- function(enrichReseult) {
  # plot gene association network
}


