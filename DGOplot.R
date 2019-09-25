if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DOSE")


BiocManager::install("clusterProfiler")



data(DO2EG)

enrichDGO <- function(gene, Gont = "MF", Dont = "DO", OrgDb =  org.Hs.eg.db,
                      pvalueCutoff = 0.05, pAdjustMethod = "BH", universe,
                      qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500,
                      readable = FALSE, pool = FALSE){
  
  # perform enrichDO and enrichGO on input of genes.
  # returns a vector of both analyses.
  
  DOanalysis <- DOSE::enrichDO(gene, Dont, pvalueCutoff, pAdjustMethod, 
                               universe, minGSSize, maxGSSize, qvalueCutoff, readable)
  
  GOanalysis <- clusterProfiler::enrichGO(gene, OrgDb, keyType = "ENTREZID", Gont, pvalueCutoff, 
                                          pAdjustMethod, universe, qvalueCutoff, minGSSize, 
                                          maxGSSize, readable, pool)
  
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


