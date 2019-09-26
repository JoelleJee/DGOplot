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

DGObarplot < function(DGOResult) {
  
  DOanalysis <- DGOResult[1]
  GOanalysis <- DGOResult[2]
  
  # plot GO and DO analysis
  DOplot <- barplot(DOanalysis)
  GOplot <- barplot(GOanalysis)
  
  # combine the data sets in the above two plots
  
  # first order the plot data by its p-value
  DOplotData <- DOplot$data[order(DOplot$data$pvalue), ]
  GOplotData <- GOplot$data[order(GOplot$data$pvalue), ]
  
  # second give each of them a column rank that numbers the data by its p-value
  DOplotData$pRank <- seq(1, length.out = nrow(DOplotData))
  GOplotData$pRank <- seq(1, length.out = nrow(GOplotData))
  
  # combind the two dataframes in alternating order
  DGOplotData <- rbind(DOplotData, GOplotData)
  DGOplotData <- DGOplotData[order(DGOplotData$pRank, DGOplotData$ID), -10] # get rid of the pRank column
  
  # 
  
}

DGOnetplot <- function(enrichReseult) {
  # plot gene association network
}


