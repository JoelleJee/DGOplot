#' Perform DO and GO enrichment analysis
#'
#' A function that performs DO and GO enrichment analysis given a list of geneIDs.
#' 
#' @param gene a vector of entrez gene id.
#' @param Gont GO categories; one of "BP", "MF", "CC", or "ALL" for all three.
#' @param Dont DO categories; one of DO or DOLite
#' @param pvalueCutoff p-valuve cutoff
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param universe background genes
#' @param qvalueCutoff q-value cutoff
#' @param minGSSize minimal size of genes annotated by NCG category for testing
#' @param maxGSSize maximal size of each geneSet for analyzing
#' @param readable whether mapping gene ID to gene Name
#' @param pool If ont='ALL', whether pool 3 GO sub-ontologies
#'
#' @return Returns a list of 2 enrichResults
#' \itemize{
#'   \item DOanalysis - an enrichResult object with DO enrichment analysis result
#'   \item GOanalysis - an enrichResult object with DO enrichment analysis result
#' }
#'
#' @examples 
#' \dontrun{
#' # load data from DOSE
#' library(DOSE)
#' data(geneList)
#' gene <- names(geneList)[abs(geneList) > 2]
#' result <- enrichDGO(gene, universe=names(geneList))
#' }
#' @export
#' @import org.Hs.eg.db
#' @import BiocManager
#' @importFrom DOSE enrichDO
#' @importFrom clusterProfiler enrichGO
#' 
enrichDGO <- function(gene, Gont = "MF", Dont = "DO",
                      pvalueCutoff = 0.05, pAdjustMethod = "BH", universe,
                      qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500,
                      readable = FALSE, pool = FALSE){
  data("DOLite2EG")
  data("EG2DOLite")
  data("DOLiteTerm")
  # perform enrichDO and enrichGO on input of genes.
  DOanalysis <- DOSE::enrichDO(gene, Dont, pvalueCutoff, pAdjustMethod, 
                               universe, minGSSize, maxGSSize, qvalueCutoff, readable)
  
  GOanalysis <- clusterProfiler::enrichGO(gene, OrgDb = org.Hs.eg.db, 
                                          keyType = "ENTREZID", Gont, 
                                          pvalueCutoff, pAdjustMethod, 
                                          universe, qvalueCutoff, 
                                          minGSSize, maxGSSize, 
                                          readable, pool)
  
  GOanalysis@result <- GOanalysis@result[ , colnames(DOanalysis@result)]
  
  DGOanalysis <- c(DOanalysis, GOanalysis)
  
  if (length(DGOanalysis) == 0) {
      stop("The input should be a vector of gene IDs with a fold change")
  }
  
  names(DGOanalysis) <- c("DO", "GO")
  
  # returns a list of both analyses.
  return(DGOanalysis)
  
}




