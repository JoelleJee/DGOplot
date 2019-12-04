#' Warning or Error generation for the number of ontology groups input.
#'
#' A helper function that checks the number of ontology groups 
#' and if showCategory exceeds them.
#' 
#' @param nDO number of DO enrichment result groups
#' @param nGO number of GO enrichment result groups
#' @param showCategory number of groups to show on the plot
#' @param type "bar" for DGObarplot and "net" for DGOnetplot
#'
#' @return warning or error depending on the input

warnOntN <- function(nDO, nGO, showCategory) {
  # warning message
  warn <- c()
  
  if (nDO == 0 & nGO == 0) {
    stop("0 ontology groups in DGOResult, or below pValueCutoff.\n")
  } else if (nDO < showCategory | nGO < showCategory){
    warn <- sprintf("Input showCategory = %d \n", showCategory)
    if (nDO < showCategory){
      warn <- sprintf("%s%d Disease Ontology groups found.Check pvalueCutoff.\n", 
                      warn, nDO)
      if (nDO == 0) {
        warn <- sprintf("%sOnly displaying GO enrichment results.\n", warn)
      } 
    } 
    if (nGO < showCategory) {
      warn <- sprintf("%s%d Gene Ontology groups found.\nCheck pvalueCutoff.\n",
                      warn, nGO)
      if (nGO == 0) {
        warn <- sprintf("%sOnly displaying GO enrichment results.\n", warn)
      }
    } 
    warning(warn)
  } else {
    return()
  }
}
