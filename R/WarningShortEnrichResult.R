warnOntN <- function(nDO, nGO, showCategory, type) {
  if (nDO == 0 & nGO == 0) {
    if (type == "bar"){
      stop("0 ontology groups in DGOResult.")
    } else if (type == "net"){
      stop("0 ontology groups in DGOResult, or below pvalueCutoff")
    }
  }
  else if (nDO == 0) {
    warn <- sprintf("0 disease ontology groups selected.\n
                    Check pvalueCutoff if using DGOnetplot().")
  } else if (nGO == 0) {
    warn <- sprintf("0 gene ontology groups selected.
                    Check pvalueCutoff if using DGOnetplot().")
  } else if (nDO < showCategory | nGO < showCategory){
    warn <- sprintf("Input showCategory = %d \n", showCategory)
    if (nDO < showCategory){
      warn <- sprintf("%sOnly %d Disease Ontology groups found.\n
                      Check pvalueCutoff if using DGOnetplot().", 
                      warn, nDO)
    } 
    if (nGO < showCategory) {
      warn <- sprintf("%sOnly %d Gene Ontology groups found.\n
                      Check pvalueCutoff if using DGOnetplot().\n",
                      warn, nGO)
    }
    warning(warn)
  } else {
    return()
  }
}
