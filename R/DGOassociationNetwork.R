#' Draw a gene association network of DO and GO enrichment analyses results.
#'
#' A function that draws a gene association network with DO and GO terms from enrichment analyses. 
#' 
#' @param DGOResult DO and GO enrichment analysis result returned from enrichDGO().
#' @param showCategory number of ontology groups to show from DO and GO each.
#' @param pvalueCutoff p-value cutoff.
#'
#' @return Returns a gene association network of DO and GO analyses.
#' \itemize{
#'   \item edges colored in gradient according to p.adjust value
#'   \item edges to DO and GO terms colored differently
#'   \item ontology term nodes size proportional to number of genes associated with it
#' }
#'
#'
#' @export
#' @import RColorBrewer
#' @import igraph
#' @import annotate
#' @import graphics

DGOnetplot <- function(DGOResult, showCategory = 6, pvalueCutoff = 0.05) {
  
  # check input for error
  checkInput(DGOResult)
  
  if (showCategory > 6) {
    stop("showCategory must be 6 or less.")
  }
  
  # combine the first showCategory number of groups from DO and GO analysis 
  DOcatN <- min(showCategory, 
                       sum(DGOResult[["DO"]]@result$p.adjust < pvalueCutoff))
  GOcatN <- min(showCategory, 
                       sum(DGOResult[["GO"]]@result$p.adjust < pvalueCutoff))
  plotDat <- rbind(DGOResult[["DO"]]@result[1:DOcatN, ], 
                      DGOResult[["GO"]]@result[1:GOcatN, ])
  
  # check if there are at least showCategory number of valid 
  # (p.adjust value less than pvalueCutoff) ontology groups
  # if not, produce warning message, or throw an error and stop if 0 groups.
  warnOntN(DOcatN, GOcatN, showCategory, type = "net")

  
  # get all the geneIDs for each ontology group
  numGroups <- nrow(plotDat)
  geneID <- c(numGroups)
  for (i in seq(numGroups)) {
    geneIDs <- strsplit(plotDat$geneID[i], "/")
    geneID[i] <- geneIDs
  }
  
  # repeat the term name and p.adjust value by the number 
  # of counts to combine with the geneIDs later 
  term <- list()
  pValue <- list()
  for (i in seq(numGroups)) {
    term[[i]] <- rep(plotDat$Description[i], plotDat$Count[i])
    pValue[[i]] <- rep(plotDat$p.adjust[i], plotDat$Count[i])
  }
  
  # combine the term and geneIDs into one data frame
  termGene <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(termGene) <- c("Description", "Gene", "p.adjust")
  for (i in seq(numGroups)){
    geneNames <- annotate::getSYMBOL(geneID[[i]], data = "org.Hs.eg")
    termGene <- rbind(termGene, 
                      data.frame("Description" = term[[i]],
                                 "Gene" = geneNames,
                                 "p.adjust" = pValue[[i]]))
  }
  
  graphNet <- igraph::graph.data.frame(termGene, directed=FALSE)
  
  
  # set the size of gene nodes to be uniform
  # number of nodes = number of ontology groups
  numGenes <- length(igraph::V(graphNet)) - numGroups         
  geneSize <- seq(3, 3, length.out = numGenes)
  
  # set the size of nodes ontology groups proportion to the count
  termSize <- plotDat$Count*0.9
  
  # set the node size
  igraph::V(graphNet)$size <- c(termSize, geneSize)
  
  # remove the labels on the terms nodes
  termLabel <- rep("", numGroups)
  geneNames <- igraph::V(graphNet)$name[(numGroups + 1):
                                          length(igraph::V(graphNet))]
  
  igraph::V(graphNet)$label <- c(termLabel, geneNames)
  
  # set color for nodes
  termCol <- RColorBrewer::brewer.pal(numGroups, "PuOr")
  geneCol <- rep("honeydew3", numGenes)
  
  igraph::V(graphNet)$color <- c(termCol, geneCol)
  
  # set color for edges: give the edges a gradient
  # according to their p.adjust value
  
  # find unique p.adjust values and give them a rank
  padjVals <- unique(plotDat$p.adjust)
  pRanks <- order(padjVals)
  
  # Now make a gradient of blue and red, each with the 
  # number of unique p.adjust values of colors blue for GO
  GOcolfunc <- grDevices::colorRampPalette(c("blue", "lightblue"))
  GOCols <- GOcolfunc(length(padjVals))
  
  # red for DO
  DOcolfunc <- grDevices::colorRampPalette(c("red", "mistyrose"))
  DOCols <- DOcolfunc(length(padjVals))
  
  # Now assign each gene with a color that's appropriate with its p.adjust value 
  
  DOEdgeCols <- c()
  
  j <- 1 # to iterate through padjustVals
  
  for (i in seq(DOcatN)){
    if (plotDat$p.adjust[i] == padjVals[j]) {
      edgeCol <- DOCols[pRanks[j]]
      j <- j + 1
    } else {
      edgeCol <- DOCols[pRanks[j - 1]]
    }
    DOEdgeCols <- c(DOEdgeCols, rep(edgeCol, plotDat$Count[i]))
  }
  
  GOEdgeCols <- c()
  
  for (i in seq(DOcatN + 1, numGroups)){
    if (plotDat$p.adjust[i] == padjVals[j]) {
      edgeCol <- GOCols[pRanks[j]]
      j <- min(length(padjVals), j + 1)
    } else {
      edgeCol <- GOCols[pRanks[j - 1]]
    }
    GOEdgeCols <- c(GOEdgeCols, rep(edgeCol, plotDat$Count[i]))
  }
  
  # set the edge colors
  igraph::E(graphNet)$color <- c(DOEdgeCols, GOEdgeCols)
  
  # Finally plot the result
  par(bg = "gray60")
  graphics::plot(graphNet,
                 vertex.label.font.cex = 1,
                 vertex.label.degree = pi/2,
                 vertex.label.dist = 0.5,
                 vertex.label.color = "black",
                 layout=layout.by.attr(graphNet, wc=1, cluster.strength = 50),
                 bg = 244)
  
  # Now add a legend for the ontology groups
  graphics::legend(x=1.2, y=0.5, plotDat$Description, pch=21, col = termCol,
                   pt.bg=termCol, cex = 0.9, pt.cex=2, bty="n", ncol=1)
  graphics::title(main = "Gene Association Network")
  
}
