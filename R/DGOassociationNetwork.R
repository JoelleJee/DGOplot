#' Draw a gene association network of DO and GO enrichment analyses results.
#'
#' A function that draws a gene association network with DO and GO terms from enrichment analyses. 
#' 
#' @param DGOResult DO and GO enrichment analysis result returned from enrichDGO().
#' @param showCategory number of ontology groups to show from DO and GO each.
#'
#' @return Returns a gene association network of DO and GO analyses.
#' \itemize{
#'   \item edges colored in gradient according to p.adjust value
#'   \item edges to DO and GO terms colored differently
#'   \item ontology term nodes size proportional to number of genes associated with it
#' }
#'
#' @examples
#' DGOnetplot(DGOResult, showCategory = 4)
#'
#' @export
#' @import RColorBrewer
#' @import ggplot2
#' @import annotate

DGOnetplot <- function(DGOResult, showCategory = 5, 
                       pvalueCutoff = 0.05,
                       foldChange = NULL, layout = "kk", 
                       colorEdge = FALSE, circular = FALSE, node_label = TRUE) {
  
  # combine the first showCategory number of groups from DO and GO analysis 
  DONumCategory <- min(showCategory, sum(DGOResult[[1]]@result$p.adjust < pvalueCutoff))
  GONumCategory <- min(showCategory, sum(DGOResult[[2]]@result$p.adjust < pvalueCutoff))
  dataToPlot <- rbind(DGOResult[[1]]@result[1:DONumCategory, ], DGOResult[[2]]@result[1:GONumCategory, ])
  
  # get all the geneIDs for each ontology group
  numGroups <- nrow(dataToPlot)
  geneID <- c(numGroups)
  for (i in seq(numGroups)) {
    geneIDs <- strsplit(dataToPlot$geneID[i], "/")
    geneID[i] <- geneIDs
  }
  
  # repeat the term name and p.adjust value by the number of counts to combine with the geneIDs later 
  term <- list()
  pValue <- list()
  for (i in seq(numGroups)) {
    term[[i]] <- rep(dataToPlot$Description[i], dataToPlot$Count[i])
    pValue[[i]] <- rep(dataToPlot$p.adjust[i], dataToPlot$Count[i])
  }
  
  # combine the term and geneIDs into one data frame
  termGeneCombined <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(termGeneCombined) <- c("Description", "Gene", "p.adjust")
  for (i in seq(numGroups)){
    geneNames <- annotate::getSYMBOL(geneID[[i]], data = "org.Hs.eg")
    termGeneCombined <- rbind(termGeneCombined, data.frame("Description" = term[[i]],
                                                           "Gene" = geneNames,
                                                           "p.adjust" = pValue[[i]]))
  }
  
  graphNet <- graph.data.frame(termGeneCombined, directed=FALSE)
  
  
  # set the size of gene nodes to be uniform
  numGenes <- length(V(graphNet)) - numGroups           # number of nodes - number of ontology groups
  geneNodeSize <- seq(3, 3, length.out = numGenes)
  
  # set the size of nodes ontology groups proportion to the count
  termNodeSize <- dataToPlot$Count*0.9
  
  # set the node size
  V(graphNet)$size <- c(termNodeSize, geneNodeSize)
  
  # remove the labels on the terms nodes
  termLabel <- rep("", numGroups)
  geneNames <- V(graphNet)$name[(numGroups + 1):length(V(graphNet))]
  
  V(graphNet)$label <- c(termLabel, geneNames)
  
  # set color for nodes
  termColor <- RColorBrewer::brewer.pal(numGroups, "PuOr")
  geneColor <- rep("honeydew3", numGenes)
  
  V(graphNet)$color <- c(termColor, geneColor)
  
  # set color for edges: give the edges a gradient according to their p.adjust value
  
  # find unique p.adjust values and give them a rank
  padjustVals <- unique(dataToPlot$p.adjust)
  pRanks <- order(padjustVals)
  
  # Now make a gradient of blue and red, each with the number of unique p.adjust values of colors
  # blue for GO
  GOcolfunc <- RColorBrewer::colorRampPalette(c("blue", "lightblue"))
  GOColors <- GOcolfunc(length(padjustVals))
  
  # red for DO
  DOcolfunc <- RColorBrewer::colorRampPalette(c("red", "mistyrose"))
  DOColors <- DOcolfunc(length(padjustVals))
  
  # Now assign each gene with a color that's appropriate with its p.adjust value 
  
  DOEdgeColors <- c()
  
  j <- 1 # to iterate through padjustVals
  
  for (i in seq(DONumCategory)){
    if (dataToPlot$p.adjust[i] == padjustVals[j]) {
      edgeColor <- DOColors[pRanks[j]]
      j <- j + 1
    } else {
      edgeColor <- DOColors[pRanks[j - 1]]
    }
    DOEdgeColors <- c(DOEdgeColors, rep(edgeColor, dataToPlot$Count[i]))
  }
  
  GOEdgeColors <- c()
  
  for (i in seq(DONumCategory + 1, numGroups)){
    if (dataToPlot$p.adjust[i] == padjustVals[j]) {
      edgeColor <- GOColors[pRanks[j]]
      j <- min(length(padjustVals), j + 1)
    } else {
      edgeColor <- GOColors[pRanks[j - 1]]
    }
    GOEdgeColors <- c(GOEdgeColors, rep(edgeColor, dataToPlot$Count[i]))
  }
  
  # set the edge colors
  E(graphNet)$color <- c(DOEdgeColors, GOEdgeColors)
  
  # png("GeneAssociationNetwork.png", 1200, 1200)
  
  # Finally plot the result and return it
  graphNetPlot <- igraph::plot(graphNet, 
                       vertex.label.font.cex = 1,
                       vertex.label.degree = pi/2,
                       vertex.label.dist = 0.5,
                       vertex.label.color = "grey0",
                       layout=layout.by.attr(graphNet, wc=1))
  
  # Now add a legend for the ontology groups
  igraph::legend(x=1.5, y=0.5, dataToPlot$Description, pch=21, col = termColor,
         pt.bg=termColor, pt.cex=2, cex=.8, bty="n", ncol=1)
  
}