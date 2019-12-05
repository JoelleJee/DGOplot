#' Draw a gene association network of DO and GO enrichment analyses results.
#'
#' A function that draws a gene association network with DO and GO terms from enrichment analyses. 
#' 
#' @param DGOResult DO and GO enrichment analysis result returned from enrichDGO().
#' @param showCategory number of ontology groups to show from DO and GO each.
#' @param pAdjustCutoff p-value cutoff.
#' @param cluster.strength clustering strength of each ontology group on the network graph
#' @param GOcol vector of colors for gradient GO edges
#' @param DOcol vector of colors for gradient DO edges
#' @param termCol name of ColorBrewer color palettes for ontology term nodes
#' @param geneCol color of gene nodes
#'
#' @return Returns a gene association network of DO and GO analyses.
#' \itemize{
#'   \item edges colored in gradient according to p.adjust value
#'   \item edges to DO and GO terms colored differently
#'   \item ontology term nodes size proportional to number of genes associated with it
#' }
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' # load data from DOSE
#' library(DOSE)
#' data(geneList)
#' gene <- names(geneList)[abs(geneList) > 2]
#' result <- enrichDGO(gene, universe=names(geneList))
#' 
#' # tamper the data little bit to see 
#' # the full functionalities of DGOplot
#' result$DO@result$p.adjust <-  result$DO@result$p.adjust / 10
#' 
#' DGOnetplot(result)
#' DGOnetplot(result,
#'            showCategory = 4,
#'            cluster.strength = 5)
#' }
#' 
#' @import RColorBrewer
#' @import igraph
#' @import annotate
#' @import graphics

DGOnetplot <- function(DGOResult, 
                       showCategory = 6, 
                       pAdjustCutoff = 0.05,
                       cluster.strength = 10,
                       GOcol = c("blue", "lightblue"),
                       DOcol = c("red", "mistyrose"),
                       termCol = "PuOr",
                       geneCol = "honeydew3") {
  
  # show Category must be > 2 for brewr.pal to work
  if (showCategory < 2) {
    stop("showCategory must be at least 2")
  }
  
  # check input for error
  checkInput(DGOResult)
  
  # check if there are at least showCategory number of valid 
  # (p.adjust value less than pAdjustCutoff) ontology groups
  # if not, produce warning message, or throw an error and stop if 0 groups.
  DOResult <- DGOResult[["DO"]]@result
  GOResult <- DGOResult[["GO"]]@result
  DOcatN <- min(showCategory, 
                sum(DOResult$p.adjust < pAdjustCutoff))
  GOcatN <- min(showCategory, 
                sum(GOResult$p.adjust < pAdjustCutoff))
  # produce warning message is insufficient number 
  # of goups found.
  warnOntN(DOcatN, GOcatN, showCategory)
  
  # combine the first showCategory number of groups from DO and GO analysis 
  plotDat <- rbind(DOResult[0:DOcatN, ], 
                   GOResult[0:GOcatN, ])
  
  # get all the geneIDs for each ontology group
  numGroups <- nrow(plotDat)
  
  # set term node color
  termCol <- try(RColorBrewer::brewer.pal(numGroups, termCol), 
                 silent = TRUE)
  
  # number of Groups cannot exceed the the max value of 
  # parameter n in brewer.pal(n, color)
  if (length(termCol) < numGroups) {
    numGroups <- length(termCol)
    plotDat <- plotDat[1:numGroups, ]
  }
  
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
  
  # set color for gene nodes
  geneCol <- rep(geneCol, numGenes)
  
  igraph::V(graphNet)$color <- c(termCol, geneCol)
  
  # set color for edges: give the edges a gradient
  # according to their p.adjust value
  
  # find unique p.adjust values and give them a rank
  padjVals <- unique(plotDat$p.adjust)
  pRanks <- order(padjVals)
  
  # Now make a gradient of blue and red, each with the 
  # number of unique p.adjust values of colors blue for GO
  GOcolfunc <- grDevices::colorRampPalette(GOcol)
  GOCols <- GOcolfunc(length(padjVals))
  
  # red for DO
  DOcolfunc <- grDevices::colorRampPalette(DOcol)
  DOCols <- DOcolfunc(length(padjVals))
  
  # Now assign each gene with a color that's appropriate with its p.adjust value 
  
  DOEdgeCols <- c()
  
  j <- 1 # to iterate through padjustVals
  if (DOcatN > 0) {
    for (i in seq(DOcatN)){
      if (DOcatN == 0) {
        break()
      }
      if (plotDat$p.adjust[i] == padjVals[j]) {
        edgeCol <- DOCols[pRanks[j]]
        j <- j + 1
      } else {
        edgeCol <- DOCols[pRanks[j - 1]]
      }
      DOEdgeCols <- c(DOEdgeCols, rep(edgeCol, plotDat$Count[i]))
    }
  }  
  GOEdgeCols <- c()
  
  if (GOcatN > 0) {
    for (i in seq(DOcatN + 1, numGroups)){
      if (plotDat$p.adjust[i] == padjVals[j]) {
        edgeCol <- GOCols[pRanks[j]]
        j <- min(length(padjVals), j + 1)
      } else {
        edgeCol <- GOCols[pRanks[j - 1]]
      }
      GOEdgeCols <- c(GOEdgeCols, rep(edgeCol, plotDat$Count[i]))
    }
  }  
  # set the edge colors
  igraph::E(graphNet)$color <- c(DOEdgeCols, GOEdgeCols)
  # make legend
  lgnd <- makeLegend(plotDat$Description)
  
  # Finally plot the result
  opar <- par(no.readonly = TRUE)
  par(bg = "gray60",
      oma = c(0,0,0,9),
      mar = c(0,0,0,0))
  coords <- netLayout(graphNet, 
                      wc=1, 
                      cluster.strength = cluster.strength)
  graphics::plot(graphNet,
                 vertex.label.font.cex = 1,
                 vertex.label.degree = pi/2,
                 vertex.label.dist = 0.5,
                 vertex.label.color = "black",
                 layout=coords)
  
  # reset margins to add legend
  par(opar)
  par(mar = c(0,0,0,0),
      oma = c(0,0,0,0.5))
  graphics::legend(x=0.75, y=0.85,lgnd, pch=21, col = termCol,
                   pt.bg=termCol, cex = 0.7, pt.cex=1.5, bty="n", ncol=1)
  par(mar = c(0,0,3,0))
  graphics::title(main = "Gene Association Network",
                  cex.main = 1.5)
  
  return(recordPlot())
  
}
  