if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

library(ggplot2)

BiocManager::install("DOSE")
library(DOSE)

BiocManager::install("clusterProfiler")
library(clusterProfiler)

install.packages("ggnewscale")
library(ggnewscale)

install.packages("igraph")
library(igraph)

library(org.Hs.eg.db)
library(annotate)

install.packages("qgraph")
library(qgraph)
 
install.packages("RColorBrewer")
library(RColorBrewer)

enrichDGO <- function(gene, Gont = "MF", Dont = "DO", OrgDb = "org.Hs.eg.db",
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

DGObarplot <- function(DGOResult, showCategory = 8) {
  
  DOanalysis <- DGOResult[[1]]
  GOanalysis <- DGOResult[[2]]
  
  # plot GO and DO analysis
  DOplot <- barplot(DOanalysis, showCategory)
  GOplot <- barplot(GOanalysis, showCategory)
  
  # combine the data sets in the above two plots
  
  # first order the plot data by its p-value
  DOplotData <- DOplot$data[order(DOplot$data$p.adjust), ]
  GOplotData <- GOplot$data[order(GOplot$data$p.adjust), ]
  
  # second give each of them a column rank that numbers the data by its p-value
  DOplotData$pRank <- seq(1, length.out = nrow(DOplotData))
  GOplotData$pRank <- seq(1, length.out = nrow(GOplotData))
  
  # combind the two dataframes in alternating order
  DGOplotData <- rbind(DOplotData, GOplotData)
  DGOplotData <- DGOplotData[order(DGOplotData$pRank, DGOplotData$ID), ]
  
  
  # Now have to make "groups" to make a double bar plot of DO and GO analysis.
  
  # ontology type
  ont <- substr(DGOplotData$ID, 1, 2)
  
  # ontology term combine by rank or p.adjust
  numDO <- nrow(DOplotData)
  numGO <- nrow(GOplotData)
  
  numGroups <- min(numDO, numGO) + abs(numDO - numGO) # to account for the difference of analysis results
  
  # denote the groups by term
  term <- character(numGroups)
  
  # indices of DGOplotData to combine GO and DO terms for grouping
  iterGroups = c(seq(1, by = 2, length.out = min(numDO, numGO)), 
                 seq((2 * min(numDO, numGO) + 1), nrow(DGOplotData)))
  
  # iterate DGOplotData to store groups in term:
  # one term for each ontology group
  iterTerm <- 1
  for (i in iterGroups) {
    if (i <= min(numDO, numGO) * 2) {
      term[c(iterTerm, iterTerm + 1)] <- paste(DGOplotData$Description[i], "\n", DGOplotData$Description[i + 1])
      iterTerm <- iterTerm + 2
      }else {
      term[iterTerm] <- paste(DGOplotData$Description[i])
      iterTerm <- iterTerm + 1
    }
  }
  
  # combine all to one data frame for plotting a double bar graph
  doubleBarplotData <- data.frame(ont, term, DGOplotData$Count, DGOplotData$p.adjust, DGOplotData$pRank)
  names(doubleBarplotData) <- c("ont", "term", "count", "p.adjust","pRank")
  
  # plot a double bar graph; group by "ont" and fill by p.adjust value
  # below is inspired by teunbrand 
  # https://stackoverflow.com/questions/57613428/grouping-scale-fill-gradient-continuous-grouped-bar-chart
  doubleBarplot <- ggplot(doubleBarplotData, aes(x=reorder(term, -pRank), y=count)) +
    geom_bar(stat="identity", aes(col=ont, group=ont, fill=p.adjust), position="dodge") +
    ylim(0, max(doubleBarplotData$count) + 0.6) + xlab("") + 
    scale_fill_continuous() + 
    coord_flip()
    
  # take coordinates of layer data of doubleBarplot and match them back to the original data.
  ldDG <- layer_data(doubleBarplot)
  ldDG <- ldDG[, c("xmin", "xmax", "ymin", "ymax")]
  
  matches <- seq(12, 1)
  
  ldDG$p.adjust <- doubleBarplotData$p.adjust[matches]
  ldDG$term <- doubleBarplotData$term[matches]
  ldDG$ont <- doubleBarplotData$ont[matches]
  
  # make a new plot with geom_rect as layers; 1 for DO and 1 for GO.
  doublePlot <- ggplot(mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)) +
    geom_rect(data=ldDG[ldDG$ont=="GO", ], aes(fill=p.adjust)) +
    scale_fill_gradient(low="blue", high="lightskyblue1",
                        limits=c(min(ldDG$p.adjust), max(ldDG$p.adjust[ldDG$ont == "GO"])),
                        name="GO p.adjust") +
    new_scale_fill() +
    geom_rect(data=ldDG[ldDG$ont=="DO", ], aes(fill=p.adjust)) +
    scale_fill_gradient(low="red", high="darksalmon",
                        limits=c(min(ldDG$p.adjust), max(ldDG$p.adjust[ldDG$ont == "DO"])),
                        name="DO p.adjust") +
    scale_x_continuous(breaks=seq_along(unique(doubleBarplotData$term)),
                       labels=rev(unique(doubleBarplotData$term))) +
    coord_flip()
  
  return (doublePlot)
  
}

# A simplified function from NetPathMiner from ahmohamed 
# https://stackoverflow.com/questions/28715736/how-to-spread-out-community-graph-made-by-using-igraph-package-in-r/28722680#28722680
 
layout.by.attr <- function(graph, wc, cluster.strength=1,layout=layout.auto) {  
  g <- graph.edgelist(get.edgelist(graph)) # create a lightweight copy of graph w/o the attributes.
  E(g)$weight <- 1
  
  attr <- cbind(id=1:vcount(g), val=wc)
  g <- g + vertices(unique(attr[,2])) + igraph::edges(unlist(t(attr)), weight=cluster.strength)
  
  l <- layout(g, weights=E(g)$weight)[1:vcount(graph),]
  return(l)
}

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
    geneNames <- getSYMBOL(geneID[[i]], data = "org.Hs.eg")
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
  termColor <- brewer.pal(numGroups, "PuOr")
  geneColor <- rep("honeydew3", numGenes)
  
  V(graphNet)$color <- c(termColor, geneColor)
  
  # set color for edges: give the edges a gradient according to their p.adjust value
  
  # find unique p.adjust values and give them a rank
  padjustVals <- unique(dataToPlot$p.adjust)
  pRanks <- order(padjustVals)
  
  # Now make a gradient of blue and red, each with the number of unique p.adjust values of colors
  # blue for GO
  GOcolfunc <- colorRampPalette(c("blue", "lightblue"))
  GOColors <- GOcolfunc(length(padjustVals))
  
  # red for DO
  DOcolfunc <- colorRampPalette(c("red", "mistyrose"))
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
  graphNetPlot <- plot(graphNet, 
                       vertex.label.font.cex = 1,
                       vertex.label.degree = pi/2,
                       vertex.label.dist = 0.5,
                       vertex.label.color = "grey0",
                       layout=layout.by.attr(graphNet, wc=1))
  
  # Now add a legend for the ontology groups
  legend(x=1.5, y=0.5, dataToPlot$Description, pch=21, col = termColor,
         pt.bg=termColor, pt.cex=2, cex=.8, bty="n", ncol=1)
  
  # dev.off()
  
  # return(graphNetPlot)
  
}



