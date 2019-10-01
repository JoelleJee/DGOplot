if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

library(ggplot2)

BiocManager::install("DOSE")


BiocManager::install("clusterProfiler")

install.packages("ggnewscale")
library(ggnewscale)


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
  return (ggplot(mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)) +
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
    coord_flip())
  
}

DGOnetplot <- function(enrichReseult) {
  # plot gene association network
}


