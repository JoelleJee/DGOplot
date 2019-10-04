#' Draw a bar graph on the DO and GO analysis.
#'
#' A function that draws a double bar graph 
#' of DO and GO analysis, ordered and colored by group and p.adjust value.
#' 
#' @param DGOResult DO and GO enrichment analysis result returned from enrichDGO().
#' @param showCategory number of ontology groups to show from DO and GO each.
#'
#' @return Returns a double bar graph of DO and GO analysis.
#' \itemize{
#'   \item different color for DO and GO results
#'   \item color gradient by p.adjust value
#'   \item ordered by p.adjust value
#' }
#'
#' @examples
#' DGOgraph <- DGObarplot(DGOResult)
#'
#' @export
#' @import ggplot2


DGObarplot <- function(DGOResult, showCategory = 8) {
  
  DOanalysis <- DGOResult["DO"]
  GOanalysis <- DGOResult["GO"]
  
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
  doubleBarplot <- ggplot2::ggplot(doubleBarplotData, aes(x=reorder(term, -pRank), y=count)) +
    geom_bar(stat="identity", aes(col=ont, group=ont, fill=p.adjust), position="dodge") +
    ylim(0, max(doubleBarplotData$count) + 0.6) + xlab("") + 
    scale_fill_continuous() + 
    coord_flip()
  
  # take coordinates of layer data of doubleBarplot and match them back to the original data.
  ldDG <- ggplot2::layer_data(doubleBarplot)
  ldDG <- ldDG[, c("xmin", "xmax", "ymin", "ymax")]
  
  matches <- seq(12, 1)
  
  ldDG$p.adjust <- doubleBarplotData$p.adjust[matches]
  ldDG$term <- doubleBarplotData$term[matches]
  ldDG$ont <- doubleBarplotData$ont[matches]
  
  # make a new plot with geom_rect as layers; 1 for DO and 1 for GO.
  doublePlot <- ggplot2::ggplot(mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)) +
    ggplot2::geom_rect(data=ldDG[ldDG$ont=="GO", ], aes(fill=p.adjust)) +
    ggplot2::scale_fill_gradient(low="blue", high="lightskyblue1",
                        limits=c(min(ldDG$p.adjust), max(ldDG$p.adjust[ldDG$ont == "GO"])),
                        name="GO p.adjust") +
    ggplot2::new_scale_fill() +
    ggplot2::geom_rect(data=ldDG[ldDG$ont=="DO", ], aes(fill=p.adjust)) +
    ggplot2::scale_fill_gradient(low="red", high="darksalmon",
                        limits=c(min(ldDG$p.adjust), max(ldDG$p.adjust[ldDG$ont == "DO"])),
                        name="DO p.adjust") +
    ggplot2::scale_x_continuous(breaks=seq_along(unique(doubleBarplotData$term)),
                       labels=rev(unique(doubleBarplotData$term))) +
    ggplot2::coord_flip()
  
  return (doublePlot)
  
}