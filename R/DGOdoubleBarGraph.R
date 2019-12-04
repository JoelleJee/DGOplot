#' Draw a bar graph on the DO and GO analysis.
#'
#' A function that draws a double bar graph 
#' of DO and GO analysis, ordered and colored by group and p.adjust value.
#' 
#' @param DGOResult DO and GO enrichment analysis result returned from enrichDGO().
#' @param showCategory number of ontology groups to show from DO and GO each.
#' @param DOcol color to use for plotting DO.
#' @param GOcol color to use for plotting GO.
#' @param pAdjustCutoff p.adjust cutoff value for ontology group to plot.
#'
#' @return Returns a double bar graph of DO and GO analysis.
#' \itemize{
#'   \item different color for DO and GO results
#'   \item color gradient by p.adjust value
#'   \item ordered by p.adjust value
#' }
#'
#'
#' @export
#' @import ggplot2
#' @import ggnewscale
#' @import graphics


DGObarplot <- function(DGOResult, showCategory = 8, 
                       DOcol = "red", GOcol = "blue",
                       pAdjustCutoff = 0.05) {
  
  checkInput(DGOResult)
  
  DOanalysis <- DGOResult[["DO"]]
  GOanalysis <- DGOResult[["GO"]]
  
  # now adjust the pvalueCutoff value in 
  # both DOanalysis and GOanalysis to user input value
  GOanalysis@pvalueCutoff <- pAdjustCutoff
  DOanalysis@pvalueCutoff <- pAdjustCutoff
  
  # warning message if there are showCategory number of resulting
  # ontology groups in DO and GOplot or throws error if there are none.
  nGO = sum(GOanalysis@result$p.adjust <= pAdjustCutoff)
  nDO = sum(DOanalysis@result$p.adjust <= pAdjustCutoff)
  warnOntN(nDO, nGO, showCategory)
  
  # plot GO and DO analysis
  # if there are no ontology groups below pAdjustCutoff
  # just return one of the above graphs
  if (nGO == 0) {
    DOplot <- graphics::barplot(DOanalysis) +
              ggplot2::ylab("Number of Genes") + 
              ggplot2::ggtitle("DO Enrichment Plot")
    return(DOplot)
  } else if (nDO == 0) {
    GOplot <- graphics::barplot(GOanalysis) +
              ggplot2::ylab("Number of Genes") + 
              ggplot2::ggtitle("GO Enrichment Plot") 
    return(GOplot)
  }
  
  
  # combine the data sets in the above two plots
  # first order the plot data by its p-value
  DOpData <- DOanalysis@result[order(DOanalysis@result$p.adjust), ]
  DOpData <- DOpData[1:min(showCategory, nDO), ]
  GOpData <- GOanalysis@result[order(GOanalysis@result$p.adjust), ]
  GOpData <- GOpData[1:min(showCategory, nGO), ]
  
  # second give each of them a column rank
  # that numbers the data by its p-value
  DOpData$pRank <- seq(1, length.out = nrow(DOpData))
  GOpData$pRank <- seq(1, length.out = nrow(GOpData))
  
  # combind the two dataframes in alternating order
  DGOpData <- rbind(DOpData, GOpData)
  DGOpData <- DGOpData[order(DGOpData$pRank, DGOpData$ID), ]
  
  
  # Now have to make "groups" to make a double bar plot of DO and GO analysis.
  
  # ontology type
  ont <- substr(DGOpData$ID, 1, 2)
  
  # ontology term combine by rank or p.adjust
  numDO <- nrow(DOpData)
  numGO <- nrow(GOpData)
  
  numGroups <- 2*min(numDO, numGO) + abs(numDO - numGO) # to account for the difference of analysis results
  
  # denote the groups by term
  term <- character(numGroups)
  
  # indices of DGOplotData to combine GO and DO terms for grouping
  iterGroups = c(seq(1, by = 2, length.out = min(numDO, numGO)))
  if (min(numDO, numGO) < showCategory ){
    iterGroups <- c(iterGroups, seq((2 * min(numDO, numGO) + 1), nrow(DGOpData)))
  }
  # iterate DGOplotData to store groups in term:
  # one term for each ontology group
  iTerm <- 1
  for (i in iterGroups) {
    if (i <= 2*min(numDO, numGO)) {
      term[c(iTerm, iTerm + 1)] <- paste(DGOpData$Description[i+1], 
                                         "\n", 
                                         DGOpData$Description[i])
      iTerm <- iTerm + 2
    }else {
      term[iTerm] <- paste(DGOpData$Description[i])
      iTerm <- iTerm + 1
    }
  }
  
  # combine all to one data frame for plotting a double bar graph
  dblBarData <- data.frame(ont, term, 
                           DGOpData$Count, 
                           DGOpData$p.adjust, 
                           DGOpData$pRank)
  names(dblBarData) <- c("ont", "term", "count", "p.adjust","pRank")
  
  # plot a double bar graph; group by "ont" and fill by p.adjust value
  # below is inspired by teunbrand 
  # https://stackoverflow.com/questions/57613428/grouping-scale-fill-gradient-continuous-grouped-bar-chart
  dblBarplot <- ggplot2::ggplot(dblBarData, 
                                ggplot2::aes(x=stats::reorder(term, -pRank), y=count)) +
                ggplot2::geom_bar(stat="identity", 
                                  ggplot2::aes(col=ont, group=ont, fill=p.adjust), 
                                  position="dodge") +
                ggplot2::ylim(0, max(dblBarData$count) + 0.6) + 
                ggplot2::xlab("") + 
                ggplot2::scale_fill_continuous() + 
                ggplot2::coord_flip()
  
  # take coordinates of layer data of doubleBarplot and match them back to the original data.
  ldDG <- ggplot2::layer_data(dblBarplot)
  ldDG <- ldDG[, c("xmin", "xmax", "ymin", "ymax")]
  
  matches <- seq(nrow(dblBarData), 1)
  
  ldDG$p.adjust <- dblBarData$p.adjust[matches]
  ldDG$term <- dblBarData$term[matches]
  ldDG$ont <- dblBarData$ont[matches]
  
  # now for the color scheme of DO and GO
  if (DOcol == "red" && GOcol == "blue") {
    colHighD <- "darksalmon"
    colLowD <- "red"
    colHighG <- "lightskyblue"
    colLowG <- "blue"
  } else {
    colHighD <- DOcol[1]
    colLowD <- DOcol[2]
    colHighG <- GOcol[1]
    colLowG <- GOcol[2]
  }
  
  
  # make a new plot with geom_rect as layers; 1 for DO and 1 for GO.
  ggplot2::ggplot(mapping=ggplot2::aes(xmin=xmin, xmax=xmax, 
                                       ymin=ymin, ymax=ymax)) + 
    ggplot2::geom_rect(data=ldDG[ldDG$ont=="GO", ], ggplot2::aes(fill=p.adjust)) +
    ggplot2::scale_fill_gradient(low=colLowG, high=colHighG,
                                 limits=c(min(ldDG$p.adjust), 
                                          max(ldDG$p.adjust[ldDG$ont == "GO"])),
                                 name="GO p.adjust") +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_rect(data=ldDG[ldDG$ont=="DO", ], ggplot2::aes(fill=p.adjust)) +
    ggplot2::scale_fill_gradient(low=colLowD, high=colHighD,
                                 limits=c(min(ldDG$p.adjust), 
                                          max(ldDG$p.adjust[ldDG$ont == "DO"])),
                                 name="DO p.adjust") +
    ggplot2::scale_x_continuous(breaks=seq_along(unique(dblBarData$term)),
                                labels=rev(unique(dblBarData$term))) +
    ggplot2::coord_flip() + 
    ggplot2::ylab("Number of Genes") + 
    ggplot2::ggtitle("DGO Enrichment Plot")
  
}
