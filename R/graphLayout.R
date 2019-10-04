#' A simplified function from NetPathMiner by ahmohamed 
#' https://stackoverflow.com/questions/28715736/how-to-spread-out-community-graph-made-by-using-igraph-package-in-r/28722680#28722680
#'
#' A function with an algorithm to spread out the nodes in an association network. 
#' 
#' @param graph garph
#' @param wc number of vertices added for algorithm
#' @param cluster.strength cluster strength
#' @param layout
#'
#' @return Returns a layout for the gene association network
#'
#' @examples
#' DGOgraph <- DGOnetplot(DGOResult)
#'
#' @import qgraph

layout.by.attr <- function(graph, wc, cluster.strength=1,layout=layout.auto) {  
  g <- graph.edgelist(get.edgelist(graph)) # create a lightweight copy of graph w/o the attributes.
  E(g)$weight <- 1
  
  attr <- cbind(id=1:vcount(g), val=wc)
  g <- g + vertices(unique(attr[,2])) + qgraph::edges(unlist(t(attr)), weight=cluster.strength)
  
  l <- layout(g, weights=E(g)$weight)[1:vcount(graph),]
  return(l)
}
