#' MAkes a legend for DGOnetplot
#'
#' A helper function that parses the ontology names to generate a legend.
#' 
#' @param descr a vector of gene ontology names
#'
#' @return Returns a legend to use for DGPnetplot
#'
makeLegend <- function(descr) {
  # construct legend (separate long names with newline)
  numLgnd <- length(descr)
  lgnd <- c(numLgnd)
  for (i in seq(descr)){
    title <- descr[i]
    # if title is too long
    if (nchar(title) > 20) {
      # get the location of white space
      space <- unlist(gregexpr(" ", title))
      # if there is just one space in title,
      # split location is the first space
      if (length(space) == 1) {
        split <- space[1]
      } else {
        # otherwise split location is the second space
        split <- space[2]}
      f <- substring(descr[i], 
                     1, split-1)
      s <- substring(descr[i],
                     split+1)
      title <- paste(f, s, sep = "\n")
    }
    # store the title into lgnd
    lgnd[i] <- toupper(title)
  }
  return(lgnd)
}
