#' Return an error if input for DGObarplot and DGOnetplot is faulty.
#'
#' A helper function that checks the input of DGObarplot and DGOnetplot.
#' 
#' @param DGOResult a list of 2 enrichResults objects, DO and GO
#'
#' @return Returns an error if input is corrupted.
#'
#'

checkInput <- function(DGOResult) {
  if (length(DGOResult) != 2) {
    stop("The input should be a list of 2 enricResult objects.")
  } 
  if (is.null(names(DGOResult)) ||
      names(DGOResult)!= c("DO", "GO")) {
    stop("The input should be a list of length 2: \n
         DO enrichment results and GO enrichment results")
  } 
  if (class(DGOResult[["DO"]]) != "enrichResult" ||
      class(DGOResult[["GO"]]) != "enrichResult") {
    stop("The input should be a list of 2 of enrichResult objects.")
  }
  return()
}
