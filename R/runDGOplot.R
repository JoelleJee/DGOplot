
#' Launch the shiny app for package DGOplot

#'

#' A function that launches the shiny app for this package.

#' The code has been placed in \code{./inst/shiny-scripts}.

#'

#' @return No return value but open up a shiny page.

#'

#' @examples

#' \dontrun{

#' runTestingPackage()

#' }

#'

#' @export

#' @importFrom shiny runApp



runDGOplot <- function() {
  
  appDir <- system.file("shiny-scripts",
                        
                        package = "DGOplot")
  
  shiny::runApp(appDir, display.mode = "normal")
  
  return()
  
}