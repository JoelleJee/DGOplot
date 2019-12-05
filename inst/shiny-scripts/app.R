#' Launch the shiny app for package DGOplot
#'
#' A function that launches the shiny app for this package.
#' The code has been placed in /code{.inst/shiny-scripts}
#' 
#' @return No return value but opens up a shiny page.
#'
#' @examples 
#' \dontrun{
#' runDGOplot()
#' }
#' 
#' @import DOSE
#' @import shiny
#' @importFrom DT DTOutput renderDT

# load geneList data
library(DOSE)
data("geneList")
data("DO2EG")
data("EG2DOLite")
data("DOLiteTerm")

# perform enrichDGO()
gene <- names(geneList)[abs(geneList) > 2]
result <- enrichDGO(gene, universe=names(geneList))
# we will change the results to show the full 
# functionalities of DGOplot.
result$DO@result$p.adjust <-  result$DO@result$p.adjust / 10

ui <- shiny::fluidPage(
  # title
  shiny::titlePanel("DGOplot Overview"),
  
  # options to see DGObarplot, DGOnetplot, and enrichDGO
  shiny::navbarPage(
    title = "Test",
    id = "navbar",
    shiny::tabPanel(title = "Barplot",
             value = "bar"),
    
    shiny::tabPanel(title = "Gene Association Network",
             value = "net"),
    
    shiny::tabPanel(title = "Compute DGO Enrichment Analysis",
             value = "enrich"),
    
    shiny::tabPanel(title = "Quit",
             value = "stop",
             icon = icon("circle-o-notch"))
  ),
  
  # buttons to select values for parameters
  shiny::sidebarLayout(
    position = "left",
    shiny::sidebarPanel(
      "change argument values:",
      
      width = 3,
      
      shiny::sliderInput(
        inputId = "pvalue",
        label = "Enrichment Analysis p-value Cutoff",
        value = 0.05, min = 0.01, max = 0.5
      ),
      
      shiny::sliderInput(
        inputId = "p.adjust",
        label = "p.adjust cutoff for plotting",
        value = 0.05, min = 0.01, max = 0.5
      ),
      
      shiny::sliderInput(
        inputId = "barShow",
        label = "Number of groups to show on bar plot",
        value = 8, min = 1, max = 15
      ),
      
      shiny::sliderInput(
        inputId = "netShow",
        label = "Number of groups to show on network plot",
        value = 4, min = 2, max = 5
      ),
      
      shiny::radioButtons(
        inputId = "DOlow",
        label = "DO: low p.adjust value",
        choices = list("purple", "blue", "red"),
        selected = "red",
        inline = TRUE
      ),
      
      shiny::radioButtons(
        inputId = "DOhigh",
        label = "DO: high p.adjust value",
        choiceValues = list("#dec7de", "#e0f0fc", "#fce0f0"),
        choiceNames = list("lilac", "lightblue", "pink"),
        selected = "#fce0f0",
        inline = TRUE
      ),
      
      shiny::radioButtons(
        inputId = "GOlow",
        label = "GO: low p.adjust value",
        choices = list("purple", "blue", "red"),
        selected = "blue",
        inline = TRUE
      ),
      
      shiny::radioButtons(
        inputId = "GOhigh",
        label = "GO: high p.adjust value",
        choiceValues = list("#dec7de", "#e0f0fc", "#fce0f0"),
        choiceNames = list("lilac", "lightblue", "pink"),
        selected = "#e0f0fc",
        inline = TRUE
      ),
      
      shiny::radioButtons(
        inputId = "GOterm",
        label = "GO term palette",
        choices = list("Set3", "PuOr", "Paired"),
        selected = "PuOr"
      ),
      
      shiny::radioButtons(
        inputId = "GOoption",
        label = "GO option",
        choices = list("BP", "MF", "CC", "ALL"),
        selected = "BP"
      ),
      
      shiny::radioButtons(
        inputId = "DOoption",
        label = "DO option",
        choices = list("DO", "DOLite"),
        selected = "DO"
      )
    ),
    
    shiny::mainPanel(
      shiny::conditionalPanel(
        condition = "input.navbar == 'enrich'",
        shiny::verticalLayout(
          DT::DTOutput("GOenrich"),
          DT::DTOutput("DOenrich")
        )
      ),
      shiny::conditionalPanel(
        condition = "input.navbar != 'enrich'",
        shiny::plotOutput("plot")
      )
    )
  )
)

server <- function(input, output) {
  
  shiny::observe({
    
    # render DGObarplot
    if (input$navbar == "bar") {
      output$plot <- shiny::renderPlot(
        
        expr = DGObarplot(result, 
                          showCategory=input$barShow, 
                          pAdjustCutoff=input$p.adjust,
                          DOcol=c(input$DOlow, input$DOhigh),
                          GOcol=c(input$GOlow, input$GOhigh)),
        
        width = 600,
        height = 600,
        res = 80
        
      )
      # render DGOnetplot
    } else if (input$navbar == "net") {
      output$plot <- shiny::renderPlot(
        
        expr = DGOnetplot(result, 
                          showCategory=input$netShow, 
                          pAdjustCutoff=input$p.adjust,
                          DOcol=c(input$DOlow, input$DOhigh),
                          GOcol=c(input$GOlow, input$GOhigh),
                          termCol=input$GOterm),
        res = 100,
        width = 1300,
        height = 700
        
      )
      # render the enrichment analyses data
    } else if (input$navbar == "enrich") {
      result <<- enrichDGO(gene = gene,
                           Gont = input$GOoption,
                           Dont = input$DOoption,
                           universe = names(geneList),
                           pvalueCutoff = input$pvalue)
      # original result doesn't show significant p.adjust value
      # for the purpose of visualizing how DGObarplot and DGOnetplot work,
      # we will change the values manually
      result$DO@result[p.adjust] <<- result$DO@result[p.adjust] / 10
      output$DOenrich <- DT::renderDT(
        DT::datatable(
          result$DO@result[1:20, 
                           c("Description", "p.adjust", "Count")],
          escape = FALSE))
      
      output$GOenrich <- DT::renderDT(
        DT::datatable(
          result$GO@result[1:20, 
                           c("Description", "p.adjust", "Count")],
          escape = FALSE))
    } else if (input$navbar == "stop") {
      shiny::stopApp()
    }
  })
}  
shiny::shinyApp(ui = ui, server = server)
