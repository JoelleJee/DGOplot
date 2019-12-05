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

source("../../R/globals.R")
library(DOSE)
data(geneList)
runiverse <<- names(geneList)
rgene <<- names(geneList)[abs(geneList) > 2]
rresult <<- enrichDGO(rgene, universe=runiverse)

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
    
    shiny::tabPanel(title = "See input Files",
             value = "file"),
    
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
      
      # Input: Select a file ----
      shiny::fileInput("gene", "Upload a CSV File of EntrezIDs of gene of interest",
                multiple = FALSE,
                accept = c("text/csv",
                           ".csv")
      ),
      
      shiny::fileInput("universe", "Upload the CSV file of EntrezIDs and its expression
              levels of the universe genes",
                multiple = FALSE,
                accept = c("text/csv",
                           ".csv")
      ),
      
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
        condition = "input.navbar == 'net' | input.navbar == 'bar'",
        shiny::plotOutput("plot")
      ),
      shiny::conditionalPanel(
        condition = "input.navbar == 'file'",
        shiny::splitLayout(
          DT::DTOutput("inputGene"),
          DT::DTOutput("universe")
        )
      )
    )
  )
)

server <- function(input, output) {
  
  shiny::observe({
    if (input$navbar == "file"){
      if (!(is.null(input$gene))){
        rgene <<- read.csv(input$gene$datapath, 
                         header = FALSE)
        geneC <- rgene
        output$inputGene <- DT::renderDT(geneC)
        rgene <<- as.character(unlist(rgene))
        names(rgene) <<- c()
      } 
      if (!(is.null(input$universe))) {
        runiverse <<- read.csv(input$universe$datapath,
                             header = FALSE)
        universeC <- runiverse
        output$universe <- DT::renderDT(universeC)
        exprVal <- unlist(runiverse[2])
        geneNames <- unlist(runiverse[1])
        runiverse <<- exprVal
        names(runiverse) <<- geneNames
      } else {
        return(NULL)
      }
        
      # rendering the double bar plot  
    } else if (input$navbar == "bar") {
      if (is.null(rresult)) return(NULL)
      output$plot <- shiny::renderPlot(
        
        expr = DGObarplot(rresult, 
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
      if (is.null(rresult)) return(NULL)
      output$plot <- shiny::renderPlot(
        
        expr = DGOnetplot(rresult, 
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
      if (is.null(rgene) || is.null(runiverse)) {
        return(NULL)
      } else{ 
        rresult <<- enrichDGO(gene = rgene,
                              Gont = input$GOoption,
                              Dont = input$DOoption,
                              universe = runiverse,
                              pvalueCutoff = input$pvalue)
        
        output$DOenrich <- DT::renderDT(
          DT::datatable(
            rresult$DO@result,
            escape = FALSE))
        
        output$GOenrich <- DT::renderDT(
          DT::datatable(
            rresult$GO@result,
            escape = FALSE))
      }
    } else if (input$navbar == "stop") {
      shiny::stopApp()
    } 
  })
}  
shiny::shinyApp(ui = ui, server = server)
