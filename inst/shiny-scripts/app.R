library(shiny)
library(DOSE)
data(geneList)

gene <- names(geneList)[abs(geneList) > 2]
result <- enrichDGO(gene, universe=names(geneList))

ui <- fluidPage(
  
  navbarPage(title = "Test", id="navbar", 
             
             tabPanel(title = "Barplot", 
                      value = "bar"),
             
             tabPanel(title = "Gene Association Network",
                      value = "net"),
             
             tabPanel(title = "Show Both Plots",
                      value = "both"),
             
             tabPanel(title = "Quit", 
                      value="stop", 
                      icon = icon("circle-o-notch"))
  ),
  
  sliderInput(inputId = "num",
              label = "Choose number of ontology groups to show",
              value = 4, min = 1, max = 8),
  
  sliderInput(inputId = "p-value",
              label = "Select p-value Cutoff",
              value = 0.05, min = 0.01, max = 0.5),
  
  plotOutput("bar")
)

server <- function(input, output) {
  
  result <- enrichDGO(gene, universe=names(geneList))
  
  output$bar <- renderPlot(
    
    expr = DGObarplot(result, showCategory = input$num),
    width = 400,
    height = 600
  )
  
  observe({
    if (input$navbar == "stop") {
      stopApp()
    }
  })
  
}

shinyApp(ui = ui, server = server)

