
library(shiny)
library("tidyverse")
library("ggplot2")
library("shinythemes")

setwd("~/Desktop/qPCR-app")

# Reference scripts & functions
source("calculation.R")

# Define UI
ui <- fluidPage(

    # Application title
    titlePanel("qPCR Analyser"),

    # Sidebar with inputs 
    sidebarLayout(
        sidebarPanel(
          fileInput("qPCRdata", "Upload your qPCR data",
                    multiple = TRUE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          textInput("housekeeping", "Enter your housekeeping gene", placeholder = "e.g. gadph"),
          textInput("goi", "Enter your gene of interest", placeholder = "e.g. cft1"),
          actionButton("calc_button", "Calculate")
        ),

        # Calculation and plotting tabs
        mainPanel(
           tabsetPanel(
             tabPanel("Calculations",
                      DT::dataTableOutput("calculations")),
             tabPanel("Plot",
                      plotOutput("plot"))
           )
        )
    )
)

# Define server logic
server <- function(input, output) {

  observeEvent(input$calc_button, {
    output$calculations <-  DT::renderDataTable({
    #validate(need(input$qPCRdata, 'Please upload your data.'))
    calculateDDT(input$housekeeping, input$goi)
    })
  })
  output$plot <-  renderPlot({
    validate(need(input$qPCRdata, 'Please upload your data.'))
     })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
