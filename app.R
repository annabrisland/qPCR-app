
library(shiny)
library("tidyverse")
library("ggplot2")
library("shinythemes")

setwd("~/Desktop/qPCR-app")

# Reference scripts & functions
source("calculation.R")
source("plot.R")

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
          textInput("housekeeping", "Enter your housekeeping gene", placeholder = "e.g. gapdh"),
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
  
  data <- reactive({read.csv(input$qPCRdata$datapath, skip = 7)
  })

  observeEvent(input$calc_button, {
    output$calculations <-  DT::renderDataTable({
      req(input$qPCRdata)
    #validate(need(input$qPCRdata, 'Please upload your data.'))
    calculateDDT(data(), input$housekeeping, input$goi)
    })
  })
  output$plot <-  renderPlot({
    req(input$qPCRdata)
    qPCRplot(data(), input$housekeeping, input$goi)
     })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
