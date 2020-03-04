library(shiny)
library(shinydashboard)
library(shinythemes)
source("waterfallPlot_module.R")
source("kaplanMeierPlot_module.R")

ui <- dashboardPage(  
  
  dashboardHeader(title = "TARGET AAML1031"),
  
  ###################### DASHBOARD SIDEBAR ######################
  dashboardSidebar(
    sidebarMenu(
      
      #---------- Gene of interest input text box -------------#
      textInput("geneInput",                                   
                label = "Enter a gene symbol", 
                placeholder = "Example: MSLN"),
      
      # --------- Plot generation tabs ------------------------#
      menuItem("Waterfall plot", tabName = "wfPlot", icon = icon("chart-bar")),
      menuItem("Kaplan-Meier curves", tabName = "kmPlot", icon = icon("notes-medical"))
      # menuItem("SNV oncoprint", tabName = "oncPrint", icon = icon("dna")),
      # menuItem("External databases", tabName = "extData", icon = icon("atlas"))
      # prescription-bottle for ADC/CAR T = prescription-bottle
    )
  ),  
  
  ###################### DASHBOARD PAGES ######################
  dashboardBody(
    
    # Using some custom CSS to change the background color of the module to white
    tags$head(tags$style(HTML(
      '/* body */
      .content-wrapper, .right-side {
      background-color: #ffffff;
      },'
    ))),
    
    tabItems(
      
      # Sourcing the waterfall plot module UI component
      tabItem(tabName = "wfPlot",

              # Calling the user interface module of the Waterfall Plot app
              wfPlotUI(id = "waterfall", label = "Waterfall plot generation")
      ),
      
      # Sourcing the waterfall plot module UI component
      tabItem(tabName = "kmPlot",
              
              # Calling the user interface module of the Waterfall Plot app
              kmPlotUI(id = "kaplanmeier", label = "Kaplan-Meier plot generation")
      ) 
      
      # tabItem(tabName = "extData",
      #         
      #         mainPanel(
      #           position = "center", 
      #           fluidRow(
      #             column(4, 
      #                    align = "left", 
      #                    uiOutput("ui_open_tab"),
      #                    tags$a(href = paste0("https://www.proteinatlas.org/search/", toupper(input$geneInput)),
      #                           "The Human Protein Atlas - Blood Atlas")
      #                    )
      #           )
      #         ))
    )
  )
)