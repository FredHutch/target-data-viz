library(shiny)
library(shinydashboard)
library(shinythemes)
library(shinyWidgets)
source("waterfallPlot_module.R")
source("kaplanMeierPlot_module.R")

ui <- dashboardPage(  
  
  dashboardHeader(title = "Meshinchi Lab Data Viz Tools"),
  
  ###################### DASHBOARD SIDEBAR ######################
  dashboardSidebar(
    sidebarMenu(id = "sdbr",
      
      #---------- Gene of interest input text box -------------#
      textInput("geneInput",                                   
                label = "Enter a gene or miRNA", 
                placeholder = "Example: MSLN"),
      
      #--------- Cohort selection ---------------------------#
      radioGroupButtons("seqDataCohort", choices = c("TARGET", "Beat AML" = "BeatAML"), 
                        status = "primary", label = "Select AML cohort", 
                        selected = "TARGET", size = "sm"),
      
      # --------- Plot generation tabs ------------------------#
      menuItem("Gene expression plots", tabName = "wfPlot", icon = icon("chart-bar")),
      menuItem("Kaplan-Meier curves", tabName = "kmPlot", icon = icon("notes-medical")),
      # menuItem("UMAP", tabName = "umap", icon = icon("spinner")),
      # menuItem("DE Genes", tabName = "deg", icon = icon("spinner")),
      menuItem("External databases", tabName = "extData", icon = icon("atlas")),
      menuItem("Reference info", tabName = "refs", icon = icon("microscope"))
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
      ), 
      
      # Building the external datasets tab that will contain links to other gene expression or protein databases
      tabItem(tabName = "extData",
              mainPanel(
                position = "center",
                fluidRow(
                  infoBoxOutput("protAtlas"),
                  infoBoxOutput("gtex")
                ),
                br(), # Centering on the page
                br(),
                br(),
                fluidRow(
                  # https://renkun-ken.github.io/formattable/ <- Really interesting package for making tables prettier
                  # https://www.displayr.com/formattable/ <- Diff vignette, same package
                  box(title = "ADC and CAR T-cell therapies", status = "info", collapsible = T, solidHeader = T, width = 12,
                      DT::dataTableOutput("therapyTable") # Scrollable table of ADC/CAR T-cell study info from clinicaltrials.gov
                  )
                )
              )
      ),
      
      # tabItem(tabName = "umap",
              # mainPanel(
                  # This works, but messes up the entire dashboard! Prob isn't designed to work with Shiny Dashboard
                  # includeHTML("www/UMAP/TARGET_AMLdx_rlps_NBM_PCAselect_selfcontained.html")
                  
                  # Part of method 1, does not work, no clue why
                  # htmlOutput("test") # I think it's able to access the file, but not display it
              # )
      # ),
      
      tabItem(tabName = "refs", 
              mainPanel(
                position = "center", 
                fluidRow(
                  box(title = "RNA-seq reference information - coming soon")
                )
              ))
    )
  )
)