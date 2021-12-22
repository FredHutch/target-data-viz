# Check the file below for global scripts & variables, 
# they were removed from the server & ui scripts to help clean them up.
source("global.R")

ui <- dashboardPage(  
  
  dashboardHeader(title = "Meshinchi Lab Data Viz Tools"),
  
  ###################### DASHBOARD SIDEBAR ######################
  dashboardSidebar(
    useShinyalert(), # This line is required to be able to use this package on the server side
    sidebarMenu(id = "sdbr",
      
      #---------- Gene of interest input text box -------------#
      textInput("geneInput",                                   
                label = "Enter a gene or miRNA", 
                placeholder = "Example: MSLN"),
      
      actionButton("check", label = "Not found? Click here!", 
                   style = 'padding:4px; font-size:60%', class = "btn-primary"),
      
      #--------- Cohort selection ---------------------------#
      radioButtons("seqDataCohort", choices = c("TARGET", "Beat AML" = "BeatAML", "SWOG", "TGCA LAML" = "TCGA"), 
                   label = "Select AML cohort", 
                   selected = "TARGET"),
      
      # !!!!!! NOTE !!!!!!! This should be made into a conditional panel for TARGET dataset only
      radioGroupButtons("seqAssembly", choices = c("GRCh38" = "grch38", "GRCh37" = "grch37"), 
                        status = "primary", label = "Select genome assembly", 
                        selected = "grch38", size = "xs"),
      
      # --------- Plot generation tabs ------------------------#
      menuItem("Gene expression plots", tabName = "wfPlot", icon = icon("chart-bar")),
      menuItem("Kaplan-Meier curves", tabName = "kmPlot", icon = icon("notes-medical")),
      menuItem("DE Genes", tabName = "deTable", icon = icon("clipboard-list")), # stream, clipboard-list
      # menuItem("UMAP", tabName = "umap", icon = icon("spinner")),
      # menuItem("Protein Paint", tabName = "protPaint", icon = icon("palette")),
      menuItem("External databases", tabName = "extData", icon = icon("atlas"))
      # menuItem("Reference info", tabName = "refs", icon = icon("dna"))
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
      
      # Sourcing the waterfall plot module UI component
      tabItem(tabName = "deTable",
              
              # Calling the user interface module of the Waterfall Plot app
             deTableUI(id = "degs", label = "Differentially expressed gene table")
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
                  box(title = "ADC and CAR T-cell therapies", collapsible = T, solidHeader = F, width = 12,
                      DT::dataTableOutput("therapyTable") # Scrollable table of ADC/CAR T-cell study info from clinicaltrials.gov
                  )
                )
              )
      )
      
      # tabItem(tabName = "umap",
      #         mainPanel(
      #             # This works, but messes up the entire dashboard! Prob isn't designed to work with Shiny Dashboard
      #             # includeHTML("www/UMAP/TARGET_AML_sg7655_blackBackground_clusters2_k31_PCAselect.html")
      # 
      #             # Part of method 1, does not work, no clue why
      #             # htmlOutput("test") # I think it's able to access the file, but not display it
      #             # Maybe try an iframe instead?
      #             # https://stackoverflow.com/questions/22177974/r-shiny-using-iframe-for-local-files
      #             # https://github.com/rstudio/shiny/issues/2535
      #         )
      # ),
      # 
      # tabItem(tabName = "protPaint",
      #         mainPanel(
      #           # This works!!!!!!
      #           includeHTML("www/Protein_Paint/embed_StJude_ProteinPaint.html")
      #         )
      # )
    )
  )
)