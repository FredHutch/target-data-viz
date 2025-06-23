# NOTES:
# https://www.htmlwidgets.org/showcase_d3heatmap.html <- Interactive heatmaps
# https://cran.r-project.org/web/packages/shinyjqui/vignettes/introduction.html <- Make items resizable!!
# https://github.com/dreamRs/fresh <- Package for building dashboard themes, if I ever decide to get fancy with that

ui <- dashboardPage(
  
  dashboardHeader(title = "Meshinchi Lab DataViz"),
  
  ###################### DASHBOARD SIDEBAR ######################
  dashboardSidebar(
    tags$style(HTML("#geneInput.form-control { color: #FFFFFF; }")), # Changes color of text in the gene text input box
    tags$style(HTML("#geneInput { font-size:13px; height:30px; }")), # Changes size of text & text entry box
    tags$style(HTML("#geneInput.form-control { background-color: #4a4a4a; }")), # Changes background color of text box
    # tags$style(HTML("#gene2.form-control { background-color: white; }")), # Changes background color of text box (this doesn't work though, not sure why)
    # tags$style(HTML("#gene2 { background-color: #white; }")), # Also doesn't work
    tags$style(HTML("#gene2.form-control { color: #2096f6; }")),
    tags$style(HTML(".help-block { color: #787878 }")), # Makes help text a bit darker & easier to read
    tags$head(tags$style(".plotdwnld { vertical-align:middle; horizontal-align:middle; height:30px; width:75px; font-size:10px; padding:2px }")),
    tags$head(tags$style(HTML(".fa{font-size: 18px;}"))), # Makes dashboard icons larger
    sidebarMenu(
                
                #---------- Gene of interest input text box -------------#
                textInput("geneInput",
                          label = "Enter a gene or miRNA",
                          placeholder = " Example: MSLN"),
                
                div(
                  style = "display: flex; gap: 60px;",
                  actionButton("ensid", label = "SHOW ENSID", 
                               style = 'padding:2px; font-size:70%', 
                               class = "btn-primary"),
                  actionButton("check", label = "Not found?", 
                               style = 'padding:2px; font-size:70%; background: #FD7370; color: #FFFFFF', 
                               class = "btn")
                ),
                
                conditionalPanel(
                  condition = "input.expDataCohort == 'TARGET' && input.aligner == 'star' && output.showSelectInput",
                  selectInput("id_select", "Duplicate Gene Names Found, Choose ENSID:", choices = NULL)
                ),
                
                # the following functions are to create checkmarks for whether the input gene is aml-restricted and transmembrane
                #-------------------------------------------------------------#
                conditionalPanel(
                  condition = "!input.geneInput", #if there is no gene inputted
                  tags$div(style = "margin-left: 13px; margin-top: 16px; color: white; font-size: 14px;", 
                           icon("question", style = "font-size: 12px;"), "   AML-restricted")
                ), 
                conditionalPanel(
                  condition = "input.geneInput && !output.gene_present", #if there is a gene inputted but it is not found in the aml-restricted list
                  tags$div(style = "margin-left: 13px; margin-top: 16px; color: #FD7370; font-size: 14px;", 
                           icon("times", style = "font-size: 12px;"), "   AML-restricted")
                ), 
                conditionalPanel(
                  condition = "input.geneInput && output.gene_present", #if there is a gene inputted and it is found in the aml-restricted list
                  tags$div(style = "margin-left: 13px; margin-top: 16px; color: #2096f6; font-size: 14px;", 
                           icon("check", style = "font-size: 12px;"), "   AML-restricted")
                ),
                conditionalPanel(
                  condition = "!input.geneInput", 
                  tags$div(style = "margin-left: 13px; margin-top: 1px; color: white; font-size: 14px;", 
                           icon("question", style = "font-size: 12px;"), "   Transmembrane")
                ), 
                conditionalPanel(
                  condition = "input.geneInput && !output.trmembrane",
                  tags$div(style = "margin-left: 13px; margin-top: 1px; color: #FD7370; font-size: 14px;", 
                           icon("times", style = "font-size: 12px;"), "   Transmembrane")
                ),
                conditionalPanel(
                  condition = "input.geneInput && output.trmembrane",
                  tags$div(style = "margin-left: 13px; margin-top: 1px; color: #2096f6; font-size: 14px;", 
                           icon("check", style = "font-size: 12px;"), "   Transmembrane")
                ),
                #-------------------------------------------------------------#
                
                shinyBS::bsTooltip("check", title = "Click here for alias suggestions",
                                   placement = "right", 
                                   trigger = "hover"),
                
                #-------- Disease selection -----------------------------#
                
                radioGroupButtons("leukemiaSelection", choices = c("AML", "BALL", "TALL", "CCLE", "PCGP"), 
                                  status = "primary", label = "Select cancer", 
                                  selected = "AML", size = "xs"),
                
                #--------- Cohort selection -----------------------------#
                
                radioButtons("expDataCohort", choices = c("TARGET", "BEAT" = "BeatAML", "SWOG", "TGCA" = "TCGA", "LEUCEGENE", "PCGP AML"), 
                             label = "Select cohort", 
                             selected = "TARGET"),
                
                conditionalPanel("input['expDataCohort'] == 'TARGET'",
                                 radioGroupButtons("aligner", choices = c("STAR" = "star", "Kallisto" = "kallisto"), 
                                                   status = "primary", label = "Select alignment method", 
                                                   selected = "star", size = "xs"),
                                 
                ),
                
                conditionalPanel("input['expDataCohort'] == 'TARGET' && input['aligner'] == 'star'", 
                                 radioGroupButtons("timepoint", choices = c("Dx" = "diagnostic", "Rem" = "remission", "Rel" = "relapse"), 
                                                   status = "primary", label = "Select timepoint", 
                                                   selected = "diagnostic", size = "xs")
                ),
                
                ## a button for controlling the genome build
                # conditionalPanel("input['expDataCohort'] == 'TARGET'",
                #                  radioGroupButtons("seqAssembly", choices = c("GRCh38" = "grch38", "GRCh37" = "grch37"), 
                #                                    status = "primary", label = "Select genome assembly", 
                #                                    selected = "grch38", size = "xs")
                #                  
                # ),
                
                br(),
                
                # --------- Plot generation tabs ------------------------#
                menuItem("Expression plots", tabName = "wfPlot", icon = icon("chart-column")),
                menuItem("Timepoints", tabName = "timePlot", icon = icon("stopwatch")),
                menuItem("Gene expressors", tabName = "geneExp", icon = icon("chart-pie")),
                menuItem("Kaplan-Meier curves", tabName = "kmPlot", icon = icon("notes-medical")),
                # menuItem("Cox models", tabName = "coxPH"),
                menuItem("Oncoprints", tabName = "oncoprint", icon = icon("stream")),
                # menuItem("Risk Classification", tabName = "Classi", icon = icon("exclamation-circle")),
                # menuItem("Heatmaps", tabName = "heatmap", icon = icon("th")),
                #menuItem("DE Genes", tabName = "deTable", icon = icon("clipboard-list")),
                menuItem("UMAP", tabName = "umap", icon = icon("spinner")),
                menuItem("External databases", tabName = "extData", icon = icon("atlas")),
                menuItem("HPA Info", tabName = "HPA", icon = icon("dna")),
                menuItem("Other Cancers", tabName = "cancertype", icon = icon("disease")),
                
                br(),
                p("DataViz v2.0.0", style = "padding-left: 15px; font-style: italic;")
    )
  ),  
  
  ###################### DASHBOARD PAGES #######################
  dashboardBody(
    
    tags$head(
      tags$style(HTML("
    /* Header Gradient */
    .main-header .logo {
      font-family: 'Tahoma', sans-serif;
      font-weight: 800;
      font-size: 17px;
      letter-spacing: 0px;
      color: #ffffff !important;
      text-align: left;
      background: linear-gradient(to right, #2096f6, #2096f6);
      padding-left: 15px;
      text-shadow: 1px 1px 2px rgba(0, 0, 0, 0.2);
      box-shadow: inset -2px 0 0 rgba(255, 255, 255, 0.1);
    }

    .main-header .navbar {
      background: linear-gradient(90deg,rgba(32, 150, 246, 1) 0%, rgba(253, 115, 112, 1) 100%);
    }
    
    .main-header .sidebar-toggle: {
      text-decoration: none !important;
    }
    
    .main-header .sidebar-toggle:hover {
      background-color: #52ADF8 !important;
      color: #2096f6 !important; /* Optional: keep the icon/text visible */
      text-decoration: none !important;
    }

    /* Sidebar Gradient: black to dark gray */
    .main-sidebar {
      background: linear-gradient(to left, #000000, #2f2f2f);
      color: #ffffff;
    }

    /* Sidebar text and icons */
    .main-sidebar .sidebar a {
      color: #ffffff;
    }
    
    .main-sidebar .sidebar .sidebar-menu a > i {
      color: #2096f6 !important;
    }

    .main-sidebar .sidebar-menu > li.active > a {
      background: linear-gradient(90deg,rgba(32, 150, 246, 1) 0%, rgba(253, 115, 112, 1) 100%);
      color: #ffffff !important;
      border-left-color: #2096f6;
    }
    
    /* Make icon white in active tab */
    .main-sidebar .sidebar-menu > li.active > a > i {
      color: #ffffff !important;
    }

    .main-sidebar .sidebar-menu > li:hover > a {
      background-color: rgba(255, 255, 255, 0.1) !important;
      color: #ffffff;
      border-left-color: #FD7370;
    }

    .sidebar-menu .treeview-menu > li > a {
      color: #cccccc;
    }

    .sidebar-menu .treeview-menu > li > a:hover {
      color: #ffffff;
    }
    
      /* Custom focus styles for buttons, nav-tabs, and sidebar items */
      button:focus,
      .nav-tabs > li > a:focus,
      .nav-tabs > li > a:focus-visible,
      .sidebar-menu a:focus {
        outline: 0px solid #2096f6 !important;
        text-decoration: none !important;
      }
      .small-box {
        transition: background-color 0.3s ease, transform 0.2s ease;
      }
      
      .small-box:hover {
        background-color: #FD7370 !important; /* red on hover */
          cursor: pointer;
      }
      
      .small-box:active {
        background-color: #b40603 !important; /* darker red on click */
      }
                      
      
      
  "))
    ),
    
    
    # Using some custom CSS to...
    tags$head(tags$style(HTML('.content-wrapper, .right-side { background-color: #ffffff; },'))), # Change the background color to white
    tags$head(tags$style(HTML('.shiny-output-error-validation { color: #2096f6; }'))), # Modify color of app error messages
    # tags$head(tags$style(HTML('.box{-webkit-box-shadow: none; -moz-box-shadow: none;box-shadow: none;}'))), # Remove border around boxes
    
    tabItems(
      
      # Sourcing the waterfall plot module UI component
      tabItem(tabName = "wfPlot",
              wfPlotUI(id = "waterfall", label = "Waterfall plot generation") 
      ),
      
      tabItem(tabName = "timePlot",
              timePlotUI(id = "timepoints", label = "Timepoint plot generation") 
      ),
      
      tabItem(tabName = "kmPlot",
              kmPlotUI(id = "kaplanmeier", label = "Kaplan-Meier plot generation")
      ),
      
      tabItem(tabName = "oncoprint",
              oncoprintUI(id = "oncoprint", label = "Oncoprint generation")
      ), 
      
      # This module is not ready for prime time yet
      # tabItem(tabName = "heatmap",
      # #Calling the user interface module of the Waterfall Plot app
      #         heatmapUI(id = "heatmap", label = "Heatmap generation")
      # ),
      
      # tabItem(tabName = "deTable",
      #         deTableUI(id = "degs", label = "Differentially expressed gene table")
      # ), 
      
      # Sourcing the gene expressor module UI component
      tabItem(tabName = "geneExp",
              geneExpUI(id = "exps", label = "Identify gene-positive cases")
      ), 
      
      tabItem(tabName = "HPA",
              HPAPlotUI(id = "hpa", label = "HPA Supporting Info")
      ), 
      
      # tabItem(tabName = "Classi",
      #         ClassiPlotUI(id = "Classi", label = "Risk Classification")
      # ), 

      tabItem(tabName = "cancertype",
              CancerPlotUI(id = "cancertype", label = "Cancer Type")
      ),
      
      # Building the external datasets tab that will contain links to other gene expression or protein databases
      tabItem(tabName = "extData",
              mainPanel(
                width = 12,
                position = "center",
                fluidRow(
                  column(3, uiOutput("protAtlas")),
                  column(3, uiOutput("gtex")),
                  column(3, uiOutput("protPaint")),
                  column(3, uiOutput("cbioportal"))
                ),
                fluidRow(
                  column(3, uiOutput("deeptmhmm")),
                  column(3, uiOutput("expasy")),
                  column(3, uiOutput("ucsc")),
                  column(3, uiOutput("vizrisk")),
                ),
                # fluidRow(
                #   column(width = 4,
                #       uiOutput("tmhmm")
                #   ),
                #   column(width = 4,
                #       verbatimTextOutput("terminal_output")
                #   ),
                #   column(width = 4,
                #      div(
                #        style = "overflow-y: scroll; text-align: center;",
                #        imageOutput("tmhmm_plot")
                #      )
                #   )
                # ), # Linebreaks to center the table on the page
                fluidRow(
                  # https://renkun-ken.github.io/formattable/ <- Really interesting package for making tables prettier
                  # https://www.displayr.com/formattable/ <- Diff vignette, same package
                  box(title = "ADC and CAR T-cell therapies", 
                      collapsible = T, 
                      solidHeader = F,
                      width = 12, 
                      DT::dataTableOutput("therapyTable") # Scrollable table of ADC/CAR T-cell study info from clinicaltrials.gov
                  )
                )
              )
      ),
      
      tabItem(tabName = "umap",
              mainPanel(
                # This works, but messes up the entire dashboard! Not sure why.
                # includeHTML("www/UMAP/TARGET_AML_sg7655_blackBackground_clusters2_k31_PCAselect.html")
                
                # Different method. This produces the UI side of the iframe
                htmlOutput("umapEmbedding")
              )
      )
      # 
      # tabItem(tabName = "protPaint",
      #         sidebarPanel(
      #           helpText("Please enter the gene or region of interest in the text box to the right."),
      #           helpText("The figure was created with the ProteinPaint visualization tool found at the St. Jude PeCan Portal: https://pecan.stjude.cloud/"),
      #           helpText("Original publication in Nature Genetics: https://www.nature.com/articles/ng.3466")
      #         ),
      #         mainPanel(
      #           # The 'includeHTML' command below works, but clips off the edges of the final embedded page
      #           # includeHTML("www/Protein_Paint/embed_StJude_ProteinPaint_writeTest.html")
      #           
      #           # Displaying the HTML in an iframe from the server side works better.
      #           # More info on embedding HTML from Protein Paint:
      #           # https://stjudecloud.github.io/docs/guides/proteinpaint/developers-guide/embedding-proteinpaint/
      #           
      #           # The problem: I can't get it to populate w/ the same gene as the user has entered in the text box...
      #           # I've tried to do that by manipulating the file on the server-side, but then the embedded page doesn't
      #           # display properly. This still needs work.
      #           htmlOutput(outputId = "htmlDisplay")
      #         
      #         )
      # )
    )
  )
)
