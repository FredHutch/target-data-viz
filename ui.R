# NOTES:
# https://www.htmlwidgets.org/showcase_d3heatmap.html <- Interactive heatmaps
# https://cran.r-project.org/web/packages/shinyjqui/vignettes/introduction.html <- Make items resizable!!
# https://github.com/dreamRs/fresh <- Package for building dashboard themes, if I ever decide to get fancy with that

ui <- dashboardPage(
  
  dashboardHeader(title = "Meshinchi Lab Data Viz Tools"),
  
  ###################### DASHBOARD SIDEBAR ######################
  dashboardSidebar(
    tags$style(HTML("#geneInput.form-control { color: #93C54B; }")), # Changes color of text in the gene text input box
    tags$style(HTML("#geneInput { font-size:13px; height:30px; }")), # Changes size of text & text entry box
    tags$style(HTML("#geneInput.form-control { background-color: #4a4a4a; }")), # Changes background color of text box
    # tags$style(HTML("#gene2.form-control { background-color: white; }")), # Changes background color of text box (this doesn't work though, not sure why)
    # tags$style(HTML("#gene2 { background-color: #white; }")), # Also doesn't work
    tags$style(HTML("#gene2.form-control { color: #93C54B; }")),
    tags$style(HTML(".help-block { color: #787878 }")), # Makes help text a bit darker & easier to read
    tags$head(tags$style(".plotdwnld { vertical-align:middle; horizontal-align:middle; height:30px; width:75px; font-size:10px; padding:2px }")),
    tags$head(tags$style(HTML(".fa{font-size: 18px;}"))), # Makes dashboard icons larger
    sidebarMenu(
                
                #---------- Gene of interest input text box -------------#
                textInput("geneInput",
                          label = "Enter a gene or miRNA",
                          placeholder = " Example: MSLN"),
                
                actionButton("check", label = "Not found?", 
                             style = 'padding:2px; font-size:70%', # Extra style CSS makes the button smaller
                             class = "btn-primary"),
                
                # the following functions are to create checkmarks for whether the input gene is aml-restricted and transmembrane
                #-------------------------------------------------------------#
                conditionalPanel(
                  condition = "!input.geneInput", #if there is no gene inputted
                  tags$div(style = "margin-left: 13px; margin-top: 16px; color: white; font-size: 14px; font-family: sans-serif;", 
                           icon("question", style = "font-size: 12px;"), "   AML-restricted")
                ), 
                conditionalPanel(
                  condition = "input.geneInput && !output.gene_present", #if there is a gene inputted but it is not found in the aml-restricted list
                  tags$div(style = "margin-left: 13px; margin-top: 16px; color: #F47174; font-size: 14px; font-family: sans-serif;", 
                           icon("times", style = "font-size: 12px;"), "   AML-restricted")
                ), 
                conditionalPanel(
                  condition = "input.geneInput && output.gene_present", #if there is a gene inputted and it is found in the aml-restricted list
                  tags$div(style = "margin-left: 13px; margin-top: 16px; color: #93C54B; font-size: 14px; font-family: sans-serif;", 
                           icon("check", style = "font-size: 12px;"), "   AML-restricted")
                ),
                conditionalPanel(
                  condition = "!input.geneInput", 
                  tags$div(style = "margin-left: 13px; margin-top: 1px; color: white; font-size: 14px; font-family: sans-serif;", 
                           icon("question", style = "font-size: 12px;"), "   Transmembrane")
                ), 
                conditionalPanel(
                  condition = "input.geneInput && !output.trmembrane",
                  tags$div(style = "margin-left: 13px; margin-top: 1px; color: #F47174; font-size: 14px; font-family: sans-serif;", 
                           icon("times", style = "font-size: 12px;"), "   Transmembrane")
                ),
                conditionalPanel(
                  condition = "input.geneInput && output.trmembrane",
                  tags$div(style = "margin-left: 13px; margin-top: 1px; color: #93C54B; font-size: 14px; font-family: sans-serif;", 
                           icon("check", style = "font-size: 12px;"), "   Transmembrane")
                ),
                #-------------------------------------------------------------#
                
                shinyBS::bsTooltip("check", title = "Click here for alias suggestions",
                                   placement = "right", 
                                   trigger = "hover"),
                
                #-------- Disease selection -----------------------------#
                
                radioGroupButtons("leukemiaSelection", choices = c("AML", "ALL", "TALL", "CCLE"), 
                                  status = "primary", label = "Select Leukemia", 
                                  selected = "AML", size = "xs"),
                
                #--------- Cohort selection -----------------------------#
                
                radioButtons("expDataCohort", choices = c("TARGET", "Beat AML" = "BeatAML", "SWOG", "TGCA LAML" = "TCGA"), 
                             label = "Select Cohort", 
                             selected = "TARGET"),
                
                conditionalPanel("input['expDataCohort'] == 'TARGET'",
                                 radioGroupButtons("seqAssembly", choices = c("GRCh38" = "grch38", "GRCh37" = "grch37"), 
                                                   status = "primary", label = "Select genome assembly", 
                                                   selected = "grch38", size = "xs")
                ),
                
                br(),
                
                # --------- Plot generation tabs ------------------------#
                menuItem("Expression plots", tabName = "wfPlot", icon = icon("chart-bar")),
                menuItem("Gene expressors", tabName = "geneExp", icon = icon("chart-pie")),
                menuItem("Kaplan-Meier curves", tabName = "kmPlot", icon = icon("notes-medical")),
                # menuItem("Cox models", tabName = "coxPH"),
                menuItem("Oncoprints", tabName = "oncoprint", icon = icon("stream")),
                menuItem("Risk Classification", tabName = "Classi", icon = icon("exclamation-circle")),
                # menuItem("Heatmaps", tabName = "heatmap", icon = icon("th")),
                menuItem("DE Genes", tabName = "deTable", icon = icon("clipboard-list")),
                menuItem("UMAP", tabName = "umap", icon = icon("spinner")),
                menuItem("External databases", tabName = "extData", icon = icon("atlas")),
                menuItem("HPA Info", tabName = "HPA", icon = icon("dna"))
    )
  ),  
  
  ###################### DASHBOARD PAGES #######################
  dashboardBody(
    
    # Using some custom CSS to...
    tags$head(tags$style(HTML('.content-wrapper, .right-side { background-color: #ffffff; },'))), # Change the background color to white
    tags$head(tags$style(HTML('.shiny-output-error-validation { color: #93C54B; }'))), # Modify color of app error messages
    # tags$head(tags$style(HTML('.box{-webkit-box-shadow: none; -moz-box-shadow: none;box-shadow: none;}'))), # Remove border around boxes
    
    tabItems(
      
      # Sourcing the waterfall plot module UI component
      tabItem(tabName = "wfPlot",
              wfPlotUI(id = "waterfall", label = "Waterfall plot generation") 
      ),
      
      # Sourcing the waterfall plot module UI component
      tabItem(tabName = "kmPlot",
              kmPlotUI(id = "kaplanmeier", label = "Kaplan-Meier plot generation")
      ),
      
      tabItem(tabName = "oncoprint",
              oncoprintUI(id = "oncoprint", label = "Oncoprint generation")
      ), 
      
      # Sourcing the waterfall plot module UI component
      # This module is not ready for prime time yet
      # tabItem(tabName = "heatmap",
      # #Calling the user interface module of the Waterfall Plot app
      #         heatmapUI(id = "heatmap", label = "Heatmap generation")
      # ),
      
      # Sourcing the waterfall plot module UI component
      tabItem(tabName = "deTable",
              deTableUI(id = "degs", label = "Differentially expressed gene table")
      ), 
      
      # Sourcing the gene expressor module UI component
      tabItem(tabName = "geneExp",
              geneExpUI(id = "exps", label = "Identify gene-positive cases")
      ), 
      
      tabItem(tabName = "HPA",
              HPAPlotUI(id = "hpa", label = "HPA Supporting Info")
      ), 
      
      tabItem(tabName = "Classi",
              ClassiPlotUI(id = "Classi", label = "Risk Classification")
      ), 
      
      # Building the external datasets tab that will contain links to other gene expression or protein databases
      tabItem(tabName = "extData",
              mainPanel(
                position = "center",
                fluidRow(
                  valueBoxOutput("protAtlas", width = 3),
                  valueBoxOutput("gtex", width = 3),
                  valueBoxOutput("protPaint", width = 3)
                ),
                br(), # Linebreaks to center the table on the page
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
