# UI function for the waterfall plot module
heatmapUI <- function(id, label = "Unsupervised clustering"){
  
  library(DT)
  library(shinyjs)
  library(shinyWidgets)
  ns <- NS(id) # Setting a unique namespace for this module
  
  fluidPage(
            
            sidebarLayout(
              position = "left", 
              
              sidebarPanel(
                helpText("Upload an Excel spreadsheet or CSV file with a single column of gene symbols. 
                         Do not include a column title or header. 
                         Do not inclue any rows/columns other than 1 single column of gene symbols.
                         The file MUST be an Excel or CSV."), # need a place to paste a gene list
                helpText("Placeholder. Might want to describe difference btwn TPM and CPM"), # need a place to paste a gene list
                br(),
                
                fileInput(ns("gene_list"), 
                          label = "Upload gene list here (.csv only)",
                          multiple = F, 
                          accept = c(".csv"),
                          buttonLabel = "Upload"),
                
                checkboxGroupInput(ns("sample_types"), 
                                   label = "Select ?", 
                                   choices = c("AML", "NBM", "CD34+ PB", "MPN", "DS-AML", "TMD"), 
                                   selected = c("AML", "NBM")),
                
                radioButtons(ns("data_type"), 
                             label = "Which form of expression data?", 
                             choices = list("TPM" = "tpm",
                                            "TMM-normalized CPM" = "cpm")),
                
               checkboxInput(ns("log_transform"), 
                             label = "Log2-transform?", 
                             value = F),
              ),
              
              mainPanel(position = "right", 
                        
                        tabsetPanel(
                          tabPanel("Figures",
                                   br(),
                                   fluidRow(
                                     column(10, offset = 0, align = "left",        
                                            plotOutput(ns("plot"), width = "100%")
                                     )
                                   )
                          )
                        )
              )
            )
  )
}

heatmap <- function(input, output, session, clinData, expData, gene, dataset) {
  
  library(tidyverse)
  library(DT)
  library(shinyWidgets)
  library(ComplexHeatmap)

  #################################################################
  #------------------------- FUNCTIONS ---------------------------#
  #################################################################
  
  readGeneList <- reactive({
    
    validate(need(input$gene_list, "Please upload a list of genes to get started."))
    # Checking to make sure the extension is correct
    file <- input$gene_list
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(
      need(ext == "csv", "Please upload a .csv file.")
      )
    
    read.csv(file$datapath, header = F, blank.lines.skip = T, strip.white = T)[,1]
  })
  
  hmMatrix <- reactive({
    gene_list <- readGeneList()
    gene_list <- unique(gene_list)
    
    # There's a bunch of cols of all NA at the far right of the dataframe, need a way to exclude these (must be miRNA-seq only samples)
    mat <- expData()[intersect(gene_list, rownames(expData())),]
    mat <- mat[,colSums(is.na(mat)) == 0] # Removing any cols with an NA
    
    if (input$log_transform == TRUE) {
      mat <- log2(mat + 1)
    }
    mat <- t(scale(t(mat), scale = T, center = T)) # Creating mean-centered, scaled z-scores
    
    return(mat)
  })
  
  plotHeatmap <- reactive({
    
    mat <- hmMatrix()
    
    print("Length of gene list is...")
    print(nrow(mat))
    
    # Only displays row labels (aka gene names) if there are a small number of them
    show_geneLabs <- ifelse(ncol(mat) < 100, TRUE, FALSE)
    show_patIDs <- ifelse(nrow(mat) < 50, TRUE, FALSE)
    
    Heatmap(mat, 
            show_column_names = show_geneLabs,
            show_row_names = show_patIDs, 
            clustering_method_columns = "ward.D2", 
            clustering_method_rows = "ward.D2")
    
  })
  
  
  #################################################################
  #-------------------- FINAL MODULE OUTPUTS ---------------------#
  #################################################################
  
  output$plot <- renderPlot({
    plotHeatmap()
  })

  
}