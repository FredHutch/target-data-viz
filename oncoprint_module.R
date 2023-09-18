# UI function for the waterfall plot module
oncoprintUI <- function(id, label = "oncoprint"){
  
  ns <- NS(id) # Setting a unique namespace for this module
  
  fluidPage(
    sidebarLayout(
      position = "left", 
      sidebarPanel(
        helpText("Upload a CSV file with a single column of gene symbols. Do not include a column title or header. 
                            Do not inclue any rows/columns other than 1 single column of gene symbols. The file MUST be a CSV. 
                            If the gene list is > 30, only the first 30 will be taken."),
        
        br(),
        
        # requests a gene_list in the form of a .csv for inputting
        fileInput(ns("gene_list"), 
                  label = "Upload gene list here (.csv only)",
                  multiple = F, 
                  accept = c(".csv"),
                  placeholder = "  No file selected",
                  buttonLabel = "Upload"),
        
        # checkboxes for fusion subsetting
        checkboxGroupInput(ns("fusions"),
                           label = "Fusions:",
                           choices = c("CBFA2T3-GLIS2", 
                                       "CBFB-MYH11",
                                       "DEK-NUP214",
                                       "ETV6-MNX1",
                                       "FUS-ERG",
                                       "KAT6A-CREBBP",
                                       "KMT2A-X",
                                       "NUP98-KDM5A",
                                       "NUP98-NSD1",
                                       "RBM15-MKL1",
                                       "RUNX1-CBFA2T3",
                                       "RUNX1-RUNX1T1",
                                       "None",
                                       "Other AML"),
                           selected = c("CBFA2T3-GLIS2")
        ),
      
        br(),
        downloadButton(ns("plot_download"), 
                      label = "plot", 
                      class = "plotdwnld"),
      
        shinyBS::bsTooltip(ns("plot_download"), 
                          title = "Click here to download a copy of the plot",
                          placement = "right", 
                          trigger = "hover"), width = 3, height = 800),
      
      #|---------------------------------------------------------------------------------------------------------------------------|#
      
      mainPanel(width = 9, height = 900, position = "right", 
        tabsetPanel(
          tabPanel("Figures",
          br(),
            fluidRow(
              column(10, offset = 0, align = "left",        
                plotOutput(ns("plot"), height = 600)
              )
            )
          )
        )
      )
    )
  )
}

oncoprint <- function(input, output, session, clinData, expData, gene, dataset) {
  
  #################################################################
  #------------------------- FUNCTIONS ---------------------------#
  #################################################################
  
  #reading in the gene list from the csv input file
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
  
  # creating a dataframe that is subsetted by the input file genelist
  oncoMatrix <- reactive({
    gene_list <- readGeneList() # reads the gene list in
    gene_list <- unique(gene_list) # makes sure they're unique
    mat <- filter(expData(), rownames(expData()) %in% gene_list) # filtering the counts file for the gene list
    return(mat)
  })
    
  # filtering the matrix by fusion group (samples)
  fusionSubset <- reactive({
    mat <- oncoMatrix() # the last function returns the dataframe so this reads it in again
    updatedclinical <- filter(clinData(), clinData()$Primary.Fusion %in% input$fusions) # filters the cde by the checked fusions
    mat <- mat[,which(colnames(mat) %in% updatedclinical$PatientID)] # filtering the dataframe for only the patients we have CDE for
    mat <- mat[,colSums(is.na(mat)) == 0] # removing NAs
    
    # this is a function for counting which patients have >5 tpm for each gene (in %)
    percentage_expressed <- function(x) {
      return(sum(x >= 5) / length(x) * 100)
    }
    
    mat$percentage <- apply(mat, 1, percentage_expressed) # applying the function to the dataframe
    mat <- mat[order(-mat$percentage),] # ordering the matrix by highest percent expression
    mat <- mat[1:(length(mat)-1)] # removing the percent expressed column
    
    # ifelse conditionals for limiting the input gene list to 30 genes max, otherwise it would not be readable
    if (length(rownames(mat)) > 30) {
      mat <- mat[c(1:30),]
    }
    else {
      mat <- mat
    }

    # the actual function for determining whether a patient "expresses" the gene at 5 TPM or more
    mat <- apply(mat, c(1,2), function(x) ifelse(x >= 5, "yes", ""))
    mat <- mat[, order(colnames(mat))] # making sure the colnames are ordered correctly
    mat <- as.matrix(mat) # turning the dataframe into a matrix for the oncoprint
    return(mat)
  })
  
  # this function is for the heatmap annotation which will react to the checked input boxes
  oncoAnnotation <- reactive({
    mat <- fusionSubset() # once again rereading the matrix
    updatedclinical <- clinData()[which(clinData()$PatientID %in% colnames(mat)),] # which patients in the CDE are in the matrix
    updatedclinical <- filter(updatedclinical, updatedclinical$Primary.Fusion %in% input$fusions) # this is the reactionary code which will subset the cde for the checked fusions
    updatedclinical <- updatedclinical[order(updatedclinical$PatientID),] # ordering the patients alphabetically to match the matrix
    
    # this is a heatmap annotation
    ha <- HeatmapAnnotation(
      column_bar = anno_oncoprint_barplot(), # keeps the frequency barplot at the top
      Fusion = updatedclinical$Primary.Fusion, # adds the fusion bar in the heatmap annotation
        col = list(Fusion = c("CBFA2T3-GLIS2" = "#e6194b", # manually assigning colors from the pals library
                            "CBFB-MYH11" = "#3cb44b", 
                            "DEK-NUP214" = "#ffe119", 
                            "ETV6-MNX1" = "#4363d8", 
                            "FUS-ERG" = "#f58231", 
                            "KAT6A-CREBBP" = "#911eb4", 
                            "KMT2A-X" = "#42d4f4",
                            "NUP98-KDM5A" = "#f032e6",
                            "NUP98-NSD1" = "#bfef45",
                            "RBM15-MKL1" = "#fabed4",
                            "RUNX1-CBFA2T3" = "#469990",
                            "RUNX1-RUNX1T1" = "#dcbeff",
                            "None" = "#9a6324",
                            "Other AML" = "#fffac8")),
      annotation_height = unit.c(unit(1.5, "cm"), unit(0.8, "cm")), # changes the annotation heights
      annotation_name_gp = gpar(fontsize = 10), # changes the fontsize for annotation label
      gp = gpar(fontsize = 8) # not sure what this does to be honest (I think sets a baseline fontsize but not sure if it needs to go or stay)
    )
  
    return(ha)
  })
  
  
  # plotting the oncoprint (finally!)
  plotOncoprint <- reactive({
    mat <- fusionSubset() # reading in the matrix
    ha <- oncoAnnotation() # reading in the heatmap annotation
     
    col = c(yes = "#2096f3") # setting a color for the bars
    
    # this sets the oncoprint dimensions and rules (how tall the bars are compared to the background, etc.)
    # the grid.rect(w*0.9, h*0.9) is important for setting the bar height
    alter_fun = list(
        yes = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["yes"], col = NA))
    )
    
    # the oncoprint plotting recipe
    d1 <- oncoPrint(mat, # the matrix
            alter_fun_is_vectorized = TRUE, # critical otherwise it'll take forever to plot
            alter_fun = alter_fun, # the function above
            col = col, # the color above
            row_names_gp = gpar(fontsize = 10), # gene fontsize
            pct_gp = gpar(fontsize = 10), # % fontsize
            column_title = NULL, # don't need a title
            top_annotation = ha, # heatmap annotation from the previous function
            remove_empty_columns = FALSE, # can turn to TRUE if you want to get rid of empty columns (harder to visualize the % though)
            heatmap_legend_param = list(title = "TPM >= 5")) # setting the legend name
    
    d1 <- draw(d1, merge_legend = TRUE) # this keeps the legends clean
    
  })
  
  
  
  #################################################################
  #-------------------- FINAL MODULE OUTPUTS ---------------------#
  #################################################################
  
  output$plot <- renderPlot({
    validate(
      need(input$fusions, 'Please select at least one fusion')) # this validation makes sure that a box is checked
      plotOncoprint()
  })
  
  # Adding a download button widget for the plot
  output$plot_download <- downloadHandler(
    filename = function() {
      paste0("oncoprint_", format(Sys.time(), "%m.%d.%Y"), ".pdf")
    }, 
    content = function(file) {
      pdf(file = file, width = 10, height = 8)
      d1 <- plotOncoprint()
      print(d1)
      dev.off()
    }
  )
  
}
