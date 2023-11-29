# UI function for the waterfall plot module, the id is a unique identifier for this module
oncoprintUI <- function(id, label = "oncoprint") {
  
  # Setting a unique namespace for this module
  ns <- NS(id) 
  
  # Creating a fluid page that will move with window adjustments
  fluidPage(
    
    # Having a sidebar positioned on the left
    sidebarLayout(
      position = "left", 
      
      # Adjusting the panel
      sidebarPanel(
        
        # Instructions for uploading a CSV file
        helpText("Upload a CSV file with a single column of gene symbols. Do not include a column title or header. 
                            Do not include any rows/columns other than 1 single column of gene symbols.
                            If the gene list is > 30, only the first 30 will be taken."),
        
        # Line breaks for readability
        br(),
        
        # File input for uploading gene list in CSV format
        fileInput(ns("gene_list"), 
                  label = "Upload gene list here (.csv only)",
                  multiple = FALSE, 
                  accept = c(".csv"),
                  placeholder = "  No file selected",
                  buttonLabel = "Upload"),
        
        # Slider input for setting the minimum TPMs, this starts at 2 by default
        sliderInput(ns("tpms"), "Minimum TPMs", min = 0, max = 50, value = 2),
        
        # Checkbox group input for selecting fusions
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
                                       "Other AML")),
        
        # Action button for (de)selecting all fusions
        actionButton(ns("selectall"), "(De)Select All", class = "btn-primary btn-sm"),
        
        br(),
        br(),
        
        # Additional annotation options
        helpText("Additional Annotation Options:"),
        
        # Checkbox inputs for various annotation options
        checkboxInput(ns("cnv"), "Primary CNV"),
        checkboxInput(ns("age"), "Age Category"),
        checkboxInput(ns("efs"), "Event Type ID"),
        checkboxInput(ns("flt3"), "FLT3-ITD"),
        checkboxInput(ns("wt1"), "WT1 Mutation"),
        
        # Select input for splitting heatmap by different options (dropdown menu)
        selectInput(ns("split"), label = "Split Heatmap By:", choices = c("No Split", "Fusion", "CNV", "Age", "EFS", "FLT3-ITD", "WT1")),
        
        br(),
        br(),
        
        # Download button for downloading the plot
        downloadButton(ns("plot_download"), 
                      label = "plot", 
                      class = "plotdwnld"),
        
        # Tooltip for download button
        shinyBS::bsTooltip(ns("plot_download"), 
                          title = "Click here to download a copy of the plot",
                          placement = "right", 
                          trigger = "hover"), 
        width = 3, height = 800),

#|---------------------------------------------------------------------------------------------------------------------------|#
      
      # Main panel for displaying the oncoprint plot
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

#|---------------------------------------------------------------------------------------------------------------------------|#

# server function for the oncoprint module
oncoprint <- function(input, output, session, clinData, expData, gene, dataset) {

  # This is a handy function that will add a select all or deselect all button
  observe({
    if(is.null(input$selectall)) return(NULL)
    else if (input$selectall%%2 == 0){
      updateCheckboxGroupInput(session, "fusions", "Fusions:", choices = c("CBFA2T3-GLIS2", 
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
                                                                            "Other AML"))
    }
    
    else {
      updateCheckboxGroupInput(session, "fusions", "Fusions:", choices = c("CBFA2T3-GLIS2", 
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
                               
                                                              selected = c("CBFA2T3-GLIS2", 
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
                                                                            "Other AML"))
      
    }
})
  
  # Reading in the gene list from the csv input file
  readGeneList <- reactive({
    validate(need(input$gene_list, "Please upload a list of genes to get started.")) # Checking to make sure the gene list is there
    file <- input$gene_list 
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext == "csv", "Please upload a .csv file.")) # Checking to make sure the extension is correct
    read.csv(file$datapath, header = F, blank.lines.skip = T, strip.white = T)[,1]
  })
  
  # Creating a dataframe that is subsetted by the input file genelist
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
    
    # this is a function for counting which patients have greater or equal to the chosen tpm cutoff for each gene (in %)
    percentage_expressed <- function(x) {
      return(sum(x >= input$tpms) / length(x) * 100)
    }
    
    # Applying the function to the dataframe
    mat$percentage <- apply(mat, 1, percentage_expressed)
    mat <- mat[order(-mat$percentage),] # ordering the matrix by highest percent expression
    mat <- mat[1:(length(mat)-1)] # removing the percent expressed column
    
    # ifelse conditionals for limiting the input gene list to 30 genes max, otherwise it would not be readable
    if (length(rownames(mat)) > 30) {
      mat <- mat[c(1:30),]
    }
    else {
      mat <- mat
    }
    
    # the actual function for determining whether a patient "expresses" the gene at the TPM cutoff and 
    # will also label the gene intracellular or transmembrane based on Ensembl data (may need updating)
    mat <- apply(mat, c(1,2), function(x) ifelse(x >= input$tpms, "yes", ""))
    mat <- mat[, order(colnames(mat))] # making sure the colnames are ordered correctly
    common_rows <- intersect(rownames(mat), transmembrane_genelist$Gene.name)
    mat[common_rows, ] <- ifelse(mat[common_rows, ] == "yes", "Transmembrane", "")
    mat[setdiff(rownames(mat), common_rows), ] <- ifelse(mat[setdiff(rownames(mat), common_rows), ] == "yes", "Intracellular", "")
    mat <- as.matrix(mat)
    return(mat)
  })
  
  # this function is for the heatmap annotation which will react to the checked input boxes
  oncoAnnotation <- reactive({
    mat <- fusionSubset() # once again rereading the matrix
    updatedclinical <- clinData()[which(clinData()$PatientID %in% colnames(mat)),] # which patients in the CDE are in the matrix
    updatedclinical <- filter(updatedclinical, updatedclinical$Primary.Fusion %in% input$fusions) # this is the reactionary code which will subset the cde for the checked fusions
    updatedclinical <- updatedclinical[order(updatedclinical$PatientID),] # ordering the patients alphabetically to match the matrix
    
    
    # custom color palettes based on colorbrewer that make the oncoprint look cleaner and also make sense with the displayed annotations:
    
    colorlist <- c("#4E79A7", "#E15759", "#EDc948", "#B07AA1", "#F28E2B", "#59A14F", "#9C755F", 
                   "#a5bcd5", "#f0abac", '#f6e4a3', "#d7bcd0", "#f8c695", "#aad3a4", "#cebaae")
    
    cnvcolorlist <- c("#4E79A7", "#E15759", "#EDc948", "#B07AA1", "#F28E2B", "#59A14F", "#9C755F", "#a5bcd5", "#f2f2f2")
    
    efscolorlist <- c("#4E79A7", "#E15759", "#EDc948", "#B07AA1", "#F28E2B", "#f2f2f2")
    
    gradcolorlist <- c("#c4d0e1", "#89a1c4", "#4E79A7", "#344967", "#161f2c", "#f2f2f2")
    
    yesnocolorlist <- c("#4E79A7", "#f0abac", "#f2f2f2")
  
  
    # creating levels in each of the clinical data elements columns so that the information is ordered nicely
    updatedclinical$Primary.Fusion <- factor(updatedclinical$Primary.Fusion, levels = c("CBFA2T3-GLIS2", "CBFB-MYH11", "DEK-NUP214", "ETV6-MNX1", "FUS-ERG", 
                                                                                        "KAT6A-CREBBP", "KMT2A-X", "NUP98-KDM5A", "NUP98-NSD1", "RBM15-MKL1",
                                                                                        "RUNX1-CBFA2T3", "RUNX1-RUNX1T1", "None", "Other AML"))
    updatedclinical$Primary.CNV <- factor(updatedclinical$Primary.CNV, levels = c("CBL deletion", "del5q", "monosomy7", "trisomy8", "trisomy21", "trisomy8/trisomy21", "-Y", "No Relevant CNV", "Unknown"))
    updatedclinical$Age.Category <- factor(updatedclinical$Age.Category, levels = c("Less than 3 years", "Between 3 and 5 years", "Between 5 and 10 years", "Between 10 and 18 years", "Greater than 18 years", "Unknown"))
    updatedclinical$`EFS event type ID` <- factor(updatedclinical$`EFS event type ID`, levels = c("Censored", "Relapse", "Death", "Death without remission", "Induction failure", "Unknown"))
    updatedclinical$FLT3.ITD <- factor(updatedclinical$FLT3.ITD, levels = c("Yes", "No", "Unknown"))
    updatedclinical$WT1.Mutation <- factor(updatedclinical$WT1.Mutation, levels = c("Yes", "No", "Unknown"))
    
    # create a color palette for each annotation, ensuring all levels are covered, flexible based on the number of levels included 
    color_palette1 <- setNames(colorlist[1:length(levels(updatedclinical$Primary.Fusion))], levels(updatedclinical$Primary.Fusion))
    color_palette2 <- setNames(cnvcolorlist[1:length(levels(updatedclinical$Primary.CNV))], levels(updatedclinical$Primary.CNV))
    color_palette3 <- setNames(gradcolorlist[1:length(levels(updatedclinical$Age.Category))], levels(updatedclinical$Age.Category))
    color_palette4 <- setNames(efscolorlist[1:length(levels(updatedclinical$`EFS event type ID`))], levels(updatedclinical$`EFS event type ID`))
    color_palette5 <- setNames(yesnocolorlist[1:length(levels(updatedclinical$FLT3.ITD))], levels(updatedclinical$FLT3.ITD))
    color_palette6 <- setNames(yesnocolorlist[1:length(levels(updatedclinical$WT1.Mutation))], levels(updatedclinical$WT1.Mutation))

    # count the number of annotations based on user input
    num_annotations <- sum(c(input$cnv, input$age, input$efs, input$flt3, input$wt1))
    
    # set up a vector of annotation heights with a default height
    annotation_heights <- c(unit(1.5, "cm"), rep(unit(0.8, "cm"), num_annotations+1))
  
    # create annotations based on user input that is flexible with the annotations added
    ha <- HeatmapAnnotation(
      column_bar = anno_oncoprint_barplot(),
      Fusion = updatedclinical$Primary.Fusion,
      CNV = if (input$cnv) updatedclinical$Primary.CNV else NULL,
      Age = if (input$age) updatedclinical$Age.Category else NULL,
      EFS = if (input$efs) updatedclinical$`EFS event type ID` else NULL,
      FLT3_ITD = if (input$flt3) updatedclinical$FLT3.ITD else NULL,
      WT1 = if (input$wt1) updatedclinical$WT1.Mutation else NULL,
      col = list(
        Fusion = color_palette1,
        CNV = if (input$cnv) color_palette2 else NULL,
        Age = if (input$age) color_palette3 else NULL,
        EFS = if (input$efs) color_palette4 else NULL,
        FLT3_ITD = if (input$flt3) color_palette5 else NULL,
        WT1 = if (input$wt1) color_palette6 else NULL
      ),
      annotation_height = annotation_heights,
      annotation_name_gp = gpar(fontsize = 10),
      gp = gpar(fontsize = 8)
    )
    
    return(ha)
    
  })
  
  # this provides the choices in the drop down menu and maps them to the clinical data elements column names
  radio_choices <- c("No Split", "Fusion", "CNV", "Age", "EFS", "FLT3-ITD", "WT1")
  column_names <- c("", "Primary.Fusion", "Primary.CNV", "Age.Category", "EFS event type ID", "FLT3.ITD", "WT1.Mutation")
  column_mapping <- setNames(column_names, radio_choices)
  
  # plotting the oncoprint (finally!)
  plotOncoprint <- reactive({
    mat <- fusionSubset() # reading in the matrix
    ha <- oncoAnnotation() # reading in the heatmap annotation
    
    updatedclinical <- clinData()[which(clinData()$PatientID %in% colnames(mat)),] # which patients in the CDE are in the matrix
    updatedclinical <- filter(updatedclinical, updatedclinical$Primary.Fusion %in% input$fusions) # this is the reactionary code which will subset the cde for the checked fusions
    updatedclinical <- updatedclinical[order(updatedclinical$PatientID),]
    
    # doing the levelling again so that on the actual oncoprint they get ordered correctly (from left to right) when split
    updatedclinical$Primary.Fusion <- factor(updatedclinical$Primary.Fusion, levels = c("CBFA2T3-GLIS2", "CBFB-MYH11", "DEK-NUP214", "ETV6-MNX1", "FUS-ERG", 
                                                                                        "KAT6A-CREBBP", "KMT2A-X", "NUP98-KDM5A", "NUP98-NSD1", "RBM15-MKL1",
                                                                                        "RUNX1-CBFA2T3", "RUNX1-RUNX1T1", "None", "Other AML"))
    updatedclinical$Primary.CNV <- factor(updatedclinical$Primary.CNV, levels = c("CBL deletion", "del5q", "monosomy7", "trisomy8", "trisomy21", "trisomy8/trisomy21", "-Y", "No Relevant CNV", "Unknown"))
    updatedclinical$Age.Category <- factor(updatedclinical$Age.Category, levels = c("Less than 3 years", "Between 3 and 5 years", "Between 5 and 10 years", "Between 10 and 18 years", "Greater than 18 years", "Unknown"))
    updatedclinical$`EFS event type ID` <- factor(updatedclinical$`EFS event type ID`, levels = c("Censored", "Relapse", "Death", "Death without remission", "Induction failure", "Unknown"))
    updatedclinical$FLT3.ITD <- factor(updatedclinical$FLT3.ITD, levels = c("Yes", "No", "Unknown"))
    updatedclinical$WT1.Mutation <- factor(updatedclinical$WT1.Mutation, levels = c("Yes", "No", "Unknown"))
    
    # setting a color for the bars
    col = c(Transmembrane = "#4E79A7", Intracellular = "#E15759")
    
    # this sets the oncoprint dimensions and rules (how tall the bars are compared to the background, etc.)
    # the grid.rect(w*0.9, h*0.9) is important for setting the bar height
    alter_fun = list(
        Transmembrane = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Transmembrane"], col = NA)),
        Intracellular = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Intracellular"], col = NA))
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
            column_split = if (!is.null(input$split) && input$split != "No Split") updatedclinical[[column_mapping[input$split]]] else NULL, # splits based on the dropdown menu for splitting
            remove_empty_columns = FALSE, # can turn to TRUE if you want to get rid of empty columns (harder to visualize the % though)
            heatmap_legend_param = list(title = "Domain"))
    
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
