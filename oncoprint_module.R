# UI function for the waterfall plot module, the id is a unique identifier for this module
oncoprintUI <- function(id, label = "oncoprint") {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    tags$head(
      tags$style(HTML("
          .sidebar-container {
            display: flex;
            height: 100vh;
          }
          .custom-sidebar {
            background-color: #f8f9fa;
            padding: 15px;
            width: 250px;
            flex-shrink: 0;
            overflow-y: auto;
          }
          .main-content {
            flex-grow: 1;
            padding: 15px;
          }
          .plotdwnld, .custom-file-upload .btn {
            background-color: #2096f6 !important;
            color: white !important;
            width: 100% !important;
            padding: 3px 7px;
            font-size: 1.2rem;
          }
          .custom-file-upload .form-control {
            font-size: 1.2rem;
            padding: 3px 7px;
          }
         .js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #FD7370}
        "))
    ),
    fluidPage(
      useShinyjs(),
      theme = shinythemes::shinytheme(theme = "paper"),
      div(class = "sidebar-container",
          
          #------------------ SIDEBAR ------------------#
          div(class = "custom-sidebar",

              radioButtons(ns("gene_input_method"),
                           label = "Gene Input Method:",
                           choices = c("Manual Entry" = "text", "Upload CSV" = "file"),
                           selected = "text",
                           inline = TRUE),
              
              conditionalPanel(
                condition = sprintf("input['%s'] == 'text'", ns("gene_input_method")),
                textAreaInput(ns("gene_text"),
                              label = "Enter gene symbols (comma or newline separated):",
                              placeholder = "e.g. MSLN, PTPN11, WT1",
                              rows = 5)
              ),
              
              conditionalPanel(
                condition = sprintf("input['%s'] == 'file'", ns("gene_input_method")),
                div(class = "custom-file-upload",
                    fileInput(ns("gene_list"), 
                              label = "Upload gene list (.csv only, one column, no header)",
                              multiple = FALSE, 
                              accept = c(".csv"),
                              placeholder = "   No file selected",
                              buttonLabel = "Upload")
                )
                
              ),
              
              
              sliderInput(ns("tpms"), "Minimum TPMs", min = 0, max = 50, value = 2),
              
              div(id = ns("fusion_group_container"),
                  checkboxGroupInput(ns("fusions"),
                                     label = "Fusions:",
                                     choices = c("CBFA2T3-GLIS2",
                                                 "CBFB-MYH11",
                                                 "DEK-NUP214",
                                                 "ETS-Family",
                                                 "KMT2A-X",
                                                 "MLLT10-X",
                                                 "NUP98-KDM5A",
                                                 "NUP98-NSD1",
                                                 "RUNX1-RUNX1T1",
                                                 "Other AML",
                                                 "None")),
                  actionButton(ns("selectall"), "(De)Select All", class = "btn-primary btn-sm w-100"),
                  br(), br()
              ),
              
              
              selectInput(ns("annot"), label = "Additional Annotation Options:", multiple = TRUE,
                          choices = c("Primary CNV", "Age Category", "Event Type ID", "FLT3-ITD", "WT1")),
              
              selectInput(ns("split"), label = "Split Heatmap By:",
                          choices = c("No Split", "Fusion", "Primary CNV", "Age Category", "Event Type ID", "FLT3-ITD", "WT1"),
                          selected = "No Split"),
              
              br(),
              
              div(style = "margin-bottom: 10px;",
                  downloadButton(ns("plot_download"),
                                 label = "Download Plot",
                                 class = "btn-primary btn-sm w-100")),
              shinyBS::bsTooltip(ns("plot_download"),
                                 title = "Click here to download a copy of the plot",
                                 placement = "right",
                                 trigger = "hover"),
          ),
          
          #------------------ MAIN PANEL ------------------#
          div(class = "main-content",
              tabsetPanel(
                tabPanel("Figures",
                         br(),
                         fluidRow(
                           column(12, align = "left",
                                  plotOutput(ns("plot"), height = "70vh")
                           )
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
oncoprint <- function(input, output, session, clinData, expData, dataset, aligner) {
  

  
  # Reactive value to store current fusion choices
  fusion_choices <- reactiveVal()
  
  observeEvent(dataset(), {
    
    if (dataset() %in% c("BeatAML", "SWOG", "TCGA", "LEUCEGENE")) {
      shinyjs::hide("fusion_group_container")
    } else {
      shinyjs::show("fusion_group_container")
    }
    
    clinical <- clinData()
    
    if (dataset() %in% c("TARGET", "PCGP AML")) {
      
      # Clean up Primary.Fusion column for PCGP
      if (dataset() == "PCGP AML") {
        clinical$Primary.Fusion <- ifelse(
          clinical$Primary.Fusion == "No Fusion", "None",
          ifelse(clinical$Primary.Fusion == "Rare Fusion", "Other AML", clinical$Primary.Fusion)
        )
        
        updateSelectInput(session = session,
                          inputId = "annot", 
                          label = "Additional Annotation Options:", 
                          choices = c("Age Category", "CEBPA", "FLT3-ITD", "NPM1"))
        
        updateSelectInput(session = session,
                          inputId = "split", 
                          label = "Split Heatmap By:", 
                          choices = c("No Split", "Fusion", "Age Category", "CEBPA", "FLT3-ITD", "NPM1"),
                          selected = "No Split")
        
      } else{
        updateSelectInput(session = session,
                          inputId = "annot", 
                          label = "Additional Annotation Options:", 
                          choices = c("Primary CNV", "Age Category", "Event Type ID", "CEBPA", "FLT3-ITD", "NPM1", "TP53", "WT1"))
        
        updateSelectInput(session = session,
                          inputId = "split", 
                          label = "Split Heatmap By:", 
                          choices = c("No Split", "Fusion", "Primary CNV", "Age Category", "Event Type ID", "CEBPA", "FLT3-ITD", "NPM1", "TP53", "WT1"),
                          selected = "No Split")
      }
      
      # Get unique fusion values, removing NA
      choices <- unique(na.omit(clinical$Primary.Fusion))
      # Remove "Other AML" and "None" to avoid duplication
      choices <- setdiff(choices, c("Other AML", "None"))
      # Sort alphabetically
      choices <- sort(choices)
      # Append "Other AML" and "None" at the end
      final_choices <- c(choices, "Other AML", "None")
      # Update reactive value with character vector (not factor)
      fusion_choices(final_choices)
      # Update the checkbox group in the UI
      updateCheckboxGroupInput(session, "fusions", "Fusions:", choices = final_choices, selected = final_choices)
      
    } else {
      
      fusion_choices <- "AML"
      updateCheckboxGroupInput(session, "fusions", "Generate Oncoprint:", choices = fusion_choices, selected = fusion_choices)
      
          if (dataset() == "BeatAML") {

            updateSelectInput(session = session,
                              inputId = "annot",
                              label = "Additional Annotation Options:",
                              choices = c("Age Category", "CEBPA", "FLT3-ITD", "NPM1", "TP53", "WT1"))

            updateSelectInput(session = session,
                              inputId = "split",
                              label = "Split Heatmap By:",
                              choices = c("No Split", "Age Category", "CEBPA", "FLT3-ITD", "NPM1", "TP53", "WT1"),
                              selected = "No Split")

          } else if (dataset() == "SWOG") {

            updateSelectInput(session = session,
                              inputId = "annot",
                              label = "Additional Annotation Options:",
                              choices = c("Age Category", "NPM1", "TP53", "WT1"))

            updateSelectInput(session = session,
                              inputId = "split",
                              label = "Split Heatmap By:",
                              choices = c("No Split", "Age Category", "NPM1", "TP53", "WT1"),
                              selected = "No Split")

          } else if (dataset() == "TCGA") {

            updateSelectInput(session = session,
                              inputId = "annot",
                              label = "Additional Annotation Options:",
                              choices = c("Age Category", "CEBPA", "FLT3-ITD", "NPM1", "TP53"))

            updateSelectInput(session = session,
                              inputId = "split",
                              label = "Split Heatmap By:",
                              choices = c("No Split", "Age Category", "CEBPA", "FLT3-ITD", "NPM1", "TP53"),
                              selected = "No Split")

          } else if (dataset() == "LEUCEGENE") {

            updateSelectInput(session = session,
                              inputId = "annot",
                              label = "There are no annotation options for LEUCEGENE",
                              choices = c(""))

            updateSelectInput(session = session,
                              inputId = "split",
                              label = "Split Heatmap By:",
                              choices = c("No Split"),
                              selected = "No Split")
          }
    }
  })

  # Select all / Deselect all logic
  observe({
    req(fusion_choices())  # ensure choices are available
    
    if (dataset() %in% c("PCGP AML", "TARGET")) {
      shinyjs::enable("selectall")
    } else {
      shinyjs::disable("selectall")
    }
    
    if (is.null(input$selectall)) return(NULL)
    
    if (input$selectall %% 2 == 0) {
      # Deselect all
      updateCheckboxGroupInput(session, "fusions", selected = character(0))
    } else {
      # Select all
      updateCheckboxGroupInput(session, "fusions", selected = fusion_choices())
    }
  })

  # Reading in the gene list from the csv input file
  readGeneList <- reactive({
    req(input$gene_input_method)
    
    if (input$gene_input_method == "file") {
      req(input$gene_list$datapath)
      
      gene_df <- tryCatch({
        read.csv(input$gene_list$datapath, stringsAsFactors = FALSE, header = FALSE)
      }, error = function(e) {
        showNotification("Error reading gene list file. Please upload a valid .csv file.", type = "error")
        return(NULL)
      })
      
      genes <- unique(trimws(gene_df[[1]]))
      genes <- genes[genes != ""]
      validate(need(length(genes) > 0, "Uploaded gene list is empty."))
      return(genes)
      
    } else {
      validate(need(input$gene_text != "", "Please enter gene names to get started."))
      
      genes <- unlist(strsplit(input$gene_text, "[,\n\r\t ]+"))
      genes <- unique(trimws(genes))
      genes <- genes[genes != ""]
      return(genes)
    }
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
    mat <- oncoMatrix()
    
    updatedclinical <- if (dataset() %in% c("TARGET", "PCGP AML")) {
      
      updatedclinical <- clinData()
      updatedclinical <- filter(updatedclinical, updatedclinical$Disease.Group == "AML")
      
      if (dataset() == "PCGP AML") {
        updatedclinical$Primary.Fusion <- ifelse(updatedclinical$Primary.Fusion == "Rare Fusion", "Other AML",
                                                 ifelse(updatedclinical$Primary.Fusion == "No Fusion", "None", updatedclinical$Primary.Fusion))
      }
      filter(updatedclinical, Primary.Fusion %in% input$fusions)
    } else {
      updatedclinical <- clinData()
      updatedclinical <- filter(updatedclinical, updatedclinical$Disease.Group == "AML")
    }
    
    # Subset mat to columns present in updatedclinical$PatientID
    mat <- mat[, colnames(mat) %in% updatedclinical$PatientID, drop = FALSE]
    
    # Filter out columns with any NAs
    mat <- mat[, colSums(is.na(mat)) == 0, drop = FALSE]
    
    # Order columns alphabetically by their names
    mat <- mat[, sort(colnames(mat)), drop = FALSE]
    
    if (ncol(mat) == 0 || nrow(mat) == 0) {
      validate(need(FALSE, "No samples meet the criteria"))
    }
    
    percentage_expressed <- function(x) sum(x >= input$tpms) / length(x) * 100
    percentages <- apply(mat, 1, percentage_expressed)
    mat <- mat[order(-percentages), , drop = FALSE]
    
    if (nrow(mat) > 30) mat <- mat[1:30, , drop = FALSE]
    
    gene_type <- ifelse(rownames(mat) %in% transmembrane_genelist$Gene.name, "Transmembrane", "Intracellular")
    mat <- ifelse(mat >= input$tpms, gene_type[row(mat)], "")

    return(as.matrix(mat))
  })
  
  # this function is for the heatmap annotation which will react to the checked input boxes
  updatedClinical <- reactive({
    mat <- fusionSubset()
    sample_ids <- colnames(mat)
    clinical <- clinData()
    
    # Filter only patients present in the matrix
    updatedclinical <- clinical[clinical$PatientID %in% sample_ids,]
    
    if (dataset() == "PCGP AML") {
      updatedclinical$Primary.Fusion <- ifelse(updatedclinical$Primary.Fusion == "Rare Fusion", "Other AML",
                                               ifelse(updatedclinical$Primary.Fusion == "No Fusion", "None", updatedclinical$Primary.Fusion))
    }
  
    # Further filter by selected fusions or cyto
    if (dataset() %in% c("TARGET", "PCGP AML")) {
      updatedclinical <- updatedclinical[updatedclinical$Primary.Fusion %in% input$fusions, ]
    }
    
    updatedclinical <- updatedclinical[order(updatedclinical$PatientID), ]
    
    if (dataset() == "TARGET" && aligner() == "kallisto") {
      
          updatedclinical$EFS.event.type.ID <- updatedclinical$`EFS event type ID`
          updatedclinical$TP53.Mutation <- updatedclinical$`TP53 mutation?`
          
          updatedclinical$Primary.CNV <- ifelse(is.na(updatedclinical$Primary.CNV), "Unknown", updatedclinical$Primary.CNV)
          updatedclinical$Age.Category <- ifelse(is.na(updatedclinical$Age.Category), "Unknown", updatedclinical$Age.Category)
          updatedclinical$EFS.event.type.ID <- ifelse(is.na(updatedclinical$EFS.event.type.ID), "Unknown", updatedclinical$EFS.event.type.ID)
          updatedclinical$FLT3.ITD <- ifelse(is.na(updatedclinical$FLT3.ITD), "No", updatedclinical$FLT3.ITD)
          updatedclinical$CEBPA.Mutation <- ifelse(is.na(updatedclinical$CEBPA.Mutation), "No", updatedclinical$CEBPA.Mutation)
          updatedclinical$NPM1.Mutation <- ifelse(is.na(updatedclinical$NPM1.Mutation), "No", updatedclinical$NPM1.Mutation)
          updatedclinical$TP53.Mutation <- ifelse(is.na(updatedclinical$TP53.Mutation), "No", updatedclinical$TP53.Mutation)

          # creating levels in each of the clinical data elements columns so that the information is ordered nicely
          # Get all unique levels except "Other AML" and "None"
          # Extract unique fusions
          existing_fusions <- unique(updatedclinical$Primary.Fusion)
          
          # Remove "Other AML" and "None" for sorting
          core_levels <- sort(setdiff(existing_fusions, c("Other AML", "None")))
          
          # Check if "Other AML" and "None" exist, and add in correct order
          optional_levels <- c()
          if ("Other AML" %in% existing_fusions) optional_levels <- c(optional_levels, "Other AML")
          if ("None" %in% existing_fusions) optional_levels <- c(optional_levels, "None")
          
          # Combine for final level order
          final_levels <- c(core_levels, optional_levels)
          
          # Apply factor level ordering
          updatedclinical$Primary.Fusion <- factor(updatedclinical$Primary.Fusion, levels = final_levels)
          
          updatedclinical$Primary.CNV <- factor(updatedclinical$Primary.CNV, levels = c("CBL deletion", "del5q", "monosomy7", "trisomy8", "trisomy21", "trisomy8/trisomy21", "-Y", "No Relevant CNV", "Unknown"))
          updatedclinical$Age.Category <- factor(updatedclinical$Age.Category, levels = c("Less than 3 years", "Between 3 and 5 years", "Between 5 and 10 years", "Between 10 and 18 years", 
                                                                                          "Greater than 18 years", "Unknown"))
          updatedclinical$EFS.event.type.ID <- updatedclinical$`EFS event type ID`
          updatedclinical$EFS.event.type.ID <- factor(updatedclinical$EFS.event.type.ID, levels = c("Censored", "Relapse", "Death", "Death without remission", "Induction failure", "Unknown"))
          updatedclinical$FLT3.ITD <- factor(updatedclinical$FLT3.ITD, levels = c("Yes", "No"))
          updatedclinical$WT1.Mutation <- factor(updatedclinical$WT1.Mutation, levels = c("Yes", "No"))
          updatedclinical$CEBPA.Mutation <- factor(updatedclinical$CEBPA.Mutation, levels = c("Yes", "No"))
          updatedclinical$NPM1.Mutation <- factor(updatedclinical$NPM1.Mutation, levels = c("Yes", "No"))
          updatedclinical$TP53.Mutation <- ifelse(!updatedclinical$`TP53 mutation` %in% c("", "No"), "Yes", "No")
          updatedclinical$TP53.Mutation <- factor(updatedclinical$TP53.Mutation, levels = c("Yes", "No"))

    } else if (dataset() == "TARGET" && aligner() == "star") {
      
      updatedclinical$Primary.CNV <- ifelse(is.na(updatedclinical$Primary.CNV), "Unknown", updatedclinical$Primary.CNV)
      updatedclinical$Age.Category <- ifelse(is.na(updatedclinical$Age.Category), "Unknown", updatedclinical$Age.Category)
      updatedclinical$EFS.event.type.ID <- ifelse(is.na(updatedclinical$EFS.event.type.ID), "Unknown", updatedclinical$EFS.event.type.ID)
      updatedclinical$FLT3.ITD <- ifelse(is.na(updatedclinical$FLT3.ITD), "No", updatedclinical$FLT3.ITD)
      updatedclinical$WT1.Mutation <- ifelse(is.na(updatedclinical$WT1.Mutation), "No", updatedclinical$WT1.Mutation)
      updatedclinical$CEBPA.Mutation <- ifelse(is.na(updatedclinical$CEBPA.Mutation), "No", updatedclinical$CEBPA.Mutation)
      updatedclinical$NPM1.Mutation <- ifelse(is.na(updatedclinical$NPM1.Mutation), "No", updatedclinical$NPM1.Mutation)
      updatedclinical$TP53.Mutation <- ifelse(is.na(updatedclinical$TP53.Mutation), "No", updatedclinical$TP53.Mutation)
      
      # Extract unique fusions
      existing_fusions <- unique(updatedclinical$Primary.Fusion)
      
      # Remove "Other AML" and "None" for sorting
      core_levels <- sort(setdiff(existing_fusions, c("Other AML", "None")))
      
      # Check if "Other AML" and "None" exist, and add in correct order
      optional_levels <- c()
      if ("Other AML" %in% existing_fusions) optional_levels <- c(optional_levels, "Other AML")
      if ("None" %in% existing_fusions) optional_levels <- c(optional_levels, "None")
      
      # Combine for final level order
      final_levels <- c(core_levels, optional_levels)
      
      # Apply factor level ordering
      updatedclinical$Primary.Fusion <- factor(updatedclinical$Primary.Fusion, levels = final_levels)
      
      
      updatedclinical$Primary.CNV <- factor(updatedclinical$Primary.CNV, levels = c("CBL deletion", "del5q", "monosomy7", "trisomy8", "trisomy21", "trisomy8/trisomy21", "minusY", "Unknown"))
      updatedclinical$Age.Category <- factor(updatedclinical$Age.Category, levels = c("Less than 3 years", "Between 3 and 5 years", "Between 5 and 10 years", "Between 10 and 18 years", "Greater than 18 years", "Unknown"))
      updatedclinical$EFS.event.type.ID <- factor(updatedclinical$EFS.event.type.ID, levels = c("Censored", "Relapse", "Death", "Death without remission", "Induction failure", "Unknown"))
      updatedclinical$FLT3.ITD <- factor(updatedclinical$FLT3.ITD, levels = c("Yes", "No"))
      updatedclinical$WT1.Mutation <- factor(updatedclinical$WT1.Mutation, levels = c("Yes", "No"))
      updatedclinical$CEBPA.Mutation <- factor(updatedclinical$CEBPA.Mutation, levels = c("Yes", "No"))
      updatedclinical$NPM1.Mutation <- factor(updatedclinical$NPM1.Mutation, levels = c("Yes", "No"))
      updatedclinical$TP53.Mutation <- factor(updatedclinical$TP53.Mutation, levels = c("Yes", "No"))

    } else if (dataset() == "BeatAML") {
      
      updatedclinical$ageAtDiagnosis <- as.numeric(updatedclinical$ageAtDiagnosis)
      updatedclinical$Age.Category <- ifelse(updatedclinical$ageAtDiagnosis < 20, "Less than 20 years",
             ifelse(updatedclinical$ageAtDiagnosis >= 20 & updatedclinical$ageAtDiagnosis <= 60, "Between 20 and 60 years",
                    ifelse(updatedclinical$ageAtDiagnosis > 60, "Greater than 60 years", NA)))
      
      updatedclinical$Age.Category <- ifelse(is.na(updatedclinical$Age.Category), "Unknown", updatedclinical$Age.Category)
      updatedclinical$FLT3.ITD <- ifelse(is.na(updatedclinical$FLT3.ITD), "No", updatedclinical$FLT3.ITD)
      updatedclinical$WT1.Mutation <- ifelse(is.na(updatedclinical$WT1.Mutation), "No", updatedclinical$WT1.Mutation)
      updatedclinical$CEBPA.Mutation <- ifelse(is.na(updatedclinical$CEBPA.Mutation), "No", updatedclinical$CEBPA.Mutation)
      updatedclinical$NPM1.Mutation <- ifelse(is.na(updatedclinical$NPM1.Mutation), "No", updatedclinical$NPM1.Mutation)
      updatedclinical$TP53.Mutation <- ifelse(is.na(updatedclinical$TP53.Mutation), "No", updatedclinical$TP53.Mutation)
      
      updatedclinical$Age.Category <- factor(updatedclinical$Age.Category, levels = c("Less than 20 years", "Between 20 and 60 years", "Greater than 60 years"))
      updatedclinical$CEBPA.Mutation <- ifelse(updatedclinical$CEBPA.Mutation != "negative", "Yes", "No")
      updatedclinical$CEBPA.Mutation <- factor(updatedclinical$CEBPA.Mutation, levels = c("Yes", "No"))
      updatedclinical$FLT3.ITD <- ifelse(updatedclinical$FLT3.ITD == "positive", "Yes", "No")
      updatedclinical$FLT3.ITD <- factor(updatedclinical$FLT3.ITD, levels = c("Yes", "No"))
      updatedclinical$NPM1.Mutation <- ifelse(updatedclinical$NPM1.Mutation == "positive", "Yes", "No")
      updatedclinical$NPM1.Mutation <- factor(updatedclinical$NPM1.Mutation, levels = c("Yes", "No"))
      updatedclinical$TP53.Mutation <- ifelse(updatedclinical$TP53.Mutation != "negative", "Yes", "No")
      updatedclinical$TP53.Mutation <- factor(updatedclinical$TP53.Mutation, levels = c("Yes", "No"))
      updatedclinical$WT1.Mutation <- ifelse(updatedclinical$WT1.Mutation != "negative", "Yes", "No")
      updatedclinical$WT1.Mutation <- factor(updatedclinical$WT1.Mutation, levels = c("Yes", "No"))
      
    } else if (dataset() == "SWOG") {
      
      updatedclinical$Age.Category <- ifelse(is.na(updatedclinical$Age.Category), "Unknown", updatedclinical$Age.Category)
      updatedclinical$WT1.Mutation <- ifelse(is.na(updatedclinical$WT1.Mutation), "No", updatedclinical$WT1.Mutation)
      updatedclinical$NPM1.Mutation <- ifelse(is.na(updatedclinical$NPM1.Mutation), "No", updatedclinical$NPM1.Mutation)
      updatedclinical$TP53.Mutation <- ifelse(is.na(updatedclinical$TP53.Mutation), "No", updatedclinical$TP53.Mutation)
      
      updatedclinical$Age.Category <- factor(updatedclinical$Age.Category, levels = c("Between 18 and 40 years", "Between 40 and 60 years", "Greater than 60 years", "Unknown"))
      updatedclinical$NPM1.Mutation <- ifelse(updatedclinical$NPM1.Mutation == "Positive", "Yes", "No")
      updatedclinical$NPM1.Mutation <- factor(updatedclinical$NPM1.Mutation, levels = c("Yes", "No"))
      updatedclinical$TP53.Mutation <- ifelse(updatedclinical$TP53.Mutation == "Positive", "Yes", "No")
      updatedclinical$TP53.Mutation <- factor(updatedclinical$TP53.Mutation, levels = c("Yes", "No"))
      updatedclinical$WT1.Mutation <- ifelse(updatedclinical$WT1.Mutation == "Positive", "Yes", "No")
      updatedclinical$WT1.Mutation <- factor(updatedclinical$WT1.Mutation, levels = c("Yes", "No"))
      
    } else if (dataset() == "TCGA") {
      
      updatedclinical$Age.Category <- ifelse(is.na(updatedclinical$Age.Category), "Unknown", updatedclinical$Age.Category)
      updatedclinical$CEBPA.Mutation <- ifelse(is.na(updatedclinical$CEBPA.Mutation), "No", updatedclinical$CEBPA.Mutation)
      updatedclinical$NPM1.Mutation <- ifelse(is.na(updatedclinical$NPM1.Mutation), "No", updatedclinical$NPM1.Mutation)
      updatedclinical$TP53.Mutation <- ifelse(is.na(updatedclinical$TP53.Mutation), "No", updatedclinical$TP53.Mutation)
      updatedclinical$Age.Category <- factor(updatedclinical$Age.Category, levels = c("Between 18 and 40 years", "Between 40 and 60 years", "Greater than 60 years"))
      updatedclinical$CEBPA.Mutation <- ifelse(updatedclinical$CEBPA.Mutation == "CEBPA mutation", "Yes", "No")
      updatedclinical$CEBPA.Mutation <- factor(updatedclinical$CEBPA.Mutation, levels = c("Yes", "No"))
      updatedclinical$NPM1.Mutation <- ifelse(updatedclinical$NPM1.Mutation == "NPM1 mutation", "Yes", "No")
      updatedclinical$NPM1.Mutation <- factor(updatedclinical$NPM1.Mutation, levels = c("Yes", "No"))
      updatedclinical$TP53.Mutation <- ifelse(updatedclinical$TP53.Mutation == "TP53 mutation", "Yes", "No")
      updatedclinical$TP53.Mutation <- factor(updatedclinical$TP53.Mutation, levels = c("Yes", "No"))
      updatedclinical$FLT3.ITD <- ifelse(grepl("in_frame_ins", updatedclinical$FLT3), "Yes", "No")
      updatedclinical$FLT3.ITD <- ifelse(is.na(updatedclinical$FLT3.ITD), "No", updatedclinical$FLT3.ITD)
      updatedclinical$FLT3.ITD <- factor(updatedclinical$FLT3.ITD, levels = c("Yes", "No"))
      
    } else if (dataset() == "PCGP AML") {
      
      updatedclinical$Age.Category <- ifelse(is.na(updatedclinical$Age.Category), "Unknown", updatedclinical$Age.Category)
      updatedclinical$CEBPA <- ifelse(is.na(updatedclinical$CEBPA), "No", updatedclinical$CEBPA)
      updatedclinical$NPM1 <- ifelse(is.na(updatedclinical$NPM1), "No", updatedclinical$NPM1)
      updatedclinical$FLT3.ITD <- ifelse(is.na(updatedclinical$FLT3.ITD), "No", updatedclinical$FLT3.ITD)
      
      # Extract unique fusions
      existing_fusions <- unique(updatedclinical$Primary.Fusion)
      
      # Remove "Other AML" and "None" for sorting
      core_levels <- sort(setdiff(existing_fusions, c("Other AML", "None")))
      
      # Check if "Other AML" and "None" exist, and add in correct order
      optional_levels <- c()
      if ("Other AML" %in% existing_fusions) optional_levels <- c(optional_levels, "Other AML")
      if ("None" %in% existing_fusions) optional_levels <- c(optional_levels, "None")
      
      # Combine for final level order
      final_levels <- c(core_levels, optional_levels)
      
      # Apply factor level ordering
      updatedclinical$Primary.Fusion <- factor(updatedclinical$Primary.Fusion, levels = final_levels)
      
      
      updatedclinical$Age.Category <- factor(updatedclinical$Age.Category, levels = c("Less than 3 years", "Between 3 and 5 years", "Between 5 and 10 years", "Between 10 and 18 years", "Greater than 18 years", "Unknown"))
      updatedclinical$CEBPA.Mutation <- updatedclinical$CEBPA
      updatedclinical$NPM1.Mutation <- updatedclinical$NPM1
      updatedclinical$CEBPA.Mutation <- factor(updatedclinical$CEBPA.Mutation, levels = c("Yes", "No"))
      updatedclinical$NPM1.Mutation <- factor(updatedclinical$NPM1.Mutation, levels = c("Yes", "No"))
      updatedclinical$FLT3.ITD <- factor(updatedclinical$FLT3.ITD, levels = c("Yes", "No"))
      
    }
    
    return(updatedclinical)
  })
    
  
  oncoAnnotation <- reactive({
    updatedclinical <- updatedClinical()
    updatedclinical <- filter(updatedclinical, updatedclinical$Disease.Group == "AML")

    # Color palettes
    colorlist <- c("#4E79A7", "#E15759", "#EDc948", "#B07AA1", "#F28E2B", "#59A14F", "#9C755F", 
                   "#a5bcd5", "#f0abac", '#f6e4a3', "#d7bcd0", "#f8c695", "#aad3a4", "#cebaae")
    cnvcolorlist <- c("#4E79A7", "#E15759", "#EDc948", "#B07AA1", "#F28E2B", "#59A14F", "#9C755F", "#a5bcd5", "#ffffff")
    efscolorlist <- c("#4E79A7", "#E15759", "#EDc948", "#B07AA1", "#F28E2B", "#ffffff")
    gradcolorlist <- c("#c4d0e1", "#89a1c4", "#4E79A7", "#344967", "#161f2c", "#ffffff")
    yesnocolorlist <- c("#161f2c", "#c4d0e1")
    
    # Initialize annotation data
    annotation_data <- list(column_bar = anno_oncoprint_barplot())
    color_list <- list()
    
    # Helper to add annotation safely if the column exists
    add_annotation <- function(name, column, palette) {
      if (column %in% colnames(updatedclinical)) {
        annotation_data[[name]] <<- updatedclinical[[column]]
        levels_ <- levels(as.factor(updatedclinical[[column]]))
        color_list[[name]] <<- setNames(palette[seq_along(levels_)], levels_)
      }
    }
    
    # Always include Fusion if TARGET/PCGP
    if (dataset() %in% c("TARGET", "PCGP AML")) {
      add_annotation("Fusion", "Primary.Fusion", colorlist)
    } else {
      updatedclinical$AML <- updatedclinical$Disease.Group
      
      add_annotation("AML", "AML", colorlist)
    }
    
    # Conditional annotations based on input$annot and existence in data
    if (any(grepl("CNV", input$annot))) add_annotation("CNV", "Primary.CNV", cnvcolorlist)
    if (any(grepl("Age", input$annot))) add_annotation("Age", "Age.Category", gradcolorlist)
    if (any(grepl("Event", input$annot))) add_annotation("EFS", "EFS.event.type.ID", efscolorlist)
    if (any(grepl("FLT3", input$annot))) add_annotation("FLT3", "FLT3.ITD", yesnocolorlist)
    if (any(grepl("WT1", input$annot))) add_annotation("WT1", "WT1.Mutation", yesnocolorlist)
    if (any(grepl("CEBPA", input$annot))) add_annotation("CEBPA", "CEBPA.Mutation", yesnocolorlist)
    if (any(grepl("NPM1", input$annot))) add_annotation("NPM1", "NPM1.Mutation", yesnocolorlist)
    if (any(grepl("TP53", input$annot))) add_annotation("TP53", "TP53.Mutation", yesnocolorlist)
    
    # Annotation heights
    if (length(annotation_data) == 1) {
      annotation_heights <- unit(1.5, "cm")
    } else {
      annotation_heights <- c(unit(1.5, "cm"), rep(unit(0.8, "cm"), length(annotation_data) - 1))
    }
    
    
    # Construct HeatmapAnnotation safely
    ha <- do.call(HeatmapAnnotation, c(
      annotation_data,
      list(
        col = color_list,
        annotation_height = annotation_heights,
        annotation_name_gp = gpar(fontsize = 10),
        gp = gpar(fontsize = 8)
      )
    ))
    
    return(ha)
  })
  
  splitList <- reactive({
    uc <- updatedClinical()
    split_input <- input$split
    
    # Mapping from user-friendly labels to column names
    column_mapping <- c(
      "Fusion" = "Primary.Fusion",
      "Primary CNV" = "Primary.CNV",
      "Age Category" = "Age.Category",
      "Event Type ID" = "EFS.event.type.ID",
      "FLT3-ITD" = "FLT3.ITD",
      "WT1" = "WT1.Mutation",
      "CEBPA" = "CEBPA.Mutation",
      "NPM1" = "NPM1.Mutation",
      "TP53" = "TP53.Mutation"
    )
    
    # Handle "No Split" or empty input
    if (is.null(split_input) || split_input == "No Split" || !split_input %in% names(column_mapping)) {
      return(NULL)
    }
    
    split_col <- column_mapping[[split_input]]
    print(split_col)
    
    # Validate column exists
    if (!split_col %in% colnames(uc)) {
      print("DOES NOT EXIST")
      return(NULL)
    }
    
    # Return column as factor for splitting
    return(as.factor(uc[[split_col]]))
  })
  
  
  
  # plotting the oncoprint (finally!)
  plotOncoprint <- reactive({
    
    validate(
      need(!dataset() %in% c("PCGP", "StJude", "GMKF", "CCLE"),
           "Oncoprint function is not currently available for this dataset")
    )
    
    mat <- fusionSubset() # reading in the matrix
    updatedclinical <- updatedClinical() # reading in the subsetted clinical data
    ha <- oncoAnnotation() # reading in the heatmap annotation
    
    validate(
      need(!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0, "No genes or samples to plot."),
      need(!all(is.na(mat) | mat == ""), "No samples meet the required TPM cutoff for the selected genes")
    )
    
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
            column_split = splitList(),
            remove_empty_columns = FALSE, # can turn to TRUE if you want to get rid of empty columns (harder to visualize the % though)
            heatmap_legend_param = list(title = "Domain"))
    
    d1 <- draw(d1, merge_legend = TRUE) # this keeps the legends clean
    
  })
  
  #################################################################
  #-------------------- FINAL MODULE OUTPUTS ---------------------#
  #################################################################
  
  output$plot <- renderPlot({
    # Only validate fusions input if dataset is PCGP AML or TARGET
    if (dataset() %in% c("PCGP AML", "TARGET")) {
      validate(
        need(length(input$fusions) > 0, "Please select at least one fusion")
      )
    }
    
    # Plot regardless for other datasets
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
