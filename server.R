server <- function(input, output, session) { 

  # ## build the temporary directory at the start
  # current_dir <- getwd()
  # temp_dir <- tempdir()
  # unlink(file.path(temp_dir, "*"), recursive = TRUE)
  # 
  # ## duplicate gene names?
  # duplicates <- FALSE
  # gene_choices <- c(NULL)
  
  # the following functions are to create checkmarks for whether the input gene is aml-restricted and transmembrane
  #-------------------------------------------------------------#
  output$gene_present <- reactive({
    tolower(input$geneInput) %in% tolower(aml_restricted_genelist$Gene)
  })
  
  outputOptions(output, "gene_present", suspendWhenHidden = FALSE)
  
  output$trmembrane <- reactive({
    tolower(input$geneInput) %in% tolower(transmembrane_genelist$Gene.name)
  })
  
  outputOptions(output, "trmembrane", suspendWhenHidden = FALSE)
  #-------------------------------------------------------------#
  
  observeEvent(input$leukemiaSelection, {
    
    choices <- switch(input$leukemiaSelection,
                      "AML" = dataset_choices$aml,
                      "BALL" = dataset_choices$all, 
                      "TALL" = dataset_choices$tall,
                      "CCLE" = dataset_choices$ccle,
                      "PCGP" = dataset_choices$pcgp) 
    
    selected <- switch(input$leukemiaSelection,
                       "AML" = "TARGET",
                       "BALL" = "StJude",
                       "TALL" = "GMKF",
                       "CCLE" = "CCLE",
                       "PCGP" = "PCGP")
    
    updateRadioButtons(
      session = session,
      inputId = "expDataCohort", 
      choices = choices,
      selected = selected
    )
  }, ignoreInit = T, ignoreNULL = T) # ignoreInit parameter is required to work!
  
  # Variable representing the *name* of the selected cohort (as a character string)
  cohort <- reactive({
      input$expDataCohort
  })
  
  expData <- reactive({
    # Select base matrix based on cohort
    matrix <- switch(input$expDataCohort,
                     "SWOG" = swog_expData,
                     "BeatAML" = beatAML_expData,
                     "TARGET" = NULL,  # Placeholder, will be assigned below
                     "TCGA" = laml_expData,
                     "StJude" = stjude_expData,
                     "GMKF" = gmkf_expData,
                     "CCLE" = ccle_expData,
                     "LEUCEGENE" = leuce_expData,
                     "PCGP AML" = pcgp_aml_expData,
                     "PCGP" = pcgp_expData)
    
    # Handle TARGET cohort separately
    if (input$expDataCohort == "TARGET") {
      if (input$aligner == "star") {
        # Default to diagnostic timepoint
        matrix <- switch(input$timepoint,
                         "diagnostic" = target_expData38_STAR_Dx, 
                         "remission" = target_expData38_STAR_Rem,
                         "relapse" = target_expData38_STAR_Rel)
        
        genenames <- matrix$name
        matrix <- matrix[,c(-1:-2)]
        rownames(matrix) <- genenames
        
      } else if (input$aligner == "kallisto") {
        # No timepoint selection for kallisto
        matrix <- target_expData38
      }
    }
    return(matrix)
  })
  
  # Reactive variable that stores the clinical data elements for the selected cohort
  studyData <- reactive({
    if (input$expDataCohort == "TARGET") {
      switch(input$aligner,
             "star" = target_mani,    # Use target_mani for STAR aligner
             "kallisto" = target_cde) # Use target_cde for Kallisto aligner
    } else {
      switch(input$expDataCohort,
             "SWOG" = swog_cde,
             "BeatAML" = beatAML_cde,
             "TCGA" = laml_cde,
             "StJude" = stjude_cde,
             "GMKF" = gmkf_cde,
             "CCLE" = ccle_cde,
             "LEUCEGENE" = leuce_mani,
             "PCGP AML" = pcgp_mani,
             "PCGP" = pcgp_total_mani)
    }
  })
  
  # Creating a variable that will be used to reactively pass the gene of interest into each module,
  # See https://tbradley1013.github.io/2018/07/20/r-shiny-modules--using-global-inputs/ for more  
  # info on passing global Shiny variables into a module
  target <- reactive({
    if (grepl("^hsa\\-mir*|mir\\-*", input$geneInput, ignore.case = T)) {
      symbol <- gsub("[Hh][Ss][Aa]-[Mm][Ii][Rr]", "hsa-miR", input$geneInput) # Casting the R to uppercase since this is all mature miR data
    } else if (grepl("^MIMAT", input$geneInput, ignore.case = T)) { # Mapping MIMAT ID back to hsa ID, the user can enter either one
      symbol <- miRmapping$hsa.ID.miRbase21[match(toupper(input$geneInput), miRmapping$MIMAT.ID)]
    } else if (grepl("orf", input$geneInput, ignore.case = T)) {
      symbol <- gsub("[Oo][Rr][F]", "orf", input$geneInput)
    } else {
      symbol <- toupper(input$geneInput)
    }
    return(symbol)
  })
  
  alignment <- reactiveVal("star")  # Default aligner is "star"
  
  # Track last known aligner for TARGET only
  targetAligner <- reactiveVal("star")
  
  # Update aligner value when aligner changes and cohort is TARGET
  observeEvent(input$aligner, {
    if (input$expDataCohort == "TARGET") {
      alignment(input$aligner)
      targetAligner(input$aligner)
    }
  })
  
  # Handle changes to cohort
  observeEvent(input$expDataCohort, {
    if (input$expDataCohort == "TARGET") {
      # Restore last used TARGET aligner or default to "star"
      updateRadioGroupButtons(session, "aligner", selected = targetAligner())
      alignment(targetAligner())
      shinyjs::enable("aligner")
    } else {
      # Force aligner to "star" internally and disable UI
      updateRadioGroupButtons(session, "aligner", selected = "star")
      alignment("star")
      shinyjs::disable("aligner")
    }
  })
  
  
  observeEvent(input$ensid, {
    req(input$expDataCohort, input$aligner, input$geneInput)
    
    if (input$expDataCohort == "TARGET" && input$aligner == "star") {
      
      matrix <- switch(input$timepoint,
                       "diagnostic" = target_expData38_STAR_Dx, 
                       "remission" = target_expData38_STAR_Rem,
                       "relapse" = target_expData38_STAR_Rel)
      
      filteredData <- matrix[matrix$name == toupper(input$geneInput), ]
    
      if (nrow(filteredData) > 0) {
        ensid_info <- paste("<strong>Matching ENSID(s):</strong><br>", 
                            paste(unique(filteredData$id), collapse = "<br>"))
      } else {
        ensid_info <- "<strong style='color: red;'>No ENSID found for this gene.</strong>"
      }
    } else {
      ensid_info <- "<strong style='color: red;'>ENSID information not available for this dataset.</strong>"
    }
    
    showModal(modalDialog(
      title = div(style = "font-size: 20px; font-weight: bold; color: #41B0FA;", "ENSID Information"),
      HTML(paste("<div style='padding: 10px; font-size: 16px;'>", ensid_info, "</div>")),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  observeEvent(input$check, {
    option <- grep(input$geneInput, rownames(expData()), value = TRUE, ignore.case = TRUE)
    
    msg <- if (length(option) == 0) {
      "No alternate names found!"
    } else if (length(option) > 20) {
      paste0(paste0(option[1:20], collapse = "  or\n"), ",\n...?")
    } else {
      paste0(paste0(option, collapse = "  or\n"), "?")
    }
    
    shinyalert(
      title = "Did you mean...",
      text = msg,
      closeOnEsc = TRUE,
      closeOnClickOutside = TRUE,
      html = FALSE,
      type = "input",
      inputType = "text",
      inputValue = if (length(option) > 0) option[1] else "",
      inputPlaceholder = "Type choice here",
      showConfirmButton = TRUE,
      showCancelButton = FALSE,
      confirmButtonText = "Retry",
      confirmButtonCol = "#41B0FA",
      timer = 0,
      imageUrl = "",
      animation = TRUE,
      callbackR = function(x){
        if (!is.null(x) && is.character(x) && nzchar(x) && x != "false") {
          updateTextInput(session, "geneInput", label = NULL, value = x)
        }
      }
    )
  })
  
  
  #--------------------- WF plot & KM plot tabs --------------------- #
  
  # Calling the waterfall plot module
  # IMPORTANT NOTE: the "target" & "cohort" variables are actually a reactive function, and would usually be called by target() & cohort(), 
  # but when passing a reactive value into a module, you *must* pull off the parentheses and pass the naked variable name as an argument.
  # Within the modules themselves, these variables are a reactive function!
  callModule(wfPlot, id = "waterfall", 
             clinData = studyData, 
             expData = expData, 
             adc_cart_targetData = adc_cart_targetData,
             gene = target, 
             aligner = alignment,
             dataset = cohort,
             parent = session) # See https://stackoverflow.com/questions/51708815/accessing-parent-namespace-inside-a-shiny-module
                               # for an explanation of the 'parent' parameter
  
  callModule(timePlot, id = "timepoints", 
             clinData = studyData, 
             expData = expData, 
             gene = target, 
             aligner = alignment,
             dataset = cohort)
  
  # Calling the Kaplan-Meier curve module
  callModule(kmPlot, id = "kaplanmeier", 
             clinData = studyData, 
             expData = expData, 
             aligner = alignment,
             dataset = cohort,
             gene = target)
  
  # This module is not ready for prime time yet
  # callModule(heatmap, id = "heatmap",
  #            clinData = studyData,
  #            expData = expData,
  #            dataset = cohort,
  #            gene = target)
  
  callModule(oncoprint, id = "oncoprint",
             clinData = studyData,
             expData = expData,
             dataset = cohort,
             aligner = alignment)
  
  # Calling the DEG table module
  # callModule(deTable, id = "degs",
  #            table = degTables37,
  #            gene = target)
  
  # Calling the DEG table module
  callModule(geneExp, id = "exps",
             clinData = studyData, 
             expData = expData, 
             gene = target,
             dataset = cohort)
  
  callModule(HPAPlot, id = "hpa",
             gene = target)

  # # Calling the HPA module
  # callModule(ClassiPlot, id = "Classi")

  callModule(CancerPlot, id = "cancertype",
             gene = target)
  
  #--------------------- External databases tab --------------------- #

  # TO DO: Add a searchable AML-restricted gene list to this tab
  
  customBox <- function(title, subtitle, link, icon_name) {
    tags$a(
      href = link, target = "_blank", style = "text-decoration: none;",
      tags$div(
        class = "small-box",
        style = paste(
          "background: #2096f6;",
          "color: #FFFFFF;",                    # dark blue text
          "min-height: 120px;",
          "border-radius: 12px;",
          "box-shadow: 0 4px 8px rgba(0,0,0,0.1);",
          "padding: 15px;"
        ),
        tags$div(
          class = "inner",
          tags$p(title, style = "font-size: 200%; font-weight: bold; margin-top: -10px;"),
          tags$h4(subtitle, style = "font-size: 120%; font-weight: normal; color: #ffffff; margin-top: -10px;")
        ),
        tags$div(class = "icon-large", icon(icon_name), style = "margin-right: 5px;")
      )
    )
  }
  
  output$protAtlas <- renderUI({
    customBox("Human Protein Atlas", "Protein Expression",
              paste0("https://www.proteinatlas.org/search/", target()), "dna")
  })
  
  output$gtex <- renderUI({
    customBox("GTEx", "Normal Tissue Expression",
              paste0("https://gtexportal.org/home/gene/", target(), "#geneExpression"), "lungs")
  })
  
  output$protPaint <- renderUI({
    customBox("ProteinPaint", "St. Jude PeCan",
              paste0("https://proteinpaint.stjude.org/?genome=hg38&gene=", target(), "&dataset=pediatric"), "brush")
  })
  
  output$cbioportal <- renderUI({
    customBox("cBioPortal", "Cancer Datasets",
              "https://www.cbioportal.org/", "database")
  })
  
  output$deeptmhmm <- renderUI({
    customBox("DeepTMHMM", "Transmembrane Prediction",
              "https://dtu.biolib.com/DeepTMHMM", "map-location-dot")
  })
  
  output$expasy <- renderUI({
    customBox("Expasy", "Nucleotide-Protein Translation",
              "https://web.expasy.org/translate/", "language")
  })
  
  output$ucsc <- renderUI({
    customBox("UCSC", "Genome Browser",
              "https://genome.ucsc.edu/index.html", "binoculars")
  })
  
  output$vizrisk <- renderUI({
    customBox("VizRisk", "Risk Classification App",
              "https://vizrisk.fredhutch.org/", "notes-medical")
  })
  
  
  
  
 ## This is the code for the DeepTMHMM Button ###################################################
  ## Here's the order of operations:
  ## 1. Press Button
  ## 2. Clear the temporary directory and the output plot.png
  ## 3. Set a "waiting for..." message
  ## 4. Translate the gene name into a fasta sequence
  ## 5. Write the fasta sequence to a fasta file
  ## 6. Run the biolib DeepTMHMM function
  ## 7. Print the cleaned terminal output 
  ## 8. Output the plot.png
  
  # # function for creating the action button for embedding DeepTMHMM ###############################
  # output$tmhmm <- renderUI({
  #   validate(
  #     need(target(), FALSE)
  #   )
  #   
  #   actionButton("start_deeptmhmm", 
  #                label = div(
  #                  "DeepTMHMM",
  #                  div("Protein localization", style = "text-transform: none; color: white; font-size: 15px; font-weight: normal; margin-top:20px; margin-bottom:20px;") # Additional white text below the label
  #                ),
  #                style = "text-transform: none; background-color: #3c8dbc; box-shadow: none; text-align: left; font-size: 21px; font-weight: bold; height: 110px; width: 100%; padding: 10px;",
  #                class = "btn-box"
  #   )
  # })
  # 
  # # creating a reactive value that will change once the output is finished
  # output_completed <- reactiveVal(FALSE)
  # 
  # poll_terminal_output <- function(myTerm) {
  #   function() {
  #     setwd(temp_dir)
  #     output <- NULL
  #     
  #     # Check if running inside RStudio
  #     if (Sys.getenv("RSTUDIO") == "1") {  
  #       full_output <- rstudioapi::terminalBuffer(myTerm)
  #       
  #       # Check if "Done" is in the terminal output
  #       if (any(grepl("Done", full_output))) {
  #         Sys.sleep(10)  # Wait for 10 seconds
  #         output_completed(TRUE)  # Mark the output as completed
  #       } else {
  #         output_completed(FALSE)
  #       }
  #     } 
  #     else {
  #       # Check if output_log.txt exists outside RStudio
  #       if (file.exists("output_log.txt")) {
  #         full_output <- readLines("output_log.txt")
  #         # Check if "Done" is in the file output
  #         if (any(grepl("Done", full_output))) {
  #           Sys.sleep(10)  # Wait for 10 seconds
  #           output_completed(TRUE)  # Mark the output as completed
  #         } else {
  #           output_completed(FALSE)
  #         }
  #       } else {
  #         # If output_log.txt doesn't exist
  #         output <- output[output != "" & !is.null(output)]
  #         output_completed(FALSE)
  #       }
  #     }
  #     return(output)  # Return the output if any
  #   }
  # }
  # 
  # 
  # # when the action button is pushed, we're starting the DeepTMHMM code ###########################
  # observeEvent(input$start_deeptmhmm, {
  #   
  #   output_completed(FALSE)
  #   
  #   output$tmhmm_plot <- renderImage({
  #     filename <- normalizePath(file.path(current_dir, 'data', 'loading.png'))
  #     
  #     list(
  #       src = filename,
  #       alt = "Loading...",
  #       width = "80%",
  #       height = "auto"
  #     )
  #   }, deleteFile = FALSE)
  #   
  #   # cleaning the temporary directory
  #   unlink(file.path(temp_dir, "*"), recursive = TRUE)
  #   
  #   # setting the temporary directory
  #   setwd(temp_dir)
  #   
  #   # this gives a little loading text before the output starts
  #   output$terminal_output <- renderText({
  #       "Waiting for DeepTMHMM to start, this may take a few seconds... \n"
  #   })
  #   
  #   # the goal of this is to take the gene name and find the longest peptide sequence
  #   ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #   gene_name <- target()
  #   sequences <- getSequence(id = gene_name, type = "hgnc_symbol", seqType = "peptide", mart = ensembl)
  #   
  #   # longest sequence function
  #   get_longest_sequence <- function(sequences) {
  #     sequence_lengths <- nchar(sequences$peptide)
  #     max_length <- max(sequence_lengths)
  #     longest_indices <- which(sequence_lengths == max_length)
  #     longest_sequence <- sequences$peptide[longest_indices[1]]
  #     return(longest_sequence)
  #   }
  #   
  #   # finding the longest sequence
  #   longest_sequence <- get_longest_sequence(sequences)
  #   
  #   # create a temporary fasta file
  #   temp_fasta <- tempfile(fileext = ".fasta")
  #   writeLines(paste0(">header\n", longest_sequence), con = temp_fasta)
  #   
  #   # now we run the python biolib package inside a terminal
  #   myTerm <- NULL
  #   if (Sys.getenv("RSTUDIO") == "1") { 
  #     myTerm <- rstudioapi::terminalCreate(show = FALSE)
  #     tmhmm <- paste("biolib run DTU/DeepTMHMM --fasta", temp_fasta)
  #     rstudioapi::terminalSend(myTerm, paste0(tmhmm, "\n"))
  #   } else {
  #     system(paste("biolib run DTU/DeepTMHMM --fasta", temp_fasta, "> output_log.txt 2>&1 &"))
  #   }
  #   
  #   # reactive function to monitor terminal output
  #   terminal_output <- reactivePoll(1000, session, checkFunc = poll_terminal_output(myTerm), valueFunc = poll_terminal_output(myTerm))
  #   
  #   observe({
  #     # Check the value of output_completed reactively
  #     if (output_completed() == TRUE) {
  #       
  #       # Set the working directory to temp_dir
  #       setwd(temp_dir)
  #       
  #       # Render the image when it's available
  #       output$tmhmm_plot <- renderImage({
  #       
  #         filename <- normalizePath(file.path(temp_dir, 'biolib_results', 'plot.png'))
  # 
  #         list(
  #           src = filename,
  #           alt = "TMHMM Plot",
  #           width = "80%",
  #           height = "auto"
  #         )
  #       }, deleteFile = TRUE)
  #       
  #     } else {
  #       # Render the placeholder/loading image
  #       output$tmhmm_plot <- renderImage({
  #         filename <- normalizePath(file.path(current_dir, 'data', 'loading.png'))
  #         
  #         list(
  #           src = filename,
  #           alt = "Loading...",
  #           width = "80%",
  #           height = "auto"
  #         )
  #       }, deleteFile = FALSE)
  #     }
  #   })
  #   
  # })
  
  
  output$therapyTable <- DT::renderDataTable({
    validate(
      need(target(), "Please enter a gene symbol in the text box.") %then%
        need(target() %in% adc_cart_targetData$`Gene target`, paste0("We do not have record of ", target(), " being targeted\n by ADC or CAR T-cell therapies.")))
    
    table <- adc_cart_targetData %>%
      filter(`Gene target` == target()) 
    
    DT::datatable(table, 
                  options = list(scrollY = "50vh",
                                 pageLength = 25,
                                 searchHighlight = TRUE), 
                  escape = F)
  })

  
  #--------------------- Protein Paint tab --------------------- #
  
  # UPDATE: Moved this entire section to the "External Databases" tab, this was extremely finicky and I had trouble getting it to work.
  
  # Can't figure out a way to reactively update the HTML embedding to
  # reflect the user-specified gene (in the Shiny app text box).
  # I would need to update the "positionbygene" parameter in the actual HTML file to do this,
  # url <- reactive({
  #   # https://stackoverflow.com/questions/28982722/shiny-iframe-reactive
  #   user_entry <- isolate(target()) # https://stackoverflow.com/questions/45886021/how-can-you-pass-a-url-to-an-iframe-via-textinput-in-r-shiny
  #   final <- isolate(paste0("https://proteinpaint.stjude.org/?genome=hg19&gene=", target(), "&dataset=pediatric"))
  #   return(final)
  # })
  
  # output$htmlDisplay <- renderUI({
  #   # validate( # This validate statement will only be needed when the embedded HTML is synced w/ the user-specified gene entry.
  #     # need(target(), "Please enter a gene symbol in the text box."))
  #   
  #   tags$iframe(style = "border-width: 0;",
  #               width = 1300,
  #               height = 800,
  #               src = url())
  #               # src = "https://proteinpaint.stjude.org/?genome=hg19&gene=MSLN&dataset=pediatric")
  #               # src = "Protein_Paint/embed_StJude_ProteinPaint.html") # src param must be a filename in the www folder,
  # })                                                                  # don't include the "www/" prefix or it won't work!
  
  # Trying to write out a modified version of the ProteinPaint HTML file whenever the target() reactive variable is changed.
  # This modified HTML file could then be supplied to the 'src' parameter in the iframe.
  # This would automatically update the HTML embedding to reflect the user-supplied gene.
  # I'm having issues getting this to work, though, so it will be commented out for now. 
  # Hopefully will have more time to work this out in the future - I think I may need to add an "Update"
  # action button
  # writeFile <- reactive({
  #   protPaint_html_mod <- gsub("(?<=positionbygene\\:\\').+(?=\\')", target(), protPaint_html, perl = T)
  #   write_file(protPaint_html_mod, "www/Protein_Paint/embed_StJude_ProteinPaint_writeTest.html")   
  # })
   # writeFile()
  
  #--------------------- UMAP tab --------------------- #
  
  # Following this post, but it doesn't work: https://stackoverflow.com/questions/24875943/display-html-file-in-shiny-app
  # This person is having the same issue I am:
  # https://stackoverflow.com/questions/56064805/displaying-html-file-using-includehtml-in-shiny-is-not-working-with-renderui
  
  output$umapEmbedding <- renderUI({

    ########### Method 1 ##############
    # includeHTML() is designed to work with HTML fragments, so a "self contained" HTML file is needed,
    # aka only the <body> section with an <html> </html> tag layer outside of it
    # includeHTML("www/UMAP/TARGET_AML_sg7655_blackBackground_clusters2_k31_PCAselect.html")

    # Passing r variables to html IF the html is included in a markdown:
    # https://stackoverflow.com/questions/61543937/pass-r-variable-to-html-in-r-markdown
    # {{ uiOutput("score_value") }} <- I think this is the syntax I need to embed in the HTML file?

    ########### Method 2 #############
    # see iframe details at https://plotly-r.com/saving.html,
    # using same parameters as the unused code above for Protein Paint.
    # To do: get a folder of UMAPs from Jenny, and allow the user to select the base plot they want to manipulate.
    # NOTE: Don't include 'www/' in filepath, see
    # https://stackoverflow.com/questions/41784631/include-link-to-local-html-file-in-datatable-in-shiny
    # for an explanation.
    tags$iframe(seamless = "seamless",
                style = "border-width: 0;",
                src = "UMAP/TARGET_AML_sg7655_blackBackground_clusters2_k31_PCAselect.html",
                height = 700, width = 1300, scrolling = "yes")
  })
  
}
