server <- function(input, output, session) { 
  
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
                     "ALL" = dataset_choices$all) 
    
    selected <- switch(input$leukemiaSelection,
                      "AML" = "TARGET",
                      "ALL" = "StJude")
    
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
  
  # Reactive variable that stores the expression data matrix for the selected cohort
  expData <- reactive({
    matrix <- switch(input$expDataCohort,
                     "SWOG" = swog_expData,
                     "BeatAML" = beatAML_expData,
                     "TARGET" = target_expData38,
                     "TCGA" = laml_expData,
                     "StJude" = stjude_expData)
    
    # For the TARGET dataset only, we have both GRCh37 & GRCh38-aligned datasets available. 
    # This will allow the user to select one of those alignments, but ONLY if the TARGET AML dataset has been selected.
    if (input$expDataCohort == "TARGET" ) {
      assembly <- switch(input$seqAssembly,
                         "grch37" = target_expData37,
                         "grch38" = target_expData38)
      return(assembly)
    } else {
      return(matrix)
    }
  })
  
  # Reactive variable that stores the clinical data elements for the selected cohort
  studyData <- reactive({
    switch(input$expDataCohort,
           "SWOG" = swog_cde,
           "BeatAML" = beatAML_cde,
           "TARGET" = target_cde,
           "TCGA" = laml_cde,
           "StJude" = stjude_cde)
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
  
  observeEvent(input$check, {
    option <- grep(input$geneInput, rownames(expData()), value = T, ignore.case = T)
    msg <- if (length(option) == 0) {
      "No alternate names found!"
    } else if (length(option) > 20) {
      paste0(paste0(option[1:20], collapse = "  or\n"), ",\n...?")
    } else {
      paste0(paste0(option, collapse = "  or\n"), "?")
    }
    
    # Check if miRNA name exists in the miRbase21 miRNA-seq data
    shinyalert(
      title = "Did you mean...",
      text = msg,
      closeOnEsc = TRUE,
      closeOnClickOutside = TRUE,
      html = FALSE,
      type = "input",
      inputType = "text",
      inputValue = option[1],
      inputPlaceholder = "Type choice here",
      showConfirmButton = TRUE,
      showCancelButton = FALSE,
      confirmButtonText = "Retry",
      confirmButtonCol = "#41B0FA",
      timer = 0,
      imageUrl = "",
      animation = TRUE,
      callbackR = function(x){
        if (!is.null(x)) updateTextInput(session, "geneInput", label = NULL, value = input$shinyalert)
      })
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
             dataset = cohort,
             parent = session) # See https://stackoverflow.com/questions/51708815/accessing-parent-namespace-inside-a-shiny-module
                               # for an explanation of the 'parent' parameter
  
  # Calling the Kaplan-Meier curve module
  callModule(kmPlot, id = "kaplanmeier", 
             clinData = studyData, 
             expData = expData, 
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
             gene = target)
  
  # Calling the DEG table module
  callModule(deTable, id = "degs",
             table = degTables37,
             gene = target)
  
  # Calling the DEG table module
  callModule(geneExp, id = "exps",
             clinData = studyData, 
             expData = expData, 
             gene = target,
             dataset = cohort)
  
  callModule(HPAPlot, id = "hpa",
             gene = target)

  # Calling the HPA module
  callModule(ClassiPlot, id = "Classi")

  callModule(CancerPlot, id = "cancertype",
             gene = target)
  
  #--------------------- External databases tab --------------------- #

  # TO DO: Add a searchable AML-restricted gene list to this tab
  
  output$protAtlas <- renderValueBox({
    validate(
      need(target(), "Please enter a gene symbol in the text box."))
    
        valueBox(value = tags$p("Human Protein\nAtlas", style = "font-size: 60%"),
                 subtitle = "Protein expression", 
                 color = "light-blue", 
                 icon = icon("prescription-bottle"),
                 href = paste0("https://www.proteinatlas.org/search/", target()))
  })
  
  output$gtex <- renderValueBox({
    validate(
      need(target(), "Please enter a gene symbol in the text box."))
    
    valueBox(value = tags$p("GTEx", style = "font-size: 60%"),
             subtitle = "Normal tissue expression",
             color = "light-blue",
             icon = icon("prescription-bottle"),
             href = paste0("https://gtexportal.org/home/gene/", target(), "#geneExpression"))
  })
  
  output$protPaint <- renderValueBox({
    validate(
      need(target(), "Please enter a gene symbol in the text box."))
    
    valueBox(value = tags$p("ProteinPaint", style = "font-size: 60%"),
             subtitle = "St. Jude PeCan visualization",
             color = "light-blue",
             icon = icon("prescription-bottle"), 
             href = paste0("https://proteinpaint.stjude.org/?genome=hg19&gene=", target(), "&dataset=pediatric"))
  })

output$tmhmm <- renderUI({
    validate(
      need(target(), "Please enter a gene symbol in the text box.")
    )
    
    actionButton("start_rselenium", 
                 label = div(
                   "DeepTMHMM",
                   div("Protein Localization, Note: Needs Firefox to Run", style = "text-transform: none; color: white; font-size: 15px; font-weight: normal; margin-top:20px; margin-bottom:20px;") # Additional white text below the label
                 ),
                 style = "text-transform: none; background-color: #3c8dbc; box-shadow: none; text-align: left; margin-left: 15px; margin-right: 20px; font-size: 21px; font-weight: bold; height: 110px; width: 33%; padding: 10px;",
                 class = "btn-box"
    )
  })
  
  observeEvent(input$start_rselenium, {
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    # Retrieve all amino acid sequences for a given gene
    gene_name <- target()
    sequences <- getSequence(id = gene_name, type = "hgnc_symbol", seqType = "peptide", mart = ensembl)
    
    # Function to find and keep only the longest sequence
    get_longest_sequence <- function(sequences) {
      # Compute lengths of all sequences
      sequence_lengths <- nchar(sequences$peptide)
      
      # Identify the maximum length
      max_length <- max(sequence_lengths)
      
      # Get indices of sequences with the maximum length
      longest_indices <- which(sequence_lengths == max_length)
      
      # Keep only the first longest sequence in case of a tie
      longest_sequence <- sequences$peptide[longest_indices[1]]
      
      return(longest_sequence)
    }
    
    # Get the longest sequence
    longest_sequence <- get_longest_sequence(sequences)
    
    # Remove the asterisk at the end if present
    remove_asterisk <- function(sequence) {
      if (substr(sequence, nchar(sequence), nchar(sequence)) == "*") {
        sequence <- substr(sequence, 1, nchar(sequence) - 1)
      }
      return(sequence)
    }
    
    cleaned_sequence <- remove_asterisk(longest_sequence)
    
    # Print the cleaned sequence
    print(cleaned_sequence)
    
    # Dynamically assign a new port for each RSelenium session
    rD <- rsDriver(browser = "firefox", port = free_port(), check = TRUE)
    remDr <- rD[["client"]]
    
    # Navigate to the DeepTMHMM website
    deep_tmhmm_url <- "https://dtu.biolib.com/DeepTMHMM"
    remDr$navigate(deep_tmhmm_url)
    
    # Wait for the page to load
    Sys.sleep(5)
    
    # Find the text box element and enter the sequence
    text_box <- remDr$findElement(using = "css selector", value = "textarea")
    text_box$sendKeysToElement(list(cleaned_sequence))
    
    # Wait for the sequence to be entered
    Sys.sleep(2)
    
    # Find and click the run button
    run_button <- remDr$findElement(using = "css selector", value = ".bp4-popover-target")
    run_button$clickElement()
    
    
  })
  
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
