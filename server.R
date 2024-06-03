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
      need(target(), FALSE))
    
    valueBox(value = tags$p("GTEx", style = "font-size: 60%"),
             subtitle = "Normal tissue expression",
             color = "light-blue",
             icon = icon("prescription-bottle"),
             href = paste0("https://gtexportal.org/home/gene/", target(), "#geneExpression"))
  })
  
  output$protPaint <- renderValueBox({
    validate(
      need(target(), FALSE))
    
    valueBox(value = tags$p("ProteinPaint", style = "font-size: 60%"),
             subtitle = "St. Jude PeCan visualization",
             color = "light-blue",
             icon = icon("prescription-bottle"), 
             href = paste0("https://proteinpaint.stjude.org/?genome=hg19&gene=", target(), "&dataset=pediatric"))
  })
  
  
  # function for creating the action button for embedding DeepTMHMM
  output$tmhmm <- renderUI({
    validate(
      need(target(), FALSE)
    )
    
    actionButton("start_deeptmhmm", 
                 label = div(
                   "DeepTMHMM",
                   div("Protein localization", style = "text-transform: none; color: white; font-size: 15px; font-weight: normal; margin-top:20px; margin-bottom:20px;") # Additional white text below the label
                 ),
                 style = "text-transform: none; background-color: #3c8dbc; box-shadow: none; text-align: left; font-size: 21px; font-weight: bold; height: 110px; width: 100%; padding: 10px;",
                 class = "btn-box"
    )
  })
  
  # creating a reactive value that will change once the output is finished
  output_completed <- reactiveVal(FALSE)
  
  # this is the function for outputting the terminal and cleaning it up
  poll_terminal_output <- function(myTerm) {
    function() {
      output <- NULL
      
      # this is a weird workaround, you can't use the system() function if it's inside of RSTUDIO for some reason
      # instead you have to use this rstudioapi package to create a new terminal
      if (Sys.getenv("RSTUDIO") == "1") { 
        full_output <- rstudioapi::terminalBuffer(myTerm)
        # extract lines starting from "Running DeepTMHMM..." and ending with "Step 4/4"
        start_idx <- grep("^Running DeepTMHMM...", full_output)
        end_idx <- grep("^Step 4/4", full_output)
        if (length(start_idx) > 0) {
          if (length(end_idx) > 0) {
            output <- full_output[start_idx:end_idx]
            output_completed(TRUE)  # update the reactive value when the task is completed
          } else {
            output <- full_output[start_idx:length(full_output)]
          }
          # filter out empty and NULL lines
          output <- output[output != "" & !is.null(output)]
        }
      }
      if (is.null(output)) {
        output <- character(0)
      }
      output
    }
  }
  
  # when the action button is pushed, we're starting the deeptmhmm code
  observeEvent(input$start_deeptmhmm, {
    
    if (dir.exists("biolib_results")) {
      unlink("biolib_results", recursive = TRUE)
    }
    
    # the goal of this is to take the gene name and find the longest peptide sequence (using that as canonical for now)
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    gene_name <- target()
    sequences <- getSequence(id = gene_name, type = "hgnc_symbol", seqType = "peptide", mart = ensembl)
    
    # this is the longest sequence function 
    get_longest_sequence <- function(sequences) {
      sequence_lengths <- nchar(sequences$peptide)
      max_length <- max(sequence_lengths)
      longest_indices <- which(sequence_lengths == max_length)
      longest_sequence <- sequences$peptide[longest_indices[1]]
      return(longest_sequence)
    }
    
    longest_sequence <- get_longest_sequence(sequences)
    
    # the asterisk indicates a stop codon, I'm pretty sure this doesn't matter whether we remove this or not
    # but for now, I'm keeping it in
    remove_asterisk <- function(sequence) {
      if (substr(sequence, nchar(sequence), nchar(sequence)) == "*") {
        sequence <- substr(sequence, 1, nchar(sequence) - 1)
      }
      return(sequence)
    }
    
    cleaned_sequence <- remove_asterisk(longest_sequence)
    
    
    # create a temporary fasta file
    temp_fasta <- tempfile(fileext = ".fasta")
    
    # construct the file
    writeLines(paste0(">header\n", cleaned_sequence), con = temp_fasta)
    
    # now we run the python biolib package inside a terminal we've opened up in R :D
    myTerm <- NULL
    if (Sys.getenv("RSTUDIO") == "1") { 
      myTerm <- rstudioapi::terminalCreate(show = FALSE)
      tmhmm <- paste("biolib run DTU/DeepTMHMM --fasta", temp_fasta)
      rstudioapi::terminalSend(myTerm, paste0(tmhmm, "\n"))
    } else {
      system(paste("biolib run DTU/DeepTMHMM --fasta", temp_fasta))
    }
    
    # this is a reactive function that will write the terminal output as it goes
    terminal_output <- reactivePoll(1000, session, checkFunc = poll_terminal_output(myTerm), valueFunc = poll_terminal_output(myTerm))
    
    # this gives a little loading text before the output starts
    output$terminal_output <- renderText({
      output <- terminal_output()
      if (length(output) == 0) {
        "Waiting for DeepTMHMM to start, this may take a few seconds... \n"
      } else {
        paste(output, collapse = "\n")
      }
    })
  })
  
  # this gets the file path for the result image and will paste it in shiny
  observeEvent(output_completed(), {
    if (output_completed()) {
      output$tmhmm_plot <- renderImage({
        filename <- normalizePath(file.path('biolib_results', 'plot.png'))
        
        if (!file.exists(filename)) {
          stop("File does not exist: ", filename)
        }
        
        list(
          src = filename,
          alt = "TMHMM Plot",
          width = "80%",
          height = "auto"
        )
      }, deleteFile = FALSE)
    }
  })
  
  
  output$therapyTable <- DT::renderDataTable({
    validate(
      need(target(), "Please enter a gene symbol in the text box.") %then%
        need(target() %in% adc_cart_targetData$`Gene target`, paste0("We do not have record of ", target(), " being targeted\n by ADC or CAR T-cell therapies.")))
    
    table <- adc_cart_targetData %>%
      filter(`Gene target` == target()) 
    
    DT::datatable(table, 
                  options = list(scrollY = "25vh",
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
