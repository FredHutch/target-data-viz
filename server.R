# Check the file below for global scripts & variables, 
# they were removed from the server & ui scripts to help clean them up.
# source("global.R")

# Explicitly setting this to null before app startup. 
# This will allow a function (sourced in global.R) to trigger 
# all source file loading at script startup ONLY.
target_expData <- NULL 

server <- function(input, output, session) { 
  
  cohort <- reactive({
    input$seqDataCohort
  })
  
  seqData <- reactive({
    matrix <- switch(input$seqDataCohort,
                     "SWOG" = swog_expData,
                     "BeatAML" = beatAML_expData,
                     "TARGET" = target_expData38,
                     "TCGA" = laml_expData)
    
    assembly <- if (input$seqDataCohort == "TARGET" ) {
      switch(input$seqAssembly,
             "grch37" = target_expData37,
             "grch38" = target_expData38)
    } else {
      matrix
    }
  })
  
  studyData <- reactive({
    switch(input$seqDataCohort,
           "SWOG" = swog_cde,
           "BeatAML" = beatAML_cde,
           "TARGET" = target_cde,
           "TCGA" = laml_cde)
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
    option <- grep(input$geneInput, rownames(seqData()), value = T, ignore.case = T)
    msg <- if (length(option) == 0) {
      "No alternate names found!"
    } else if (length(option) > 20) {
      paste0(paste0(option[1:20], collapse = "  or\n"), ",\n...?")
    } else {
      paste0(paste0(option, collapse = "  or\n"), "?")
    }
    
    # Check if miRNA name exists in the miRbase21 miRNA-seq data
    # shinyalert(
    #   title = "Did you mean...",
    #   text = msg,
    #   # size = "xs",
    #   closeOnEsc = TRUE,
    #   closeOnClickOutside = FALSE,
    #   html = FALSE,
    #   type = "input",
    #   inputType = "text",
    #   inputValue = option[1],
    #   inputPlaceholder = "Type choice here",
    #   showConfirmButton = TRUE,
    #   showCancelButton = FALSE,
    #   confirmButtonText = "Retry",
    #   confirmButtonCol = "#41B0FA",
    #   timer = 0,
    #   imageUrl = "",
    #   animation = TRUE,
    #   callbackR = function(x){
    #     if (!is.null(x)) updateTextInput(session, "geneInput", label = NULL, value = input$shinyalert)
    #   })
    
    modalDialog(
      msg,
      title = "Did you mean...",
      size = "s",
      easyClose = F,
      placeholder = "Type choice here",
    )
  })
  
  #--------------------- WF plot & KM plot tabs --------------------- #
  
  # Calling the waterfall plot module
  # IMPORTANT NOTE: the "target" & "cohort" variables are actually a reactive function, and would usually be called by target() & cohort(), 
  # but when passing a reactive value into a module, you *must* pull off the parentheses and pass the naked variable name as an argument.
  # Within the modules themselves, these variables are a reactive function!
  callModule(wfPlot, id = "waterfall", 
             clinData = studyData, 
             expData = seqData, 
             adc_cart_targetData = adc_cart_targetData,
             gene = target, 
             dataset = cohort,
             parent = session) # See https://stackoverflow.com/questions/51708815/accessing-parent-namespace-inside-a-shiny-module
                               # for an explanation of the 'parent' parameter
  
  # Calling the Kaplan-Meier curve module
  callModule(kmPlot, id = "kaplanmeier", 
             clinData = studyData, 
             expData = seqData, 
             dataset = cohort,
             gene = target)
  
  
  # Calling the DEG table module
  callModule(deTable, id = "degs",
             table = degTables37,
             gene = target)
  
  #--------------------- External databases tab --------------------- #

  # TO DO: Add a searchable AML-restricted gene list to this tab
  
  output$protAtlas <- renderInfoBox({
    validate(
      need(target(), "Please enter a gene symbol in the text box."))
    
        infoBox(value = "Human Protein Atlas", 
                 title = "Protein expression",
                 color = "red", 
                 icon = icon("prescription-bottle"), href = paste0("https://www.proteinatlas.org/search/", target()))
  })
  
  output$gtex <- renderInfoBox({
    validate(
      need(target(), "Please enter a gene symbol in the text box."))
    
    infoBox(value = "GTEX Gene \nExpression", 
            title = "Normal tissue expression",
            color = "orange", fill = F,
            icon = icon("prescription-bottle"), href = paste0("https://gtexportal.org/home/gene/", target(), "#geneExpression"))
  })
  
  output$protPaint <- renderInfoBox({
    validate(
      need(target(), "Please enter a gene symbol in the text box."))
    
    infoBox(value = "St. Jude \nProteinPaint", 
            title = "PeCan visualization",
            color = "blue", fill = F,
            icon = icon("prescription-bottle"), href = paste0("https://proteinpaint.stjude.org/?genome=hg19&gene=", target(), "&dataset=pediatric"))
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
    tags$iframe(seamless = "seamless",
                style = "border-width: 0;",
                src = "UMAP/TARGET_AML_sg7655_blackBackground_clusters2_k31_PCAselect.html",
                height = 800, width = 1300, scrolling = "yes")
  })
  
}
