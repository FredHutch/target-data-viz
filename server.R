library(shiny) 
library(shinyalert)
library(tidyverse)
library(DT)
library(dplyr)
source("waterfallPlot_module.R")
source("kaplanMeierPlot_module.R")
`%then%` <- shiny:::`%OR%` # See https://shiny.rstudio.com/articles/validation.html for details on the %then% operator

#################### Loading external data ########################################
# PLEASE NOTE: Large expression datasets required for this app to function are *not* stored in the Github repo,
# as the filesizes are >1 GB. To access the expression data & run this app on your local machine, 
# the appropriate "data" directory will need to be copied down from our AWS S3 bucket to replace the "data" folder in the local repo.
# When pushing/pulling back to Github, the "data" folder in this repo will be ignored! 
# To modify the data used by this app, the Shiny app "data" object in S3 must be modified. Changes to the
# local "data" folder accessed by the Shiny scripts will NOT affect the web version of the app.

readData <- function(target_cde, target_expData, beatAML_cde, beatAML_expData, adc_cart_targetData) {
  # Creating a progress bar to let the user know when the expression data is done loading.
  # Code adapted from http://www.mazsoft.com/blog/post/2018/01/01/show-progress-bar-when-pre-loading-data-in-shiny-app
  progress <- Progress$new()
  progress$set(value = 0.0, message = 'Please wait, loading data...')
  progress$set(value = 0.10, message = 'Loading mRNA expression data...')
  target_expData <<- readRDS("data/mRNA/TARGET_RBD_Dx_AML_ExpressionData_TPM_filt4dupGenes_with_cellLines_CD34posNBM_DSAML_MPN_10.16.2020_FinalforShiny.RDS")
  progress$set(value = 0.25, message = 'Loading mRNA expression data...')
  beatAML_expData <<- readRDS("data/mRNA/BeatAML_Supplementary_Tables_TPM_Linear_Scale.RDS") %>%
    column_to_rownames("geneSymbol")
  progress$set(value = 0.50, message = 'Loading mature miRNA data...')
  load("data/miRNA/TARGET_AML_AAML1031_expn_matrix_mimat_norm_miRNA_RPM_01.07.2019_FinalforShiny.RData", .GlobalEnv)
  miRmapping <<- read.csv("data/miRNA/hsa_gff3_IDMap.csv")
  progress$set(value = 0.75, message = 'Loading clinical data...')
  load("data/Clinical/Beat_AML_Supplementary_ClinicalData_FinalforShiny.RData", .GlobalEnv)
  load("data/Clinical/TARGET_AML_0531_1031_merged_CDEs_Shareable_9.18.20_FinalforShiny.RData", .GlobalEnv)
  load("data/ADC_and_CARTcell_Targets_Database_ADCReview_clinicaltrialsGov_FinalforShiny.RData", .GlobalEnv)
  progress$set(value = 1, message = 'Done loading!')
  progress$close()
}

##############################################################################

target_expData <- NULL # Setting this to null so it will only be read one time

server <- function(input, output, session) { 
  
  # Reading in the expression & clinical data one time, immediately after the app starts up
  if (is.null(target_expData)) {
    readData(target_cde, target_expData, beatAML_cde, beatAML_expData, adc_cart_targetData)
    # print("Testing mode - data already in environment")
  }
  
  cohort <- reactive({
    input$seqDataCohort
  })
  
  seqData <- reactive({
    switch(input$seqDataCohort,
           "BeatAML" = beatAML_expData,
           "TARGET" = target_expData)
  })
  
  studyData <- reactive({
    switch(input$seqDataCohort,
           "BeatAML" = beatAML_cde,
           "TARGET" = target_cde)
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
    # EDGE-CASES: What if nothing is found? What will the message in the popup box be?
    # What if 100+ values are found? The box is so overloaded you can't exit out of it,
    # needs a max length of values that it can return
    # Also, when you click OUTSIDE the box, it just returns "false"
    shinyalert(
      title = "Did you mean...",
      text = msg,
      # size = "xs",
      closeOnEsc = TRUE,
      closeOnClickOutside = FALSE,
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
  
  #--------------------- External databases tab --------------------- #

  output$protAtlas <- renderInfoBox({
    
    validate(
      need(target(), "Please enter a gene symbol or miRNA in the text box."))
    
        infoBox(value = "Human Protein Atlas", 
                 title = "Protein expression",
                 color = "red", 
                 icon = icon("prescription-bottle"), href = paste0("https://www.proteinatlas.org/search/", target()))
  })
  
  output$gtex <- renderInfoBox({
    
    validate(
      need(target(), "Please enter a gene symbol or miRNA in the text box."))
    
    infoBox(value = "GTEX Gene \nExpression", 
            title = "Normal tissue expression",
            color = "orange", fill = F,
            icon = icon("prescription-bottle"), href = paste0("https://gtexportal.org/home/gene/", target(), "#geneExpression"))
  })
  
  output$therapyTable <- DT::renderDataTable({
    
    validate(
      need(target(), "Please enter a gene symbol in the text box.") %then%
        need(target() %in% adc_cart_targetData$`Gene target`, paste0("We do not have record of ", target(), "being targeted\n for ADC or CAR T-cell therapies.")))
    
    table <- adc_cart_targetData %>%
      filter(`Gene target` == target()) 
    
    DT::datatable(table, 
                  options = list(scrollY = "50vh"), 
                  escape = F)
  })

  
  # Following this post, but it doesn't work: https://stackoverflow.com/questions/24875943/display-html-file-in-shiny-app
  # This person is having the same issue I am:
  # https://stackoverflow.com/questions/56064805/displaying-html-file-using-includehtml-in-shiny-is-not-working-with-renderui
  
  # output$test <- renderUI({

    # Method 1
    # includeHTML() is designed to work with HTML fragments, so a "self contained" HTML file is needed
    # includeHTML("www/UMAP/TARGET_AMLdx_rlps_NBM_PCAselect_selfcontained.html") # Doesn't work, produces a blank page
    # HTML(knitr::knit2html("About.Rmd", fragment.only = TRUE))

    # Method 2, see iframe details at https://plotly-r.com/saving.html
    # tags$iframe(seamless="seamless",
                # src="www/UMAP/plotly/TARGET_AMLdx_rlps_NBM_PCAselect_not_selfcontained_bodyOnly.html",
                # height=600, width=800, scrolling="yes")
    # })
  
}
