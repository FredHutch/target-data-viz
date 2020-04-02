library(shiny)
library(tidyverse)
library(DT)
source("waterfallPlot_module.R")
source("kaplanMeierPlot_module.R")
`%then%` <- shiny:::`%OR%` # See https://shiny.rstudio.com/articles/validation.html for details on the %then% operator

################ Loading external data ########################################
# Creating a progress bar to let the user know when the expression data is done loading.
# Code adapted from http://www.mazsoft.com/blog/post/2018/01/01/show-progress-bar-when-pre-loading-data-in-shiny-app
# NOTE: GitHub has a max allowed filesize of 100 MB, so the counts data had to be broken into multiple smaller files in order to be uploadable
readData <- function(clinData, countsData) {
  progress <- Progress$new()
  progress$set(value = 0.0, message = 'Loading expression data...')
  progress$set(value = 0.10, message = 'Loading expression data part 1...')
  p1 <- readRDS("data/TARGET_AAML1031_0531_RBD_Dx_Relapse_TPMs_filt4dupGenes_with_cellLines_CD34posNBM_forShinyApp_PART1_12.27.2019.RDS")
  progress$set(value = 0.25, message = 'Loading expression data part 2...')
  p2 <- readRDS("data/TARGET_AAML1031_0531_RBD_Dx_Relapse_TPMs_filt4dupGenes_with_cellLines_CD34posNBM_forShinyApp_PART2_12.27.2019.RDS")
  progress$set(value = 0.45, message = 'Loading expression data part 3...')
  p3 <- readRDS("data/TARGET_AAML1031_0531_RBD_Dx_Relapse_TPMs_filt4dupGenes_with_cellLines_CD34posNBM_forShinyApp_PART3_12.27.2019.RDS")
  progress$set(value = 0.85, message = 'Loading expression data part 4...')
  p4 <- readRDS("data/TARGET_AAML1031_0531_RBD_Dx_Relapse_TPMs_filt4dupGenes_with_cellLines_CD34posNBM_forShinyApp_PART4_12.27.2019.RDS")
  countsData <<- rbind(p1, p2, p3, p4)
  progress$set(value = 0.95, message = 'Loading clinical data...')
  clinData <<- readRDS("data/AAML1031_TARGET_CDEs_with_HiAR_PrimaryCyto_and_FusionCalls_9.3.19.RDS") %>% mutate_at(vars(Reg.), ~as.character(.))
  mutData <<- readxl::read_excel("data/TARGET_AML_1031_merged_CDEs_mutationsOnly_01.08.20.xlsx")
  clinData <<- left_join(clinData, mutData, by = c("USI", "Reg."))
  adc_cart_targetData <<- readxl::read_excel("data/ADC_and_CARTcell_Targets_Database_ADCReview_clinicaltrialsGov.xlsx")
  progress$set(value = 1, message = 'Done loading!')
  progress$close()
}
##############################################################################

countsData <- NULL
clinData <- NULL

server <- function(input, output, session) { 
  
  # Reading in the counts & clinical data one time, immediately after the app starts up
  if (is.null(countsData)) {
    readData(clinData, countsData)
  }
  
  # Turning clinical trial identifier into a URL to use for hyperlinking: https://clinicaltrials.gov/ct2/about-site/link-to
  adc_cart_targetData <- adc_cart_targetData %>%
    dplyr::select(`Treatment type`, `Gene symbol of target (Final)`, `Associated cancer(s)`, 
                  `If currently in clinical trials, drug/trial ID number`, `Development sponsor`, `Development status`) %>%
    mutate_at(vars(`If currently in clinical trials, drug/trial ID number`), ~ifelse(is.na(.), ., paste0("<a href='", "https://clinicaltrials.gov/show/", ., "'>", ., "</a>"))) %>%
    rename(`Gene target` = `Gene symbol of target (Final)`)
  
  # Creating a variable that will be used to reactively pass the gene of interest into each module,
  # see https://tbradley1013.github.io/2018/07/20/r-shiny-modules--using-global-inputs/ for more  
  # info on passing global Shiny variables into a module
  target <- reactive({
    if (startsWith(input$geneInput, "hsa-mir")) {
      newGene <- input$geneInput
    } else {
      newGene <- toupper(input$geneInput)
    }
    return(newGene)
  })
  
  # IMPORTANT NOTE: the "target" variable is actually a reactive function, and would usually be called by target(), 
  # but when passing it into a module, you *must* pull off the parentheses and pass the naked variable name as an argument.
  # Within the modules themselves, this variable is called gene() and is a reactive function!
  
  # Calling the waterfall plot module
  callModule(wfPlot, id = "waterfall", 
             clinData = clinData, 
             countsData = countsData, 
             adc_cart_targetData = adc_cart_targetData,
             gene = target, 
             parent = session) # See https://stackoverflow.com/questions/51708815/accessing-parent-namespace-inside-a-shiny-module
                               # for an explanation of the 'parent' parameter
  
  # Calling the Kaplan-Meier curve module
  callModule(kmPlot, id = "kaplanmeier", clinData = clinData, countsData = countsData, gene = target)
  
  #--------------------- External databases tab --------------------- #

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
  
  output$therapyTable <- DT::renderDataTable({
    
    validate(
      need(target(), "Please enter a gene symbol in the text box.") %then%
        need(target() %in% adc_cart_targetData$`Gene target`, "We do not have record of that gene being targeted\n for ADC development or CAR T-cell therapies."))
    
    table <- adc_cart_targetData %>%
      filter(`Gene target` == target()) 
    
    DT::datatable(table, 
                  options = list(scrollY = "50vh"), 
                  escape = F)
  })

  
  
}