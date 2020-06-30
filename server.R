library(shiny)
library(tidyverse)
library(DT)
library(dplyr)
source("waterfallPlot_module.R")
source("kaplanMeierPlot_module.R")
`%then%` <- shiny:::`%OR%` # See https://shiny.rstudio.com/articles/validation.html for details on the %then% operator

################ Loading external data ########################################
# Creating a progress bar to let the user know when the expression data is done loading.
# Code adapted from http://www.mazsoft.com/blog/post/2018/01/01/show-progress-bar-when-pre-loading-data-in-shiny-app
# NOTE: GitHub has a max allowed filesize of 100 MB, so the expression data had to be broken into multiple smaller files in order to be uploadable
readData <- function(target_cde, target_expData) {
  
  # NOTE: Could potentially use the switch() function to load different datasets, depending on which dataset is selected in the dashboard
  progress <- Progress$new()
  progress$set(value = 0.0, message = 'Loading mRNA expression data...')
  target_expData <<- readRDS("data/mRNA/TARGET_AAML1031_0531_RBD_Dx_Relapse_TPMs_filt4dupGenes_with_cellLines_CD34posNBM_forShinyApp_11.20.2019.RDS")
  beatAML_expData <<- readRDS("data/mRNA/BeatAML_Supplementary_Tables_TPM_Linear_Scale.RDS") %>%
    column_to_rownames("geneSymbol")
  progress$set(value = 0.85, message = 'Loading mature miRNA data...')
  miRdata <- read.csv("data/miRNA/TARGET_AML_AAML1031_expn_matrix_mimat_norm_miRNA_RPM_01.07.2019.csv", stringsAsFactors = F) %>%
    rename_at(vars(-mir), ~gsub("\\..*", "", .)) %>% # Cutting off the aliquot ID and sample timepoint info (they're all Dx)
    separate(mir, into = c("hsa ID", "Accession Number"), sep = "\\.") 
  
  miRmapping <<- miRdata[,c("hsa ID", "Accession Number")]
  
  miRdata <<- miRdata %>% 
    dplyr::select(-"hsa ID") %>%
    column_to_rownames("Accession Number")
  
  progress$set(value = 0.95, message = 'Loading clinical data...')
  beatAML_cde <<- readxl::read_excel("data/Clinical/Beat_AML_Supplementary_ClinicalData.xlsx") %>%
    filter(isDenovo == TRUE) %>%
    rename(Primary.Fusion = finalFusion) %>%
    rename(USI = LabId) %>%
    rename(Event.Type = causeOfDeath) %>%
    rename(Cyto.Risk = ELN2017) %>%
    rename(FLT3.ITD = `FLT3-ITD`) %>%
    mutate(Group = "AML",
           MLL.Fusion = ifelse(str_detect(pattern = "KMT2A", string = Primary.Fusion), Primary.Fusion, "Other AML"), 
           Age.Category = case_when(ageAtDiagnosis < 10 ~ "Less than 10 years", 
                                    ageAtDiagnosis >= 10 & ageAtDiagnosis < 18 ~ "Between 10 and 18 years", 
                                    ageAtDiagnosis >= 18 & ageAtDiagnosis < 40 ~ "Between 18 and 40 years", 
                                    ageAtDiagnosis >= 40 & ageAtDiagnosis < 60 ~ "Between 40 and 60 years", 
                                    ageAtDiagnosis >= 60 ~ "Greater than 60 years"), 
           `Recoded OS ID` = ifelse(vitalStatus == "Dead", 1, 0), 
           `OS time (days)` = overallSurvival) %>%
    mutate_at(vars(Cyto.Risk), ~case_when(. == "FavorableOrIntermediate" ~ "Favorable Or Intermediate",
                                           . == "IntermediateOrAdverse" ~ "Intermediate Or Adverse",
                                           TRUE ~ .))
  
  target_cde <- readRDS("data/Clinical/AAML1031_TARGET_CDEs_with_HiAR_PrimaryCyto_and_FusionCalls_9.3.19.RDS") %>%
    mutate_at(vars(Reg.), ~as.character(.)) %>%
    dplyr::rename(Primary.Cyto = `Primary Cytogenetic Code`) %>%
    dplyr::rename(Cyto.Risk = `Cyto/Fusion/Molecular Risk`) %>%
    dplyr::rename(CNS.Disease.Category = `CNS disease at on-study`) %>%
    dplyr::rename(Event.Type = `EFS event type ID`) %>%
    mutate_at(vars(Cytogenetic.Category.2), ~gsub("\\.|\\_", " ", .)) %>%
    mutate_at(vars(Rare.Fusions), ~gsub("\\.", "\\-", .)) %>%
    mutate_at(vars(SNVs), ~gsub("\\.", " ", .)) %>%
    mutate_at(vars(SNVs), ~gsub("mutation", "mut", .)) %>%
    mutate_at(vars(SNVs), ~gsub("with", "w\\/", .)) %>%
    mutate_at(vars(Cytogenetic.Category.2, Rare.Fusions, SNVs), ~gsub("OtherAML", "Other AML", .)) %>%

    # Adding columns to use for re-categorizing MLL and other fusions (if they occur in less than 10 patients)
    mutate(Primary.Fusion = ifelse(str_detect(pattern = "KMT2A-", string = `Primary Fusion/CNV`), "KMT2A-X", `Primary Fusion/CNV`)) %>%
    mutate(MLL.Fusion = ifelse(str_detect(pattern = "KMT2A", string = `Primary Fusion/CNV`), `Primary Fusion/CNV`, "Other AML")) %>%
    group_by(MLL.Fusion) %>%
    mutate(Fusion.Count = n()) %>%
    ungroup() %>%
    mutate_at(vars(MLL.Fusion), ~ifelse(Fusion.Count < 10, "Other MLL", .)) %>%
    group_by(Primary.Fusion) %>%
    mutate(Fusion.Count = n()) %>%
    ungroup() %>%
    mutate_at(vars(Primary.Fusion), ~ifelse(Fusion.Count < 10, "Other AML", .)) %>%
    
    # Adding rows for the cell lines & other controls, which are not typically included in the clinical data
    tibble::add_row(USI = c("K562.01", "K562.02", "ME1", "MO7E", "NOMO1", "Kasumi.D1", "MV4.11.D1"), Group = NA) %>%
    tibble::add_row(USI = grep("^RO|^BM", colnames(target_expData), value = T), Group = "NBM") %>%
    tibble::add_row(USI = grep("34POS", colnames(target_expData), value = T), Group = "CD34+ PB")
  
  mutData <- readxl::read_excel("data/Clinical/TARGET_AML_1031_merged_CDEs_mutationsOnly_01.08.20.xlsx")
  target_cde <<- left_join(target_cde, mutData, by = c("USI", "Reg.")) %>%
    mutate_at(vars(RAS.Gene), ~factor(., levels = c("KRAS", "NRAS", "Unknown", "None")))
  adc_cart_targetData <<- readxl::read_excel("data/ADC_and_CARTcell_Targets_Database_ADCReview_clinicaltrialsGov.xlsx") %>%
    dplyr::select(`Treatment type`, `Gene symbol of target (Final)`, `Associated cancer(s)`, 
                  `If currently in clinical trials, drug/trial ID number`, `Development sponsor`, `Development status`) %>%
    mutate_at(vars(`If currently in clinical trials, drug/trial ID number`), ~ifelse(is.na(.), ., paste0("<a href='", "https://clinicaltrials.gov/show/", ., "'>", ., "</a>"))) %>%
    rename(`Gene target` = `Gene symbol of target (Final)`)
  
  progress$set(value = 1, message = 'Done loading!')
  progress$close()
}
##############################################################################

target_expData <- NULL

server <- function(input, output, session) { 
  
  # Reading in the expression & clinical data one time, immediately after the app starts up
  if (is.null(target_expData)) {
    readData(target_cde, target_expData)
  }
  
  # Creating a variable that will be used to reactively pass the gene of interest into each module,
  # see https://tbradley1013.github.io/2018/07/20/r-shiny-modules--using-global-inputs/ for more  
  # info on passing global Shiny variables into a module
  target <- reactive({
    if (grepl(x = input$geneInput, "$hsa\\-mir*|mir\\-*", ignore.case = T)) {
      symbol <- miRmapping$`Accession Number`[match(tolower(input$geneInput), miRmapping$`hsa ID`)] 
    } else {
      symbol <- toupper(input$geneInput)
    }
    return(symbol)
  })
  
  cohort <- reactive({
    input$seqDataCohort
  })
  
  seqData <- reactive({
    if (input$seqDataCohort == "BeatAML") beatAML_expData else
      if (input$seqDataCohort == "TARGET" && input$seqDataType == "mRNA") target_expData else
        if (input$seqDataCohort == "TARGET" && input$seqDataType == "miRNA") miRdata
  })
  
  studyData <- reactive({
    switch(input$seqDataCohort,
           "BeatAML" = beatAML_cde,
           "TARGET" = target_cde)
  })
  
  #--------------------- WF plot & KM plot tabs --------------------- #
  
  # Calling the waterfall plot module - this is broken, doesn't recognize the gene/target
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
        need(target() %in% adc_cart_targetData$`Gene target`, "We do not have record of that gene being targeted\n for ADC or CAR T-cell therapies."))
    
    table <- adc_cart_targetData %>%
      filter(`Gene target` == target()) 
    
    DT::datatable(table, 
                  options = list(scrollY = "50vh"), 
                  escape = F)
  })

  
  
}