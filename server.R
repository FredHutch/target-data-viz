library(shiny)
library(tidyverse)
source("waterfallPlot_module.R")
source("kaplanMeierPlot_module.R")

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
  
  # Calling the waterfall plot module
  callModule(wfPlot, id = "waterfall", clinData = clinData, countsData = countsData, gene = target)
  
  # Calling the Kaplan-Meier curve module
  callModule(kmPlot, id = "kaplanmeier", clinData = clinData, countsData = countsData, gene = target)
}