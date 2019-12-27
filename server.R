library(shiny)
library(tidyverse)
source("waterfallPlot_module.R")
source("kaplanMeierPlot_module.R")

################ Loading external data ########################################
# Creating a progress bar to let the user know when the expression data is done loading.
# Code adapted from http://www.mazsoft.com/blog/post/2018/01/01/show-progress-bar-when-pre-loading-data-in-shiny-app
# NOTE: GitHub has a max allowed filesize of 100 MB, so the counts data had to be broken into 3 files in order to be uploadable
readData <- function(clinData, countsData){
  progress <- Progress$new()
  progress$set(value = 0.0, message = 'Loading expression data part 1...')
  progress$set(value = 0.15, message = 'Loading expression data part 1...')
  c1 <- readRDS("data/TARGET_AAML1031_0531_RBD_Dx_Relapse_TPMs_filt4dupGenes_with_cellLines_CD34posNBM_forShinyApp_PART1_12.27.2019.RDS")
  progress$set(value = 0.35, message = 'Loading expression data part 2...')
  c2 <- readRDS("data/TARGET_AAML1031_0531_RBD_Dx_Relapse_TPMs_filt4dupGenes_with_cellLines_CD34posNBM_forShinyApp_PART2_12.27.2019.RDS")
  progress$set(value = 0.55, message = 'Loading expression data part 3...')
  c3 <- readRDS("data/TARGET_AAML1031_0531_RBD_Dx_Relapse_TPMs_filt4dupGenes_with_cellLines_CD34posNBM_forShinyApp_PART3_12.27.2019.RDS")
  progress$set(value = 0.75, message = 'Loading full expression dataset...')
  countsData <<- rbind(c1, c2, c3)
  progress$set(value = 0.95, message = 'Loading clinical data...')
  clinData <<- readRDS("data/AAML1031_TARGET_CDEs_with_HiAR_PrimaryCyto_and_FusionCalls_9.3.19.RDS")
  progress$set(value = 1, message = 'Done loading!')
  progress$close()
}
##############################################################################

countsData <- NULL
clinData <- NULL

server <- function(input, output, session) { 
  
  # Reading in the counts & clinical data one time, immediately when the app starts up
  if(is.null(countsData)){
    readData(clinData, countsData)
  }
  
  # Creating a variable that will be used to reactively pass the gene of interest into each module
  # See https://tbradley1013.github.io/2018/07/20/r-shiny-modules--using-global-inputs/ for more  
  # info on passing global Shiny variables into a module
  target <- reactive(input$geneInput)

  # Calling the waterfall plot module
  callModule(wfPlot, id = "waterfall", clinData = clinData, countsData = countsData, gene = target)
  
  # Calling the Kaplan-Meier curve generator module
  callModule(kmPlot, id = "kaplanmeier", clinData = clinData, countsData = countsData, gene = target)
}