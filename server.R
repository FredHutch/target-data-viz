# SHINY DAHSBOARD SERVER

############  TO PUSH UPDATES TO THE SHINYAPPS.IO SERVER ######################
# library(rsconnect)
# rsconnect::deployApp('/Volumes/homes/Shiny/kaplanMeierPlot_app')

library(shiny)
library(ggplot2)
library(tidyverse)
library(data.table)
library(DT)
source("waterfallPlot_module.R")
source("kaplanMeierPlot_module.R")

# Function to read in counts % clinical data, if they haven't been used yet
readData <- function(clinData, countsData){
  progress <- Progress$new()
  progress$set(value = 0.0, message = 'Loading expression data...')
  progress$set(value = 0.25, message = 'Loading expression data...')
  progress$set(value = 0.45, message = 'Loading expression data...')
  countsData <<- readRDS("data/TARGET_AAML1031_0531_RBD_Dx_Rlps_TPMs_filt_dupGenes_w_cellLines_CD34posNBM.RDS")
  progress$set(value = 0.75, message = 'Loading expression data...')
  progress$set(value = 0.90, message = 'Loading expression data...')
  clinData <<- readRDS("data/AAML1031_TARGET_CDEs_with_HiAR_PrimaryCyto_and_FusionCalls_9.3.19.RDS")
  progress$set(value = 1, message = 'Done loading!')
  progress$close()
}

################ Loading external data ########################################
# Creating a progress bar to let the user know when the expression data is done loading.
# Code adapted from http://www.mazsoft.com/blog/post/2018/01/01/show-progress-bar-when-pre-loading-data-in-shiny-app
countsData <- NULL
clinData <- NULL


##############################################################################

server <- function(input, output, session) { 
  
  if(is.null(countsData)){
    readData(clinData, countsData)
  }
  
  # See https://tbradley1013.github.io/2018/07/20/r-shiny-modules--using-global-inputs/ for more info 
  # on passing global Shiny variables into a module
  target <- reactive(input$geneInput)
  
  callModule(wfPlot, id = "waterfall", clinData = clinData, countsData = countsData, gene = target)
  
  callModule(kmPlot, id = "kaplanmeier", clinData = clinData, countsData = countsData, gene = target)
}