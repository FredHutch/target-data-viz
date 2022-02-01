######### Loading packages
library(tidyverse)
library(shiny)
library(shinyalert)
library(shinydashboard)
library(shinythemes)
library(shinyWidgets)
library(shinyjs)
library(shinyBS)
library(data.table)
library(DT)

# https://rstudio.github.io/bslib/articles/bslib.html#bootswatch New themes package to check out

######### Loading modules
source("waterfallPlot_module.R")
source("kaplanMeierPlot_module.R")
source("degTable_module.R")
source("geneExpressors_module.R")
source("heatmap_module.R")

######### Loading external data
# PLEASE NOTE: Large expression datasets required for this app to function are *not* stored in the Github repo,
# as the filesizes are >1 GB. To access the expression data & run this app on your local machine, 
# the appropriate "data" directory will need to be copied down from our AWS S3 bucket to replace the "data" folder in the local repo.
# When pushing/pulling back to Github, the "data" folder in this repo will be ignored! 
# To modify the data used by this app, the Shiny app "data" object in S3 must be modified. Changes to the
# local "data" folder accessed by the Shiny scripts will NOT affect the web version of the app.

readData <- function(x) {
  
  # Gene-level expression matrices
  target_expData38 <<- readRDS("data/mRNA/TARGET_RBD_Dx_AML_ExpressionData_TPM_filt4dupGenes_with_cellLines_CD34posNBM_DSAML_MPN_GRCh38_12.17.2020_FinalforShiny.RDS")
  target_expData37 <<- readRDS("data/mRNA/TARGET_RBD_Dx_AML_ExpressionData_TPM_filt4dupGenes_with_cellLines_CD34posNBM_DSAML_MPN_GRCh37_10.16.2020_FinalforShiny.RDS")
  beatAML_expData <<- readRDS("data/mRNA/BeatAML_Supplementary_Tables_TPM_Linear_Scale.RDS") %>%
    column_to_rownames("geneSymbol")
  swog_expData <<- readRDS("data/mRNA/SWOG_AML_ExpressionData_TPM_GRCh38_FinalforShiny.RDS")
  laml_expData <<- readRDS("data/mRNA/TCGA_LAML_ExpressionData_TPM_FinalforShiny.RDS")

  # miRNA expression matrices (for TARGET dataset only)
  load("data/miRNA/TARGET_AML_AAML1031_expn_matrix_mimat_norm_miRNA_RPM_01.07.2019_FinalforShiny.RData", .GlobalEnv)
  miRmapping <<- read.csv("data/miRNA/hsa_gff3_IDMap.csv")

  # Clinical data
  load("data/Clinical/Beat_AML_Supplementary_ClinicalData_FinalforShiny.RData", .GlobalEnv)
  load("data/Clinical/TARGET_AML_merged_CDEs_Shareable_FinalforShiny.RData", .GlobalEnv)
  load("data/Clinical/SWOG_AML_Merged_CDEs_FinalforShiny.RData", .GlobalEnv)
  load("data/Clinical/TCGA_LAML_ClinicalData_FinalforShiny.RData", .GlobalEnv)

  # Misc accessory files
  load("data/ADC_and_CARTcell_Targets_Database_ADCReview_clinicaltrialsGov_FinalforShiny.RData", .GlobalEnv)
  load("data/DEGs/TARGET_AML_vs_NBM_and_Others_Ribodepleted_DEGs_per_Group_GRCh37_12.18.2020_FinalforShiny.RData", .GlobalEnv)
  deColKey <<- read.csv("data/Limma_Column_Descriptions.csv")
  colMapping <<- read.csv("data/Dataset_Column_Mapping_File.csv", check.names = F, na.strings = "")
}

testing <- TRUE

if (testing == TRUE) {
  print("Testing mode - data already in environment")
} else {
  readData()
}

####### Misc global functions & files
`%then%` <- function(a, b) {
  if (is.null(a)) b else a
}

bs <- 16 # Base font size for figures