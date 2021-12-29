# Loading packages
library(tidyverse)
library(shiny)
library(shinyalert)
library(shinydashboard)
library(shinythemes)
library(shinyWidgets)

# Loading modules
source("waterfallPlot_module.R")
source("kaplanMeierPlot_module.R")
source("degTable_module.R")

#################### Loading external data ########################################
# PLEASE NOTE: Large expression datasets required for this app to function are *not* stored in the Github repo,
# as the filesizes are >1 GB. To access the expression data & run this app on your local machine, 
# the appropriate "data" directory will need to be copied down from our AWS S3 bucket to replace the "data" folder in the local repo.
# When pushing/pulling back to Github, the "data" folder in this repo will be ignored! 
# To modify the data used by this app, the Shiny app "data" object in S3 must be modified. Changes to the
# local "data" folder accessed by the Shiny scripts will NOT affect the web version of the app.

# May be able to retire this function - the use of a `start.R` script negates need for data loading progress bar.
# Will keep code here for a bit longer until I am certain this can be removed.
# readData <- function() {
#   # Creating a progress bar to let the user know when the expression data is done loading.
#   # Code adapted from http://www.mazsoft.com/blog/post/2018/01/01/show-progress-bar-when-pre-loading-data-in-shiny-app
#   progress <- Progress$new()
#   progress$set(value = 0.0, message = 'Please wait, loading data...')
#   progress$set(value = 0.10, message = 'Loading mRNA expression data...')
#   target_expData38 <<- readRDS("data/mRNA/TARGET_RBD_Dx_AML_ExpressionData_TPM_filt4dupGenes_with_cellLines_CD34posNBM_DSAML_MPN_GRCh38_12.17.2020_FinalforShiny.RDS")
#   target_expData37 <<- readRDS("data/mRNA/TARGET_RBD_Dx_AML_ExpressionData_TPM_filt4dupGenes_with_cellLines_CD34posNBM_DSAML_MPN_GRCh37_10.16.2020_FinalforShiny.RDS")
#   progress$set(value = 0.25, message = 'Loading mRNA expression data...')
#   beatAML_expData <<- readRDS("data/mRNA/BeatAML_Supplementary_Tables_TPM_Linear_Scale.RDS") %>%
#     column_to_rownames("geneSymbol")
#   swog_expData <<- readRDS("data/mRNA/SWOG_AML_ExpressionData_TPM_GRCh38_FinalforShiny.RDS")
#   laml_expData <<- readRDS("data/mRNA/TCGA_LAML_ExpressionData_TPM_FinalforShiny.RDS")
#   
#   progress$set(value = 0.50, message = 'Loading mature miRNA data...')
#   load("data/miRNA/TARGET_AML_AAML1031_expn_matrix_mimat_norm_miRNA_RPM_01.07.2019_FinalforShiny.RData", .GlobalEnv)
#   miRmapping <<- read.csv("data/miRNA/hsa_gff3_IDMap.csv")
#   
#   progress$set(value = 0.75, message = 'Loading clinical data...')
#   load("data/Clinical/Beat_AML_Supplementary_ClinicalData_FinalforShiny.RData", .GlobalEnv)
#   load("data/Clinical/TARGET_AML_merged_CDEs_Shareable_FinalforShiny.RData", .GlobalEnv)
#   load("data/Clinical/SWOG_AML_Merged_CDEs_FinalforShiny.RData", .GlobalEnv)
#   load("data/Clinical/TCGA_LAML_ClinicalData_FinalforShiny.RData", .GlobalEnv)
#   
#   load("data/ADC_and_CARTcell_Targets_Database_ADCReview_clinicaltrialsGov_FinalforShiny.RData", .GlobalEnv)
#   load("data/DEGs/TARGET_AML_vs_NBM_and_Others_Ribodepleted_DEGs_per_Group_GRCh37_12.18.2020_FinalforShiny.RData", .GlobalEnv)
#   deColKey <<- read.csv("data/Limma_Column_Descriptions.csv")
#   progress$set(value = 1, message = 'Done loading!')
#   progress$close()
# }

# readData <- function(x) {
  # Reading in required data for app to function
  target_expData38 <<- readRDS("data/mRNA/TARGET_RBD_Dx_AML_ExpressionData_TPM_filt4dupGenes_with_cellLines_CD34posNBM_DSAML_MPN_GRCh38_12.17.2020_FinalforShiny.RDS")
  target_expData37 <<- readRDS("data/mRNA/TARGET_RBD_Dx_AML_ExpressionData_TPM_filt4dupGenes_with_cellLines_CD34posNBM_DSAML_MPN_GRCh37_10.16.2020_FinalforShiny.RDS")
  beatAML_expData <<- readRDS("data/mRNA/BeatAML_Supplementary_Tables_TPM_Linear_Scale.RDS") %>%
    column_to_rownames("geneSymbol")
  swog_expData <<- readRDS("data/mRNA/SWOG_AML_ExpressionData_TPM_GRCh38_FinalforShiny.RDS")
  laml_expData <<- readRDS("data/mRNA/TCGA_LAML_ExpressionData_TPM_FinalforShiny.RDS")
  
  load("data/miRNA/TARGET_AML_AAML1031_expn_matrix_mimat_norm_miRNA_RPM_01.07.2019_FinalforShiny.RData", .GlobalEnv)
  miRmapping <<- read.csv("data/miRNA/hsa_gff3_IDMap.csv")
  
  load("data/Clinical/Beat_AML_Supplementary_ClinicalData_FinalforShiny.RData", .GlobalEnv)
  load("data/Clinical/TARGET_AML_merged_CDEs_Shareable_FinalforShiny.RData", .GlobalEnv)
  load("data/Clinical/SWOG_AML_Merged_CDEs_FinalforShiny.RData", .GlobalEnv)
  load("data/Clinical/TCGA_LAML_ClinicalData_FinalforShiny.RData", .GlobalEnv)
  
  load("data/ADC_and_CARTcell_Targets_Database_ADCReview_clinicaltrialsGov_FinalforShiny.RData", .GlobalEnv)
  load("data/DEGs/TARGET_AML_vs_NBM_and_Others_Ribodepleted_DEGs_per_Group_GRCh37_12.18.2020_FinalforShiny.RData", .GlobalEnv)
  deColKey <<- read.csv("data/Limma_Column_Descriptions.csv")
  colMapping <<- read.csv("data/Dataset_Column_Mapping_File.csv", check.names = F, na.strings = "")
  protPaint_html <<- includeHTML("www/Protein_Paint/embed_StJude_ProteinPaint_test.html")
# }

# if (is.null(target_expData38)) {
  # readData()
  # print("Testing mode - data already in environment")
# }

# Misc global functions & files
`%then%` <- function(a, b) {
  if (is.null(a)) b else a
}
