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
library(viridis)
library(viridisLite)
library(ggplot2)
library(plotly)
library(fst)
library(ComplexHeatmap)
library(survival)
library(survminer)
library(ggsurvfit)
library(biomaRt)
#library(RSelenium)
#library(netstat)

# https://rstudio.github.io/bslib/articles/bslib.html#bootswatch New themes package to check out

######### Loading modules
source("waterfallPlot_module.R")
source("time_module.R")
source("kaplanMeierPlot_module.R")
#source("degTable_module.R")
source("geneExpressors_module.R")
#source("heatmap_module.R")
source("HPA_module.R")
source("oncoprint_module.R")
#source("classification_module.R")
source("cancerModule.R")

######### Loading external data
# PLEASE NOTE: Large expression datasets required for this app to function are *not* stored in the Github repo,
# as the filesizes are >1 GB. To access the expression data & run this app on your local machine, 
# the appropriate "data" directory will need to be copied down from our AWS S3 bucket to replace the "data" folder in the local repo.
# When pushing/pulling back to Github, the "data" folder in this repo will be ignored! 
# To modify the data used by this app, the Shiny app "data" object in S3 must be modified. Changes to the
# local "data" folder accessed by the Shiny scripts will NOT affect the web version of the app.

readData <- function(x) {
  
  # Gene-level expression matrices
  target_expData38_STAR_Dx <<- readRDS("data/mRNA/target_tpm_star_grch38_Dx_032425.RDS")
  target_expData38_STAR_Rem <<- readRDS("data/mRNA/target_tpm_star_grch38_Rem_032425.RDS")
  target_expData38_STAR_Rel <<- readRDS("data/mRNA/target_tpm_star_grch38_Rel_032425.RDS")
  
  target_expData38 <<- readRDS("data/mRNA/TARGET_RBD_Dx_AML_ExpressionData_TPM_filt4dupGenes_with_cellLines_CD34posNBM_DSAML_MPN_GRCh38_12.17.2020_FinalforShiny.RDS")
  target_expData37 <<- readRDS("data/mRNA/TARGET_RBD_Dx_AML_ExpressionData_TPM_filt4dupGenes_with_cellLines_CD34posNBM_DSAML_MPN_GRCh37_10.16.2020_FinalforShiny.RDS")
  beatAML_expData <<- readRDS("data/mRNA/BeatAML_Supplementary_Tables_TPM_Linear_Scale_v2.RDS")
  swog_expData <<- readRDS("data/mRNA/SWOG_AML_ExpressionData_TPM_GRCh38_FinalforShiny_v2.RDS")
  laml_expData <<- readRDS("data/mRNA/TCGA_LAML_ExpressionData_TPM_FinalforShiny_v2.RDS")
  stjude_expData <<- readRDS("data/mRNA/St_Jude_Expression_Data_TPM_filt4dupGenes_FinalforShiny_1.RDS")
  gmkf_expData <<- readRDS("data/mRNA/GMKF_TALL_TPM_Expression.RDS")
  ccle_expData <<- readRDS("data/mRNA/CCLE_TPM_Expression.RDS")
  leuce_expData <<- readRDS("data/mRNA/Leucegene_expression_data_TPM_v2.RDS")
  pcgp_aml_expData <<- readRDS("data/mRNA/stjude_aml_pcgp_v2.rds")
  pcgp_expData <<- readRDS("data/mRNA/pcgp_total_data.RDS")
  
  # miRNA expression matrices (for TARGET dataset only)
  load("data/miRNA/TARGET_AML_AAML1031_expn_matrix_mimat_norm_miRNA_RPM_01.07.2019_FinalforShiny.RData", .GlobalEnv)
  miRmapping <<- read.csv("data/miRNA/hsa_gff3_IDMap.csv")

  # Clinical data
  target_mani <<- read.csv("data/Clinical/complete_manifest_052825_exp4.csv")
  beatAML_cde <<- read.csv("data/Clinical/beatAML_cde_v3.csv")
  swog_cde <<- read.csv("data/Clinical/swog_cde_v3.csv")
  laml_cde <<- read.csv("data/Clinical/laml_cde_v3.csv")
  leuce_mani <<- read.csv("data/Clinical/Leucegene_manifest_v4.csv")
  pcgp_mani <<- read.csv("data/Clinical/stjude_aml_pcgp_manifest_v3.csv")
  pcgp_total_mani <<- read.csv("data/Clinical/pcgp_total_mani_v2.csv")
  
  #load("data/Clinical/Beat_AML_Supplementary_ClinicalData_FinalforShiny.RData", .GlobalEnv)
  load("data/Clinical/TARGET_AML_merged_CDEs_Shareable_FinalforShiny_Updated_12_17_24.RData", .GlobalEnv)
  #load("data/Clinical/SWOG_AML_Merged_CDEs_FinalforShiny.RData", .GlobalEnv)
  #load("data/Clinical/TCGA_LAML_ClinicalData_FinalforShiny.RData", .GlobalEnv)
  load("data/Clinical/StJude_ALL_ClinicalData_FinalforShiny.RData", .GlobalEnv)
  load("data/Clinical/GMKF_TALL_Clinical.RData", .GlobalEnv)
  load("data/Clinical/CCLE_Clinical_Data.RData", .GlobalEnv)
  

  
  # Misc accessory files
  load("data/ADC_and_CARTcell_Targets_Database_ADCReview_clinicaltrialsGov_FinalforShiny.RData", .GlobalEnv)
  load("data/DEGs/TARGET_AML_vs_NBM_and_Others_Ribodepleted_DEGs_per_Group_GRCh37_12.18.2020_FinalforShiny.RData", .GlobalEnv)
  deColKey <<- read.csv("data/Limma_Column_Descriptions.csv")
  colMapping <<- read.csv("data/Dataset_Column_Mapping_File.csv", check.names = F, na.strings = "")
  aml_restricted_genelist <<- read.csv("data/aml_restricted_genelist.csv")
  transmembrane_genelist <<- read.csv("data/transmembraneprot.csv")
  
  immdata <<- read.fst("data/hpa/rna_immune_cell_sample.fst") #the HPA data
  all_genes <<- readRDS("data/hpa/all_genes.RDS") #a list of genes from both of the datasets for autocorrection
  subloc <<- read.fst("data/hpa/subcellular_location.fst") 
  protein <<- read.fst("data/hpa/tissue_data.fst")
  
  # for the classification module
  # classification <<- read.csv("data/classification.csv")
  # km_cde <<- read.csv("data/km_updated_1_29_24.csv")

  # tcga_cancer <<- readRDS("data/concat_gtex_tcga_data.RDS")
  # tcga_csv <<- read.csv("data/listforcancers_bothlocat_10in90_final.csv")
  # gtex_tcga_combined <<- readRDS("data/concatenated_for_comparison_tcga_gtex.RDS")
  #tcga_newcsv <<- readRDS("data/tcga_expression_matrix.rds")
  tcga_newcsv <<- readRDS("data/tcga_concatenated_full.rds")
  tcga_manifest <<- read.csv("data/tcga_manifest.csv")
  
  gtex_csv <<- readRDS("data/gtex_concatenated_full.rds")
  gtex_manifest <<- read.csv("data/gtex_manifest.csv")
  gtex_manifest$Tissue <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", gtex_manifest$Tissue, perl=TRUE)
  
  
}

testing <- FALSE

if (identical(testing, TRUE)) {
  print("Testing mode - Data already in environment")
} else if (exists("target_expData38")) {
  print("Data already exists in the environment")
} else {
  readData()
}

####### Misc global functions & files
`%then%` <- function(a, b) {
  if (is.null(a)) b else a
}

bs <- 16 # Base font size for figures
dataset_choices <- list(
  aml = c("TARGET", "BEAT" = "BeatAML", "SWOG", "TGCA" = "TCGA", "LEUCEGENE", "PCGP AML"),
  all = c("St. Jude" = "StJude"),
  tall = c("GMKF" = "GMKF"),
  ccle = c("CCLE" = "CCLE"),
  pcgp = c("PCGP" = "PCGP")
)
