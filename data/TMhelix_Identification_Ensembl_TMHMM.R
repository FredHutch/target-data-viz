library("tidyverse")
library("biomaRt")

# Reading in the counts data used by the Shiny app for plotting
p1 <- readRDS("TARGET_AAML1031_0531_RBD_Dx_Relapse_TPMs_filt4dupGenes_with_cellLines_CD34posNBM_forShinyApp_PART1_12.27.2019.RDS")
p2 <- readRDS("TARGET_AAML1031_0531_RBD_Dx_Relapse_TPMs_filt4dupGenes_with_cellLines_CD34posNBM_forShinyApp_PART2_12.27.2019.RDS")
p3 <- readRDS("TARGET_AAML1031_0531_RBD_Dx_Relapse_TPMs_filt4dupGenes_with_cellLines_CD34posNBM_forShinyApp_PART3_12.27.2019.RDS")
p4 <- readRDS("TARGET_AAML1031_0531_RBD_Dx_Relapse_TPMs_filt4dupGenes_with_cellLines_CD34posNBM_forShinyApp_PART4_12.27.2019.RDS")
countsData <<- rbind(p1, p2, p3, p4)

ids2symbols_conversions <- read.csv("/Volumes/homes/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/analysis/0000.00.02_Reference_GeneInfo/GeneSymbol_Ensembl_ID_Conversion_GRCh37.69_FromBCCA.csv", header = TRUE, stringsAsFactors = F) %>%
  filter(geneSymbol %in% rownames(countsData))


# Uses vector of gene stable IDs to query Ensembl for further info about the gene & associated transcripts, proteins etc.
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "uswest.ensembl.org")
results <- getBM(values = ids2symbols_conversions$gene_id, filter = "ensembl_gene_id",
                 attributes = c("external_gene_name", "ensembl_gene_id", "tmhmm"),
                 mart = mart)

write.csv(results, "TMhelix_Prediction_Results_TMHMM_Ensembl_biomaRt_04.02.2020.csv", row.names = F)
