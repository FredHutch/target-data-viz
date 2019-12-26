setwd("/Volumes/homes/Shiny/waterfallPlot_app/data/")

# Loading "Relapse" batch RNA-seq data
load("TARGET_AML1031_Ribodepleted_Relapse_Counts_filteredForDupGenes_07.17.2019.RData")

# Identifying column names for cell lines
cellLines <- c("CSGH", "K562", "ME1", "MO7E", "NOMO1")
cellLines <- unlist(lapply(cellLines, function(x) grep(x, colnames(counts_relapse$TPMs), value = T)))

# Identifying column names for CD34+ NBM samples
cd34_NBM <- grep("R0|BM+.34", colnames(counts_relapse$TPMs), value = T)

# Casting to a data frame because it's currently a tibble, and subsetting it will drop all the rownames
relapse <- as.data.frame(counts_relapse$TPMs, drop = FALSE, stringsAsFactors = F)

# Subsetting the full relapse dataset to only contain the CD34+ NBMs and cell lines
relapse <- relapse[,c(cd34_NBM, cellLines), drop=FALSE]

# Making sure the rownames are still there
head(rownames(relapse))
head(rownames(counts_relapse$TPMs))

# Cleaning the colnames to remove aliquot and timepoint info
colnames(relapse) <- gsub("TARGET\\.[0-2]{2}\\.", "", colnames(relapse))
colnames(relapse) <- gsub("\\.[0-9]{2}A\\.01R", "", colnames(relapse))
colnames(relapse)

# Reading in ribodepleted sequencing batch 1 and 2 data
batch1 <- readRDS("TARGET_AAML1031_0531_Ribodepleted_TPMs_filtered_for_dupGenes_with_cellLines_11.20.2019.RDS")
batch1 <- as.data.frame(batch1, drop = F, stringsAsFactors = F)

# Removing the azacytidine experiment cell line samples
batch1 <- batch1[,!colnames(batch1) %in% grep("Kasumi.AZA|MV4.11.AZA", colnames(batch1), value = T)]

identical(rownames(relapse), rownames(batch1)) # Making sure these are the same so I can just bind the whole thing together into 1 dataframe
head(rownames(batch1))
head(rownames(relapse))

full_countsData <- cbind(batch1, relapse)

saveRDS(full_countsData, "TARGET_AAML1031_0531_Ribodepleted_Diagnostic_Relapse_TPMs_filtered_for_dupGenes_with_cellLines_CD34posNBM_11.20.2019.RDS")
