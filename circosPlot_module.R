#=======================================================================#
# Title: Circos Plots Module                                            #
# Created: 01-03-2025                                                   #  
# Last Update: 03-17-2025                                               #
#                                                                       #
# Authors:                                                              #
#   - Lena Allen, UO BGMP Student (lenarallen19@gmail.com)              #
#   - Lauren Williams, UO BGMP Student (laurenrwilliams.12@gmail.com)   #
#========================================================================


#################################################################
##                UI FUNCTION FOR CIRCOS MODULE                ##
#################################################################
circosPlotUI <- function(id, label = "Circos plot parameters"){
  
  #set namespace for module
  ns <- NS(id)
  
  #Create list of fusion types 
  all_fusions <- list("RUNX1-RUNX1T1", "MLLT10-DDX3X", "NUP98-NSD1", "DEUP1-EYA1", "CADM2-MGAM", "CSNK1G1-OR6J1", "DEUP1-RALYL", "GNPDA1-ABCC2", "PTPRR-COL6A6", "RUNX1-KCNB2", "RUNX1-LRRC69")
  names(all_fusions) <- list("RUNX1-RUNX1T1", "MLLT10-DDX3X", "NUP98-NSD1", "DEUP1-EYA1", "CADM2-MGAM", "CSNK1G1-OR6J1", "DEUP1-RALYL", "GNPDA1-ABCC2", "PTPRR-COL6A6", "RUNX1-KCNB2", "RUNX1-LRRC69")
  
  #Create list of chromosome choices
  chromosome_choices <- list("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11","chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21","chr22", "chrX", "chrY")
  names(chromosome_choices) <- list("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11","chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21","chr22", "chrX", "chrY")
  
  ##-------------------------------------------------------------##
  ##                          Page Setup                         ##
  ##-------------------------------------------------------------##
  
  tagList(
    useShinyjs(),
  fluidPage(
    theme = shinythemes::shinytheme(theme = "paper"),
    
    ##-------------------------------------------------------------##
    ##                            Sidebar                          ##
    ##-------------------------------------------------------------##
    
    # Placing the sidebar on the left hand side of the screen
    # sidebarLayout(
    #   position = "left",
    sidebarPanel(
      
      # Creating dropdown menu to select the fusion class for circos plot generation
      selectInput(ns("fusion_group"),
                  label = "Select a fusion to view circos plots of patients with relevant chromosomal rearrangements:",
                  choices = all_fusions),
      
      # Buttons to select whether circos plots generated will show all chromosomes 
      # or two chromosomes where specific translocation is occurring
      radioButtons(ns("all_chroms"), 
                   label = "View all chromosomes?", 
                   choices = list("All chromosomes" = "all", 
                                  "Two chromosomes" = "two", 
                                  "One chromosome" = "one")),
      
      # Conditional panel for if the user chooses to view only two chromosomes
      conditionalPanel(
        condition = paste0("input['", ns("all_chroms"), "'] == 'two'"),
        
        # Dropdown input for one of the selected chromosomes
        selectInput(ns("chromA"), label = "Select first chromosome:", choices = chromosome_choices),
        
        # Dropdown input for the second of the selected chromosomes -- GOING TO NEED TO MAKE SURE THIS CANNOT BE THE SAME AS CHROM_A
        selectInput(ns("chromB"), label = "Select second chromosome:", choices = chromosome_choices)),
      
      conditionalPanel(
        condition = paste0("input['", ns("all_chroms"), "'] == 'one'"),
        
        # Dropdown input for one of the selected chromosomes
        selectInput(ns("chromC"), label = "Select chromosome:", choices = chromosome_choices)),
      
      
      actionButton(ns("updateplot"), 
                   label = "Update Plot"),
      #Text on sidebar
      helpText(" "),
      helpText(" "),
      helpText(" "),
      helpText(" "),
      helpText("Translocations indicative of canonical pediatric AML fusions (NUP98-NSD1, RUNX1-RUNX1T1, and MLLT10-DDX3X) are depicted in blue."),
      helpText(" "),
      helpText(" "),
      helpText(" "),
      helpText(" "),
      helpText("KEY:"),
      helpText("Purple Histogram: Insertions"),
      helpText("Green Histogram: Deletions"),
      helpText("Green Boxes/Bars: Inversions"),
      helpText("Green Links: Duplications"),
      helpText("Red Links/Blue Links: Translocations"),
    ),
    
    ##-----------------------------------------------------------##
    ##                       Main Plot Panel                     ##
    ##-----------------------------------------------------------##
    
    
    #Creating space for plot to go!
    #Using SlickR to make carousel if there is more than one patients per fusion
    mainPanel(
      add_busy_spinner(spin = "fading-circle"),
      slickROutput(ns("plot"))
      
    )))
}


#################################################################
##              SERVER FUNCTION FOR CIRCOS MODULE              ##
#################################################################
circosPlot <- function(input, output, session){
  
  ##-------------------------------------------------------------##
  ##                          Data Prep                          ##
  ##-------------------------------------------------------------##
  
  #read in SV data RDS object
  svdata <- readRDS("./data/sv_circos_data.rds")
  
  #read in canonical fusion
  all_canonical_fusions <- read_delim("./data/full_canonical_fusions_list.tsv", delim = "\t")
  possible_new_fusions <- read_delim("./data/possible_new_fusions_list.tsv", delim = "\t")
  
  #since we're no longer working with VCF files, fix the names to match the patient_fusion_dt
  ptnames <- names(svdata)
  ptnames <- lapply(ptnames, function(x)strsplit(x, split = "_")[[1]])
  ptnames <- lapply(ptnames, function(l) l[[1]])
  names(svdata) <- ptnames
  ptnames
  
  ##-------------------------------------------------------------##
  ##               Find Canonical Fusions in Dataset             ##
  ##-------------------------------------------------------------##
  
  # making an empty data frame to store the results
  found_fusions <- data.frame(Patient = character(),
                                    Type = character(),
                                    Gene = character(),
                                    FusionMate = character(),
                                    VCF_CHROM = character(),
                                    VCF_POS = integer(),
                                    BND_MATE_CHROM = character(),
                                    BND_MATE_POS = integer(),
                                    stringsAsFactors = FALSE)
  
  # loop through each vcf in circos
  for (patient_name in names(svdata)) {
    # access the vcf for the current patient
    vcf_data <- svdata[[patient_name]]$vcf
    
    # grep for the BNDs (filter rows where ID contains "BND") idk if this needs to be as specific as the python code, if it does, can change
    vcf_data <- filter(vcf_data, grepl("BND", vcf_data$ID))
    
    # extract mate chromosome and position and add to the dataframe
    vcf_data <- vcf_data %>%
      mutate(
        BND_MATE_CHROM = str_extract(ALT, "\\w+(?=:)"),
        BND_MATE_POS = as.numeric(str_extract(ALT, "(?<=:)\\d+"))
      )
    
    # extract only the chrom, start pos, then the mate chrom and start pos
    vcf_data <- vcf_data[, c(1, 2, 9, 10)]
    
    # loop through each row in the vcf_data
    for (i in 1:nrow(vcf_data)) {
      # Extract the CHROM, POS, BND_MATE_CHROM, and BND_MATE_POS for the current row
      vcf_chrom <- vcf_data$CHROM[i]
      vcf_pos <- vcf_data$POS[i]
      bnd_mate_chrom <- vcf_data$BND_MATE_CHROM[i]
      bnd_mate_pos <- vcf_data$BND_MATE_POS[i]
      
      # find matching rows in canonical_fusions based on chrom and pos
      first_fusion <- all_canonical_fusions %>%
        filter(Chromosome == vcf_chrom,
               Start <= vcf_pos & vcf_pos <= End) # Check if pos is within Start and End
      
      # do it again with the mate chrom and pos
      second_fusion <- all_canonical_fusions %>%
        filter(Chromosome == bnd_mate_chrom,
               Start <= bnd_mate_pos & bnd_mate_pos <= End) # Check if BND_MATE_POS is within Start and End
      
      # find matching rows in canonical_fusions based on chrom and pos
      first_newfusion <- possible_new_fusions %>%
        filter(Chromosome == vcf_chrom,
               Start <= vcf_pos & vcf_pos <= End) # Check if pos is within Start and End
      
      # do it again with the mate chrom and pos
      second_newfusion <- possible_new_fusions %>%
        filter(Chromosome == bnd_mate_chrom,
               Start <= bnd_mate_pos & bnd_mate_pos <= End) # Check if BND_MATE_POS is within Start and End
      
      # if both first_fusion partner and second_fusion partner are found
      if (nrow(first_fusion) > 0 && nrow(second_fusion) > 0) {
        # get the Gene from the first match
        gene <- first_fusion$Gene[1]
        mate <- second_fusion$Gene[1]
        
        # store the result and add the patient file name
        found_fusions <- rbind(found_fusions, data.frame(
          Patient = patient_name, # save the patient name (filename)
          Type = "Canonical",
          Gene = gene,
          FusionMate = mate,
          VCF_CHROM = vcf_chrom,
          VCF_POS = vcf_pos,
          BND_MATE_CHROM = bnd_mate_chrom,
          BND_MATE_POS = bnd_mate_pos
        ))
      } else if (nrow(first_newfusion) > 0 && nrow(second_newfusion) > 0) {
        # get the Gene from the first match
        gene <- first_newfusion$Gene[1]
        mate <- second_newfusion$Gene[1]
        
        # store the result and add the patient file name
        found_fusions <- rbind(found_fusions, data.frame(
          Patient = patient_name, # save the patient name (filename)
          Type = "New",
          Gene = gene,
          FusionMate = mate,
          VCF_CHROM = vcf_chrom,
          VCF_POS = vcf_pos,
          BND_MATE_CHROM = bnd_mate_chrom,
          BND_MATE_POS = bnd_mate_pos
        ))
      }
    }
  }
  
  found_canon_fusions <- found_fusions %>% subset(Type == "Canonical")
  ptlist <- unique(found_fusions$Patient)
  
  ##-------------------------------------------------------------##
  ##                            CIRCOS                           ##
  ##-------------------------------------------------------------##
  
  #reactive expression to get list of patients with the user-selected fusion
  cp_objs <- eventReactive(input$updateplot, {
    circos.clear()
    fusion <- input$fusion_group
    circostype <- input$all_chroms
    fusion <- fusion[1]
    pts_with_fusion <- c()
    
    #pull out 
    fusion_genes <- lapply(fusion, function(x)strsplit(x, split = "-")[[1]])
    fusion_gene <- fusion_genes[[1]][1]
    fusion_gene_mate <- fusion_genes[[1]][2]
    
    for (i in 1:nrow(found_fusions)) { 
      found_gene <- found_fusions$Gene[i]
      found_gene_mate <- found_fusions$FusionMate[i]
      pat <- found_fusions$Patient[i]

      if (found_gene == fusion_gene & found_gene_mate == fusion_gene_mate){
        pts_with_fusion <- append(pts_with_fusion, pat)
      } else if (found_gene == fusion_gene_mate & found_gene_mate == fusion_gene){
        pts_with_fusion <- append(pts_with_fusion, pat)
      }
    }
    
    pts_with_fusion <- unique(pts_with_fusion)
    
    ##-----------------------------------------------------------------------##
    ##                   User Chooses to View ALL Chromosomes                ##
    ##-----------------------------------------------------------------------##
    if (length(pts_with_fusion) > 0 & circostype == "all" ) {
      circos.clear()
      length_pts <- length(pts_with_fusion)
      for(i in 1:length(pts_with_fusion)) {
        
        pat <- pts_with_fusion[i]
        
        #get vcf data for patient
        vcf <- svdata[pat][[1]]
        
        #convert to BED
        vcf_bed <- suppressWarnings(vcf2bed(vcf, other = c("ID", "INFO", "ALT")))
        vcf_bed <- vcf_bed %>% filter(chr!="chrM") %>% mutate(value = 1)
        vcf_bed <- data.frame(vcf_bed)
        
        ##--------------------------------------------------##
        ##                   Initialize Files               ##
        ##--------------------------------------------------##
        
        #initializing BED files for each variant
        bed_INS <- vcf_bed %>% filter(grepl("INS", ID)) %>% distinct(chr, start, .keep_all = TRUE)
        bed_DEL <- vcf_bed %>% filter(grepl("DEL", ID)) %>% distinct(chr, start, .keep_all = TRUE)
        bed_INV <- vcf_bed %>% filter(grepl("INV", ID)) %>% distinct(chr, start, .keep_all = TRUE)
        bed_DUP <- vcf_bed %>% filter(grepl("DUP", ID)) %>% distinct(chr, start, .keep_all = TRUE)
        bed_BND <- vcf_bed %>% filter(grepl("BND", ID)) %>% distinct(chr, start, .keep_all = TRUE)
        
        #initializing bed_BND_MATE
        mates <- bed_BND$ALT
        mates <- strsplit(mates, ":")
        chrom_mate <- sapply(mates, `[[`,1)
        chrom_mate <- str_match(chrom_mate, "chr\\w+")
        pos_mate <- sapply(mates, `[[`, 2)
        pos_mate <- str_match(pos_mate, "\\d+")
        pos_mate <- as.numeric(pos_mate)
        
        chr <- c(chrom_mate)
        start <- c(pos_mate)
        end <- c(pos_mate)
        value <- c(1)
        bed_BND_MATES <- data.frame(chr, start, end, value)
        
        #fixing bed_DEL end values
        del_end <- data.frame(matrix(nrow=0,ncol=1))
        del_end <- del_end %>% rename(end = matrix.nrow...0..ncol...1.)
        info <- bed_DEL$INFO
        info <- strsplit(info, "END=")
        NEWEND <- sapply(info, `[[`, 2)
        NEWEND <- str_match(NEWEND, "\\d+")
        NEWEND <- as.numeric(NEWEND)
        
        bed_DEL <- bed_DEL %>% mutate(NEWEND = NEWEND)
        bed_DEL <- bed_DEL %>% rename(OLDend = end)
        bed_DEL <- bed_DEL %>% mutate(OLDend = NULL)
        bed_DEL <- bed_DEL %>% rename(end = NEWEND)
        bed_DEL <- bed_DEL[, c("chr", "start", "end", "ID", "INFO", "value")]
        
        #fixing bed_INV end values
        inv_end <- data.frame(matrix(nrow=0,ncol=1))
        inv_end <- inv_end %>% rename(end = matrix.nrow...0..ncol...1.)
        info <- bed_INV$INFO
        info <- strsplit(info, "END=")
        NEWEND <- sapply(info, `[[`, 2)
        NEWEND <- str_match(NEWEND, "\\d+")
        NEWEND <- as.numeric(NEWEND)
        
        bed_INV <- bed_INV %>% mutate(NEWEND = NEWEND)
        bed_INV <- bed_INV %>% rename(OLDend = end)
        bed_INV <- bed_INV %>% mutate(OLDend = NULL)
        bed_INV <- bed_INV %>% rename(end = NEWEND)
        bed_INV <- bed_INV[, c("chr", "start", "end", "ID", "INFO", "value")]
        
        #fixing bed_DUP end values
        dup_end <- data.frame(matrix(nrow=0,ncol=1))
        dup_end <- dup_end %>% rename(end = matrix.nrow...0..ncol...1.)
        info <- bed_DUP$INFO
        info <- strsplit(info, "END=")
        NEWEND <- sapply(info, `[[`, 2)
        NEWEND <- str_match(NEWEND, "\\d+")
        NEWEND <- as.numeric(NEWEND)
        
        bed_DUP <- bed_DUP %>% mutate(NEWEND = NEWEND)
        bed_DUP <- bed_DUP %>% rename(OLDend = end)
        bed_DUP <- bed_DUP %>% mutate(OLDend = NULL)
        bed_DUP <- bed_DUP %>% rename(end = NEWEND)
        bed_DUP <- bed_DUP[, c("chr", "start", "end", "ID", "INFO", "value")]
        
        #marking canonical fusions
        if (pat %in% found_canon_fusions$Patient){
          
          #subset just the rows of the patient
          patient_canon_fusions <- found_canon_fusions %>% subset(Patient == pat) 
          
          #make a dataframe named bed_BND_CANON with VCF_CHROM, VCF_POS, VCF_POS, Gene, value (1) 
          bed_BND_CANON <- data.frame(
            chr = patient_canon_fusions$VCF_CHROM, 
            pos = patient_canon_fusions$VCF_POS, 
            pos = patient_canon_fusions$VCF_POS, 
            gene = patient_canon_fusions$Gene, 
            value = 1)
          
          #make a dataframe named bed_BND_MATES_CANON with BND_MATE_CHROM, BND_MATE_POS, BND_MATE_POS, FusionMate, value (1)
          bed_BND_MATES_CANON <- data.frame(
            chr_mate = patient_canon_fusions$BND_MATE_CHROM, 
            pos_mate = patient_canon_fusions$BND_MATE_POS, 
            pos_mate = patient_canon_fusions$BND_MATE_POS, 
            gene_mate = patient_canon_fusions$FusionMate, 
            value = 1)
        }
        
        ##--------------------------------------------------##
        ##                   Make CIRCOS Plot               ##
        ##--------------------------------------------------##
        # initialize with ideogram (all chromosomes)
        circos.clear()
        png(paste(i, ".png", sep=""), units="px", width=800, height=800, res=160)
        circos.par(track.height = 0.08, start.degree=90)
        circos.initializeWithIdeogram(species = "hg38")
        
        # INSERTION TRACK
        circos.trackHist(bed_INS$chr, x = bed_INS$start, bin.size = 1000000, col = "purple", border = "purple")
        
        # DELETION TRACK
        circos.trackHist(bed_DEL$chr, x = bed_DEL$start, bin.size = 1000000, col = "green3", border = "green3")
        
        # INVERSION TRACK
        circos.genomicTrack(bed_INV, ylim = c(0,1), panel.fun = function(region, value, ...) {circos.genomicRect(region, value, col = "green3", border = "green3")})
        
        # DUPLICATION TRACK
        circos.genomicLink(bed_DUP, bed_DUP, col = "green3", border = "green3")
        
        # TRANSLOCATION TRACK
        circos.genomicLink(bed_BND, bed_BND_MATES, col = "red2", border = "red2")
        
        # CANONICAL TRANSLOCATIONS TRACK
        if (pat %in% found_fusions$Patient){
          circos.genomicLink(bed_BND_CANON, bed_BND_MATES_CANON, col = "blue", border = "blue")
        }
        dev.off()
        
      } 
      
      rl <- sprintf("%i.png", 1:length_pts)
      rl} 
    
    ##-----------------------------------------------------------------------##
    ##                   User Chooses to View TWO Chromosomes                ##
    ##-----------------------------------------------------------------------##
    else if (length(pts_with_fusion) > 0 & circostype == "two" ) {
      circos.clear()
      chra <- input$chromA
      chrb <- input$chromB
      
      if (chra == chrb) {
        validate("Please select two different chromosomes.")
        NULL}
      
      else {
      length_pts <- length(pts_with_fusion)
      for(i in 1:length(pts_with_fusion)) {
        
        pat <- pts_with_fusion[i]
        
        
        #get vcf data for patient
        vcf <- svdata[pat][[1]]
        
        #convert to BED
        vcf_bed <- suppressWarnings(vcf2bed(vcf, other = c("ID", "INFO", "ALT")))
        vcf_bed <- vcf_bed %>% filter(chr!="chrM") %>% mutate(value = 1)
        vcf_bed <- data.frame(vcf_bed)
        
        ##--------------------------------------------------##
        ##                   Initialize Files               ##
        ##--------------------------------------------------##
        
        #initializing BED files for each variant
        bed_INS <- vcf_bed %>% filter(grepl("INS", ID)) %>% distinct(chr, start, .keep_all = TRUE)
        bed_DEL <- vcf_bed %>% filter(grepl("DEL", ID)) %>% distinct(chr, start, .keep_all = TRUE)
        bed_INV <- vcf_bed %>% filter(grepl("INV", ID)) %>% distinct(chr, start, .keep_all = TRUE)
        bed_DUP <- vcf_bed %>% filter(grepl("DUP", ID)) %>% distinct(chr, start, .keep_all = TRUE)
        bed_BND <- vcf_bed %>% filter(grepl("BND", ID)) %>% distinct(chr, start, .keep_all = TRUE) %>% filter(grepl(paste(chra, "$|", chrb, "$", sep =""), chr)) %>% filter(grepl(paste(chra, ":|", chrb, ":", sep =""), ALT))
        bed_BND_CANON <- NULL
        bed_BND_MATES_CANON <- NULL
        
        if(length(bed_BND$chr) != 0) {
        #initializing bed_BND_MATE
        mates <- bed_BND$ALT
        mates <- strsplit(mates, ":")
        chrom_mate <- sapply(mates, `[[`,1)
        chrom_mate <- str_match(chrom_mate, "chr\\w+")
        pos_mate <- sapply(mates, `[[`, 2)
        pos_mate <- str_match(pos_mate, "\\d+")
        pos_mate <- as.numeric(pos_mate)

        chr <- c(chrom_mate)
        start <- c(pos_mate)
        end <- c(pos_mate)
        value <- c(1)
        bed_BND_MATES <- data.frame(chr, start, end, value)}

        #fixing bed_DEL end values
        del_end <- data.frame(matrix(nrow=0,ncol=1))
        del_end <- del_end %>% rename(end = matrix.nrow...0..ncol...1.)
        info <- bed_DEL$INFO
        info <- strsplit(info, "END=")
        NEWEND <- sapply(info, `[[`, 2)
        NEWEND <- str_match(NEWEND, "\\d+")
        NEWEND <- as.numeric(NEWEND)

        bed_DEL <- bed_DEL %>% mutate(NEWEND = NEWEND)
        bed_DEL <- bed_DEL %>% rename(OLDend = end)
        bed_DEL <- bed_DEL %>% mutate(OLDend = NULL)
        bed_DEL <- bed_DEL %>% rename(end = NEWEND)
        bed_DEL <- bed_DEL[, c("chr", "start", "end", "ID", "INFO", "value")]

        #fixing bed_INV end values
        inv_end <- data.frame(matrix(nrow=0,ncol=1))
        inv_end <- inv_end %>% rename(end = matrix.nrow...0..ncol...1.)
        info <- bed_INV$INFO
        info <- strsplit(info, "END=")
        NEWEND <- sapply(info, `[[`, 2)
        NEWEND <- str_match(NEWEND, "\\d+")
        NEWEND <- as.numeric(NEWEND)

        bed_INV <- bed_INV %>% mutate(NEWEND = NEWEND)
        bed_INV <- bed_INV %>% rename(OLDend = end)
        bed_INV <- bed_INV %>% mutate(OLDend = NULL)
        bed_INV <- bed_INV %>% rename(end = NEWEND)
        bed_INV <- bed_INV[, c("chr", "start", "end", "ID", "INFO", "value")]

        #fixing bed_DUP end values
        dup_end <- data.frame(matrix(nrow=0,ncol=1))
        dup_end <- dup_end %>% rename(end = matrix.nrow...0..ncol...1.)
        info <- bed_DUP$INFO
        info <- strsplit(info, "END=")
        NEWEND <- sapply(info, `[[`, 2)
        NEWEND <- str_match(NEWEND, "\\d+")
        NEWEND <- as.numeric(NEWEND)

        bed_DUP <- bed_DUP %>% mutate(NEWEND = NEWEND)
        bed_DUP <- bed_DUP %>% rename(OLDend = end)
        bed_DUP <- bed_DUP %>% mutate(OLDend = NULL)
        bed_DUP <- bed_DUP %>% rename(end = NEWEND)
        bed_DUP <- bed_DUP[, c("chr", "start", "end", "ID", "INFO", "value")]

        #marking canonical fusions
        if (pat %in% found_canon_fusions$Patient) {
          
          #subset just the rows of the patient
          patient_canon_fusions <- found_canon_fusions %>% subset(Patient == pat)
          
          #make a dataframe named bed_BND_CANON with VCF_CHROM, VCF_POS, VCF_POS, Gene, value (1) 
          bed_BND_CANON <- data.frame(
            chr = patient_canon_fusions$VCF_CHROM,
            chr_mate = patient_canon_fusions$BND_MATE_CHROM,
            pos = patient_canon_fusions$VCF_POS, 
            pos = patient_canon_fusions$VCF_POS, 
            gene = patient_canon_fusions$Gene, 
            value = 1)
          
          #make a dataframe named bed_BND_MATES_CANON with BND_MATE_CHROM, BND_MATE_POS, BND_MATE_POS, FusionMate, value (1)
          bed_BND_MATES_CANON <- data.frame(
            chr_mate = patient_canon_fusions$BND_MATE_CHROM, 
            chr = patient_canon_fusions$VCF_CHROM,
            pos_mate = patient_canon_fusions$BND_MATE_POS, 
            pos_mate = patient_canon_fusions$BND_MATE_POS, 
            gene_mate = patient_canon_fusions$FusionMate, 
            value = 1)
        }

        ##--------------------------------------------------##
        ##                  Make CIRCOS Plot                ##
        ##--------------------------------------------------##
        # initialize with ideogram (two chromosomes)
        circos.clear()
        png(paste(i, ".png", sep=""), units="px", width=800, height=800, res=160)
        circos.par(track.height = 0.08, start.degree =90)
        circos.initializeWithIdeogram(species = "hg38", chromosome.index = c(chra, chrb))

        #initialize with subset of chromosomes
        bed_INS_two <- bed_INS %>% filter(grepl(paste(chra, "$|", chrb, "$", sep =""), chr))
        bed_DEL_two <- bed_DEL %>% filter(grepl(paste(chra, "$|", chrb, "$", sep =""), chr))
        bed_INV_two <- bed_INV %>% filter(grepl(paste(chra, "$|", chrb, "$", sep =""), chr))
        bed_DUP_two <- bed_DUP %>% filter(grepl(paste(chra, "$|", chrb, "$", sep =""), chr))
        
        if(length(bed_BND$chr) != 0){
        bed_BND_two <- bed_BND
        bed_BND_MATES_two <- bed_BND_MATES
        
        if(is.null(bed_BND_CANON) != TRUE & is.null(bed_BND_MATES_CANON) != TRUE) {

        bed_BND_CANON_two <- bed_BND_CANON %>% filter(chr == chra & chr_mate == chrb | chr == chrb & chr_mate == chra) %>% dplyr::select(-chr_mate)
        bed_BND_MATES_CANON_two <- bed_BND_MATES_CANON %>% filter(chr == chra & chr_mate == chrb | chr == chrb & chr_mate == chra) %>% dplyr::select(-chr)}}

        # INSERTION TRACK
        circos.trackHist(bed_INS_two$chr, x = bed_INS_two$start, bin.size = 1000000, col = "purple", border = "purple")

        # DELETION TRACK
        circos.trackHist(bed_DEL_two$chr, x = bed_DEL_two$start, bin.size = 1000000, col = "green3", border = "green3")

        # INVERSION TRACK
        circos.genomicTrack(bed_INV_two, ylim = c(0,1), panel.fun = function(region, value, ...) {circos.genomicRect(region, value, col = "green3", border = "green3")})

        # DUPLICATION TRACK
        circos.genomicLink(bed_DUP_two, bed_DUP_two, col = "green3", border = "green3")

        # TRANSLOCATION TRACK ------ this can sometimes throw an error, but sometimes not
        if(length(bed_BND$chr) != 0){
        circos.genomicLink(bed_BND_two, bed_BND_MATES_two, col = "red2", border = "red2")

          if (pat %in% found_fusions$Patient){
          circos.genomicLink(bed_BND_CANON_two, bed_BND_MATES_CANON_two, col = "blue", border = "blue")
        }}

        # #CANONICAL TRANSLOCATIONS --------- this can also sometimes throw an error, but sometimes not
        # if (pat %in% found_canon_fusions$Patient){
        # circos.genomicLink(bed_BND_CANON_two, bed_BND_MATES_CANON_two, col = "blue", border = "blue")
        # }

        dev.off()
      }
      
      al <- sprintf("%i.png", 1:length_pts)
      al
      }}

    ##-------------------------------------------------------------------------
    ##                   User Chooses to View One Chromosome               --
    ##-------------------------------------------------------------------------
    else if (length(pts_with_fusion) > 0 & circostype == "one" ) {
      circos.clear()
      chrc <- input$chromC #fix this
      
      length_pts <- length(pts_with_fusion)
      for(i in 1:length(pts_with_fusion)) {
        
        pat <- pts_with_fusion[i]
        
        #get vcf data for patient
        vcf <- svdata[pat][[1]]
        
        #convert to BED
        vcf_bed <- suppressWarnings(vcf2bed(vcf, other = c("ID", "INFO", "ALT")))
        vcf_bed <- vcf_bed %>% filter(chr!="chrM") %>% mutate(value = 1)
        vcf_bed <- data.frame(vcf_bed)
        
        ##--------------------------------------------------##
        ##                  Initialize Files                ##
        ##--------------------------------------------------##
        
        #initializing BED files for each variant
        bed_INS <- vcf_bed %>% filter(grepl("INS", ID)) %>% distinct(chr, start, .keep_all = TRUE)
        bed_DEL <- vcf_bed %>% filter(grepl("DEL", ID)) %>% distinct(chr, start, .keep_all = TRUE)
        bed_INV <- vcf_bed %>% filter(grepl("INV", ID)) %>% distinct(chr, start, .keep_all = TRUE)
        bed_DUP <- vcf_bed %>% filter(grepl("DUP", ID)) %>% distinct(chr, start, .keep_all = TRUE)

        #fixing bed_DEL end values
        del_end <- data.frame(matrix(nrow=0,ncol=1))
        del_end <- del_end %>% rename(end = matrix.nrow...0..ncol...1.)
        info <- bed_DEL$INFO
        info <- strsplit(info, "END=")
        NEWEND <- sapply(info, `[[`, 2)
        NEWEND <- str_match(NEWEND, "\\d+")
        NEWEND <- as.numeric(NEWEND)
        
        bed_DEL <- bed_DEL %>% mutate(NEWEND = NEWEND)
        bed_DEL <- bed_DEL %>% rename(OLDend = end)
        bed_DEL <- bed_DEL %>% mutate(OLDend = NULL)
        bed_DEL <- bed_DEL %>% rename(end = NEWEND)
        bed_DEL <- bed_DEL[, c("chr", "start", "end", "ID", "INFO", "value")]
        
        #fixing bed_INV end values
        inv_end <- data.frame(matrix(nrow=0,ncol=1))
        inv_end <- inv_end %>% rename(end = matrix.nrow...0..ncol...1.)
        info <- bed_INV$INFO
        info <- strsplit(info, "END=")
        NEWEND <- sapply(info, `[[`, 2)
        NEWEND <- str_match(NEWEND, "\\d+")
        NEWEND <- as.numeric(NEWEND)
        
        bed_INV <- bed_INV %>% mutate(NEWEND = NEWEND)
        bed_INV <- bed_INV %>% rename(OLDend = end)
        bed_INV <- bed_INV %>% mutate(OLDend = NULL)
        bed_INV <- bed_INV %>% rename(end = NEWEND)
        bed_INV <- bed_INV[, c("chr", "start", "end", "ID", "INFO", "value")]
        
        #fixing bed_DUP end values
        dup_end <- data.frame(matrix(nrow=0,ncol=1))
        dup_end <- dup_end %>% rename(end = matrix.nrow...0..ncol...1.)
        info <- bed_DUP$INFO
        info <- strsplit(info, "END=")
        NEWEND <- sapply(info, `[[`, 2)
        NEWEND <- str_match(NEWEND, "\\d+")
        NEWEND <- as.numeric(NEWEND)
        
        bed_DUP <- bed_DUP %>% mutate(NEWEND = NEWEND)
        bed_DUP <- bed_DUP %>% rename(OLDend = end)
        bed_DUP <- bed_DUP %>% mutate(OLDend = NULL)
        bed_DUP <- bed_DUP %>% rename(end = NEWEND)
        bed_DUP <- bed_DUP[, c("chr", "start", "end", "ID", "INFO", "value")]
        
        
        ##--------------------------------------------------##
        ##                  Make CIRCOS Plot                ##
        ##--------------------------------------------------##
        # Initialize with ideogram (one chromosome)
        circos.clear()
        png(paste(i, ".png", sep=""), units="px", width=800, height=800, res=160)
        circos.par(track.height = 0.08, start.degree =90)
        circos.initializeWithIdeogram(species = "hg38", chromosome.index = c(chrc))
        
        # Initialize with one chromosome
        bed_INS_one <- bed_INS %>% filter(grepl(paste(chrc, "$", sep=""), chr))
        bed_DEL_one <- bed_DEL %>% filter(grepl(paste(chrc, "$", sep=""), chr))
        bed_INV_one <- bed_INV %>% filter(grepl(paste(chrc, "$", sep=""), chr))
        bed_DUP_one <- bed_DUP %>% filter(grepl(paste(chrc, "$", sep=""), chr))
        
        # INSERTION TRACK
        circos.trackHist(bed_INS_one$chr, x = bed_INS_one$start, bin.size = 1000000, col = "purple", border = "purple")
        
        # DELETION TRACK
        circos.trackHist(bed_DEL_one$chr, x = bed_DEL_one$start, bin.size = 1000000, col = "green3", border = "green3")
        
        # INVERSION TRACK
        circos.genomicTrack(bed_INV_one, ylim = c(0,1), panel.fun = function(region, value, ...) {circos.genomicRect(region, value, col = "green3", border = "green3")})
        
        # DUPLICATION TRACK
        circos.genomicLink(bed_DUP_one, bed_DUP_one, col = "green3", border = "green3")
        
        # NO TRANSLOCATION TRACK FOR ONE CHROMOSOME VISUALIZATION
        # NO CANONICAL TRANSLOCATIONS FOR ONE CHROMOSOME VISUALIZATION

        dev.off()
      }
      
      pl <- sprintf("%i.png", 1:length_pts)
      pl
    }
    
  })
  
  ##-------------------------------------------------------------##
  ##                       RENDER CIRCOS PLOTS                   ##
  ##-------------------------------------------------------------##  
  output$plot <- renderSlickR({
    
    pts <- cp_objs()
    if(is.null(pts) != TRUE) {
      slickR(pts)
    } 
  })
}

