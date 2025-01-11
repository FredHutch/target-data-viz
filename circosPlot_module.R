#################################################################
##                UI FUNCTION FOR CIRCOS MODULE                ##
#################################################################
circosPlotUI <- function(id, label = "Circos plot parameters"){
  
  #set namespace for module
  ns <- NS(id)
  
  #Create list of fusion types 
  canonical_fusions <- list("CBFA2T3-GLIS2", "CBFB-MYH11", "DEK-NUP214", "KAT6A-CREBBP", "RBM15-MKL1", "RUNX1-RUNX1T1", "NPM1-MLF1", "ERG-X", "ETV6-X", "FEV-X", "FLI1-X", "KMT2A-X", "MECOM-X", "MLLT10-X", "NUP98-X")
  names(canonical_fusions) <- list("CBFA2T3-GLIS2", "CBFB-MYH11", "DEK-NUP214", "KAT6A-CREBBP", "RBM15-MKL1", "RUNX1-RUNX1T1", "NPM1-MLF1", "ERG-X", "ETV6-X", "FEV-X", "FLI1-X", "KMT2A-X", "MECOM-X", "MLLT10-X", "NUP98-X")
  
  #Create list of chromosome choices
  chromosome_choices <- list("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11","chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21","chr22", "chrX", "chrY", "chrM")
  names(chromosome_choices) <- list("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11","chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21","chr22", "chrX", "chrY", "chrM")
  
  ###################################################################
  ##                          Page Setup                           ##
  ###################################################################

  tagList(
  useShinyjs(),
  fluidPage(
  theme = shinythemes::shinytheme(theme = "paper"),


  ###################################################################
  ##                          Sidebar                              ##
  ###################################################################
  
  # Placing the sidebar on the left hand side of the screen
  sidebarLayout(
    position = "left",
    sidebarPanel(
      
      # Creating dropdown menu to select the fusion class for circos plot generation
      selectInput(ns("fusion_group"),
                  label = "Select a fusion type",
                  choices = canonical_fusions),
      
      # Buttons to select whether circos plots generated will show all chromosomes 
      # or two chromosomes where specific translocation is occurring
      radioButtons(ns("all_chroms"), 
                   label = "View all chromosomes?", 
                   choices = list("All chromosomes" = "all", 
                                  "Two chromosomes" = "two"))),
    
      # Temporary fix for not being able to hover over circos plots and get more info
      # Narrowing down to two chromosomes has similar impact

      # Conditional panel for if the user chooses to view only two chromosomes
      conditionalPanel(
        condition = paste0("input['", ns("all_chroms"), "'] == 'two'"),
        
        # Dropdown input for one of the selected chromosomes
        selectInput(ns("chromA"), label = "Select first chromosome:", choices = chromosome_choices),
        
        # Dropdown input for the second of the selected chromosomes -- GOING TO NEED TO MAKE SURE THIS CANNOT BE THE SAME AS CHROM_A
        selectInput(ns("chromB"), label = "Select second chromosome:", choices = chromosome_choices))

        ),
  
        # Txt on sidebar
        helpText("Circos plots will be generated to show canonical fusions in blue."),
  
      # Creating download button so Circos plots can be downloaded locally
      downloadButton(ns("plot_download"), 
                 label = "plot", 
                 class = "plotdwnld"),
  
      shinyBS::bsTooltip(ns("plot_download"), 
                     title = "Click here to download a copy of the plot",
                     placement = "right", 
                     trigger = "hover"),

  
  ###############################################################
  #----------------------- MAIN PLOT PANEL ---------------------#
  ###############################################################


    #Creating space for plot to go!
    #Not sure how we will go about showing multiple circos plots at once - patchwork-style?
    mainPanel(
      position = "right",
      tabsetPanel(

        #Circos plot tab panel -  TITLE OF TAB PANEL, NOT THE PLOT OBJECT!
        tabPanel("Plot",
                 # Linebreaks to help center the plot on the page -
                 # they used this in some of the existing modules,
                 # unsure if we'll need it yet but keeping it for now :)
                 br(),
                 br(),

                 # not sure what fluidRow() is doing here
                 # this is taken from their waterfall plot module, but seems necessary
                 # their comments say:
                 # "This will be a reactive object that is linked to an item in the
                 # output list, created in the "server" script"
                 fluidRow(
                   column(12, offset = 0, align = "left",
                          plotOutput(ns("plot"), width = "600px")
                   )))))))
  }


#################################################################
##              SERVER FUNCTION FOR CIRCOS MODULE              ##
#################################################################
circosPlot <- function(input, output, session){
  
  
  ##---------------------------------------------------------------
  ##                          Data Prep                          --
  ##---------------------------------------------------------------

  # Set up dropdown choices - same as above?
  # canonical_fusions <- as.list(c("CBFA2T3-GLIS2", "CBFB-MYH11", "DEK-NUP214", "KAT6A-CREBBP", "RBM15-MKL1", "RUNX1-RUNX1T1", "NPM1-MLF1", "ERG-X", "ETV6-X", "FEV-X", "FLI1-X", "KMT2A-X", "MECOM-X", "MLLT10-X", "NUP98-X"))
  # names(canonical_fusions) <- c("CBFA2T3-GLIS2", "CBFB-MYH11", "DEK-NUP214", "KAT6A-CREBBP", "RBM15-MKL1", "RUNX1-RUNX1T1", "NPM1-MLF1", "ERG-X", "ETV6-X", "FEV-X", "FLI1-X", "KMT2A-X", "MECOM-X", "MLLT10-X", "NUP98-X")
  
  #prepare dataframe of patient ids and which fusions are observed in each
  #columns are fusion types, rownames are patient IDs
  patientids <- read.table("./data/ptlist.txt")
  patientids <- patientids$V1
  fusionnames <- c("CBFA2T3-GLIS2", "CBFB-MYH11", "DEK-NUP214", "KAT6A-CREBBP", "RBM15-MKL1", "RUNX1-RUNX1T1", "NPM1-MLF1", "ERG-X", "ETV6-X", "FEV-X", "FLI1-X", "KMT2A-X", "MECOM-X", "MLLT10-X", "NUP98-X")
  
  patient_fusion_dt <- data.frame(matrix(ncol = 15, nrow = 61))
  colnames(patient_fusion_dt) <- fusionnames
  rownames(patient_fusion_dt) <- patientids
  
  patient_fusion_dt["TARGET-20-PAURDN-03A-01D", "NUP98-X"] <- "yes"
  patient_fusion_dt["TARGET-20-PAUPIY-03A-01D", "NUP98-X"] <- "yes"
  patient_fusion_dt["TARGET-20-PAVBIH-09A-02D", "MLLT10-X"] <- "yes"
  patient_fusion_dt["TARGET-20-PAUZRY-09A-02D", "RUNX1-RUNX1T1"] <- "yes"
  patient_fusion_dt["TARGET-20-PAUNVN-09A-01D", "ETV6-X"] <- "yes"
  
  #read in SV data RDS object
  svdata <- readRDS("./data/sv_circos_data.rds")
  
  #since we're no longer working with VCF files, fix the names to match the patient_fusion_dt
  ptnames <- names(svdata)
  ptnames <- lapply(ptnames, function(x)strsplit(x, split = "_")[[1]])
  ptnames = lapply(ptnames, function(l) l[[1]])
  names(svdata) <- ptnames
  
  
  ##---------------------------------------------------------------
  ##                            CIRCOS                           --
  ##---------------------------------------------------------------
  
  #subset dataframe to patients with the fusion selected in the dropdown
  pts_with_fusion <- reactive(patient_fusion_dt %>% filter(input$fusion_group == "yes") %>% rownames())
  
  circos_objects <- c()

  #loop through all patients with fusion,
  for (pt in pts_with_fusion) {

    #get vcf data for that patient
    vcf = svdata$pt
    
    #convert to BED
    vcf_bed <- vcf2bed(vcf, other = c("ID", "INFO", "ALT"))
    vcf_bed <- vcf_bed %>% filter(chr!="chrM") %>% mutate(value = 1) 
    vcf_bed <- data.frame(vcf_bed)
    
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
    
    #MARKING CANONICAL FUSIONS
    #TARGET-20-PAVBIH-09A-02D_consensus -> MLLT10 gene on chromosome 10 with chromosome X
    if (patient_name == "TARGET-20-PAVBIH-09A-02D"){
      bed_BND_CANON <- data.frame(
        chr = c("chr10", "chr10", "chr10"), 
        pos = c(21729167, 21729137, 21729166),
        pos = c(21729167, 21729137, 21729166),
        gene = c("MLLT10", "MLLT10", "MLLT10"),
        value = c(1))
      bed_BND_MATES_CANON <- data.frame(
        chr_mate = c("chrX", "chrX", "chrX"),  
        pos_mate = c(41341864, 41341865, 41341865), 
        pos_mate = c(41341864, 41341865, 41341865), 
        gene_mate = c("", "", ""),
        value = c(1))
    }
    
    #MARKING CANONICAL FUSIONS
    #TARGET-20-PAUZRY-09A-02D_consensus -> RUNX1 gene on chromosome 21 with RUNX1T1 gene on chromosome 8
    if (patient_name == "TARGET-20-PAUZRY-09A-02D"){
      bed_BND_CANON <- data.frame("chr21", 34838337, 34838337, "RUNX1", 1)
      bed_BND_MATES_CANON <- data.frame("chr8", 92058671, 92058671, "RUNX1T1", 1)
    }
    
    #MARKING CANONICAL FUSIONS
    #TARGET-20-PAUNVN-09A-01D_Tumor_consensus -> ETV6 gene on chromosome 12 with chromosome 4 
    if (patient_name == "TARGET-20-PAUNVN-09A-01D"){
      chr <- c("chr12", "chr12")
      chr_mate <- c("chr4", "chr4")
      pos <- c(117749121, 117749135)
      pos_mate <- c(59084883, 59081772)
      value <- c(1)
      gene <- c("ETV6", "ETV6")
      gene_mate <- c("", "")
      bed_BND_CANON <- data.frame(chr, pos, pos, gene, value)
      bed_BND_MATES_CANON <- data.frame(chr_mate, pos_mate, pos_mate, gene_mate, value)
    }
    
    #MARKING CANONICAL FUSIONS
    #TARGET-20-PAURDN-03A-01D_consensus -> NUP98 gene on chromosome 11 with chromosome 5
    if (patient_name == "TARGET-20-PAURDN-03A-01D"){
      bed_BND_CANON <- data.frame("chr11", 3743985, 3743985, "NUP98", 1)
      bed_BND_MATES_CANON <- data.frame("chr5", 177233312, 177233312, " ", 1)
    }
    
    #MARKING CANONICAL FUSIONS
    #PAVESI-03A-01D_Tumor_consensus -> NUP98 gene on chromosome 11 with chromosome 5
    if (patient_name == "PAVESI-03A-01D_Tumor"){
      chr <- c("chr11", "chr11")
      chr_mate <- c("chr5", "chr5")
      pos <- c(3743776, 3743777)
      pos_mate <- c(177234100, 177234102)
      value <- c(1)
      gene <- c("NUP98", "NUP98")
      gene_mate <- c("", "")
      bed_BND_CANON <- data.frame(chr, pos, pos, gene, value)
      bed_BND_MATES_CANON <- data.frame(chr_mate, pos_mate, pos_mate, gene_mate, value)
    }
    
    #MARKING CANONICAL FUSIONS
    #TARGET-20-PAUPIY-03A-01D_consensus -> NUP98 gene on chromosome 11 with chromosome 5
    if (patient_name == "TARGET-20-PAUPIY-03A-01D"){
      chr <- c("chr11", "chr11", "chr11")
      chr_mate <- c("chr5", "chr5", "chr5")
      pos <- c(3739358, 3739361, 3680098)
      pos_mate <- c(177216203, 177216209, 177190111)
      value <- c(1)
      gene <- c("NUP98", "NUP98", "NUP98")
      gene_mate <- c("", "", "")
      bed_BND_CANON <- data.frame(chr, pos, pos, gene, value)
      bed_BND_MATES_CANON <- data.frame(chr_mate, pos_mate, pos_mate, gene_mate, value)
    }
    
    #make CIRCOS plots - need to save plot to a variable so it can be added to a list for plotting (in case of more than one patient)
    
    # INITIALIZE WITH IDEOGRAM (all chromosomes)
    circos.par(track.height = 0.08)
    circos.initializeWithIdeogram(species = "hg38")
    
    # ADD PATIENT NAME
    text(1, 1, pt, cex = 1)
    
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
    circos.genomicLink(bed_BND_CANON, bed_BND_MATES_CANON, col = "blue", border = "blue")
    
    #add circos objects to list
    
    
  }
  
  #if more than one patient, patchwork to make plot all circos plots?
  
  
  #if only one patient, just plot
  
  
  
  
  
  
  
  
}
  