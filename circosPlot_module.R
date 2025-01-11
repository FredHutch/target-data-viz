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
  
  ##----------------------------------------------------------------
  ##                          Page Setup                           --
  ##----------------------------------------------------------------

  tagList(
  useShinyjs(),
  fluidPage(
  theme = shinythemes::shinytheme(theme = "paper"),


  sidebarLayout(
    #sidebar on left screen
    position = "left",
    sidebarPanel(

      #dropdown to select fusion class you want to see circos plots for
      selectInput(ns("fusion_group"),
                  label = "Select a fusion type",
                  choices = canonical_fusions),

      #buttons to select whether circos plots generated will show all chromosomes
      #or only two chromosomes
      radioButtons(ns("all_chroms"),
                   label = "View all chromosomes?",
                   choices = list("All chromosomes" = "all",
                                  "Two chromosomes" = "two")),

      #conditional panel for if the user chooses to view only two chromosomes
      conditionalPanel(
        condition = paste0("input['", ns("all_chroms"), "'] == 'two'"),
        #dropdown input for one of the selected chromosomes
        selectInput(ns("chromA"), label = "Select first chromosome:", choices = chromosome_choices),
        #dropdown input for the second of the selected chromosomes -- GOING TO NEED TO MAKE SURE THIS CANNOT BE THE SAME AS CHROM_A
        selectInput(ns("chromB"), label = "Select second chromosome:", choices = chromosome_choices))

        ),



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
                   ))))))))
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

    #make CIRCOS plots - need to save plot to a variable so it can be added to a list for plotting (in case of more than one patient)


    #add circos objects to list
    
    
  }
  
  #if more than one patient, patchwork to make plot all circos plots?
  
  
  #if only one patient, just plot
  
  
  
  
  
  
  
  
}
  