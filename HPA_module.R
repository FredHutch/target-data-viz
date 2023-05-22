#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

HPAPlotUI <- function(id, label = "Human Protein Atlas Supporting"){

############################ ----- PACKAGES ----- ##########################################################################

library(shiny)
library(shinythemes)
library(viridis)
library(viridisLite)
library(tidyverse)
library(ggplot2)
library(dqshiny)
library(thematic)
library(shinyalert)
library(plotly)
library(DT)
library(fst)

ns <- NS(id)
  
############################ ----- DATA ----- ##############################################################################

#read in the data
data <- read.fst("data/rna_immune_cell_sample.fst") #the HPA data
all_genes <- readRDS("data/all_genes.RDS") #a list of genes from both of the datasets for autocorrection
sub <- read.fst("data/subcellular_location.fst")
protein <- read.fst("data/tissue_data.fst")




############################ ----- UI ----- ################################################################################

# Define UI for application that draws a histogram
  tagList(
    fluidPage(
      theme = shinythemes::shinytheme(theme = "paper"), #this is the theme Amanda chose
      
      # Sidebar layout
      sidebarLayout(
        sidebarPanel(
          
          #text for information
          p("HPA Immune Cell scData was pulled from the Human Protein Atlas. It displays transcript expression 
        levels for immune cells in normalized transcripts per million (TPM) for a gene inputted above. 
        The data was extracted from 109 samples."),
          p("Tissue type data was pulled from the Human Protein Atlas from tissue micro arrays. Protein expression scores 'low, medium, high' are literature-based annotations."),
          p(em("WARNING: If the plot is blank it means the gene does not
        exist in the dataset.")),
          hr(),
          p(strong("Cell Type")),
          p("Uhlen M et al., A genome-wide transcriptomic analysis of protein-coding genes in human blood cells. Science. (2019) 366(6472) PubMed: 31857451 DOI: 10.1126/science.aax9198"),
          p("Data available from v22.0.proteinatlas.org: 19. rna_immune_cell_sample.tsv.zip"),
          p(strong("Tissue Type")),
          p("Uhlén M et al., Tissue-based map of the human proteome. Science (2015) PubMed: 25613900 DOI: 10.1126/science.1260419"),
          p("Data available from v22.0.proteinatlas.org: 1. normal_tissue.tsv.zip"),
          p(em("Accessed on: 05-09-2023")),
          #p(em("Header graphic created with BioRender.com")),
          width = 3, height = 900),
        
        
        mainPanel(width = 9, height = 900,
                  #creates tabs
                  tabsetPanel(
                    #outputs the plots and sets the tab to graph ratio
                    tabPanel("Cell Types", br(), plotlyOutput(ns("Plot1"))),
                    tabPanel("Tissue Types", br(), plotlyOutput(ns("Plot2"))),
                  ),
                  
                  br(),
                  br(),
                  br(),
                  br(),
                  
                  #creates a row below the plots for the table
                  fluidRow(
                    br(),
                    br(),
                    br(),
                    br(),
                    column(12, br(), br(), DT::dataTableOutput(ns("mytable")))
                  ),
                  br(),
                  br(),
                  br(),
                  br(),
                  
        )
      )
    )
)
}
    ############################ ----- SERVER ----- ##########################################################################
    
    # Define server logic required to draw a histogram
    HPAPlot <- function(input, output, session, gene) {
      
      
      #makes the shiny app reactive to the gene input and plotting
      toListen <- reactive({
        list(input$geneInput, input$Plot, input$table)
      })
      
    
      
      observeEvent(toListen(),{
        
        filtData <- reactive({
          
          data <- data %>%
            filter(gene == gene()) %>%
            mutate(group = factor(group, levels <- c("basophil", "eosinophil", "neutrophil", "classical monocyte",
                                                         "non-classical monocyte", "intermediate monocyte",
                                                         "T-reg", "gdT-cell", "MAIT T-cell",
                                                         "memory CD4 T-cell", "naive CD4 T-cell",
                                                         "memory CD8 T-cell", "naive CD8 T-cell",
                                                         "memory B-cell", "naive B-cell", "plasmacytoid DC",
                                                         "myeloid DC", "NK-cell", "total PBMC")))
          data <- data %>%
            filter(gene == gene()) %>%
            mutate(celltype = factor(celltype, levels <- c("Granulocytes", "Monocytes", 
                                                           "T-cells", "B-cells", "Dendritic cells",
                                                           "NK-cells", "Total PBMC")))
          
        })
        
      filtProt <- reactive({
        
          protein <- protein %>%
            filter(`Gene name` == gene()) %>%
            mutate(Level = factor(Level, levels <- c("Low", "Medium", "High")))
          
        })
        
      filtSub <- reactive({
        sub <- sub[1:6] %>%
          filter(`Gene name` == gene())
        })
        
        ############################ ----- Plots ----- ###########################################################      
      plotHist <- reactive({
        #bar chart for the HPA data
        
        HPAplot <- filtData() %>% ggplot(aes(x=group, y=TPM, fill=celltype)) + 
          scale_fill_viridis(option="viridis", discrete = TRUE, name = "Cell Type") +
          geom_bar(stat='identity') +
          theme_classic(base_size = 12) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_line(color = "lightgray", linewidth = 0.5, linetype = 1)) +
          ggtitle("Immune Cell-Type Expression in Transcripts Per Million") +
          xlab("") +
          ylab("TPM")
        
      })
        
        
      plotProt <- reactive({ 
        
        Protplot <- filtProt() %>% ggplot(aes(x=Tissue, y=Level, fill=Tissue, text = Level)) +
          scale_fill_viridis(option="viridis", discrete = TRUE) +
          geom_bar(stat='identity') +
          theme_classic(base_size = 12) +
          theme(
            legend.position="none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_line(color = "lightgray", linewidth = 0.5, linetype = 1)) +
          ggtitle("Protein Expression Per Tissue Type") +
          xlab("") +
          ylab("Protein Expression Score")
        
      }) 
      
      
        #creates the dropdown tab for switching between plots    
        output$Plot1 <- renderPlotly({
          ggplotly(plotHist(), tooltip = "y", height = 600)})
        
        output$Plot2 <- renderPlotly({
          ggplotly(plotProt(), tooltip = "text", height = 600)})
        
        output$mytable <- DT::renderDataTable(
          filtSub(), 
          selection = 'none',
          escape = FALSE,
          filter = "none",
          class = "row-border",
          rownames = FALSE,
          options = list(
            paging = FALSE,
            autowidth = TRUE,
            searching = FALSE,
            dom = "t",
            ordering = FALSE)
        )
      })
    }
    
  ############################ ----- RUN ----- ###############################################################################
  
  #Run the application 
  #shinyApp(ui = ui, server = server)
  
  ############################################################################################################################