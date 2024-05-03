#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
CancerPlotUI <- function(id, label = "TCGA Expression by Cancer Type") {
  ns <- NS(id)
  
  ############################ ----- UI ----- ################################################################################
  
  # Define UI for application that draws a histogram
  tagList(
    fluidPage(
      theme = shinythemes::shinytheme(theme = "paper"), #this is the theme Amanda chose
      
      tabsetPanel(
        tabPanel("Gene Expression", 
                 fluidRow(
                   br(),
                   column(6, radioButtons(ns("scaleType"), label = "Select Transformation:",
                                          choices = c("TPM", "Log2(TPM + 1)"),
                                          selected = "TPM")
                   )
                 ),
                 fluidRow(height = 900, plotOutput(ns("cancerPlot")))
        ),
        
        tabPanel("With GTEx", # Add a new tab called "Total Genes"
                column(3, 
                       br(),
                       fluidRow(height = 250, radioButtons(ns("scaleType1"), label = "Select Transformation:",
                                                           choices = c("TPM", "Log2(TPM + 1)"),
                                                           selected = "TPM")),
                       fluidRow(height = 650, checkboxGroupInput(ns("locations"), label = "Pick a Cancer Type", choices = unique(gtex_tcga_combined$Location), selected = NULL))),
                column(9, plotOutput(ns("canGtexPlot")))
                ),
        
        
        tabPanel("Total Genes", # Add a new tab called "Total Genes"
                 div(style = "overflow-y: scroll; height: 750px;",
                     DT::dataTableOutput(ns("tcgacsv")))
        )
      ),
      div(
        style = "font-size: 24px; text-align: center;",
        uiOutput(ns("message"))
      )
    )
  )
}


############################ ----- SERVER ----- ##########################################################################

CancerPlot <- function(input, output, session, gene) {
  
  filtData <- reactive({
    if (gene() == "" | !gene() %in% tcga_cancer$sample) {
      return(NULL)
    } else {
      dataframe <- subset(tcga_cancer, sample == gene())
      return(dataframe)
    }
  })
  
  output$cancerPlot <- renderPlot({
    if (is.null(filtData())) {
      return(NULL)
    } else {
      
      if (input$scaleType == "TPM") {
      
          p <- ggplot(filtData(), aes(x = cancer, y = as.numeric(TPM), color = cancer, fill = cancer)) +
            geom_jitter(size = 0.75) +
            scale_y_continuous(labels = scales::number_format()) +
            theme_classic(base_size = 16) +
            theme(
              legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_line(color = "lightgray", linewidth = 0.5, linetype = 1)) +
            ggtitle("Gene Expression Per Cancer Type") +
            xlab("") +
            ylab("TPM") +
            stat_summary(fun.y=median, geom="crossbar", size=0.5, width=0.8, color="black")
      
      } else{
        p <- ggplot(filtData(), aes(x = cancer, y = log2(as.numeric(TPM)+1), color = cancer, fill = cancer)) +
          geom_jitter(size = 0.75) +
          scale_y_continuous(labels = scales::number_format()) +
          theme_classic(base_size = 16) +
          theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_line(color = "lightgray", linewidth = 0.5, linetype = 1)) +
          ggtitle("Gene Expression Per Cancer Type") +
          xlab("") +
          ylab("Log2(TPM + 1)") +
          stat_summary(fun.y=median, geom="crossbar", size=0.5, width=0.8, color="black")
      }
    
      print(p)
    }

  }, height = 750)
  
  filtered_df <- reactive({
      if (gene() == "" | !gene() %in% gtex_tcga_combined$Genes) {
        return(NULL)
      } else {
        dataframe <- subset(gtex_tcga_combined, Genes == gene())
            if (length(input$locations) > 0) {
              filtered <- dataframe %>% filter(Location %in% input$locations)
            } else {
              filtered <- NULL
            }
      }
    return(filtered)
  })
  
  output$canGtexPlot <- renderPlot({
    if (is.null(filtered_df())) {
      return(NULL)  # Return nothing if filtered data is NULL
    }else{
      
      if (input$scaleType1 == "TPM") {
      
        p <- ggplot(filtered_df(), aes(x = reorder(Cancer_or_Tissue, -startsWith(Cancer_or_Tissue, "TCGA")), y = as.numeric(TPM), color = Source)) +
          geom_jitter(size = 0.75) +
          facet_grid(. ~ Location, scales = "free_x", space = "free_x") +
          scale_colour_manual(values = c("TCGA" = "#FF4242","GTEx" = "#47abd8")) +
          scale_y_continuous(labels = scales::number_format()) + # This line changes the y-axis formatting
          theme_bw(base_size = 16) +
          theme(
            legend.position = "right",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_line(color = "lightgray", linewidth = 0.5, linetype = 1)
          ) +
          ggtitle("Gene Expression Per Cancer Type") +
          xlab("Location") +
          ylab("TPM") +
          stat_summary(fun.y=median, geom="crossbar", size=0.5, width=0.8, color="black") +
          guides(color = guide_legend(override.aes = list(size = 3)))
      }else{
        
        p <- ggplot(filtered_df(), aes(x = reorder(Cancer_or_Tissue, -startsWith(Cancer_or_Tissue, "TCGA")), y = log2(as.numeric(TPM)+1), color = Source)) +
          geom_jitter(size = 0.75) +
          facet_grid(. ~ Location, scales = "free_x", space = "free_x") +
          scale_colour_manual(values = c("TCGA" = "#FF4242","GTEx" = "#47abd8")) +
          scale_y_continuous(labels = scales::number_format()) + # This line changes the y-axis formatting
          theme_bw(base_size = 16) +
          theme(
            legend.position = "right",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_line(color = "lightgray", linewidth = 0.5, linetype = 1)
          ) +
          ggtitle("Gene Expression Per Cancer Type") +
          xlab("Location") +
          ylab("Log2(TPM+1)") +
          stat_summary(fun.y=median, geom="crossbar", size=0.5, width=0.8, color="black") +
          guides(color = guide_legend(override.aes = list(size = 3)))
      
      }
    
    print(p)
    }
  }, height = 750)
  
  
  output$message <- renderUI({
    
    data <- filtData()
    
    if (gene() == "") {
      return("Type a gene into the text box to the upper left.")
    } else if (!gene() %in% data$sample) {
      return("This gene is not in the dataset.")
    } else {
      return(NULL)
    }
  })
  
  
  output$tcgacsv <- DT::renderDataTable({
    
    tcga_csv <- tcga_csv[,c(1:3)]

    DT::datatable(tcga_csv, 
                  class = "compact nowrap hover row-border order-column", # Defines the CSS formatting of the final table, can string multiple options together
                  options = list(
                    dom = 'frtip',  # 'f' for filter/search, 'r' for processing info, 't' for table
                    scrollY = TRUE,
                    searchHighlight = TRUE,
                    pageLength = 50), 
                  escape = FALSE) %>%
      DT::formatStyle(columns = c(1, 2, 3), fontSize = "100%")  # Apply formatting to all three columns
  })
  
}

############################ ----- RUN ----- ###############################################################################

#Run the application 
#shinyApp(ui = ui, server = server)

############################################################################################################################