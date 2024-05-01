#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(data.table)

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