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
        tabPanel("TCGA Expression", 
                 fluidRow(
                   br(),
                   column(6, radioButtons(ns("scaleType_t"), label = "Select Transformation:",
                                          choices = c("TPM", "Log2(TPM + 1)"),
                                          selected = "TPM")
                   )
                 ),
                 fluidRow(height = 900, plotOutput(ns("tcgaPlot")))
        ),
        
        tabPanel("TCGA Summary Table", # Add a new tab called "Total Genes"
                 br(),
                 div(style = "overflow-y: scroll; height: 750px;",
                     DT::dataTableOutput(ns("tcga_table")))
        ),
        
        tabPanel("GTEx Expression", 
                 fluidRow(
                   br(),
                   column(6, radioButtons(ns("scaleType_g"), label = "Select Transformation:",
                                          choices = c("TPM", "Log2(TPM + 1)"),
                                          selected = "TPM")
                   )
                 ),
                 fluidRow(height = 900, plotOutput(ns("gtexPlot")))
        ),
        
        tabPanel("GTEx Summary Table", # Add a new tab called "Total Genes"
                 br(),
                 div(style = "overflow-y: scroll; height: 750px;",
                     DT::dataTableOutput(ns("gtex_table")))
        ),
        
        tabPanel("With GTEx", # Add a new tab called "Total Genes"
                 column(3, 
                        br(),
                        fluidRow(height = 250, radioButtons(ns("scaleType_tg"), label = "Select Transformation:",
                                                            choices = c("TPM", "Log2(TPM + 1)"),
                                                            selected = "TPM")),
                        fluidRow(height = 650, checkboxGroupInput(ns("locations"), label = "Pick a Cancer Type", choices = unique(sort(tcga_manifest$Location)), selected = NULL))),
                 column(9, 
                        br(),
                        plotOutput(ns("canGtexPlot")))
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
  
  ## filtering the tcga dataset
  filtData <- reactive({
    if (gene() == "" | !gene() %in% tcga_newcsv$sample) {
      return(NULL)
    } else {
      dataframe <- subset(tcga_newcsv, sample == gene())
      
      dataframe <- dataframe %>%
        pivot_longer(cols = -sample,    # All columns except PRAME
                     names_to = "ID",  # Create a new "Sample" column for the sample names
                     values_to = "TPM")
      
      dataframe <- merge(dataframe, tcga_manifest, by.x = "ID", by.y = "Sample")
      
      return(dataframe)
    }
  })
  
  ## filtering the gtex dataset
  filtData_g <- reactive({
    if (gene() == "" | !gene() %in% gtex_csv$Description) {
      return(NULL)
    } else {
      dataframe <- subset(gtex_csv, Description == gene())
      
      dataframe <- dataframe %>%
        pivot_longer(cols = -Description,    # All columns except PRAME
                     names_to = "ID",  # Create a new "Sample" column for the sample names
                     values_to = "TPM")
      
      dataframe <- merge(dataframe, gtex_manifest, by.x = "ID", by.y = "Sample")
      
      dataframe$Tissue <- gsub("_", " ", dataframe$Tissue)
      
      dataframe$TPM <- as.numeric(dataframe$TPM)
      dataframe <- na.omit(dataframe)
      
      return(dataframe)
    }
  })
  
  ## organizing the gtex/tcga combined dataset for comparison
  filtData_tg <- reactive({
    if (gene() == "" | !gene() %in% tcga_newcsv$sample) {
      return(NULL)
    } else {
      dataframe_gtex <- subset(gtex_csv, Description == gene())
      dataframe_gtex <- dataframe_gtex %>%
        pivot_longer(cols = -Description,    # All columns except PRAME
                     names_to = "ID",  # Create a new "Sample" column for the sample names
                     values_to = "TPM")
      dataframe_gtex$Source <- "GTEx"
      
      dataframe_tcga <- subset(tcga_newcsv, sample == gene())
      dataframe_tcga <- dataframe_tcga %>%
        pivot_longer(cols = -sample,    # All columns except PRAME
                     names_to = "ID",  # Create a new "Sample" column for the sample names
                     values_to = "TPM")
      dataframe_tcga$Source <- "TCGA"
    
      dataframe_gtex <- merge(dataframe_gtex, gtex_manifest, by.x = "ID", by.y = "Sample")
      dataframe_tcga <- merge(dataframe_tcga, tcga_manifest, by.x = "ID", by.y = "Sample")
      
      dataframe_gtex <- dataframe_gtex[,-6]
      dataframe_tcga <- dataframe_tcga[,-6]
      
      colnames(dataframe_gtex)[2] <- "Sample"
      colnames(dataframe_tcga)[2] <- "Sample"
      colnames(dataframe_gtex)[5] <- "Cancer_or_Tissue"
      colnames(dataframe_tcga)[5] <- "Cancer_or_Tissue"
      
      dataframe_gtex$Cancer_or_Tissue <- paste0("GTEx ", dataframe_gtex$Cancer_or_Tissue)
      dataframe_tcga$Cancer_or_Tissue <- paste0("TCGA ", dataframe_tcga$Cancer_or_Tissue)
      dataframe_gtex$Cancer_or_Tissue <- gsub("_", " ", dataframe_gtex$Cancer_or_Tissue)
      
      dataframe <- rbind(dataframe_gtex, dataframe_tcga)
      
      if (length(input$locations) > 0) {
        filtered <- dataframe %>% filter(Location %in% input$locations)
        return(filtered)
      } else {
        return(NULL)
      }
    }
  })
  
  
  output$tcgaPlot <- renderPlot({
    if (is.null(filtData())) {
      return(NULL)
    } else {
      
      if (input$scaleType_t == "TPM") {
        
        p <- ggplot(filtData(), aes(x = Cancer, y = as.numeric(TPM), color = Cancer, fill = Cancer)) +
          geom_jitter(size = 0.75, height = 0) +
          scale_y_continuous(labels = scales::number_format()) +
          theme_classic(base_size = 16) +
          theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_line(color = "lightgray", linewidth = 0.5, linetype = 1),
            plot.margin = margin(t = 10, r = 30, b = 0, l = 40)) +
          ggtitle(paste(gene(), "Expression Per Cancer Type")) +
          xlab("") +
          ylab("TPM") +
          stat_summary(fun.y=median, geom="crossbar", size=0.5, width=0.8, color="black")
        
      } else{
        p <- ggplot(filtData(), aes(x = Cancer, y = log2(as.numeric(TPM)+1), color = Cancer, fill = Cancer)) +
          geom_jitter(size = 0.75, height = 0) +
          scale_y_continuous(labels = scales::number_format()) +
          theme_classic(base_size = 16) +
          theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_line(color = "lightgray", linewidth = 0.5, linetype = 1),
            plot.margin = margin(t = 10, r = 30, b = 0, l = 40)) +
          ggtitle(paste(gene(), "Expression Per Cancer Type")) +
          xlab("") +
          ylab("Log2(TPM + 1)") +
          stat_summary(fun.y=median, geom="crossbar", size=0.5, width=0.8, color="black")
      }
      
      print(p)
    }
    
  }, height = 750)
  
  output$gtexPlot <- renderPlot({
    if (is.null(filtData_g())) {
      return(NULL)
    } else {
      
      if (input$scaleType_g == "TPM") {
        
        p <- ggplot(filtData_g(), aes(x = Tissue, y = as.numeric(TPM), color = Tissue, fill = Tissue)) +
          geom_jitter(size = 0.75, height = 0) +
          scale_y_continuous(labels = scales::number_format()) +
          theme_classic(base_size = 16) +
          theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_line(color = "lightgray", linewidth = 0.5, linetype = 1),
            plot.margin = margin(t = 10, r = 30, b = 0, l = 40)) +
          ggtitle(paste(gene(), "Expression Per Tissue Type")) +
          xlab("") +
          ylab("TPM") +
          stat_summary(fun.y=median, geom="crossbar", size=0.5, width=0.8, color="black")
        
      } else{
        p <- ggplot(filtData_g(), aes(x = Tissue, y = log2(as.numeric(TPM)+1), color = Tissue, fill = Tissue)) +
          geom_jitter(size = 0.75, height = 0) +
          scale_y_continuous(labels = scales::number_format()) +
          theme_classic(base_size = 16) +
          theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_line(color = "lightgray", linewidth = 0.5, linetype = 1),
            plot.margin = margin(t = 10, r = 30, b = 0, l = 40)) +
          ggtitle(paste(gene(), "Expression Per Tissue Type")) +
          xlab("") +
          ylab("Log2(TPM + 1)") +
          stat_summary(fun.y=median, geom="crossbar", size=0.5, width=0.8, color="black")
      }
      
      print(p)
    }
    
  }, height = 750)
  
  output$canGtexPlot <- renderPlot({
    if (is.null(filtData_tg())) {
      return(NULL)  # Return nothing if filtered data is NULL
    }else{
      
      if (input$scaleType_tg == "TPM") {
        
        print(filtData_tg())
        
        p <- ggplot(filtData_tg(), aes(x = reorder(Cancer_or_Tissue, -startsWith(Cancer_or_Tissue, "TCGA")), y = as.numeric(TPM), color = Source)) +
          geom_jitter(size = 0.75) +
          facet_grid(. ~ Location, scales = "free_x", space = "free_x") +
          scale_colour_manual(values = c("TCGA" = "#FF4242","GTEx" = "#47abd8")) +
          scale_y_continuous(labels = scales::number_format()) + # This line changes the y-axis formatting
          theme_bw(base_size = 16) +
          theme(
            legend.position = "right",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_line(color = "lightgray", linewidth = 0.5, linetype = 1),
            strip.background=element_rect(colour="black", fill="#3c8dbd"),
            strip.text = element_text(face = "bold")
          ) +
          ggtitle(paste(gene(), "Expression in TCGA Cancer and GTEx Normal Tissue")) +
          xlab("Location") +
          ylab("TPM") +
          stat_summary(fun.y=median, geom="crossbar", size=0.5, width=0.8, color="black") +
          guides(color = guide_legend(override.aes = list(size = 3)))
      }else{
        
        p <- ggplot(filtData_tg(), aes(x = reorder(Cancer_or_Tissue, -startsWith(Cancer_or_Tissue, "TCGA")), y = log2(as.numeric(TPM)+1), color = Source)) +
          geom_jitter(size = 0.75) +
          facet_grid(. ~ Location, scales = "free_x", space = "free_x") +
          scale_colour_manual(values = c("TCGA" = "#FF4242","GTEx" = "#47abd8")) +
          scale_y_continuous(labels = scales::number_format()) + # This line changes the y-axis formatting
          theme_bw(base_size = 16) +
          theme(
            legend.position = "right",
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_line(color = "lightgray", linewidth = 0.5, linetype = 1),
            strip.background=element_rect(colour="black", fill="#3c8dbd"),
            strip.text = element_text(face = "bold")
          ) +
          ggtitle(paste(gene(), "Expression in TCGA Cancer and GTEx Normal Tissue")) +
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
  
  
  ### This is the table for showing the data by cancer
  output$tcga_table <- DT::renderDataTable({
    if (gene() == "") {
      return(NULL)
    }
    else{
      data <- filtData()
      colnames(data)[2] <- "Sample"
      
      data <- data %>%
        group_by(Cancer) %>%
        summarize(
          N = n(),  # Count of samples for each cancer type
          `Median TPM` = round(median(TPM, na.rm = TRUE), 2),  # Median TPM expression
          `Mean TPM` = round(mean(TPM, na.rm = TRUE), 2),  # Mean TPM expression
          `Min TPM` = round(min(TPM, na.rm = TRUE), 2),
          `Max TPM` = round(max(TPM, na.rm = TRUE), 2),
          `Range TPM` = round(as.numeric(max(TPM, na.rm = TRUE) - min(TPM, na.rm = TRUE)), 2)
        ) %>%
        ungroup()
      
      DT::datatable(data, 
                    class = "compact nowrap hover row-border order-column",
                    extensions = 'Buttons', # Defines the CSS formatting of the final table
                    options = list(
                      dom = 'Bfrtip',  # 'f' for filter/search, 'r' for processing info, 't' for table
                      scrollY = TRUE,
                      searchHighlight = TRUE,
                      pageLength = 100,
                      buttons = c('csv', 'excel')), 
                    escape = FALSE) %>%
        DT::formatStyle(columns = c(1, 2, 3, 4, 5, 6, 7), fontSize = "100%") %>%  # Apply formatting to all columns
        DT::formatStyle(
          columns = 3,  # Specify the columns to style
          background = styleColorBar(range(data[, 3], na.rm = TRUE), "#b0d2e6"),  # Corrected this line
          backgroundSize = '90% 88%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'right')
    }
  })
  
  ### This is the table for showing the data by tissue
  output$gtex_table <- DT::renderDataTable({
    if (gene() == "") {
      return(NULL)
    }
    else{
      data <- filtData_g()
      colnames(data)[2] <- "Sample"
      
      data <- data %>%
        group_by(Tissue) %>%
        summarize(
          N = n(),  # Count of samples for each cancer type
          `Median TPM` = round(median(TPM, na.rm = TRUE), 2),  # Median TPM expression
          `Mean TPM` = round(mean(TPM, na.rm = TRUE), 2),  # Mean TPM expression
          `Min TPM` = round(min(TPM, na.rm = TRUE), 2),
          `Max TPM` = round(max(TPM, na.rm = TRUE), 2),
          `Range TPM` = round(as.numeric(max(TPM, na.rm = TRUE) - min(TPM, na.rm = TRUE)), 2)
        ) %>%
        ungroup()
      
      DT::datatable(data, 
                    class = "compact nowrap hover row-border order-column",
                    extensions = 'Buttons',# Defines the CSS formatting of the final table
                    options = list(
                      dom = 'Bfrtip',  # 'f' for filter/search, 'r' for processing info, 't' for table
                      scrollY = TRUE,
                      searchHighlight = TRUE,
                      pageLength = 100,
                      buttons = c('csv', 'excel')), 
                    escape = FALSE) %>%
        DT::formatStyle(columns = c(1, 2, 3, 4, 5, 6, 7), fontSize = "100%") %>%  # Apply formatting to all columns
        DT::formatStyle(
          columns = 3,  # Specify the columns to style
          background = styleColorBar(range(data[, 3], na.rm = TRUE), "#b0d2e6"),  # Corrected this line
          backgroundSize = '90% 88%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'right')
    }
  })
}

############################ ----- RUN ----- ###############################################################################

#Run the application 
#shinyApp(ui = ui, server = server)

############################################################################################################################
