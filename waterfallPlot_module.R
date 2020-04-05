# UI function for the waterfall plot module
wfPlotUI <- function(id, label = "Gene expression plot parameters"){
  
  library(DT)
  library(shinyjs)
  ns <- NS(id) # Setting a unique namespace for this module
  
  tagList(
    useShinyjs(),
            
            ###############################################################
            #----------------------- SIDEBAR -----------------------------#
            ###############################################################
            
            sidebarLayout(
              
              # Placing the sidebar on the left side of the screen
              position = "left",
              
              sidebarPanel(
                
                # Dropdown menu to select variable to use for arranging/grouping patients in waterfall plot
                selectInput(ns("grouping_var"), 
                            label = "Select a grouping variable",             # The name of each list item is what is shown in the box;
                            choices = list("AML vs. normal" = "Group",        # the value corresponds to a column of the CDEs
                                           "Cell lines" = "Cell.Lines", 
                                           "Primary cytogenetics" = "Primary.Cytogenetic.Code", 
                                           "Primary fusions" = "Primary.Fusion",
                                           "Rare fusions" = "Rare.Fusions", 
                                           "KMT2A/MLL fusions" = "MLL.Fusion",
                                           "SNVs" = "SNVs", 
                                           "Age category" = "Age.Category", 
                                           "Cyto/molecular risk" = "Cyto.Risk", 
                                           "Event type" = "Event.Type",
                                           "CNS disease" = "CNS.Disease.Category")),
                
                radioButtons(ns("plot_type"), 
                             label = "Select a type of plot to generate", 
                              choices = list("Waterfall plot" = "wf", 
                                             "Box/violin plots" = "bx")),
                
                conditionalPanel(
                  condition = paste0("input['", ns("plot_type"), "'] == 'bx'"),
                  checkboxInput(ns("log"),                                                  
                            label = "Log2 transform the data",
                            value = FALSE)),
                
                conditionalPanel(
                  condition = paste0("input['", ns("plot_type"), "'] == 'bx'"),
                  checkboxInput(ns("test"),                                                  
                                label = "Perform significance tests",
                                value = FALSE)),
                
                conditionalPanel(
                  condition = paste0("input['", ns("test"), "'] == 1"),
                  checkboxGroupInput(ns("comparisons"),                                                  
                                label = "Select 2 groups to compare",
                                choices = c("A", "B", "C", "D"))), # These are just placeholders and will be replaced
                                                                   # on the server side with the appropriate categories
                
                helpText("The grouping variable will be used to arrange patients along the x axis (for waterfall plots) 
                          or to group patients together (for box and violin plots) based on a common clinical characteristic 
                          to help highlight expression patterns within the groups." ),
                
                helpText("NOTE: If cell lines is selected, please reference the 'Summary Stats' tab for expression data, 
                         as only one sample is available for most cell lines."),
                
                br(),
                
                hidden(
                  actionButton(ns("adc_flag"), "See ADC or CAR T-cell therapies")
                )
                
              ),
              
              
              ###############################################################
              #----------------------- MAIN PLOT PANEL ---------------------#
              ###############################################################
              
              mainPanel(
                
                position = "right",
                tags$head(tags$style(HTML(".small-box {height: 40px}"))),
                tags$style(HTML(".fa{font-size: 20px;}")),
                
                tabsetPanel(
                  
                  #-------------------- Waterfall plot -----------------------#
                  tabPanel("Waterfall plot", # This is the title of the tab panel, NOT the name of the plot object!
                           br(),   
                           br(),             # Linebreaks to help center the plot on the page
                           fluidRow(
                             column(10, offset = 0, align = "left",                   # This will be a reactive object that is linked to an item in the 
                                    plotOutput(ns("plot"), width = "100%")),          # output list, created in the "server" script
                             column(2, offset = 0, align = "right",                   
                                    downloadButton(ns("plot_download"), 
                                                   label = "Download plot")),
                             column(2, offset = 0, align = "right", 
                                    downloadButton(ns("ggplot_download"), 
                                                   label = "ggplot2 object", 
                                                   style = 'padding:5px; font-size:70%; margin-top:10px',
                                                   class = "btn-info"))
                           )
                  ),
                  
                  #-------------------- Summary table -----------------------#
                  tabPanel("Summary stats", 
                           br(),
                           br(),
                           fluidRow(
                             column(10, offset = 0, align = "left", 
                                    DT::dataTableOutput(ns("table"))), # Table of summary stats for plot
                             column(1, offset = 0, align = "right", 
                                    downloadButton(ns("table_download"), 
                                                   label = "Download table"))
                           )
                  )
                )
              ) 
            )
  )
}



# Server function for the waterfall plot module
wfPlot <- function(input, output, session, clinData, countsData, adc_cart_targetData, gene, parent){
  
  library(tidyverse)
  library(DT)
  `%then%` <- shiny:::`%OR%` # See https://shiny.rstudio.com/articles/validation.html for details on the %then% operator
  
  #-------------------- Data preparation -----------------------#
  
  # Selecting the ADC & CAR T-cell therapy data for only the gene of interest
  therapyData <- reactive({
    filter(adc_cart_targetData, `Gene target` == gene())
  })
  
  # Cleaning up the clinical data grouping columns so they're displayed nicely in the plots
  clinData <- clinData %>%
    rename(Primary.Cytogenetic.Code = `Primary Cytogenetic Code`) %>%
    rename(Cyto.Risk = `Cyto/Fusion/Molecular Risk`) %>%
    rename(CNS.Disease.Category = `CNS disease at on-study`) %>%
    rename(Event.Type = `EFS event type ID`) %>%
    mutate_at(vars(Cytogenetic.Category.2), ~gsub("\\.|\\_", " ", .)) %>%
    mutate_at(vars(Rare.Fusions), ~gsub("\\.", "\\-", .)) %>%
    mutate_at(vars(SNVs), ~gsub("\\.", " ", .)) %>%
    mutate_at(vars(SNVs), ~gsub("mutation", "mut", .)) %>%
    mutate_at(vars(SNVs), ~gsub("with", "w\\/", .)) %>%
    mutate_at(vars(Cytogenetic.Category.2, Rare.Fusions, SNVs), ~gsub("OtherAML", "Other AML", .)) %>%
    
    # Adding columns to use for re-categorizing MLL and other fusions if they occur in less than 10 patients 
    mutate(Primary.Fusion = ifelse(str_detect(pattern = "KMT2A-", string = `Primary Fusion/CNV`), "KMT2A-X", `Primary Fusion/CNV`)) %>%
    mutate(MLL.Fusion = ifelse(str_detect(pattern = "KMT2A", string = `Primary Fusion/CNV`), `Primary Fusion/CNV`, "Other AML")) %>%
    group_by(MLL.Fusion) %>%
    mutate(Fusion.Count = n()) %>%
    ungroup() %>%
    mutate_at(vars(MLL.Fusion), ~ifelse(Fusion.Count < 10, "Other MLL", .)) %>%
    group_by(Primary.Fusion) %>%
    mutate(Fusion.Count = n()) %>%
    ungroup() %>%
    mutate_at(vars(Primary.Fusion), ~ifelse(Fusion.Count < 10, "Other AML", .)) %>%
    
    # Adding rows for the cell lines, which are not typically included in the clinical data
    tibble::add_row(USI = c("K562.01", "K562.02", "ME1", "MO7E", "NOMO1", "Kasumi.D1", "MV4.11.D1"), Group = NA) %>%
    tibble::add_row(USI = grep("^RO|^BM", colnames(countsData), value = T), Group = "NBM") %>%
    tibble::add_row(USI = grep("34POS", colnames(countsData), value = T), Group = "CD34+ PB")
  
  # Filtering the counts data to only retain the gene of interest & throwing
  # errors if a non-existent gene is provided
  expData <- reactive({
    validate(
      need(gene(), "Please enter a gene symbol in the text box.") %then%
        need(gene() %in% rownames(countsData), "That gene symbol does not exist in the counts data! \nDouble-check the symbol, or try an alias/synonym."))
    
    countsData[gene(), intersect(colnames(countsData), clinData$USI)]
  })
  
  
  # Transforming the counts into a long-format dataframe to use with ggplot
  plotData <- reactive({
    
    plotData <- expData() %>%
      gather(USI, Expression) %>%
      mutate_at(vars(Expression), ~as.numeric(.)) %>%
      left_join(., clinData, by = "USI") %>%
      mutate(Log2 = log2(Expression + 1),
             Cell.Lines = ifelse(USI %in% c("K562.01", "K562.02", "ME1", "MO7E", "NOMO1", "Kasumi.D1", "MV4.11.D1"), USI, "AML samples")) %>%
      
      # Modifying the chosen grouping variable to keep the NBMs and PBs from being categorized as NA
      mutate_at(vars(!!input$grouping_var), ~case_when(Group == "NBM" ~ "NBM", 
                                                       Group == "CD34+ PB" ~ "CD34+ PB",
                                                       TRUE ~ .)) %>%
      
      mutate_at(vars(MLL.Fusion, Rare.Fusions, Primary.Fusion, SNVs), ~forcats::fct_relevel(., "Other AML", after = Inf)) %>%
      mutate_at(vars(Cell.Lines), ~forcats::fct_relevel(., "AML samples", after = Inf)) %>%
      mutate_at(vars(!!input$grouping_var), ~forcats::fct_relevel(., "CD34+ PB", after = Inf)) %>%
      mutate_at(vars(!!input$grouping_var), ~forcats::fct_relevel(., "NBM", after = Inf)) %>%
      group_by(!!input$grouping_var) %>%                                                              
      arrange_(input$grouping_var, "Expression")  # Reordering patients so that the specified groups are 
                                                  # grouped together and ordered by increasing expression
    
    # Setting the patient order using factor levels, so they won't be rearranged 
    # alphabetically by ggplot (this step is required for ggplot waterfall plots)
    plotData$USI <- factor(plotData$USI, levels = plotData$USI)
    
    plotData
  })
  
  # Updating the options for significance testing to reflect the
  # grouping variable chosen by the user
  observe({
    x <- unique(plotData()[[input$grouping_var]])
    x <- x[!is.na(x)]
    
    updateCheckboxGroupInput(session, 
                             inputId = "comparisons", 
                             label = "Select 2 of the following to compare",
                             choices = x)
  })
  
  # Making a function that will generate the waterfall plot and 
  # can be called from multiple places in the script
  plotFun <- reactive({ 
    
    if (input$plot_type == "bx") {   # Modifying the axis labels and columns used 
      if (input$log == TRUE) {       # if the selected plot type is a boxplot 
        expCol <- "Log2"
        yaxlab <- paste0(gene(), " Expression (log2 TPM + 1)\n")
      } else {
        expCol <- "Expression"
        yaxlab <- paste0(gene(), " Expression (TPM)\n")
      }
      
      p <- plotData() %>% 
        drop_na(input$grouping_var) %>%
        ggplot(aes_string(x = input$grouping_var, y = expCol, fill = input$grouping_var)) +
        theme_classic() +
        labs(x = "\nCategories", y = yaxlab, fill = gsub("\\.", " ", input$grouping_var)) +
        theme(axis.text.x = element_blank(),
              plot.title = element_text(size = 15, hjust = 0.5),
              axis.ticks = element_blank(),
              legend.position = "bottom") +
        geom_violin(scale = "width", aes_string(color = input$grouping_var)) +
        geom_boxplot(width = 0.1, fill = "white") +
        guides(color = FALSE)
      
    } else { # Generating a waterfall plot if boxplot/violin plot is not selected
      
      p <- plotData() %>% 
        drop_na(input$grouping_var) %>%
        ggplot(aes_string(x = "USI", y = "Expression", fill = input$grouping_var)) +
        theme_classic() +
        labs(x = "\nPatients", y = paste0(gene(), " Expression (TPM)\n"), fill = gsub("\\.", " ", input$grouping_var)) +
        theme(axis.text.x = element_blank(),
              plot.title = element_text(size = 15, hjust = 0.5),
              axis.ticks = element_blank(),
              legend.position = "bottom") +
        geom_bar(stat = "identity", width = 1, position = position_dodge(width = 0.4))
      p
    }
    
    if (length(input$comparisons) > 1) {
      
      validate(
        need(length(input$comparisons > 1), "Please select 2 groups to compare."))
      
      c <- p + ggpubr::stat_compare_means(method = "wilcox.test", comparisons = list(input$comparisons))
      c
    } else {
      p
    }
  })
  
  # Function to generate an expression summary table from the plot data
  tableFun <- reactive({
    plotData() %>%
      drop_na(input$grouping_var) %>%
      group_by_(input$grouping_var) %>%
      # group_by(!!input$grouping_var) %>%
      summarize(N = n(), 
                `Mean (TPM)` = round(mean(Expression, na.rm = T), 2), 
                `Median (TPM)` = round(median(Expression, na.rm = T), 2), 
                `Range (TPM)` = paste0(round(min(Expression), 2), " - ", round(max(Expression), 2)))
  })
  
  #-------------------- Waterfall plot tab -----------------------#
  
  # Saving the plot to the output list object so it can be run & saved reactively
  output$plot <- renderPlot({
    plotFun()
  })
  
  # Adding a download button widget for the plot
  output$plot_download <- downloadHandler(
    filename = function() {
      paste0("TARGET_AAML1031_", gene(), "_", input$grouping_var, "_", input$plot_type, "_generated_", format(Sys.time(), "%m.%d.%Y"), ".png")
    }, 
    content = function(file) {
      ggsave(filename = file, plot = plotFun(), width = 7, height = 5, device = "png", dpi = 150)
    }
  )
  
  output$ggplot_download <- downloadHandler(
    filename = function() {
      paste0("TARGET_AAML1031_", gene(), "_", input$grouping_var, "_ggplotObject_generated_", format(Sys.time(), "%m.%d.%Y"), ".RDS")
    },
    content = function(file) {
      withProgress(message = "Saving RDS file", detail = "This may take a while...", value = 0, {
        for (i in 1:50) {
          incProgress(1/70)
          Sys.sleep(0.25)
        }
      })
      saveRDS(object = plotFun(), file = file, compress = F)
    }
  )
  
  #-------------------- Summary table tab -----------------------#
  
  output$table <- DT::renderDataTable({
    # DT::datatable(tableFun(), options = list(dom = "t", paging = FALSE, scrollY = "500px"), autoHideNavigation = T, rownames = F)
    DT::datatable(tableFun(), 
                  callback = JS("$('table.dataTable.no-footer').css('border-bottom', 'none');"),
                  options = list(dom = "t", paging = FALSE, scrollY = "600px"), 
                  autoHideNavigation = T, 
                  rownames = F)
  })
  
  # Adding a download button widget for the table
  output$table_download <- downloadHandler(
    filename = function(){
      paste0("TARGET_AAML1031_", gene(), "_Summary_Table_generated_", format(Sys.time(), "%m.%d.%Y"), ".xlsx")
    }, 
    content = function(file){
      write.xlsx(file = file, x = tableFun())
    }
  )

  # This will hide the ADC/CAR T action buttons if the gene isn't the target of any clinical trials
  observeEvent(gene(), {
    if (gene() %in% therapyData()$`Gene target`) {
      shinyjs::show("adc_flag")
    }
  })
  
  # The server side of the ADC/CAR T action button that will switch the user 
  # when clicked
  observeEvent(input$adc_flag, {
    newval <- "extData"
    updateTabItems(session = parent, "sdbr", selected = newval)
  })
  
  
}