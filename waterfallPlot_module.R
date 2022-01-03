# UI function for the waterfall plot module
wfPlotUI <- function(id, label = "Gene expression plot parameters"){
  
  library(DT)
  library(shinyjs)
  library(shinyWidgets)
  ns <- NS(id) # Setting a unique namespace for this module
  
  # Creating a list of dropdown choices for the plot type selection
  choices <- as.list(filter(colMapping, Module_Code != "Mutation" & !is.na(Final_Column_Label))$Final_Column_Name)
  names(choices) <- filter(colMapping, Module_Code != "Mutation" & !is.na(Final_Column_Label))$Final_Column_Label
  # https://github.com/dreamRs/shinyWidgets <- Give this a shot for some of the new checkbox options
  # https://www.htmlwidgets.org/showcase_d3heatmap.html <- Interactive heatmaps
  # https://cran.r-project.org/web/packages/shinyjqui/vignettes/introduction.html <- Make items resizable!!
  
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
                pickerInput(ns("grouping_var"), 
                            label = "Select a grouping variable",             # The name of each list item is what is shown in the box;
                            choices = choices),                               # the value corresponds to a column of the CDEs
                          
                radioButtons(ns("plot_type"), 
                             label = "Select a type of plot to generate", 
                              choices = list("Waterfall plot" = "wf", 
                                             "Box/violin plots" = "bx",
                                             "Strip plot" = "str")),
                
                conditionalPanel(
                  condition = paste0("input['", ns("plot_type"), "'] == 'bx' || input['", ns("plot_type"), "'] == 'str'"),
                  checkboxInput(ns("labels"),                                                  
                                label = "Add x-axis labels",
                                value = FALSE)),
                
                conditionalPanel(
                  condition = paste0("input['", ns("plot_type"), "'] == 'bx' || input['", ns("plot_type"), "'] == 'str'"),
                  checkboxInput(ns("log"),                                                  
                            label = "Log2 transform the data",
                            value = FALSE)),
                
                conditionalPanel(
                  condition = paste0("input['", ns("plot_type"), "'] == 'bx' || input['", ns("plot_type"), "'] == 'str'"),
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
                          to help highlight expression patterns within the groups."),
                
                helpText("NOTE: If cell lines is selected, please reference the 'Summary Stats' tab for expression data, 
                         as only one sample is available for most cell lines."),
                
                br(),
                
                # http://timelyportfolio.github.io/buildingwidgets/week25/sweetalert_examples.html <- Alert to notify
                # user that they're switching tabs
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
                  tabPanel("Plot",           # This is the title of the tab panel, NOT the name of the plot object!
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
wfPlot <- function(input, output, session, clinData, expData, adc_cart_targetData, gene, dataset, parent){
  
  library(tidyverse)
  library(DT)
  library(shinyWidgets)
  
  bs <- 17 # Base font size for figures
  
  #################################################################
  #------------------------- FUNCTIONS ---------------------------#
  #################################################################
  
  # See https://shiny.rstudio.com/articles/validation.html for details on the %then% operator
  # `%then%` <- shiny:::`%OR%` 
  
  # Version above seems to be deprecated, use below code instead:
  `%then%` <- function(a, b) {
    if (is.null(a)) b else a
  }
  
  reLevel_cols <- function(col) {
    new_col <- col
    if (any(grepl("Other AML", col))) {
      new_col <- forcats::fct_relevel(new_col, "Other AML", after = Inf)
    }
    if (any(grepl("No Relevant CNV", col))) {
      new_col <- forcats::fct_relevel(col, "No Relevant CNV", after = Inf)
    }
    return(new_col)
  }
  
  #----------------- Plot generation function -------------------#
  plotFun <- reactive({ 
    
    # Selecting units to display on the y-axis
    if (input$log == TRUE) {
      expCol <- "Log2"
      yaxLab <- paste0("\n", gene(), " Expression (log2 TPM + 1)") # The extra newline is to keep x-axis labels 
    } else {                                                       # from running off the side of the plot
      expCol <- "Expression"
      yaxLab <- paste0("\n", gene(), " Expression (TPM)")
    }
    
    # Customizing the x-axis labels based on user input
    xaxLabs <- if (input$labels == TRUE) {
      element_text(hjust = 1, vjust = 1, angle = 20) 
    } else {
      element_blank()
    }
    
    # Specifying location of the plot legend
    plotLegend <- ifelse(input$labels == TRUE, "none", "bottom")
    
    if (input$plot_type == "bx") { # Generating box plots
      p <- plotData() %>% 
        drop_na(input$grouping_var) %>%
        ggplot(aes_string(x = input$grouping_var, y = expCol, fill = input$grouping_var)) +
        theme_classic(base_size = bs) +
        labs(x = NULL, y = yaxLab, fill = gsub("\\.", " ", input$grouping_var)) +
        theme(axis.text.x = xaxLabs,
              plot.title = element_text(size = bs + 2, hjust = 0.5),
              legend.position = plotLegend,
              legend.text = element_text(size = bs - 6),
              legend.title = element_blank()) +
        geom_violin(scale = "width", aes_string(color = input$grouping_var)) +
        geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.4) +
        guides(color = "none")
      p
      
    } else if (input$plot_type == "str") { # Generating strip plots w/ jittered raw data points
      p <- plotData() %>% 
        drop_na(input$grouping_var) %>%
        ggplot(aes_string(x = input$grouping_var, y = expCol, fill = input$grouping_var, color = input$grouping_var)) +
        theme_classic(base_size = bs) +
        labs(x = NULL, y = yaxLab, fill = gsub("\\.", " ", input$grouping_var)) +
        theme(axis.text.x = xaxLabs,
              plot.title = element_text(size = bs + 2, hjust = 0.5),
              legend.position = plotLegend,
              legend.text = element_text(size = bs - 6),
              legend.title = element_blank()) +
        guides(color = "none") +
        geom_jitter(width = 0.3, size = 0.7) +
        stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black")
      p
      
    } else if (input$plot_type == "wf") { # Generating a waterfall plot
      p <- plotData() %>% 
        drop_na(input$grouping_var) %>%
        ggplot(aes_string(x = "PatientID", y = "Expression", fill = input$grouping_var)) +
        theme_classic(base_size = bs) +
        labs(x = "Patients", y = yaxLab, fill = gsub("\\.", " ", input$grouping_var)) +
        theme(axis.text.x = element_blank(),
              plot.title = element_text(size = bs + 2, hjust = 0.5),
              axis.ticks = element_blank(),
              legend.position = "bottom",
              legend.text = element_text(size = bs - 6),
              legend.title = element_blank()) +
        geom_bar(stat = "identity", width = 1, position = position_dodge(width = 0.4))
      p
    }
    
    # Performing hypothesis tests 
    if (length(input$comparisons) > 1) {
      validate(
        need(length(input$comparisons > 1), "Please select 2 groups to compare.")) # !!!!! Need to add another error message for the cell lines specifically !!!!!!
      c <- p + ggpubr::stat_compare_means(method = "wilcox.test", comparisons = list(input$comparisons))
      c # Return plot generated above + additional geom layer created by stat_compare_means
    } else {
      p # Return the plot as-is with no additional geom layers
    }
  })
  
  #----------------- Summary table function -------------------#
  # Function to generate an expression summary table from the plot data
  tableFun <- reactive({
    plotData() %>%
      drop_na(input$grouping_var) %>%
      group_by(!!as.name(input$grouping_var)) %>%
      dplyr::summarize(N = n(), 
                       Gene = gene(),
                       `Mean (TPM)` = round(mean(Expression, na.rm = T), 2), 
                       `Median (TPM)` = round(median(Expression, na.rm = T), 2), 
                       `Range (TPM)` = paste0(round(min(Expression), 2), " - ", round(max(Expression), 2)), 
                       .groups = "keep")
  })
  
  
  
  #################################################################
  #-------------------- DATA PREPARATION -------------------------#
  #################################################################

  # Setting up a list of grouping variables that are available for each dataset.
  dropdown_choices <- filter(colMapping, Module_Code != "Mutation" & !is.na(Final_Column_Label))$Final_Column_Name
  names(dropdown_choices) <- filter(colMapping, Module_Code != "Mutation" & !is.na(Final_Column_Label))$Final_Column_Label
  
  # Some dropdown choices are not available for all datasets - this function will filter the options
  # depending on which dataset the user has selected. 
  # NOTE: This isn't working the way I had hoped, but it'll do for now (doesn't change for each dataset)
  disabled_choices <- reactive({
    x <- filter(colMapping, Module_Code != "Mutation" & is.na(!!sym(dataset())))$Final_Column_Name
    names(x) <- filter(colMapping, Module_Code != "Mutation" & is.na(!!sym(dataset())))$Final_Column_Label
    return(x)
  })
  
  # Updating the options for significance testing to reflect the
  # grouping variable chosen by the user.
  observe({
    x <- unique(plotData()[[input$grouping_var]])
    x <- x[!is.na(x)]
    updateCheckboxGroupInput(session, 
                             inputId = "comparisons", 
                             label = "Select 2 of the following to compare",
                             choices = x)
  })
  
  # Updating plot dropdown options based on the dataset selected by the user.
  observeEvent(dataset(), {
    if (dataset() == "TARGET") {
      updatePickerInput(
        session = session, 
        inputId = "grouping_var",
        choices = dropdown_choices)
    } else {
      updatePickerInput(
        session = session, 
        inputId = "grouping_var",
        choices = dropdown_choices[!dropdown_choices %in% disabled_choices()])
    }
  }, ignoreInit = T)
  
  # Filtering the ADC & CAR T-cell therapy data to select therapies targeting the gene of interest.
  # This will be used for the conditional button that links to the therapeutic database tab.
  therapyData <- reactive({
    filter(adc_cart_targetData, `Gene target` == gene())
  })
  
  # Filtering the counts data to only retain the gene of interest.
  # An error will be thrown if non-existent or unrecognized gene is provided.
  geneData <- reactive({
    validate(
      need(gene(), "Please enter a gene symbol or miRNA in the text box to the left.") %then%
        need(gene() %in% rownames(expData()), paste0(gene(), " does not exist in the counts data!\nDouble-check the symbol or ID, or try an alias/synonym."))
      )
    
    df <- expData() %>%
      rownames_to_column("Gene") %>%
      filter(Gene == gene()) %>%
      dplyr::select(Gene, any_of(intersect(clinData()$PatientID, colnames(expData())))) %>%
      column_to_rownames("Gene")
    
    return(df)
  })
  
  # Transforming the counts into a long-format dataframe (to use with ggplot).
  plotData <- reactive({
    
    validate(
      need(!((dataset() %in% c("BeatAML", "SWOG", "TCGA")) && (input$grouping_var %in% disabled_choices())), "That grouping option is not available for this dataset.\nPlease select another option."))
    
    plotDF <- geneData() %>%
      gather(PatientID, Expression) %>%
      mutate_at(vars(Expression), ~as.numeric(.)) %>%
      left_join(., clinData(), by = "PatientID") %>%
      mutate(Log2 = log2(Expression + 1))
    
    # Modifying the chosen grouping variable to keep the NBMs and PBs from being categorized as NA.
    # This keeps them on the plot - if they're recategorized as NA, they will be removed from the final plot.
    plotDF <- plotDF %>%  
      mutate_at(vars(any_of(!!input$grouping_var), -one_of("Age.Category")), ~case_when(Disease.Group == "NBM" ~ "NBM", 
                                                                                        Disease.Group == "CD34+ PB" ~ "CD34+ PB",
                                                                                        TRUE ~ .)) %>%
      mutate(across(any_of(c("MLL.Fusion", "Rare.Fusion", "Primary.Fusion", "SNVs", "Primary.CNV")), ~reLevel_cols(.)))
    
    if (dataset() == "TARGET") {
      plotDF <- plotDF %>% 
        mutate_at(vars(Cell.Line), ~forcats::fct_relevel(., "AML", after = Inf)) %>%
        mutate_at(vars(any_of(!!input$grouping_var), -one_of("Age.Category")), ~forcats::fct_relevel(., "CD34+ PB", after = Inf)) %>%
        mutate_at(vars(any_of(!!input$grouping_var), -one_of("Age.Category")), ~forcats::fct_relevel(., "NBM", after = Inf))
    }
    
    plotDF <- plotDF %>%
      drop_na(Expression) %>% # Removing samples without expression data from the plot
      group_by(!!input$grouping_var) %>%                                                              
      arrange_(input$grouping_var, "Expression")  # Reordering patients so that the specified groups are 
                                                  # grouped together and ordered by increasing expression

    # Setting the patient order using factor levels, so they won't be rearranged 
    # alphabetically by ggplot (this step is required for waterfall plots made w/ ggplot)
    plotDF$PatientID <- factor(plotDF$PatientID, levels = plotDF$PatientID)
    
    return(plotDF)
  })

  
  
  #################################################################
  #-------------------- FINAL MODULE OUTPUTS ---------------------#
  #################################################################
  
  #-------------------- Plot tab -----------------------#
  
  # Saving the plot to the output list object so it can be run & saved reactively
  output$plot <- renderPlot({
    plotFun()
  })
  
  # Adding a download button widget for the plot
  output$plot_download <- downloadHandler(
    filename = function() {
      paste0(dataset(), "_AML_", gene(), "_", input$grouping_var, "_", input$plot_type, "_generated_", format(Sys.time(), "%m.%d.%Y"), ".png")
    }, 
    content = function(file) {
      ggsave(filename = file, plot = plotFun(), width = 5, height = 4.5, device = "png", dpi = 250)
    }
  )
  
  output$ggplot_download <- downloadHandler(
    filename = function() {
      paste0(dataset(), "_AML_", gene(), "_", input$grouping_var, "_ggplotObject_generated_", format(Sys.time(), "%m.%d.%Y"), ".RDS")
    },
    content = function(file) {
      withProgress(message = "Preparing RDS file", detail = "This may take a while...", value = 0, {
        for (i in 1:50) {
          incProgress(1/70)
          Sys.sleep(0.25)
        }
      })
      saveRDS(object = plotFun(), file = file, compress = F)
    }
  )
  
  #-------------------- Data tab -----------------------#
  
  # https://glin.github.io/reactable/articles/examples.html#conditional-styling
  output$table <- DT::renderDataTable({
    DT::datatable(tableFun(), 
                  callback = JS("$('table.dataTable.no-footer').css('border-bottom', 'none');"),
                  options = list(dom = "t", paging = FALSE, scrollY = "600px"), 
                  rownames = F)
  })
  
  # Adding a download button widget for the table
  output$table_download <- downloadHandler(
    filename = function(){
      paste0(dataset(), "_AML_", gene(), "_Summary_Table_generated_", format(Sys.time(), "%m.%d.%Y"), ".xlsx")
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