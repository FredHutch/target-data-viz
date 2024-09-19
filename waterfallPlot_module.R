# UI function for the waterfall plot module
wfPlotUI <- function(id, label = "Gene expression plot parameters"){

  ns <- NS(id) # Setting a unique namespace for this module

  # Creating a list of dropdown choices for the plot type selection
  choices <- as.list(filter(colMapping, Module_Code != "Mutation" & !is.na(Final_Column_Label))$Final_Column_Name)
  names(choices) <- filter(colMapping, Module_Code != "Mutation" & !is.na(Final_Column_Label))$Final_Column_Label

  # Using tagList() instead of fluidPage() to allow for the ADC/CAR T-cell therapy button to change the tab
  tagList(
    useShinyjs(),
    fluidPage(
      theme = shinythemes::shinytheme(theme = "paper"),


            ###############################################################
            #----------------------- SIDEBAR -----------------------------#
            ###############################################################

            sidebarLayout(
              position = "left", # Placing the sidebar on the left side of the screen
              sidebarPanel(

                # Dropdown menu to select variable to use for arranging/grouping patients in waterfall plot
                selectInput(ns("grouping_var"),
                            label = "Select a grouping variable",             # The name of each list item is what is shown in the box;
                            choices = choices),                               # the value corresponds to a column of the CDEs

                radioButtons(ns("plot_type"),
                             label = "Select a type of plot to generate",
                              choices = list("Waterfall plot" = "wf",
                                             "Box/violin plots" = "bx",
                                             "Strip plot" = "str",
                                             "Scatter plot" = "sctr")),

                conditionalPanel(
                  condition = paste0("input['", ns("plot_type"), "'] == 'sctr'"),
                  textInput(ns("gene2"),
                            label = "2nd gene for comparison")),

                conditionalPanel(
                  condition = paste0("input['", ns("plot_type"), "'] == 'bx' || input['", ns("plot_type"), "'] == 'str'"),
                  checkboxInput(ns("labels"),
                                label = "Add x-axis labels",
                                value = FALSE)),
                
                conditionalPanel(
                  condition = paste0("input['", ns("plot_type"), "'] == 'bx' || input['", ns("plot_type"), "'] == 'str' || input['", ns("plot_type"), "'] == 'sctr'"),
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
                shinyjs::hidden(
                  actionButton(ns("adc_flag"),
                               label = "See targeted therapies",
                               class = "btn-primary")
                               # width = "50%")
                ),
                br(),
                br(),
                downloadButton(ns("plot_download"),
                               label = "plot",
                               class = "plotdwnld"),
                
                shinyBS::bsTooltip(ns("plot_download"),
                                   title = "Click here to download a copy of the plot",
                                   placement = "right",
                                   trigger = "hover"),
              
              br(),
              br(),
              
              # Adds a help button to identify acronyms for disease types in cohorts
              actionButton(ns("key_button"),
                           label = "Sample Key",
                           icon = icon("info-circle"),
                           class = "btn-primary",
                           disabled = FALSE),
              
              shinyBS::bsTooltip(ns("key_button"),
                                 title = "Click for a Sample Type Key",
                                 placement = "right",
                                 trigger = "hover")
              
              ),
              
              ###############################################################
              #----------------------- MAIN PLOT PANEL ---------------------#
              ###############################################################

              mainPanel(
                position = "right",
                tabsetPanel(

                  #-------------------- Waterfall plot -----------------------#
                  tabPanel("Plot",           # This is the title of the tab panel, NOT the name of the plot object!
                           br(),
                           br(),             # Linebreaks to help center the plot on the page
                           fluidRow(
                             column(11, offset = 0, align = "left",                   # This will be a reactive object that is linked to an item in the
                                    plotlyOutput(ns("plot"), height = "100%")           # output list, created in the "server" script
                                    )
                           )
                  ),

                  #-------------------- Summary table -----------------------#
                  tabPanel("Summary stats",
                           br(),
                           br(),
                           fluidRow(
                             column(12, offset = 0, align = "left",
                                    DT::dataTableOutput(ns("table")))
                           )
                  )
                )
              )
            )
  )
  )
}



# Server function for the waterfall plot module
wfPlot <- function(input, output, session, clinData, expData, adc_cart_targetData, gene, dataset, parent){

  library(ggpubr)

  # bs <- 17 # Base font size for figures

  
  
  #################################################################
  #-------------------- DATA PREPARATION -------------------------#
  #################################################################

  # Setting up a list of grouping variables that are available for each dataset.
  dropdown_choices <- filter(colMapping, Module_Code != "Mutation" & !is.na(Final_Column_Label))$Final_Column_Name
  names(dropdown_choices) <- filter(colMapping, Module_Code != "Mutation" & !is.na(Final_Column_Label))$Final_Column_Label

  # Some dropdown choices are not available for all datasets - this function will filter the options
  # depending on which dataset the user has selected.
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

  # Updating plot dropdown options based on the dataset selected by the user
  observeEvent(dataset(), {
    updateSelectInput(
      session = session,
      inputId = "grouping_var",
      choices = dropdown_choices[!dropdown_choices %in% disabled_choices()])
  }, ignoreInit = T)
  
  observeEvent(input$key_button, {
    showModal(
      modalDialog(
        title = "Sample Type Key",
        HTML(
          paste(
          "AML: Acute Myeloid Leukemia",
          "CB: Cord Blood",
          "CD34+ PB: CD34+ Peripheral Blood",
          "DS: Down Syndrome AML",
          "MPN: Myeloproliferative Neoplasm",
          "NBM: Normal Bone Marrow",
          "NBM Lymph Neg: Myeloid Sorted Normal Bone Marrow",
          "NBM Lymph Pos: Lymphoid Sorted Normal Bone Marrow",
          "TMD: Transient Myeloproliferative Disorder",
          sep = "<br>")
        ),
          easyClose = TRUE)
    )
  })

  #################################################################
  #------------------------- FUNCTIONS ---------------------------#
  #################################################################

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
    
    # Requests entry of another gene symbol ONLY when the input plot type is a scatter plot
    if (input$plot_type == "sctr" && input$gene2 == "") {
      validate("Please enter a 2nd gene symbol or miRNA in the new text box.")
    }
    
    if (input$gene2 == gene()) {
      validate("Please enter a different 2nd gene symbol.")
    }
    
    # The default value for an empty text entry box = ""
    genes2keep <- if (input$plot_type == "sctr" & input$gene2 != "") {
      c(gene(), input$gene2)
    } else if (input$plot_type != "sctr") {
      gene()
    }
    
    df <- expData() %>%
      rownames_to_column("Gene") %>%
      filter(Gene %in% genes2keep) %>%
      dplyr::select(Gene, any_of(intersect(clinData()$PatientID, colnames(expData()))))
    
        return(df)
  })

  # Transforming the counts into a long-format dataframe (to use with ggplot).
  plotData <- reactive({
    
    # Should prevent the user from selecting "Malignancy" or "Tissue" as the grouping variable for the initialized TARGET dataset 
    # Needed to add in some additional checks as it is the case that when the user has a selected filter and then switches datasets
    # it might not be available for the new dataset
    validate(
      need(!((dataset() %in% c("TARGET", "BeatAML", "SWOG", "TCGA", "StJude", "GMKF", "CCLE")) && (input$grouping_var %in% disabled_choices())), "That grouping option is not available for this dataset.\nPlease select another option."))
    
    plotDF <- geneData() %>%
      pivot_longer(names_to = "PatientID", values_to = "Expression", -Gene) %>%
      drop_na(Expression) %>% # Removing samples without expression data from the dataset
      mutate(across(Expression, ~as.numeric(.))) %>%
      left_join(., clinData(), by = "PatientID") %>%
      mutate(Log2 = log2(Expression + 1))

    # Modifying the chosen grouping variable to keep the NBMs and PBs from being categorized as NA.
    # This keeps them on the plot - if they're recategorized as NA, they will be removed from the final plot.
    plotDF <- plotDF %>%
      mutate_at(vars(any_of(!!input$grouping_var), -one_of("Age.Category")), ~case_when(Disease.Group == "NBM" ~ "NBM",
                                                                                        Disease.Group == "NBM, Lymph Neg" ~ "NBM, Lymph Neg",
                                                                                        Disease.Group == "NBM, Lymph Pos" ~ "NBM, Lymph Pos",
                                                                                        Disease.Group == "CB" ~ "CB",
                                                                                        Disease.Group == "CD34+ PB" ~ "CD34+ PB",
                                                                                        TRUE ~ .)) %>%
      mutate(across(any_of(c("MLL.Fusion", "Rare.Fusion", "Primary.Fusion", "SNVs", "Primary.CNV")), ~reLevel_cols(.)))

    # Define the order of the grouping variable levels for the TARGET cohort
    if (dataset() == "TARGET") {
      plotDF <- plotDF %>%
        mutate_at(vars(Cell.Line), ~forcats::fct_relevel(., "AML", after = Inf)) %>%
        mutate_at(vars(any_of(!!input$grouping_var), -one_of("Age.Category")), ~forcats::fct_relevel(., "MPN", after = Inf)) %>%
        mutate_at(vars(any_of(!!input$grouping_var), -one_of("Age.Category")), ~forcats::fct_relevel(., "DS", after = Inf)) %>%
        mutate_at(vars(any_of(!!input$grouping_var), -one_of("Age.Category")), ~forcats::fct_relevel(., "TMD", after = Inf)) %>%
        mutate_at(vars(any_of(!!input$grouping_var), -one_of("Age.Category")), ~forcats::fct_relevel(., "Cell line", after = Inf)) %>%
        mutate_at(vars(any_of(!!input$grouping_var), -one_of("Age.Category")), ~forcats::fct_relevel(., "CB", after = Inf)) %>%
        mutate_at(vars(any_of(!!input$grouping_var), -one_of("Age.Category")), ~forcats::fct_relevel(., "CD34+ PB", after = Inf)) %>%
        mutate_at(vars(any_of(!!input$grouping_var), -one_of("Age.Category")), ~forcats::fct_relevel(., "NBM", after = Inf)) %>%
        mutate_at(vars(any_of(!!input$grouping_var), -one_of("Age.Category")), ~forcats::fct_relevel(., "NBM, Lymph Neg", after = Inf)) %>%
        mutate_at(vars(any_of(!!input$grouping_var), -one_of("Age.Category")), ~forcats::fct_relevel(., "NBM, Lymph Pos", after = Inf))
      }

    plotDF <- plotDF %>%
      # group_by(!!input$grouping_var) %>%    # Don't think this is needed....
      arrange_(input$grouping_var, "Expression")      # Reordering patients so that the specified groups are
                                                                     # grouped together and ordered by increasing expression.
    if (input$plot_type == "wf") {
      # Setting the patient order using factor levels, so they won't be rearranged
      # alphabetically by ggplot (this step is required for waterfall plots made w/ ggplot)
      plotDF$PatientID <- factor(plotDF$PatientID, levels = plotDF$PatientID)
    }

    return(plotDF)
  })


  #----------------- Plot generation function -------------------#
  plotFun <- reactive({

    if (any(grepl(gene(), c(miRmapping$Alias, miRmapping$hsa.ID.miRbase21)))) {
      units <- "RPM"
    } else {
      units <- "TPM"
    }

    # Selecting units to display on the y-axis
    if (input$log == TRUE) {
      expCol <- "Log2"
      yaxLab <- paste0("\n", gene(), " Expression (log2 ", units, " + 1)\n") # The extra newline is to keep x-axis labels
    } else {                                                                 # from running off the side of the plot
      expCol <- "Expression"
      yaxLab <- paste0("\n", gene(), " Expression (", units, ")\n")
    }

    # Customizing the x-axis labels based on user input
    xaxLabs <- if (input$labels == TRUE) {
      element_text(hjust = 1, vjust = 1, angle = 20)
    } else {
      element_blank()
    }

    # Specifying location of the plot legend
    plotLegend <- ifelse(input$labels == TRUE, "none", "right")

    if (input$plot_type == "bx") { # Generating box plots
      p <- plotData() %>%
        drop_na(input$grouping_var) %>%
        ggplot(aes_string(x = input$grouping_var, y = expCol, fill = input$grouping_var)) +
        theme_classic(base_size = bs) +
        labs(x = NULL, y = yaxLab, fill = gsub("\\.", " ", input$grouping_var)) +
        theme(axis.text.x = xaxLabs,
              axis.text.y = element_text(size = bs + 1),
              plot.title = element_text(size = bs + 2, hjust = 0.5),
              legend.position = plotLegend,
              legend.text = element_text(size = bs - 5),
              legend.title = element_blank()) +
        geom_violin(scale = "width", aes_string(color = input$grouping_var), alpha = 0.75) +
        geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white", color = "black") +
        guides(color = "none")
        
      # Try to convert to a plotly plot with interactive tooltips
      p <- ggplotly(p, tooltip = c("y", "color"))
      
      p

    } else if (input$plot_type == "str") { # Generating strip plots w/ jittered raw data points
      p <- plotData() %>%
        drop_na(input$grouping_var) %>%
        ggplot(aes_string(x = input$grouping_var, y = expCol, fill = input$grouping_var, color = input$grouping_var)) +
        theme_classic(base_size = bs) +
        labs(x = NULL, y = yaxLab, fill = gsub("\\.", " ", input$grouping_var)) +
        theme(axis.text.x = xaxLabs,
              axis.text.y = element_text(size = bs + 1),
              plot.title = element_text(size = bs + 2, hjust = 0.5),
              legend.position = plotLegend,
              legend.text = element_text(size = bs - 5),
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
              axis.text.y = element_text(size = bs),
              plot.title = element_text(size = bs + 2, hjust = 0.5),
              axis.ticks.x = element_blank(),
              legend.position = "bottom",
              legend.text = element_text(size = bs - 5),
              legend.title = element_blank()) +
        geom_bar(stat = "identity", width = 1, position = position_dodge(width = 0.4))
      p

    } else if (input$plot_type == "sctr") { # Generating a scatter plot

      p <- plotData() %>%
        drop_na(input$grouping_var) %>%
        filter(Disease.Group == c("AML")) %>%
        ungroup() %>%
        dplyr::select(PatientID, Gene, Expression, !!sym(input$grouping_var)) %>%
        pivot_wider(names_from = "Gene", values_from = "Expression")

      if (input$log == TRUE) {
        p <- p %>%
          mutate(expCol_1 = log2(!!sym(gene()) + 1),
                 expCol_2 = log2(!!sym(input$gene2) + 1))
        xaxLab <- paste0("\n", gene(), " (log2 ", units, " + 1)")
        yaxLab <- paste0(input$gene2, " (log2 ", units, " + 1)\n")
      } else {
        xaxLab <- paste0("\n", gene(), " (", units, ")")
        yaxLab <- paste0(input$gene2, " (", units, ")\n")
        p <- p %>%
          rename(expCol_1 = !!sym(gene()),
                 expCol_2 = !!sym(input$gene2))
      }

      p <- p %>%
        ggpubr::ggscatter(.,  x = "expCol_1", y = "expCol_2",
                  cor.coef = TRUE,
                  cor.coef.size = 7,
                  cor.method = "spearman",
                  xlab = xaxLab,
                  ylab = yaxLab,
                  size = 2,
                  add = "reg.line",
                  conf.int = TRUE,
                  color = input$grouping_var,
                  add.params = list(color = "black", fill = "grey90", size = 1)) +
        guides(color = guide_legend(override.aes = list(size = 5, shape = 16))) +
        theme(axis.title = element_text(size = bs),
              axis.text = element_text(size = bs),
              legend.text = element_text(size = bs - 5),
              legend.title = element_blank(),
              legend.position = "bottom")
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
                       `N >= 5 (TPM)` = sum(Expression >= 5, na.rm = T),
                       `% >= 5 (TPM)` = round(sum(Expression >= 5, na.rm = T) / n() * 100, 2),
                       .groups = "keep")
  })

  #################################################################
  #-------------------- FINAL MODULE OUTPUTS ---------------------#
  #################################################################

  #-------------------- Plot tab -----------------------#

  # Saving the plot to the output list object so it can be run & saved reactively
  output$plot <- renderPlotly({
    plotFun()
  })

  # Adding a download button widget for the plot
  output$plot_download <- downloadHandler(
    filename = function() {
      paste0(dataset(), "_AML_", gene(), "_", input$grouping_var, "_", input$plot_type, "_generated_", format(Sys.time(), "%m.%d.%Y"), ".png")
    },
    content = function(file) {
      ggsave(filename = file, plot = plotFun(), width = 6, height = 4.5, device = "png", dpi = 250)
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
                  class = "compact nowrap hover row-border order-column", # Defines the CSS formatting of the final table, can string multiple options together
                  callback = DT::JS("$('table.dataTable.no-footer').css('border-bottom', 'none');"),
                  extensions = 'Buttons', # See https://rstudio.github.io/DT/extensions.html for more extensions & features
                  options = list(scrollY = "70vh",
                                 dom = 'Bfrtip',
                                 buttons = list(
                                   list(extend = 'excel', filename = paste0(dataset(), "_AML_", gene(), "_Summary_Table_generated_", format(Sys.time(), "%m.%d.%Y")))),
                                 scrollX = TRUE,
                                 # fixedColumns = list(leftColumns = 1),
                                 searchHighlight = TRUE,
                                 pageLength = 50),
                  escape = F) %>%
      DT::formatStyle(columns = c(1,2,4), fontSize = "100%")
  })

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
