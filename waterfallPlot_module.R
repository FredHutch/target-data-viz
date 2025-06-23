# UI function for the waterfall plot module
wfPlotUI <- function(id, label = "Gene expression plot parameters"){
  
  ns <- NS(id) # Setting a unique namespace for this module
  
  # Creating a list of dropdown choices for the plot type selection
  choices <- as.list(filter(colMapping, Module_Code != "Mutation" & !is.na(Final_Column_Label))$Final_Column_Name)
  names(choices) <- filter(colMapping, Module_Code != "Mutation" & !is.na(Final_Column_Label))$Final_Column_Label
  
  # Using tagList() instead of fluidPage() to allow for the ADC/CAR T-cell therapy button to change the tab
  
  tagList(
    useShinyjs(),
    tags$head(
      tags$style(HTML("
        .sidebar-container {
          display: flex;
          height: 100vh;
        }
        .custom-sidebar {
          background-color: #f8f9fa;
          padding: 15px;
          width: 250px;
          flex-shrink: 0;
        }
        .main-content {
          flex-grow: 1;
          padding: 15px;
        }
        /* Hovered or keyboard-selected option */
        .selectize-dropdown-content .option.active {
          background-color: #2096f6 !important; /* red */
          color: white !important;
        }
      
        /* Currently selected option when dropdown opens */
        .selectize-dropdown-content .option.selected {
          background-color: #2096f6 !important; /* red */
          color: white !important;
        }
      "))
    ),
    fluidPage(
      theme = shinythemes::shinytheme(theme = "paper"),
      div(class = "sidebar-container",
          div(class = "custom-sidebar",
              selectInput(ns("grouping_var"),
                          label = "Select a grouping variable",
                          choices = choices),
          
          radioButtons(ns("plot_type"),
                       label = "Select a type of plot to generate",
                       choices = list("Waterfall plot" = "wf",
                                      "Box plot" = "bx",
                                      "Strip plot" = "str",
                                      "Scatter plot" = "sctr")),
          
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'sctr'"),
            textInput(ns("gene2"),
                      label = "2nd gene for comparison")
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'bx' || input['", ns("plot_type"), "'] == 'str'"),
            checkboxInput(ns("labels"), "Add x-axis labels", FALSE)
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("plot_type"), "'] == 'bx' || input['", ns("plot_type"), "'] == 'str' || input['", ns("plot_type"), "'] == 'sctr'"),
            checkboxInput(ns("log"), "Log2 transform the data", FALSE)
          ),
          
          conditionalPanel(
            condition = paste0("input['", ns("test"), "'] == 1"),
            checkboxGroupInput(ns("comparisons"),
                               label = "Select 2 groups to compare",
                               choices = c("A", "B", "C", "D"))
          ),
          
          tags$hr(style = "margin: 15px 0;"),
          
          helpText("The grouping variable will be used to arrange patients along the x axis (for waterfall plots)
            or to group patients together (for box and violin plots)..."),
          
          helpText("NOTE: If cell lines is selected, please reference the 'Summary Stats' tab..."),
          
          tags$hr(style = "margin: 15px 0;"),
          
          # Targeted therapies button
          shinyjs::hidden(
            div(style = "margin-bottom: 10px;",
                actionButton(ns("adc_flag"),
                             label = "See targeted therapies",
                             class = "btn-primary btn-sm w-100"))
          ),
          
          # Plot download button
          div(style = "margin-bottom: 10px;",
              downloadButton(ns("plot_download"),
                             label = "Download Plot",
                             class = "btn-primary btn-sm w-100")),
          shinyBS::bsTooltip(ns("plot_download"),
                             title = "Click here to download the plot",
                             placement = "right",
                             trigger = "hover"),
          
          # Conditional Sample key button
          conditionalPanel(
            condition = paste0("input['", ns("grouping_var"), "'] == 'Disease.Group'"),
            div(style = "margin-bottom: 10px;",
                actionButton(ns("key_button"),
                             label = "Sample Key",
                             icon = icon("info-circle"),
                             class = "btn btn-sm w-100",
                             style = "background: #FD7370; color: #FFFFFF;")),
            shinyBS::bsTooltip(ns("key_button"),
                               title = "Click for a sample type key",
                               placement = "right",
                               trigger = "hover")
          )),
        
        
        ###############################################################
        #----------------------- MAIN PLOT PANEL ---------------------#
        ###############################################################
        
        mainPanel(
          position = "right",
          width = 10,
          tabsetPanel(
            
            #-------------------- Waterfall plot -----------------------#
            tabPanel("Plot",           # This is the title of the tab panel, NOT the name of the plot object!
                     br(),
                     br(),             # Linebreaks to help center the plot on the page
                     fluidRow(
                       column(12, offset = 0, align = "left",                   # This will be a reactive object that is linked to an item in the
                              plotlyOutput(ns("plot"), height = "70vh")           # output list, created in the "server" script
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
wfPlot <- function(input, output, session, clinData, expData, adc_cart_targetData, gene, aligner, dataset, parent){
  
  library(ggpubr)
  
  # bs <- 17 # Base font size for figures
  #print(target_id())
  #print(gene())
  
  # Making the gene2 input non-case sensitive
  observeEvent(input$gene2, {
    newValue <- toupper(input$gene2)
    updateTextInput(session, "gene2", value = newValue)
  })
  
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
  
  #Updating the options for significance testing to reflect the
  #grouping variable chosen by the user.
  # observe({
  #   x <- unique(plotData()[[input$grouping_var]])
  #   x <- x[!is.na(x)]
  #   updateCheckboxGroupInput(session,
  #                            inputId = "comparisons",
  #                            label = "Select 2 of the following to compare",
  #                            choices = x)
  # })
  
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
            "APL: Acute Promyelocytic Leukemia",
            "JMML: Juvenile Myelomonocytic Leukemia",
            "TMD: Transient Myeloproliferative Disorder",
            "Cell line: Cell line",
            "CB: Cord Blood",
            "CD34+ PB: CD34+ Peripheral Blood",
            "DS: Down Syndrome AML",
            "MPN: Myeloproliferative Neoplasm",
            "NBM: Normal Bone Marrow",
            "MSNBM: Myeloid Sorted Normal Bone Marrow",
            "LSNBM: Lymphoid Sorted Normal Bone Marrow",
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
  
  # # transforming the counts data if using the star aligner method
  # waterData <- reactive({
  #   if (dataset() == "TARGET" && aligner() == "star"){
  #     genenames <- expData()$name
  #     exp_matrix <- expData()[,c(-1:-2)]
  #     rownames(exp_matrix) <- genenames
  #   } else {
  #     exp_matrix <- expData()
  #   }
  #   return(exp_matrix)
  # })
  
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
    
    ## If the alignment method is kallisto or another cohort
    genes2keep <- if (input$plot_type == "sctr" & input$gene2 != "") { ## The default value for an empty text entry box = ""
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
      need(!((dataset() %in% c("TARGET", "BeatAML", "SWOG", "TCGA", "StJude", "GMKF", "CCLE", "LEUCEGENE", "PCGP AML", "PCGP")) && (input$grouping_var %in% disabled_choices())), "That grouping option is not available for this dataset.\nPlease select another option."))
    
    plotDF <- geneData() %>%
      pivot_longer(names_to = "PatientID", values_to = "Expression", -Gene) %>%
      drop_na(Expression) %>% # Removing samples without expression data from the dataset
      mutate(across(Expression, ~as.numeric(.))) %>%
      left_join(., clinData(), by = "PatientID") %>%
      mutate(Log2 = log2(Expression + 1))
    
    # Modifying the chosen grouping variable to keep the NBMs and PBs from being categorized as NA.
    # This keeps them on the plot - if they're recategorized as NA, they will be removed from the final plot.
    plotDF <- plotDF %>%
      mutate_at(vars(any_of(!!input$grouping_var), -one_of("Age.Category", "Cell.Line")), ~case_when(Disease.Group == "NBM" ~ "NBM",
                                                                                        Disease.Group == "MSNBM" ~ "MSNBM",
                                                                                        Disease.Group == "LSNBM" ~ "LSNBM",
                                                                                        Disease.Group == "CD34+ PB" ~ "CD34+ PB",
                                                                                        TRUE ~ .)) %>%
      mutate(across(any_of(c("MLL.Fusion", "Rare.Fusion", "Primary.Fusion", "SNVs", "Primary.CNV")), ~reLevel_cols(.)))
  
    normals <- c("CD34+ PB", "NBM", "MSNBM", "LSNBM")
    
    if (!(input$grouping_var %in% c("Cell.Line", "Age.Category", "Risk"))) {
      plotDF[[input$grouping_var]] <- forcats::fct_relevel(
        as.factor(plotDF[[input$grouping_var]]),
        normals,
        after = Inf
      )
    }
    
    if (dataset() == "TARGET") {
      plotDF$Age.Category <- forcats::fct_relevel(plotDF$Age.Category, c("Less than 3 years", "Between 3 and 5 years", "Between 5 and 10 years", "Between 10 and 18 years", "Greater than 18 years", "Unknown"))
      plotDF$Risk <- forcats::fct_relevel(plotDF$Risk, c("High", "Standard", "Low", normals))
    }
    if (dataset() == "BeatAML") {
      plotDF$Age.Category <- forcats::fct_relevel(plotDF$Age.Category, c("Less than 10 years", "Between 10 and 18 years","Between 18 and 40 years", "Between 40 and 60 years", "Greater than 60 years"))
      plotDF$Risk <- forcats::fct_relevel(plotDF$Risk, c("Adverse", "Intermediate Or Adverse", "Intermediate", "Favorable Or Intermediate", "Favorable", normals))
    }
    if (dataset() == "SWOG") {
      plotDF$Age.Category <- forcats::fct_relevel(plotDF$Age.Category, c("Between 18 and 40 years", "Between 40 and 60 years", "Greater than 60 years", "Unknown"))
      plotDF$Cytogenetic.Category <- forcats::fct_relevel(plotDF$Cytogenetic.Category, c("FAV", "UNF", "UNK", normals))
    }
    if (dataset() == "TCGA") {
      print(table(plotDF$Age.Category))
      plotDF$Age.Category <- forcats::fct_relevel(plotDF$Age.Category, c("Between 18 and 40 years", "Between 40 and 60 years", "Greater than 60 years"))
      plotDF$Risk <- forcats::fct_relevel(plotDF$Risk, c("Poor", "Intermediate", "Good", "N.D.", normals))
    }
    if (dataset() == "PCGP AML") {
      plotDF$Age.Category <- forcats::fct_relevel(plotDF$Age.Category, c("Less than 3 years", "Between 3 and 5 years", "Between 5 and 10 years", "Between 10 and 18 years", "Greater than 18 years", "Unknown"))
    }
    if (dataset() == "PCGP") {
      plotDF$Age.Category <- forcats::fct_relevel(plotDF$Age.Category, c("Less than 3 years", "Between 3 and 5 years", "Between 5 and 10 years", "Between 10 and 18 years", "Greater than 18 years", "Unknown"))
    }

      # Setting the patient order using factor levels, so they won't be rearranged
      # alphabetically by ggplot (this step is required for waterfall plots made w/ ggplot)
    if (input$grouping_var == "SNVs") {
        # Mutation columns of interest
        mutation_cols <- c("WT1.Mutation", "NPM1.Mutation", "FLT3.ITD", "CEBPA.Mutation")
        
        # Remove existing SNVs column if present
        if ("SNVs" %in% colnames(plotDF)) {
          plotDF <- plotDF %>% dplyr::select(-SNVs)
        }
        
        # Reshape data to long format for mutations
        long_df <- plotDF %>%
          pivot_longer(cols = all_of(mutation_cols), names_to = "SNVs", values_to = "Value") %>%
          filter(Value == "Yes") %>%
          dplyr::select(-Value)
        
        # Identify patients with mutations
        patients_with_mutations <- unique(long_df$PatientID)
        
        # Identify control samples to always include
        control_samples <- plotDF %>%
          filter(Disease.Group %in% c("CD34+ PB", "NBM", "MSNBM", "LSNBM"))
        
        # Create entries for control samples that should have "No Mutation"
        control_samples_no_mutation <- control_samples %>%
          filter(!PatientID %in% patients_with_mutations) %>%
          mutate(SNVs = Disease.Group) %>%
          dplyr::select(PatientID, SNVs, Expression, Log2)
        
        # Create entries for other patients without mutations
        other_no_mutation_df <- plotDF %>%
          filter(!PatientID %in% patients_with_mutations,
                 !PatientID %in% control_samples$PatientID) %>%
          mutate(SNVs = NA) %>%
          dplyr::select(PatientID, SNVs, Expression, Log2)
        
        # Combine all data: mutation data, control samples without mutation, and other patients without mutation
        plotDF <- bind_rows(long_df, control_samples_no_mutation, other_no_mutation_df)
        
        # Clean up SNVs column to remove ".Mutation" from mutation column names
        plotDF$SNVs <- sub("\\.Mutation$", "", plotDF$SNVs)
        plotDF$SNVs <- fct_relevel(plotDF$SNVs, c("CEBPA", "FLT3.ITD", "NPM1", "WT1", "CD34+ PB", "NBM", "MSNBM", "LSNBM"))
      }

  if (input$plot_type == "wf") {
    plotDF <- plotDF %>%
      arrange(!!sym(input$grouping_var), Expression) %>%
      mutate(SampleID = row_number())  # unique ID per row
    
    # Now use SampleID for factoring instead of PatientID
    plotDF$SampleID <- factor(plotDF$SampleID, levels = plotDF$SampleID)
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
      element_text(hjust = 1, vjust = 1, angle = 90)
    } else {
      element_blank()
    }
    
    # Removes the plot legend if you turn on x-axis labels
    plotLegend <- ifelse(input$labels == TRUE, FALSE, TRUE)
    
    if (input$plot_type == "bx") { # Generating box plots
      # Check data before plotting
      data_check <- plotData() %>%
        drop_na(input$grouping_var)  # Make sure grouping var and expCol are valid
      
      # Check if grouping and expCol have sufficient variation
      print(summary(data_check[[input$grouping_var]]))
      print(summary(data_check[[expCol]]))  # Ensure you reference expCol correctly here as a column
      
      # Dynamically refer to the expCol column using tidy eval (!!sym())
      p <- data_check %>%
        ggplot(aes_string(x = input$grouping_var, y = expCol, fill = input$grouping_var)) +
        theme_classic(base_size = bs, base_family = "Helvetica") +
        labs(x = NULL, y = yaxLab, fill = gsub("\\.", " ", input$grouping_var)) +
        theme(axis.title.y = element_text(size = bs),
              axis.text.x = xaxLabs,
              axis.text.y = element_text(size = bs),
              plot.title = element_text(size = bs + 2, hjust = 0.5),
              legend.position = "bottom",
              legend.text = element_text(size = bs - 4),
              legend.title = element_blank()) +
        geom_boxplot(width = 0.2, outlier.shape = NA, aes_string(fill = input$grouping_var), color = "black") + 
        guides(color = "none")
      
      # Check the plot without plotly conversion
      print(p)
      
      # Try to convert to a plotly plot with interactive tooltips
      p <- ggplotly(p, tooltip = c("y", "color"), dynamicTicks = TRUE)
      
      # Adjust legend position in plotly
      p <- layout(p,
                  showlegend = plotLegend,
                  legend = list(orientation = "h", 
                                y = -0.1, 
                                x = 0.5, 
                                xanchor = "center"
                  ),
                  xaxis = list(tickfont = list(size = bs)))
      
      p
    } else if (input$plot_type == "str") {
      plot_df <- plotData() %>%
        drop_na(!!sym(input$grouping_var)) %>%
        mutate(
          grouping_var = as.factor(.data[[input$grouping_var]]),
          tooltip = paste0(
            "PatientID: ", PatientID, "\n",
            "Expression: ", .data[[expCol]], " TPM\n",
            "Group: ", .data[[input$grouping_var]]
          )
        )
      
      p <- ggplot(plot_df, aes(
        x = grouping_var,
        y = !!sym(expCol),
        color = grouping_var,
        text = tooltip
      )) +
        geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.75) +
        
        # Add one black median line per group
        stat_summary(
          fun = median,
          geom = "crossbar",
          width = 0.5,
          color = "black",
          fatten = 0.8,
          aes(group = grouping_var)
        ) +
        
        labs(x = NULL, y = yaxLab, color = gsub("\\.", " ", input$grouping_var)) +
        theme_classic(base_size = bs, base_family = "Helvetica") +
        theme(
          axis.title.y = element_text(size = bs),
          axis.text.x = xaxLabs,
          axis.text.y = element_text(size = bs),
          plot.title = element_text(size = bs + 2, hjust = 0.5),
          legend.position = "bottom",
          legend.text = element_text(size = bs - 4),
          legend.title = element_blank()
        )
      
      p <- ggplotly(p, tooltip = "text")
      
      # Loop over traces and update only those with a legend entry
      for (i in seq_along(p$x$data)) {
        if (!is.null(p$x$data[[i]]$name) && p$x$data[[i]]$mode == "markers") {
          p$x$data[[i]]$marker$size <- 6  # Adjust this number to increase legend dot size
          p$x$data[[i]]$showlegend <- TRUE
        }
      }
      
      # Adjust the legend layout if needed
      p <- layout(
        p,
        showlegend = plotLegend,
        legend = list(
          orientation = "h",
          y = -0.1,
          x = 0.5,
          xanchor = "center",
          itemsizing = "constant"
        )
      )
      
      p
      
    } else if (input$plot_type == "wf") { # Generating a waterfall plot
      
          if (dataset() == "PCGP") {
            p <- plotData() %>%
              drop_na(input$grouping_var) %>%
              ggplot(aes_string(
                x = "SampleID",
                y = "Expression",
                fill = input$grouping_var
              )) +
              geom_bar(
                stat = "identity",
                width = 1,
                position = position_dodge(width = 0.4),
                inherit.aes = TRUE
              ) +
              aes(text = paste0("PatientID: ", PatientID,
                                "<br>Expression: ", round(Expression, 3), " TPM",
                                "<br>", gsub("\\.", "", input$grouping_var), ": ", get(input$grouping_var),
                                "<br>Subtype: ", Subtype)
              ) +
              theme_classic(base_size = bs, base_family = "Helvetica") +
              labs(
                x = NULL,
                y = yaxLab,
                fill = gsub("\\.", " ", input$grouping_var)
              ) +
              theme(
                axis.text.x = element_blank(),
                axis.text.y = element_text(size = bs),
                plot.title = element_text(size = bs + 2, hjust = 0.5),
                axis.ticks.x = element_blank(),
                legend.position = "bottom",
                legend.text = element_text(size = bs - 4),
                legend.title = element_blank()
              )
          } else {
      
                p <- plotData() %>%
                  drop_na(input$grouping_var) %>%
                  ggplot(aes_string(
                    x = "SampleID",
                    y = "Expression",
                    fill = input$grouping_var
                  )) +
                  geom_bar(
                    stat = "identity",
                    width = 1,
                    position = position_dodge(width = 0.4),
                    inherit.aes = TRUE
                  ) +
                  aes(text = paste0("PatientID: ", PatientID,
                                    "<br>Expression: ", round(Expression, 3), " TPM",
                                    "<br>", gsub("\\.", "", input$grouping_var), ": ", get(input$grouping_var)
                  )) +
                  theme_classic(base_size = bs, base_family = "Helvetica") +
                  labs(
                    x = NULL,
                    y = yaxLab,
                    fill = gsub("\\.", " ", input$grouping_var)
                  ) +
                  theme(
                    axis.text.x = element_blank(),
                    axis.text.y = element_text(size = bs),
                    plot.title = element_text(size = bs + 2, hjust = 0.5),
                    axis.ticks.x = element_blank(),
                    legend.position = "bottom",
                    legend.text = element_text(size = bs - 4),
                    legend.title = element_blank()
                  )
          }
      
      # Convert to plotly with custom tooltip
      p <- ggplotly(p, tooltip = "text", dynamicTicks = TRUE)
      # Update layout and styling using R-style plotly syntax
      p <- layout(p,
                  legend = list(
                    orientation = "h",
                    y = -0.1,
                    x = 0.5,
                    xanchor = "center"
                  ),
                  bargap = 0
      )
      
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
      
      # Calculate correlation and p-value
      correlation <- cor.test(p$expCol_1, p$expCol_2, method = "pearson")
      p_value <- round(correlation$p.value, 6)
      
      p <- p %>%
        ggpubr::ggscatter(.,  x = "expCol_1", y = "expCol_2",
                          cor.coef = FALSE,
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
              legend.text = element_text(size = bs - 4),
              legend.title = element_blank(),
              legend.position = "bottom")
      
      p <- ggplotly(p, dynamicTicks = TRUE)
      
      # Add p-value annotation
      p <- p %>% 
        layout(annotations = list(
          x = 1,
          y = 1,
          text = paste("p-value =", p_value),
          xref = "paper",
          yref = "paper",
          showarrow = FALSE,
          font = list(size = 22)
        ))
      
      # Adjust legend position in plotly
      p <- layout(p, 
                  legend = list(orientation = "h", 
                                y = -0.2, 
                                x = 0.5, 
                                xanchor = "center"
                  ))
      p
      
    }
    
    ### Perform Statistical Comparison ###
  #   if (length(input$comparisons) > 1) {
  #     validate(
  #       need(length(input$comparisons > 1), "Please select 2 groups to compare."))
  #     # c <- ggpubr::stat_compare_means(method = "wilcox.test", comparisons = list(input$comparisons)) ### This no longer works when using ggplotly (fails to recognize added geom)
  #     
  #     ### Manually generate the comparison and add to plot ###
  #     
  #     # Sort plot data for the comparitors of interest
  #     df <- plotData() %>% 
  #       drop_na(input$grouping_var) %>%
  #       filter(!!as.name(input$grouping_var) %in% input$comparisons) %>%
  #       ungroup()
  #     
  #     # Create the formula
  #     grouping_var <- input$grouping_var
  #     formula <- as.formula(paste0("Expression ~ ", grouping_var))
  #     test <- compare_means(formula = formula, data = df, method = "wilcox.test")
  #     
  #     p_value <- test$p
  #     
  #     # Add the result from the statistical test as an annotation to the plot
  #     y = ifelse(input$log == TRUE, max(plotData()$Log2) + 1, max(plotData()$Expression) + 1)
  #     
  #     p <- p %>% add_annotations(
  #       text = paste("p-value:", format(p_value, digits = 3)),
  #       x = 1.5, y = y,
  #       showarrow = FALSE
  #     )
  #     
  #   } else {
  #     p # Return the plot as-is with no additional geom layers
  #   }
  })
  
  #----------------- Summary table function -------------------#
  # Function to generate an expression summary table from the plot data
  
  ### Needs to be treated a little differently for cell line data, we don't want to group by Disease.Group
  ### Instead want to plot all cell lines individually for sorting. E.g., don't use Disease.Group as a grouping variable.
  
  # Function to generate an expression summary table from the plot data
  tableFun <- reactive({
    data <- plotData() %>%
      drop_na(input$grouping_var)
    
    grouped_data <- if (dataset() == "CCLE" && input$grouping_var == "Disease.Group") {
      group_by(data, Name)
    } else {
      group_by(data, !!as.name(input$grouping_var))
    }
    
    summarized_data <- grouped_data %>%
      dplyr::summarize(
        N = n(),
        Gene = gene(),
        `Mean (TPM)` = round(mean(Expression, na.rm = TRUE), 2),
        `Median (TPM)` = round(median(Expression, na.rm = TRUE), 2),
        `Range (TPM)` = paste0(round(min(Expression), 2), " - ", round(max(Expression), 2)),
        `N >= 5 (TPM)` = sum(Expression >= 5, na.rm = TRUE),
        `% >= 5 (TPM)` = round(sum(Expression >= 5, na.rm = TRUE) / n() * 100, 2),
        .groups = "keep"
      )
    
    if (dataset() == "PCGP" && "Subtype" %in% colnames(data)) {
      summarized_data <- left_join(summarized_data, data %>%
                                     group_by(!!as.name(input$grouping_var)) %>%
                                     summarize(Subtype = paste(unique(Subtype), collapse = ", "), .groups = "drop"), 
                                   by = input$grouping_var)
    }
    
    summarized_data
  })
  
  
  
  
  #################################################################
  #-------------------- FINAL MODULE OUTPUTS ---------------------#
  #################################################################
  
  #-------------------- Plot tab -----------------------#
  
  # Saving the plot to the output list object so it can be run & saved reactively
  output$plot <- renderPlotly({
    
    # Create the plot
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
                                 searchHighlight = TRUE,
                                 pageLength = 5000
                  ),
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
