# UI for the Kaplan-Meier plot module
kmPlotUI <- function(id, label = "Kaplan-Meier plot parameters") {
  ns <- NS(id)
  
  mut_choices <- as.list(filter(colMapping, Module_Code == "Mutation" & !is.na(Final_Column_Label))$Final_Column_Name)
  names(mut_choices) <- filter(colMapping, Module_Code == "Mutation" & !is.na(Final_Column_Label))$Final_Column_Label
  
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
          width: 280px;
          flex-shrink: 0;
        }
        .main-content {
          flex-grow: 1;
          padding: 15px;
          min-width: 0;    
          max-width: 100vw;  
          overflow-x: auto; 
        }
        .selectize-dropdown-content .option.active {
          background-color: #2096f6 !important;
          color: white !important;
        }
        .selectize-dropdown-content .option.selected {
          background-color: #2096f6 !important;
          color: white !important;
        }
      "))
    ),
    fluidPage(
      theme = shinythemes::shinytheme("paper"),
      div(class = "sidebar-container",
          
          # Sidebar Panel
          div(class = "custom-sidebar",
              checkboxGroupInput(ns("test_type"), 
                                 label = "Select a survival type to model", 
                                 choices = c("Event-Free Survival (EFS)" = "EFS",
                                                "Overall Survival (OS)" = "OS", 
                                                "Disease-Free Survival (DFS)" = "DFS", 
                                                "Relapse Risk (RR)" = "RR")),
              
              radioButtons(ns("time_type"), 
                           label = "Select a time interval to display", 
                           choices = list("Years" = "years", "Days" = "days")),
              
              radioButtons(ns("filter_cohort"), 
                           label = "Limit to a single subgroup?", 
                           choices = list("No - view all patients" = "no", "Yes" = "yes")),
              
              conditionalPanel(
                condition = paste0("input['", ns("filter_cohort"), "'] == 'yes'"),
                selectInput(ns("select_subgroup"), label = "Which one?",
                            choices = list(
                              "KMT2A/MLL rearranged" = "MLL|KMT2A-",
                              "inv(16)" = "inv\\(16\\)|CBFB-MYH11",
                              "t(8;21)" = "t\\(8\\;21\\)|RUNX1-RUNX1T1",
                              "DEK-NUP214" = "DEK-NUP214",
                              "ETS-Family" = "ETS-Family",
                              "MLLT10-X" = "MLLT10-X",
                              "FLT3-ITD" = "FLT3-ITD",
                              "KMT2A-PTD" = "KMT2A-PTD",
                              "WT1" = "WT1",
                              "NPM1" = "NPM1(?!\\-)",
                              "CEBPA" = "CEBPA",
                              "CBFA2T3-GLIS2" = "CBFA2T3\\-GLIS2",
                              "NUP98 fusions" = "NUP98-",
                              "NUP98-NSD1" = "NUP98-NSD1",
                              "NUP98-KDM5A" = "NUP98-KDM5A",
                              "Monosomy 7" = "[Mm]onosomy7|Monosomy 7",
                              "del5q/del7q" = "del5q|del7q",
                              "Trisomy 8" = "[Tt]risomy8",
                              "Normal karyotype" = "Normal"))
              ),
              
              selectInput(ns("strata_var"), 
                          label = "Select a method of grouping the patients", 
                          choices = list("By median" = "median", 
                                         "By quartile" = "quartile", 
                                         "By percentile" = "percentile", 
                                         "By TPM cutoff" = "tpm_value",
                                         "By mutation status" = "mutation")),
              
              conditionalPanel(
                condition = paste0("input['", ns("strata_var"), "'] == 'percentile'"),
                textInput(ns("perc_cutoff"), label = "Cutoff percentile", placeholder = "Example: 75")
              ),
              conditionalPanel(
                condition = paste0("input['", ns("strata_var"), "'] == 'tpm_value'"),
                textInput(ns("tpm_cutoff"), label = "TPM cutoff value", placeholder = "Example: 5")
              ),
              conditionalPanel(
                condition = paste0("input['", ns("strata_var"), "'] == 'mutation'"),
                selectInput(ns("mutCol"), label = "Which mutation?", choices = mut_choices)
              ),
              
              tags$hr(style = "margin: 15px 0;"),
              
              helpText("The patients are sorted by expression of the gene of interest. Kaplan-Meier curves will be generated:"),
              helpText("- for each half of patients (median),"),
              helpText("- for each quartile (quartile),"),
              helpText("- or for patients above/below a percentile (percentile)."),
              
              tags$hr(style = "margin: 15px 0;"),
              
              downloadButton(ns("plot_download"), "Download Plot", class = "btn-primary w-100"),
              
              shinyBS::bsTooltip(ns("plot_download"), 
                                 title = "Click here to download a copy of the plot",
                                 placement = "right", 
                                 trigger = "hover")
          ),
          
          # Main Panel
          div(class = "main-content",
              tabsetPanel(
                tabPanel("Kaplan-Meier curves",
                         br(),
                         div(style = 'overflow-y: auto;',
                             plotOutput(ns("plot"), height = "80vh", width = "100%"))
                ),
                tabPanel("Patient survival data",
                         br(),
                         DT::dataTableOutput(ns("table"))
                )
              )
          )
      )
    )
  )
}


# Server function for the Kaplan-Meier plot module
kmPlot <- function(input, output, session, dataset, clinData, expData, gene, aligner){
  
  # https://www.mailman.columbia.edu/research/population-health-methods/competing-risk-analysis <- Info on how to handle RR outcome data
  # http://www.math.ucsd.edu/~rxu/math284/CompRisk.pdf
  # https://www.google.com/search?q=competing+events+survival+analysis&oq=competing+events+&aqs=chrome.1.69i57j0l7.4635j1j7&sourceid=chrome&ie=UTF-8
  
  library(survminer)
  library(survival)
  library(cowplot)
  library(cmprsk)
  library(gtools)
  library(openxlsx)
  
  #################################################################
  #------------------------- FUNCTIONS ---------------------------#
  #################################################################
  
  # This is turning into an unwieldy, clunky megafunction and should be broken into some smaller functions: 
  # 1. Fitting the survival object
  # 2. Generating the actual plot
  # These could even be separated into 2 functions each, one for the KM function and the other for cumul. inc analysis
  KMplot <- function(testType) {
    
    # Identifying which event column is needed, 
    # depending on which test type is selected
    validate(
      need(grepl("EFS\\.ID|OS\\.ID|DFS\\.ID|RR\\.ID", colnames(plotData())), "This type of survival data is not available in this dataset.")
    )
    time <- grep(paste0(testType, "\\.Time"), colnames(plotData()), value = T)
    event <- grep(paste0(testType, "\\.ID"), colnames(plotData()), value = T)
    time <- plotData()[,time]
    event <- plotData()[,event]
    
    # Converting the time from days -> years, if requested by the user.
    if (input$time_type == "years") {
      time <- as.numeric(time)/365
      x_max <- 5
      x_breaks <- seq(0, x_max, by = 1)
    } else{
      time <- as.numeric(time)
      x_max <- 365 * 5  # 5 years in days
      x_breaks <- seq(0, x_max, by = 365)  # yearly intervals
      }
    
    # Creating the survival objects for EFS, DFS, and OS.
    if (testType != "RR" ) { 
      surv.obj <- Surv(time = as.numeric(time),
                       event = as.numeric(event))
    }
    
    # Fitting the Kaplan-Meier curves using the recoded survival data.
    # If test type is RR, the cumulative incidence function will be used instead.
    if (testType == "RR") {

      
      event <- ifelse(grepl("Unknown", event), NA, event)
      
      group <- ifelse(input$strata_var == "mutation", input$mutCol, input$strata_var)
      group <- as.vector(plotData()[,group])
      group <- gsub(" ", ".", group)
      fit <- cmprsk::cuminc(ftime = time, 
                            fstatus = event, 
                            cencode = "Censored", # This is the censor code for the TARGET data
                            group = group)
      
      # relapse_names <- grep("Primary event$", names(fit), value = TRUE)
      # 
      # # Check if we found any matching names
      # if (length(relapse_names) == 0) {
      #   stop("No 'Primary event' curves found in the cuminc object.")
      # }
      # 
      # # Subset to only Primary event (relapse) curves
      # fit <- fit[relapse_names]
      # 
      # # Preserve attributes and class
      # #attributes(fit) <- attributes(fit)
      # class(fit) <- "cuminc"
      
    } else {
      formula <- if (input$strata_var == "mutation") { 
        as.formula(paste0("surv.obj ~ ", input$mutCol))
      } else {
        as.formula(paste0("surv.obj ~ ", input$strata_var))
      }
      
      fit <- survminer::surv_fit(formula = formula, data = plotData())
      names(fit$strata) <- gsub(".+=", "", names(fit$strata)) # Fixing the names of the strata to make them more concise & descriptive
      
      # Adding patient count (n) to the strata name (not sure how to do this for RR quite yet)
      for (x in seq(length(names(fit$strata)))) {
        names(fit$strata)[x] <- paste0(names(fit$strata[x]), " (n = ", fit$n[[x]], ")")
      }
    }
    
    # Creating the plot title using the group & test type information
    title <- paste0("KM curves for ", testType, " by ", gene())
    
    nstrata <- if (input$strata_var == "mutation") {
      length(levels(as.factor(plotData()[[input$mutCol]])))
    } else {
      length(levels(as.factor(plotData()[[input$strata_var]])))
    }
    
    # Changing colors of ggsurvplot & ggcomprisk objects:
    # https://stackoverflow.com/questions/51387396/different-color-type-and-line-type-for-multiple-groups-in-survival-curve
    # https://stackoverflow.com/questions/64186565/customizing-a-competing-risks-plot-in-r-with-package-cmprsk
    if (testType == "RR") {
      plot <- ggcompetingrisks(fit, 
                               data = plotData(),
                               xlim = c(0,x_max),
                               multiple_panels = FALSE,
                               palette = c("#fa766d", "#06bfc1", "#7faf1a", "#c47eff"),
                               ggtheme = theme_classic(base_size = bs) +
                                 theme(plot.title = element_text(hjust = 0.5, size = bs + 2),
                                       axis.text = element_text(size = bs),
                                       axis.title = element_text(size = bs),
                                       legend.position = "bottom")) +
        guides(fill = guide_legend(nrow = 2), color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
      
      plot$mapping <- aes(x = time, y = est, colour = group, linetype = event)
      plot <- plot +
        geom_line(size = 1) +
        labs(x = paste0("Time (", input$time_type, ")"),
             title = paste0("RR curves by ", gene()),
             linetype = "Event",
             colour = "Group") +
        scale_x_continuous(limits = c(0, x_max), breaks = x_breaks) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
      
    } else {
      plot <- ggsurvplot(fit, data = plotData(),
                         pval = TRUE,
                         xlim = c(0,x_max),
                         ggtheme = theme_classic(base_size = bs) +
                           theme(plot.title = element_text(hjust = 0.5, size = bs + 2),
                                 axis.text = element_text(size = bs),
                                 axis.title = element_text(size = bs),
                                 legend.position = "bottom"),
                         legend = "bottom",
                         title = title)
      
      plot$plot <- plot$plot +
        guides(fill = guide_legend(title = NULL, nrow = nstrata),
               color = guide_legend(title = NULL, nrow = nstrata)) +
        labs(x = paste0("Time (", input$time_type, ")")) +
        scale_x_continuous(limits = c(0, x_max), breaks = x_breaks) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
      
    }
  }
  
  # Function to generate a table with outcome data used for the figures
  tableFun <- reactive({
    
    validate(
      need(input$test_type, "Please select at least one survival metric to analyze.")
    )
    
    type <- ifelse(input$strata_var == "mutation", input$mutCol, 
                   grep(input$strata_var, colnames(plotData()), value = T))
    
    if (length(input$test_type) == 1) {
      time <- grep(paste0(input$test_type, "\\.time\\.\\.days"), colnames(plotData()), value = T)
      event <- grep(paste0(input$test_type, "\\.ID"), colnames(plotData()), value = T)
    } else if (length(input$test_type) > 1) {
      time <- unlist(lapply(input$test_type, function(x) {
        grep(paste0("^", x, "\\.time\\.\\.days"), colnames(plotData()), value = T)
      }))
      event <- unlist(lapply(input$test_type, function(x) {
        grep(paste0(x, "\\.ID"), colnames(plotData()), value = T)
      }))
    }
  
    plotData() %>%
      dplyr::select(any_of(c("PatientID", !!time, !!event, !!type, "Age.Category", "Primary.Fusion", "SNVs", "Risk", "Filter.Code"))) %>%
      rename(Expression.category = !!type)
  })
  
  # Used the info here to figure out how to display multiple plots:
  # https://stackoverflow.com/questions/51302112/create-plots-based-on-check-box-selection-in-r-shiny
  finalPlot <- reactive({
    
    validate(
      need(dataset() != "StJude", "Survival data is not currently available for this cohort.")
    )
    
    validate(
      need(dataset() != "GMKF", "Survival data is not currently available for this cohort.")
    )
    
    validate(
      need(dataset() != "CCLE", "Survival data is not available for cell line data.")
    )    
    
    validate(
      need(dataset() != "LEUCEGENE", "Survival data is not currently available for this cohort.")
    )   
    
    validate(
      need(dataset() != "PCGP", "Survival data is not currently available for this cohort.")
    )   
    
    validate(
      need(input$test_type, "Please select at least one survival metric to analyze.")
    )
    if (length(input$test_type) == 1) {
      list(KMplot(input$test_type))                   # This creates a gg & ggplot object.
    } else if (length(input$test_type) > 1) {         # Generating multiple plots if the user selects more than 1 type
      lapply(input$test_type, function(x) KMplot(x)) 
    }
  })
  
  #################################################################
  #-------------------- DATA PREPARATION -------------------------#
  #################################################################
  
  observeEvent(dataset(), {
      
      if (dataset() ==  "SWOG") {
        # All columns w/ mutation data follow the same naming convention, so they can be easily selected w/ grep.
        mut_choices <- as.list(grep("mutation", colnames(swog_cde), value = T))
        names(mut_choices) <- gsub("_mutation", "", mut_choices)
        
        updateSelectInput(session = session,
                          inputId = "select_subgroup",
                          label = "Which one?",
                          choices = list(
                            "NPM1" = "NPM1",
                            "TP53" = "TP53",
                            "WT1" = "WT1",
                            "Normal karyotype" = "Normal"))
        
      } else {
        mut_choices <- as.list(filter(colMapping, Module_Code == "Mutation" & !is.na(!!sym(dataset())))$Final_Column_Name)
        names(mut_choices) <- filter(colMapping, Module_Code == "Mutation" & !is.na(!!sym(dataset())))$Final_Column_Label
      }
      
      survival_choices <- filter(colMapping, grepl("\\.ID|\\.Time", Final_Column_Name)) %>%
        filter(!is.na(!!sym(dataset()))) %>% 
        pull(Final_Column_Name) %>%
        gsub("\\..+", "", .) %>%
        unique()
      
      names(survival_choices) <- case_when(grepl("OS", survival_choices) ~ "Overall Survival (OS)",
                                           grepl("EFS", survival_choices) ~ "Event-Free Survival (EFS)",
                                           grepl("DFS", survival_choices) ~ "Disease-Free Survival (DFS)",
                                           grepl("RR", survival_choices) ~ "Relapse Risk (RR)")

      updateCheckboxGroupInput(session = session, 
                               inputId = "test_type",
                               label = "Select a survival type to model",
                               choices = survival_choices)
      
      updateSelectInput(session = session,
                        inputId = "mutCol",
                        label = "Which mutation?",
                        choices = mut_choices)
    })
  
  observeEvent(dataset(), {

    if (dataset() == "BeatAML") {

      updateSelectInput(session = session,
                        inputId = "select_subgroup",
                        label = "Which one?",
                        choices = list(
                          "KMT2A/MLL rearranged" = "MLL|KMT2A-",
                          "inv(16)" = "inv\\(16\\)|CBFB-MYH11",
                          "t(8;21)" = "t\\(8\\;21\\)|RUNX1-RUNX1T1",
                          "FLT3-ITD" = "FLT3-ITD",
                          "NPM1" = "NPM1(?!\\-)",
                          "CEBPA" = "CEBPA",
                          "WT1" = "WT1",
                          "TP53" = "TP53",
                          "PML-RARA" = "PML-RARA",
                          "Normal karyotype" = "Normal"))
      
    } else if (dataset() == "TCGA") {
      
      updateSelectInput(session = session,
                        inputId = "select_subgroup",
                        label = "Which one?",
                        choices = list(
                          "FLT3-ITD" = "FLT3-ITD",
                          "NPM1" = "NPM1(?!\\-)",
                          "CEBPA" = "CEBPA",
                          "WT1" = "WT1",
                          "TP53" = "TP53",
                          "Normal karyotype" = "Normal"))

    } else if (dataset() %in% c("LEUCEGENE", "PCGP")) {
      
      updateSelectInput(session = session,
                        inputId = "select_subgroup",
                        label = "Which one?",
                        choices = "")
    }
  })
  
  # The reactive function below contains multiple components that
  # generate the dataframe that will be used for generating Kaplan-Meier estimates
  # and plotting the survival curves
  plotData <- reactive({
    validate(
      need(gene(), "Please enter a gene symbol or miRNA in the text box to the left.") %then%
        need(gene() %in% rownames(expData()), paste0(gene(), " does not exist in the data! \nDouble-check the symbol, or try an alias."))
    )
    
    # Subsetting exp data to only retain gene of interest & merging with the 
    # clinical metadata. This will be used to categorize patients for the survival analysis
    mergedDF <- expData() %>%
      rownames_to_column("Gene") %>%
      filter(Gene == gene()) %>%
      dplyr::select(Gene, any_of(intersect(clinData()$PatientID, colnames(expData())))) %>%
      column_to_rownames("Gene") %>%
      gather(PatientID, Expression) %>% # Making the dataframe a long-format table, for use with ggplot
      mutate_at(vars(Expression), ~as.numeric(.)) %>%
      left_join(., clinData(), by = "PatientID") %>%
      filter(AML.Sample == "AML") # Removing NBMs, as they don't have outcome data
    
    # Interactively selecting columns that contain fusion or cytogenetic information 
    # (the columns differ between the available datasets).
    # This will be used to restrict the dataset to specific subgroups, if the user opts to do so.
    if (input$filter_cohort == "yes") {
      
      print(colnames(mergedDF))
      
      test <- any(grepl(input$select_subgroup, mergedDF$Filter.Code, perl = T))
      validate(
        need(test, "There are either no records in the dataset that contain the selected alteration,\nor there is currently insufficient data to identify them.")
      )
      mergedDF <- filter(mergedDF, grepl(input$select_subgroup, Filter.Code, perl = T))
    }
    
    # Adding column onto the merged dataset with user-selected plot grouping info,
    # based on the 'strata_var' variable
    if (input$strata_var == "median") {
      mergedDF <- mergedDF %>%
        mutate(median = ifelse(Expression > median(Expression, na.rm = T), "Above median", "Below median"))
      
    } else if (input$strata_var == "quartile") {
      # Throws an error message if the expression data can't be evenly divided into 4 quartiles (in which case tryCatch returns FALSE), 
      # can be checked using quantile(expData$Expression, probs = seq(0, 1, by = 0.25)) to see if Q1-Q3 have a TPM of 0
      nlevels <- try(length(levels(gtools::quantcut(mergedDF$Expression, q = 4, labels = c("Q1","Q2", "Q3","Q4")))))
      
      test <- ifelse(nlevels == 4, TRUE, FALSE)
      validate(
        need(test, "The expression data cannot be evenly divided into 4 quartiles, \nlikely because expression of the gene is very low in the majority of patients. \n\nPlease select 'By median' instead.")
      )
      mergedDF <- mergedDF %>%
        mutate(quartile = gtools::quantcut(Expression, q = 4, labels = c("Q1 - lowest quartile","Q2", "Q3","Q4 - highest quartile")))
      
    } else if (input$strata_var == "percentile") {
      validate(
        need(input$perc_cutoff, "Please enter a percentile to use as a cutoff.") %then%
          need(input$perc_cutoff %in% as.character(seq(1,100)), "The percentile must be between 0 and 100."))
      
      mergedDF <- mergedDF %>%
        mutate(ranks = percent_rank(Expression)) %>%
        mutate(percentile = case_when(ranks*100 < as.numeric(input$perc_cutoff) ~ paste0("Below ", input$perc_cutoff, "%"), 
                                      ranks*100 >= as.numeric(input$perc_cutoff) ~ paste0("Above ", input$perc_cutoff, "%")))
      
    } else if (input$strata_var == "tpm_value") {
      
      # Checking to see if the minimum TPM value is 0 - will add 1 if it's 0, otherwise it will be left as-is
      min_val <- ifelse(ceiling(min(mergedDF$Expression, na.rm = T)) == 0, 1, ceiling(min(mergedDF$Expression, na.rm = T)))
      exp_range <- seq(min_val, floor(max(mergedDF$Expression, na.rm = T)))
      validate(
        need(input$tpm_cutoff, "Please enter a TPM value to use as a cutoff.") %then%
          need(input$tpm_cutoff %in% as.character(exp_range), paste0("Must be a number between ", min_val, " and ", floor(max(mergedDF$Expression, na.rm = T)), "."))
      )
      
      mergedDF <- mergedDF %>%
        mutate(tpm_value = case_when(Expression < as.numeric(input$tpm_cutoff) ~ paste0("Below ", input$tpm_cutoff, " TPM"), 
                                     Expression >= as.numeric(input$tpm_cutoff) ~ paste0("Above ", input$tpm_cutoff, " TPM")))
      
    } else if (input$strata_var == "mutation") {
      validate(
        need(input$mutCol, "Please select a mutation of interest.")
      )
    }
    
    print(mergedDF)
    return(mergedDF)
  })

  
  
  #################################################################
  #-------------------- FINAL MODULE OUTPUTS ---------------------#
  #################################################################
  
  
  #-------------------- Plot tab -----------------------#
  
  output$plot <- renderPlot({
    cowplot::plot_grid(plotlist = finalPlot(), ncol = 2, nrow = 2)
  })
  
  # Adding a download button widget for the plot.
  # Will need to add a prompt for entering username & password prior to download.
  output$plot_download <- downloadHandler(
    filename = function(){
      paste0(as.character(dataset()), "_AML_", paste(input$test_type, collapse = "_"), "_Kaplan-Meier_curves_by_", gene(), "_generated_", format(Sys.time(), "%m.%d.%Y"), ".pdf")
    },
    content = function(file){
      # plots <- cowplot::plot_grid(plotlist = finalPlot(), ncol = 1, nrow = length(input$test_type))
      # ggsave(filename = file, plot = plots, width = plotWidth(), height = plotHeight(), device = "png", dpi = 150)
      pdf(file = file, width = 5, height = 5)
      invisible(lapply(finalPlot(), print))
      dev.off()
    }
  )
  
  # Adding a download button widget for the ggsurvplot object - Does not work yet!
  output$ggsurvplot_download <- downloadHandler(
    filename = function() {
      paste0(as.character(dataset()), "_AML_", paste(input$test_type, collapse = "_"), "_Kaplan-Meier_curves_by_", gene(), "_ggsurvplotObject_generated_", format(Sys.time(), "%m.%d.%Y"), ".RDS")
    },
    content = function(file) {
      withProgress(message = "Saving RDS file", detail = "This may take a while...", value = 0, {
        for (i in 1:70) {
          incProgress(1/70)
          Sys.sleep(0.25)
        }
      })
      plots <- lapply(input$test_type, function(x) KMplot(x)) %>%
        set_names(input$test_type)
      saveRDS(plots, file = file, compress = F)
    }
  )
  
  #-------------------- Outcome data tab -----------------------#
  
  output$table <- DT::renderDataTable({
    DT::datatable(tableFun(), 
                  class = "compact nowrap hover row-border order-column", # Defines the CSS formatting of the final table, can string multiple options together
                  extensions = 'Buttons', # See https://rstudio.github.io/DT/extensions.html for more extensions & features
                  options = list(scrollY = "70vh",
                                 dom = 'Bfrtip',
                                 buttons = list(
                                   list(extend = 'excel', 
                                        filename = paste0(as.character(dataset()), "_Outcome_Summary_Table_basedOn_", gene(), "_Expression_generated_", format(Sys.time(), "%m.%d.%Y")))),
                                 scrollX = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = 50), 
                  escape = F) %>%
      DT::formatStyle(columns = c(1,2,4), fontSize = "100%")
  })
}
