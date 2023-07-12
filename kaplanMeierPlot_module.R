# UI for the Kaplan-Meier plot module
kmPlotUI <- function(id, label = "Kaplan-Meier plot parameters"){
  ns <- NS(id) # Setting a unique namespace for this module
  
  mut_choices <- as.list(filter(colMapping, Module_Code == "Mutation" & !is.na(Final_Column_Label))$Final_Column_Name)
  names(mut_choices) <- filter(colMapping, Module_Code == "Mutation" & !is.na(Final_Column_Label))$Final_Column_Label
  
  # Using tagList() instead of fluidPage() to allow for the ADC/CAR T-cell therapy button to change the tab
  tagList(
    fluidPage(
    
              ###############################################################
              #----------------------- SIDEBAR -----------------------------#
              ###############################################################
              
              sidebarLayout(
                position = "left", 
                sidebarPanel(
                  
                  # Checkboxes to select the type of test & format of time data
                  checkboxGroupInput(ns("test_type"), 
                                     label = "Select a survival type to model", 
                                     choices = list("Event-Free Survival (EFS)" = "EFS",
                                                    "Overall Survival (OS)" = "OS", 
                                                    "Disease-Free Survival (DFS)" = "DFS", 
                                                    "Relapse Risk (RR)" = "RR")),
                  radioButtons(ns("time_type"), 
                               label = "Select a time interval to display", 
                               choices = list("Years" = "years", 
                                              "Days" = "days")), 
                  
                  radioButtons(ns("filter_cohort"), 
                               label = "Limit to a single subgroup?", 
                               choices = list("No - view all patients" = "no",
                                              "Yes" = "yes")),
                  conditionalPanel(
                    condition = paste0("input['", ns("filter_cohort"), "'] == 'yes'"),
                    selectInput(ns("select_subgroup"),                                                  
                                label = "Which one?",
                                choices = list("KMT2A/MLL rearranged" = "MLL|KMT2A-",       # The name of each list item is the option the
                                               "inv(16)" = "inv\\(16\\)",                   # user will see in a drop-down menu;
                                               "t(8;21)" = "t\\(8\\;21\\)",                 # the value is the regex that will be used to 
                                               "FLT3-ITD" = "FLT3-ITD",                     # filter the dataset for the specified alteration.
                                               "KMT2A-PTD" = "KMT2A-PTD",
                                               "WT1" = "WT1",
                                               "NPM1" = "NPM1(?!\\-)",
                                               "CEBPA" = "CEBPA",
                                               "CBFA2T3-GLIS2" = "CBFA2T3\\-GLIS2", 
                                               "NUP98 fusions" = "NUP98-",
                                               "Monosomy 7" = "[Mm]onosomy7",
                                               "del5q/del7q" = "del5q|del7q",
                                               "Trisomy 8" = "[Tt]risomy8",
                                               "Normal karyotype" = "Normal"))
                  ), 
                  
                  # Dropdown menu to select grouping variable for patients, based on expression
                  selectInput(ns("strata_var"), 
                              label = "Select a method of grouping the patients", 
                              choices = list("By median" = "median", 
                                             "By quartile" = "quartile", 
                                             "By percentile" = "percentile", 
                                             "By TPM cutoff" = "tpm_value",
                                             "By mutation status" = "mutation")),
                  
                  # Adding 2 text boxes that will only appear if "Percentile" is selected for the strata
                  # but this doesn't seem to work with Shiny Dashboard? 
                  # https://github.com/rstudio/shiny/issues/1586 
                  conditionalPanel(
                    condition = paste0("input['", ns("strata_var"), "'] == 'percentile'"),
                    textInput(ns("perc_cutoff"),                                                  
                              label = "Cutoff percentile",
                              placeholder = "Example: 75")),
                  
                  conditionalPanel(
                    condition = paste0("input['", ns("strata_var"), "'] == 'tpm_value'"),
                    textInput(ns("tpm_cutoff"),                                                  
                              label = "TPM cutoff value",
                              placeholder = "Example: 5")),
                  
                  conditionalPanel(
                    condition = paste0("input['", ns("strata_var"), "'] == 'mutation'"),
                    selectInput(ns("mutCol"),                                                  
                              label = "Which mutation?",
                              choices = mut_choices)),
                  
                  br(),
                  
                  helpText("The patients are sorted by expression of the gene of interest.
                             Kaplan-Meier curves will be generated for each half of patients (median), 
                             for each quartile of patients (quartile), 
                             or for patients that fall either above or below a percentile cutoff point specified by the user (percentile)."),
                br(),
                br(),
                downloadButton(ns("plot_download"), 
                               label = "plot", 
                               class = "plotdwnld"),
                
                shinyBS::bsTooltip(ns("plot_download"), 
                                   title = "Click here to download a copy of the plot",
                                   placement = "right", 
                                   trigger = "hover")
                ),
              
                ###############################################################
                #----------------------- MAIN PLOT PANEL ---------------------#
                ###############################################################
                
                mainPanel(
                  position = "right", 
                  tabsetPanel(
                    
                    #-------------------- Survival analysis tabs -----------------------#
                    
                    tabPanel("Kaplan-Meier curves", 
                             br(),
                             br(),
                             fluidRow(
                               style = 'height:40vh',
                               div(style = 'max-height: 900px; overflow-y: scroll; position: relative',
                                   plotOutput(ns("plot"), height = "1400px", width = "50%")),
                             )
                    ), 
                    
                    #-------------------- Patient data tab -----------------------#
                    
                    tabPanel("Patient survival data", 
                             br(),
                             br(),
                             fluidRow(
                               DT::dataTableOutput(ns("table")) # Table of summary stats for plot + outcome data
                             ),  
                    )
                  )
                )
              )
    )
  )
}

# Server function for the Kaplan-Meier plot module
kmPlot <- function(input, output, session, dataset, clinData, expData, gene){
  
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
      need(grepl("EFS\\.ID|OS\\.ID|DFS\\.ID|RR\\.ID", colnames(plotData())), "This type of survival data is not available in this dataset. --LOGAN --10AM --7/12")
    )
    time <- grep(paste0(testType, "\\.Time"), colnames(plotData()), value = T)
    event <- grep(paste0(testType, "\\.ID"), colnames(plotData()), value = T)
    time <- plotData()[,time]
    event <- plotData()[,event]
    
    # Converting the time from days -> years, if requested by the user.
    if (input$time_type == "years") {
      time <- time/365
    }
    
    # Creating the survival objects for EFS, DFS, and OS.
    if (testType != "RR" ) { 
      surv.obj <- Surv(time = as.numeric(time),
                       event = as.numeric(event))
    }
    
    # Fitting the Kaplan-Meier curves using the recoded survival data.
    # If test type is RR, the cumulative incidence function will be used instead.
    if (testType == "RR") {
      group <- ifelse(input$strata_var == "mutation", input$mutCol, input$strata_var)
      group <- as.vector(plotData()[,group])
      group <- gsub(" ", ".", group)
      fit <- cmprsk::cuminc(ftime = time, 
                            fstatus = event, 
                            cencode = "Censored", # This is the censor code for the TARGET data
                            group = group)
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
    title <- paste0("Kaplan-Meier curves for ", testType)
    
    nstrata <- if (input$strata_var == "mutation") {
      length(levels(as.factor(plotData()[[input$mutCol]])))
    } else {
      length(levels(as.factor(plotData()[[input$strata_var]])))
    }
    
    # Changing colors of ggsurvplot & ggcomprisk objects:
    # https://stackoverflow.com/questions/51387396/different-color-type-and-line-type-for-multiple-groups-in-survival-curve
    # https://stackoverflow.com/questions/64186565/customizing-a-competing-risks-plot-in-r-with-package-cmprsk
    if (testType == "RR") {
      # Relapse risk plots need to be generated w/ the cumulative incidence of risk calculation, not the Kaplan-Meier method,
      # so this will need a different plotting method to accommodate the cuminc object ('fit').
      plot <- ggcompetingrisks(fit, 
                               data = plotData(),
                               multiple_panels = F,
                               palette = RColorBrewer::brewer.pal(9, "Set1")[c(1:5, 7:9)],
                               ggtheme = theme_classic(base_size = bs) + 
                                 theme(plot.title = element_text(hjust = 0.5, size = bs + 2),
                                       legend.position = "bottom")) +
        ggplot2::annotate("text", y = 1, x = (max(time, na.rm = T) - max(time, na.rm = T)/5),
                          label = gene(), hjust = 0.5, vjust = 1, size = 8, fontface = "bold") +
        guides(fill = guide_legend(nrow = 2), color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
    
      # Changing default strata schema, rather than using a distinct linetype for each group (aka the default)
      # https://stackoverflow.com/questions/64186565/customizing-a-competing-risks-plot-in-r-with-package-cmprsk
      plot$mapping <- aes(x = time, y = est, colour = group, linetype = event)
      plot <- plot + 
        geom_line(size = 1) + # Making line wider than default, to match the KM plot for OS, EFS, etc.
        labs(x = paste0("Time (", input$time_type, ")"),
             title = "Cumulative incidence curves for RR",
             linetype = "Event",
             colour = "Group")
    } else {
      # Creates a plot for OS, EFS, and DFS, which use Kaplan-Meier curves and the logrank test
      plot <- ggsurvplot(fit, data = plotData(), # Creating the Kaplan-Meier plot
                         pval = TRUE,
                         ggtheme = theme_classic(base_size = bs) + theme(plot.title = element_text(hjust = 0.5, size = bs + 2),
                                                                         axis.text.x = element_text(size = bs + 1)),
                         legend = "bottom",
                         title = title) +
        guides(fill = guide_legend(title = NULL, nrow = nstrata), 
               color = guide_legend(title = NULL, nrow = nstrata)) +
        labs(x = paste0("Time (", input$time_type, ")"))
      
      # Adding gene name as a text annotation layer (to the ggsurvplot object), 
      # to be displayed in the top right corner of the plot
      # https://stackoverflow.com/questions/10747307/legend-placement-ggplot-relative-to-plotting-region
      plot$plot <- plot$plot +
        ggplot2::annotate("text", y = 1, x = (max(time, na.rm = T) - max(time, na.rm = T)/5), 
                          label = gene(), hjust = 0.5, vjust = 1, size = 8, fontface = "bold")
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
      need(dataset() != "StJude", "Survival data is not currently available for this cohort. Please try again later.")
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
    
    return(mergedDF)
  })

  
  
  #################################################################
  #-------------------- FINAL MODULE OUTPUTS ---------------------#
  #################################################################
  
  
  #-------------------- Plot tab -----------------------#
  
  output$plot <- renderPlot({
    cowplot::plot_grid(plotlist = finalPlot(), ncol = 1, nrow = 4)
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