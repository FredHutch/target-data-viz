# UI for the Kaplan-Meier plot module
kmPlotUI <- function(id, label = "Kaplan-Meier plot parameters"){
  
  library(DT)
  ns <- NS(id)
  
  tagList( 
    # Throwing a fluidPage() inside of the taglist so I can use a Bootstrap theme... it does work, not sure how kosher this is though?
    fluidPage(theme = shinytheme("lumen"),
              tags$head(tags$style(HTML('.shiny-output-error-validation {
                                                                         color: #93C54B;
                                                                         }'))), # Custom CSS to modify the app error messages
    
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
                  
                  # Dropdown menu to select grouping variable for patients, based on expression
                  selectInput(ns("strata_var"), 
                              label = "Select a method of grouping the patients", 
                              choices = list("By median" = "median", 
                                             "By quartile" = "quartile", 
                                             "By percentile" = "percentile", 
                                             "By mutation status" = "mutation")),
                  
                  # Adding 2 text boxes that will only appear if "Percentile" is selected for the strata
                  # but this doesn't seem to work with Shiny Dashboard? 
                  # https://github.com/rstudio/shiny/issues/1586 
                  conditionalPanel(
                    condition = paste0("input['", ns("strata_var"), "'] == 'percentile'"),
                    textInput(ns("cutoff"),                                                  
                              label = "Cutoff percentile",
                              placeholder = "Example: 75")),
                  
                  conditionalPanel(
                    condition = paste0("input['", ns("strata_var"), "'] == 'mutation'"),
                    selectInput(ns("mutCol"),                                                  
                              label = "Which mutation?",
                              choices = list("NPM1" = "NPM.mutation.", 
                                             "CEBPA" = "CEBPA.mutation.", 
                                             "WT1" = "WT1.mutation.", 
                                             "cKit (exon 8)" = "c.Kit.Mutation.Exon.8", 
                                             "cKit (exon 17)" = "c.Kit.Mutation.Exon.17", 
                                             "RAS mutation" = "RAS.Mutation", 
                                             "RAS gene" = "RAS.Gene", 
                                             "CBL" = "CBL.Mutation", 
                                             "FLT3-ITD" = "FLT3.ITD.positive.", 
                                             "FLT3 point mutation" = "FLT3.PM.category"))),
                  
                  helpText("The patients are sorted by expression of the gene of interest.
                             Kaplan-Meier curves will be generated for each half of patients (median), 
                             for each quartile of patients (quartile), 
                             or for patients that fall either above or below a percentile cutoff point specified by the user (percentile).")),
                
                ###############################################################
                #----------------------- MAIN PLOT PANEL ---------------------#
                ###############################################################
                
                mainPanel(
                  position = "right", 
                  tabsetPanel(
                    
                    #-------------------- KM plot tab -----------------------#
                    
                    tabPanel("Kaplan-Meier curves", 
                             br(),
                             br(),
                             fluidRow(
                               style = 'height:40vh',
                               column(9, offset = 0, align = "left", 
                                      div(style = 'max-height: 500px; overflow-y: scroll; position: relative',
                                          plotOutput(ns("plot"), height = "1500px"))),
                               column(2, offset = 0, align = "right", 
                                      downloadButton(ns("plot_download"), 
                                                     label = "Download plot")),
                               column(2, offset = 0, align = "right", 
                                      downloadButton(ns("ggsurvplot_download"), 
                                                     label = "ggsurvplot object", 
                                                     style = 'padding:5px; font-size:70%; margin-top:10px',
                                                     class = "btn-info"))
                             )
                    ), 
                    
                    #-------------------- Patient data tab -----------------------#
                    
                    tabPanel("Patient survival data", 
                             br(),
                             br(),
                             fluidRow(
                               column(11, offset = 0, align = "left", 
                                      DT:: dataTableOutput(ns("table"))),  # Table of summary stats for plot + outcome data
                               column(1, offset = 0, align = "right", 
                                      downloadButton(ns("table_download"), 
                                                     label = "Download table", class = NULL))
                             )
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
  library(tidyverse)
  library(cowplot)
  library(gtools)
  library(data.table)
  library(DT)
  library(openxlsx)
  `%then%` <- shiny:::`%OR%` # See https://shiny.rstudio.com/articles/validation.html for details on the %then% operator
  
  #-------------------- Data preparation -----------------------#
  
  observeEvent(dataset(), {
    if(dataset() == "BeatAML") {
      updateCheckboxGroupInput(
        session = session, 
        inputId = "test_type",
        label = "Select a survival type to model",
        choices = c("Overall Survival (OS)" = "OS"))
      
      updateSelectInput(session = session,
                        inputId = "mutCol",
                        label = "Which mutation?",
                        choices = list("NPM1", "CEBPA", "WT1", "KIT", "JAK2", "NRAS", "FLT3-ITD" = "FLT3.ITD", "TP53"))
      
    } else if (dataset() == "TARGET") {
      updateCheckboxGroupInput(
        session = session, 
        inputId = "test_type",
        label = "Select a survival type to model",
        choices = c("Event-Free Survival (EFS)" = "EFS",
                       "Overall Survival (OS)" = "OS", 
                       "Disease-Free Survival (DFS)" = "DFS", 
                       "Relapse Risk (RR)" = "RR"))
      
      updateSelectInput("mutCol",   
                        session = session,
                  label = "Which mutation?",
                  choices = list("NPM1" = "NPM.mutation.", "CEBPA" = "CEBPA.mutation.", "WT1" = "WT1.mutation.", "cKit (exon 8)" = "c.Kit.Mutation.Exon.8", 
                                 "cKit (exon 17)" = "c.Kit.Mutation.Exon.17", "RAS mutation" = "RAS.Mutation", "RAS gene" = "RAS.Gene", "CBL" = "CBL.Mutation", 
                                 "FLT3-ITD" = "FLT3.ITD.positive.", "FLT3 point mutation" = "FLT3.PM.category"))
    }
  }, ignoreInit = T)
  
  # Interactively filtering the exp data to only retain the gene of interest
  plotData <- reactive({
    
    # Outputting error messages if the gene symbol doesn't exist in our data
    validate(
      need(gene(), "Please enter a gene symbol or miRNA in the text box to the left.") %then%
        need(gene() %in% rownames(expData()), "That gene symbol does not exist in the counts data! \nDouble-check the symbol, or try an alias.")
    )
    
    # Subsetting exp data to only retain gene of interest & adding counts data onto the CDEs to use for the survival analysis
    mergedDF <- expData() %>%
      rownames_to_column("Gene") %>%
      filter(Gene == gene()) %>%
      dplyr::select(Gene, any_of(intersect(clinData()$USI, colnames(expData())))) %>%
      column_to_rownames("Gene") %>%
      gather(USI, Expression) %>%
      mutate_at(vars(Expression), ~as.numeric(.)) %>%
      left_join(., clinData(), by = "USI") %>%
      mutate_at(vars(Group), ~ifelse(grepl("^RO|^BM", USI), "NBM", "AML"))
    
    # Adding column onto the merged dataset with grouping info based on the user selected strata_var
    if (input$strata_var == "median") {
      mergedDF <- mergedDF %>%
        mutate(median = ifelse(Expression > median(Expression, na.rm = T), "Above median", "Below median"))
      
    } else if (input$strata_var == "quartile") {
      # Throws an error message if the expression data can't be evenly divided into 4 quartiles, 
      # check with quantile(expData$Expression, probs = seq(0, 1, by = 0.25)) to see if Q1-Q3 have a TPM of 0
      validate(
        need(length(levels(gtools::quantcut(mergedDF$Expression, q = 4))) == 4, "The expression data cannot be evenly divided into quartiles, \nlikely because expression of the gene is very low in the majority of patients. \n\nPlease select 'By median' instead.")
        )
      mergedDF <- mergedDF %>%
        mutate(quartile = gtools::quantcut(Expression, q = 4, labels = c("Q1 - lowest quartile", 
                                                                 "Q2", "Q3", 
                                                                 "Q4 - highest quartile")))
    } else if (input$strata_var == "percentile") {
      validate(
        need(input$cutoff, "Please enter a percentile to use as a cutoff."))
          # need(input$cutoff > 0 && input$cutoff < 100, "The percentile must be between 0 and 100.")) # This doesn't work for some reason??
      
      mergedDF <- mergedDF %>%
        mutate(ranks = percent_rank(Expression)) %>%
        mutate(percentile = case_when(ranks*100 < as.numeric(input$cutoff) ~ paste0("Below ", input$cutoff, "%"), 
                                      ranks*100 >= as.numeric(input$cutoff) ~ paste0("Above ", input$cutoff, "%")))
      
    } else if (input$strata_var == "mutation") {
      validate(
        need(input$mutCol, "Please select a mutation of interest.")
      )
    }
    
    return(mergedDF)
  })
  
  
  
  #------------------- Functions --------------------------#
  
  KMplot <- function(testType) {
    
    # Identifying which event column is needed, 
    # depending on whether or not the test type is EFS, OS, DFS, or RR
    if (testType == "DFS") {
      # Recoding an "event" as 1 (censored data is 0)
      event <- ifelse(grepl("[Ee]vent", plotData()[,"DFS from end of induction 1 for patients who are CR at EOI1 indicator"]), 1, 0) 
      time <- plotData()[,"Days to DFS from end of induction 1 for patients who are CR at EOI1"]
    } else if (testType == "RR") {
      event <- ifelse(grepl("[Ee]vent", plotData()[,"RR from CR (end of course 1) indicator"]), 1, 0)
      time <- plotData()[,"Days to RR from CR (end of course 1)"]
    } else if (testType %in% c("EFS", "OS")) {
      
      validate(
        need(paste0("Recoded ", testType, " ID") %in% colnames(plotData()), "This type of survival data is not available in this dataset.")
      )
      
      event <- plotData()[,paste0("Recoded ", testType, " ID")]
      time <- plotData()[,paste0(testType, " time (days)")]
    }
    
    # Converting the time from days -> years, if requested by the user
    if (input$time_type == "years") {
      time <- time/365
    }
    
    # Creating the survival objects
    surv.obj <- Surv(time = as.numeric(time),
                     event = as.numeric(event))
    
    formula <- if (input$strata_var == "mutation") {
      as.formula(paste0("surv.obj ~ ", input$mutCol))
    } else {
      as.formula(paste0("surv.obj ~ ", input$strata_var))
    }
    
    # Fitting the Kaplan-Meier curves using the recoded survival data
    fit <- survminer::surv_fit(formula = formula, data = plotData())
    
    # Fixing the names of the strata to make them more concise & descriptive
    names(fit$strata) <- gsub(".+=", "", names(fit$strata))
    
    # Adding patient count (n) to the strata name
    for (x in seq(length(names(fit$strata)))) {
      names(fit$strata)[x] <- paste0(names(fit$strata[x]), " (n = ", fit$n[[x]], ")")
    }
    
    # Creating the plot title using the group & test type information
    title <- paste0("Kaplan-Meier curves for ", testType)
    
    strataCol <- if (input$strata_var == "mutation") { 
      length(levels(as.factor(plotData()[[input$mutCol]])))
    } else {
      length(levels(as.factor(plotData()[[input$strata_var]])))
    }
    
    # Creating the Kaplan-Meier plot
    plot <- ggsurvplot(fit, data = plotData(),
                       pval = TRUE,
                       ggtheme = theme_bw() + theme(plot.title = element_text(hjust = 0.5, 
                                                                              size = 12)),
                       legend = "right",
                       title = title) +
      guides(
        fill = guide_legend(title = NULL, 
                            nrow = strataCol), # This is problematic for the mutation strata type, need to figure out a workaround
        color = guide_legend(title = NULL, 
                             nrow = strataCol)) +
      labs(x = paste0("Time (in ", input$time_type, ")"))
    
    # Adding gene name as a text annotation layer, to be displayed in the top right corner of the plot
    # https://stackoverflow.com/questions/10747307/legend-placement-ggplot-relative-to-plotting-region
    plot$plot <- plot$plot +
      ggplot2::annotate("text", y = 1, x = max(time, na.rm = T), label = gene(), hjust = 0.5, vjust = 1, size = 8, fontface = "bold")
    
    return(plot$plot) # Final returned object from this function is the "plot" component of the ggsurvplot object
  }
  
  # Making a function to generate a summary table with outcome data
  tableFun <- reactive({
    
    type <- ifelse(input$strata_var == "mutation", input$mutCol, 
                   grep(input$strata_var, colnames(plotData()), value = T))
    
    if (length(input$test_type) == 1) {
      time <- grep(paste0("^", input$test_type, " time"), colnames(plotData()), value = T)
      event <- grep(paste0("^", input$test_type, ".*ID"), colnames(plotData()), value = T)
    } else if (length(input$test_type) > 1) {
      time <- unlist(lapply(input$test_type, function(x) {
        grep(paste0("^", x, " time"), colnames(plotData()), value = T)
      }))
      event <- unlist(lapply(input$test_type, function(x) {
        grep(paste0("^", x, ".*ID"), colnames(plotData()), value = T)
      }))
    }
    
    plotData() %>%
      dplyr::select(any_of(c("Reg.", "USI", !!time, !!event, !!type))) %>% # Using !! before the variables to evaluate them to colnames in the select() call
      dplyr::select(any_of(c("Reg.", "USI")), contains("EFS"), contains("OS"), !!type) %>%
      rename(`Expression category` = !!type) ### PERCENTILE TABLE DOESN'T WORK, FIX!!!!!!!!
  })
  
  # https://stackoverflow.com/questions/51302112/create-plots-based-on-check-box-selection-in-r-shiny
  finalPlot <- reactive({
    # Error message that appears if the user hasn't checked either test type box yet
    validate(
      need(input$test_type, "Please select at least one survival metric to analyze.")
      )
    if (length(input$test_type) == 1) {
      list(KMplot(input$test_type))                   # This creates a gg & ggplot object
    } else if (length(input$test_type) > 1) {         # Generating multiple plots if the user selects more than 1 type
      lapply(input$test_type, function(x) KMplot(x)) 
    }
  })
  
  # Setting the dimensions of the final, downloaded plot (this does NOT apply to the plot displayed on the page)
  plotHeight <- reactive({
    4 * length(input$test_type)
  })
  
  plotWidth <- reactive({
    if (any(nchar(unique(clinData[[input$mutCol]]) > 27))) {
      6.5
    } else {
      5.5
    }
  })
  
  #-------------------- KM plot tab -----------------------------#
  
  output$plot <- renderPlot({
    cowplot::plot_grid(plotlist = finalPlot(), ncol = 1, nrow = 4)
  })
  
  # Adding a download button widget for the plot
  output$plot_download <- downloadHandler(
    filename = function(){
      paste0("TARGET_AAML1031_", gene(), "_KaplanMeier_Curves_", paste(input$test_type, collapse = "_"), "_generated_", format(Sys.time(), "%m.%d.%Y"), ".png")
    },
    content = function(file){
      plots <- cowplot::plot_grid(plotlist = finalPlot(), ncol = 1, nrow = length(input$test_type))
      ggsave(filename = file, plot = plots, width = plotWidth(), height = plotHeight(), device = "png", dpi = 150)
    }
  )
  
  # Adding a download button widget for the ggsurvplot object - THIS STILL NEEDS WORK
  output$ggsurvplot_download <- downloadHandler(
    filename = function() {
      paste0("TARGET_AAML1031_", gene(), "_KaplanMeier_Curves_", paste(input$test_type, collapse = "_"), "_ggsurvplotObject_generated_", format(Sys.time(), "%m.%d.%Y"), ".RDS")
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
  
  #-------------------- Patient data tab -----------------------#
  
  output$table <- DT::renderDataTable({
    DT::datatable(tableFun(), options = list(paging = FALSE, scrollY = "500px"), rownames = F)
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
}