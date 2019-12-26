
# UI for KAPLAN-MEIER app
kmPlotUI <- function(id, label = "Kaplan-Meier plot parameters"){
  
  ns <- NS(id)
  
  tagList( 
    
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
                  
                  # Checkboxes to select the type of test & formate of time data
                  checkboxGroupInput(ns("test_type"), 
                                     label = "Select a survival type to model", 
                                     choices = list("Event-Free Survival (EFS)" = "EFS",
                                                    "Overall Survival (OS)" = "OS", 
                                                    "Disease-Free Survival (DFS)" = "DFS", 
                                                    "Relapse-Free Survival (RFS)")),
                  radioButtons(ns("time_type"), 
                               label = "Select a time interval to display", 
                               choices = list("Years" = "years", 
                                              "Days" = "days")), 
                  
                  # Dropdown menu to select grouping variable for patients, based on expression
                  selectInput(ns("strata_var"), 
                              label = "Select a method of grouping the patients", 
                              choices = list("By median" = "median", 
                                             "By quartile" = "quartile", 
                                             "By percentile" = "percentile")),
                  
                  # Adding 2 text boxes that will only appear if "Percentile" is selected for the strata
                  # but this doesn't seem to work with Shiny Dashboard? 
                  # https://github.com/rstudio/shiny/issues/1586 
                  conditionalPanel(
                    condition = paste0("input['", ns("strata_var"), "'] == 'percentile'"),
                    textInput(ns("cutoff"),                                                  
                              label = "Cutoff percentile",
                              placeholder = "Example: 75")
                  ),
                  
                  # Additional help text
                  helpText("The patients are sorted by expression of the gene of interest.
                             Kaplan-Meier curves will be generated for each half of patients (median), 
                             for each quartile of patients (quartile), 
                             or for patients that fall either above or below a percentile cutoff point specified by the user (percentile).")
                ),
                
                
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
                               
                               column(10, offset = 0, align = "left", 
                                      div(style = 'max-height: 700px; overflow-y: scroll; position: relative',
                                          plotOutput(ns("plot"), height = "900px") 
                                      )
                               ),
                               column(1, offset = 0, align = "right", 
                                      downloadButton(ns("plot_download"), 
                                                     label = "Download plot", class = NULL)
                               )
                             )
                    ), 
                    
                    #-------------------- Patient data tab -----------------------#
                    
                    tabPanel("Patient survival data", 
                             br(),
                             br(),
                             fluidRow(
                               column(10, offset = 0, align = "left", 
                                      # Table of summary stats for plot + outcome data
                                      dataTableOutput(ns("table"))
                               ), 
                               column(1, offset = 0, align = "right", 
                                      downloadButton(ns("table_download"), 
                                                     label = "Download table", class = NULL)
                               )
                             )
                    )
                  )
                )
              )
    )
  )
}

kmPlot <- function(input, output, session, clinData, countsData, gene){
  
  library(survminer)
  library(survival)
  library(tidyverse)
  require(eply)
  require(GGally)
  library(gridExtra)
  library(grid)
  library(gtools)
  library(data.table)
  library(openxlsx)
  library(DT)
  `%then%` <- shiny:::`%OR%` # See https://shiny.rstudio.com/articles/validation.html for details on the %then% operator
  
  #-------------------- Data preparation -----------------------#
  
  # Cleaning up the clinical data grouping columns
  clinData <- clinData %>%
    rename(Primary.Cytogenetic.Code = `Primary Cytogenetic Code`) %>%
    mutate_at(vars(Cytogenetic.Category.2), ~gsub("\\.|\\_", " ", .)) %>%
    mutate_at(vars(Rare.Fusions), ~gsub("\\.", "\\-", .)) %>%
    mutate_at(vars(SNVs), ~gsub("\\.", " ", .)) %>%
    mutate_at(vars(Cytogenetic.Category.2, Rare.Fusions, SNVs), ~gsub("OtherAML", "Other AML", .))
  
  # Interactively filtering the exp data to only retain the gene of interest
  clinData2 <- reactive({
    
    # Outputting error messages if the gene symbol doesn't exist in our data
    validate(
      need(gene(), "Please enter a gene symbol in the text box to the left.") %then%
        need(toupper(gene()) %in% rownames(countsData), "That gene symbol does not exist in the counts data! \nDouble-check the symbol, or try an alias.")
    )
    
    countsData <- countsData[toupper(gene()),] # Subsetting exp data to only retain gene of interest
    
    # Adding counts data onto the CDEs to use for the survival analysis
    countsData <- countsData %>%
      gather(USI, Expression) %>%
      mutate_at(vars(Expression), ~as.numeric(.)) %>%
      left_join(., clinData, by = "USI") %>%
      mutate_at(vars(Group), ~ifelse(grepl("^RO|^BM", USI), "NBM", "AML"))
    
    # Adding column onto the merged dataset with grouping info based on the user selected strata_var
    if(input$strata_var == "median"){
      countsData <- countsData %>%
        mutate(median = ifelse(Expression > median(Expression, na.rm = T), "Above median", "Below median"))
      
    }else if(input$strata_var == "quartile"){
      countsData <- countsData %>%
        mutate(quartile = gtools::quantcut(Expression, q = 4, labels = c("Q1 - lowest quartile", 
                                                                 "Q2", "Q3", 
                                                                 "Q4 - highest quartile")))
    }else if(input$strata_var == "percentile"){
      validate(
        need(input$cutoff, "Please enter a percentile to use as a cutoff.")
      )
      
      countsData <- countsData %>%
        mutate(ranks = percent_rank(Expression)) %>%
        mutate(percentile = case_when(ranks*100 < as.numeric(input$cutoff) ~ paste0("Below ", input$cutoff, "%"), 
                                      ranks*100 >= as.numeric(input$cutoff) ~ paste0("Above ", input$cutoff, "%")))
    }
    
    countsData
  })
  
  
  
  #------------------- Functions --------------------------#
  
  
  
  KMplot <- function(testType){
    
    # Identifying which columns are needed, depending on whether or not the test type is EFS or OS
    event <- clinData2()[,paste0("Recoded ", testType, " ID")]
    
    if(input$time_type == "days"){
      time <- clinData2()[,paste0(testType, " time (days)")]
    }else{
      time <- clinData2()[,paste0(testType, " time (days)")]/365
    }
    
    # Creating the survival objects
    surv.obj <- Surv(time = as.numeric(time),
                     event = as.numeric(event))
    
    formula <- as.formula(paste0("surv.obj ~ ", input$strata_var))
    
    # Fitting the Kaplan-Meier curves using the recoded survival data
    fit <- survminer::surv_fit(formula = formula, data = clinData2())
    
    # Fixing the names of the strata to make them more concise & descriptive
    names(fit$strata) <- gsub(".+=", "", names(fit$strata))
    
    # Adding patient count (n) to the strata name
    for(x in seq(length(names(fit$strata)))){
      names(fit$strata)[x] <- paste0(names(fit$strata[x]), " (n = ", fit$n[[x]], ")")
    }
    
    # Creating the plot title using the group & test type information
    title <- paste0("Kaplan-Meier curves for ", testType)
    
    # Creating the Kaplan-Meier plot
    plot <- ggsurvplot(fit, data = clinData2(),
                       pval = TRUE,
                       ggtheme = theme_bw() + theme(plot.title = element_text(hjust = 0.5, 
                                                                              size = 12)),
                       legend = "right",
                       title = title) +
      guides(
        fill = guide_legend(title = NULL, 
                            nrow = length(levels(as.factor(clinData2()[[input$strata_var]])))),
        color = guide_legend(title = NULL, 
                             nrow = length(levels(as.factor(clinData2()[[input$strata_var]]))))) +
      labs(x = paste0("Time (in ", input$time_type, ")"))
    
    return(plot)
  }
  
  blankPlot <- ggplot() + theme_void()
  
  
  # Making a function to generate a summary table with outcome data
  tableFun <- reactive({
    
    type <- grep(input$strata_var, colnames(clinData2()), value = T)
    
    if(length(input$test_type) == 1){
      time <- grep(paste0("^", input$test_type, " time"), colnames(clinData2()), value = T)
      event <- grep(paste0("^", input$test_type, ".*ID"), colnames(clinData2()), value = T)
      
    }else if(length(input$test_type) > 1){
      
      time <- unlist(lapply(input$test_type, function(x) {
        grep(paste0("^", x, " time"), colnames(clinData2()), value = T)
      }))
      
      event <- unlist(lapply(input$test_type, function(x) {
        grep(paste0("^", x, ".*ID"), colnames(clinData2()), value = T)
      }))
    }
    
    clinData2() %>%
      dplyr::select(Reg., USI, !!time, !!event, !!type) %>% # Using !! before the variables to evaluate them in the select() call
      dplyr::select(Reg., USI, contains("EFS"), contains("OS"), !!type) %>%
      rename(`Expression category` = !!type) ### PERCENTILE TABLE DOESN'T WORK, FIX!!!!!!!!
  })
  
  
  # https://stackoverflow.com/questions/51302112/create-plots-based-on-check-box-selection-in-r-shiny
  finalPlot <- reactive({
    
    # Error message that appears if the user hasn't checked either test type box yet
    validate(need(input$test_type, "Please select at least one survival metric to analyze."))
    
    if(length(input$test_type) == 1){
      
      survp <- KMplot(input$test_type)
      cowplot::plot_grid(survp$plot, nrow = 2)
      
    # Generating multiple plots if the input$test_type object contains both EFS & OS
    }else if(length(input$test_type) > 1){
      
      plots <- lapply(input$test_type, function(x) KMplot(x)) 
      arrange_ggsurvplots(plots, ncol = 1, nrow = 2)
    }
  })
  
  # Setting the dimensions of the final, downloaded plot (this does NOT apply to the plot displayed on the page)
  plotDim <- reactive({
    if(length(input$test_type) == 1){
      6
    }else{
      14
    }
  })
  
  #-------------------- KM plot tab -----------------------------#
  
  output$plot <- renderPlot({
    finalPlot()
  })
  
  # Adding a download button widget for the plot
  output$plot_download <- downloadHandler(
    
    filename = function(){
      paste0("TARGET_AAML1031_", toupper(gene()), "_KaplanMeier_Curves_generated_", format(Sys.time(), "%m.%d.%Y"), ".png")
    },
    content = function(file){
      ggsave(filename = file, plot = finalPlot(), width = plotDim(), device = "png")
    }
  )
  
  #-------------------- Patient data tab -----------------------#
  
  output$table <- renderDataTable({
    datatable(tableFun(), options = list(paging = FALSE, scrollY = "500px"), rownames = F)
  })
  
  # Adding a download button widget for the table
  output$table_download <- downloadHandler(
    
    filename = function(){
      paste0("TARGET_AAML1031_", toupper(gene()), "_Summary_Table_generated_", format(Sys.time(), "%m.%d.%Y"), ".xlsx")
    }, 
    content = function(file){
      write.xlsx(file = file, x = tableFun())
    }
  )
}