# UI function for the waterfall plot module
geneExpUI <- function(id, label = "Identifying expressors"){
  
  library(DT)
  library(shinyjs)
  library(shinyWidgets)
  ns <- NS(id) # Setting a unique namespace for this module
  
  fluidPage(theme = shinytheme("lumen"),
            tags$head(tags$style(HTML('.shiny-output-error-validation { color: #93C54B; }'))), # Custom CSS to modify the app error messages
            tags$head(tags$style(HTML('.box{-webkit-box-shadow: none; -moz-box-shadow: none;box-shadow: none;}'))), # Removes border around boxes
            
            sidebarLayout(
              position = "left", 
              
              sidebarPanel(
                helpText("In the absence of protein expression data, it can be useful to identify expressors of a gene (or miRNA species) using available RNA-seq data. However, in order to determine samples that are positive or negative for the gene of interest, the application of arbitrary cutoff points to the continuous transcript expression data is required."),
                helpText("Recommended cutoff values include 5, or 10 TPM, which are higher than typical sequencing noise, or 1 TPM if the majority of cases express the gene/miRNA species near 0. 
                \nHowever, these values are arbitrary and should be adjusted based on summary statistics and/or the distribution of expression values seen in the histogram."), 
                helpText("The colored line represents the selected expression cutoff."), 
                br(),
                
                # Should this be numericInput? Answer = NO, not for now. numericInput forces a "starting" value and does not allow placeholder suggestions.
                textInput(ns("tpm_cutoff"),                                                  
                          label = "Use TPM cutoff value",
                          placeholder = "Example: 5"),
                
                checkboxInput(ns("other_cutoff"), 
                              label = "Use different selection method", 
                              value = FALSE),
                
                conditionalPanel(
                  condition = paste0("input['", ns("other_cutoff"), "'] == 1"), # No clue why this needs to be 1, since the docs say it's TRUE if checked!
                  radioButtons(ns("select_cutoff"), 
                               label = "Select one of the following options:",
                               choices = list("Use median" = "median",
                                              "Use top quartile" = "quartile"))
                ),
                
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
                helpText("NOTE: Please keep in mind that transcript expression does not always correlate with protein expression.")),
              
              mainPanel(position = "right", 
                        
                        tabsetPanel(
                          tabPanel("Figures",
                                   br(),
                                   fluidRow(
                                     column(10, offset = 0, align = "left",        
                                            plotOutput(ns("histogram"), width = "100%"),
                                            br(),
                                            plotOutput(ns("pieChart"), width = "100%")),  
                                   )
                                   ),
                          tabPanel("Patient list",
                                   fluidRow(
                                     box(width = 12, div(style = 'overflow-x: scroll',
                                                         DT::dataTableOutput(ns("rankedTable")))),
                                   ))
                        )
              )
            )
  )
}

geneExp <- function(input, output, session, clinData, expData, gene, dataset) {
  
  library(tidyverse)
  library(DT)
  library(shinyWidgets)
  
  #################################################################
  #------------------------- FUNCTIONS ---------------------------#
  #################################################################
  
  # Simplified syntax for selecting a cutoff method based on user-provided input
  chooseMethod <- reactive({
    case_when(input$other_cutoff == 1 & input$select_cutoff == "quartile" ~ "quartile",
              input$other_cutoff == 1 & input$select_cutoff == "median" ~ "median",
              TRUE ~ "manual_selection")
  })
  
  # Making a function to generate a long-format table of expression
  makeTable <- reactive({
    
    ############## Messages to the user about required input information ###################
    validate(
      need(gene(), "Please enter a gene symbol or miRNA in the text box to the left.")
    )
    
    if (input$tpm_cutoff == "" & input$other_cutoff == FALSE) {
      validate("No cutoff criteria has been selected. Please choose a cutoff value OR choose a diff selection type.")
    }
    ########################################################################################
    
    # Transforming expression matrix into long-form table
    expTable <- expData() %>%
      rownames_to_column("Gene") %>%
      filter(Gene == gene()) %>%
      pivot_longer(names_to = "Sample.ID", values_to = "Expression", -Gene) %>%
      filter(!is.na(Expression)) %>%
      mutate(Alterations = clinData()$Filter.Code[match(Sample.ID, clinData()$PatientID)],
             AML.Sample = clinData()$AML.Sample[match(Sample.ID, clinData()$PatientID)]) %>%
      filter(AML.Sample == "AML")
    
    # Filtering table if the user wants to restrict the analysis to a single AML subset
    if (input$filter_cohort == "yes") {
      
      # Throws an error message if the selected dataset + filter combo doesn't have any entries
      test <- filter(expTable, grepl(input$select_subgroup, Alterations, perl = T)) %>% nrow()
      validate(
        need(test > 0, "There are either no records in the dataset that contain the selected alteration,\nor there is currently insufficient data to identify them.")
        )
      
      expTable <- filter(expTable, grepl(input$select_subgroup, Alterations, perl = T))
    }
    
    if (chooseMethod() == "quartile") {
      
      nlevels <- try(length(levels(gtools::quantcut(expTable$Expression, q = 4, labels = c("Q1","Q2", "Q3","Q4")))))
      
      # Throws an error message if the expression data can't be evenly divided into 4 quartiles (in which case tryCatch returns FALSE), 
      # can be checked using quantile(expData$Expression, probs = seq(0, 1, by = 0.25)) to see if Q1-Q3 have a TPM of 0
      test <- ifelse(nlevels == 4, TRUE, FALSE)
      validate(
        need(test, "The expression data cannot be evenly divided into 4 quartiles, 
             \nlikely because expression of the gene is very low in the majority of patients. 
             \n\nPlease select a different option.")
      )
      
      expTable <- expTable %>%
        mutate(Quartile = gtools::quantcut(Expression, q = 4, labels = c("Q1","Q2", "Q3","Q4"))) %>%
        mutate(Filter.Category = case_when(Quartiles == "Q4" ~ paste0(gene(), "+"), 
                                           TRUE ~ paste0(gene(), "-")))
      
    } else if (chooseMethod() == "median") {
      
      med <<- median(expTable$Expression, na.rm = T) # Assigning this globally so it can be accessed by the function that determines the histogram x-intercept
      expTable <- mutate(expTable, Filter.Category = case_when(Expression >= med ~ paste0(gene(), "+"), # Text box entry values are stored as characters,
                                                               Expression < med ~ paste0(gene(), "-"),  # so they need to be cast to numeric for this to work
                                                               TRUE ~ NA_character_))
    } else {
      
      # Checking to see if the minimum TPM value is 0 - will add 1 if it's 0, otherwise it will be left as-is
      min_val <- ifelse(ceiling(min(expTable$Expression, na.rm = T)) == 0, 1, ceiling(min(expTable$Expression, na.rm = T)))
      exp_range <- seq(min_val, floor(max(expTable$Expression, na.rm = T)))
      
      validate(
        need(input$tpm_cutoff, "Please provide a TPM cutoff.") %then%
          need(input$tpm_cutoff %in% as.character(exp_range), paste0("Must be a number between ", min_val, " and ", floor(max(expTable$Expression, na.rm = T)), "."))
      )
      
      expTable <- mutate(expTable, Filter.Category = case_when(Expression >= as.numeric(input$tpm_cutoff) ~ paste0(gene(), "+"), # Text box entry values are stored as characters,
                                                               Expression < as.numeric(input$tpm_cutoff) ~ paste0(gene(), "-"),  # so they need to be cast to numeric for this to work
                                                               TRUE ~ NA_character_))
    }
    
    expTable <-  expTable %>%
      arrange(desc(Expression))
    
    return(expTable)
    
  })
  
  
  xint <- reactive({
    
    finalTable <- makeTable()
    
    if (chooseMethod() == "quartile") {
      x <- filter(finalTable, grepl("\\+", Filter.Category)) %>%
        pull(Expression) %>%
        min(., na.rm = T)
      
      print("Other cutoff selected, will use quartiles")
      print(paste0("Q4 starts at ", round(x, 2), " TPM"))
      print("")
      
    } else if (chooseMethod() == "median") {
      x <- med
      
      print("Other cutoff selected, will use median")
      print(paste0("Median is ", round(x, 2), " TPM"))
      print("")
    } else {
      # This needs to be cast as a numeric, since the textInput() widget creates a character string
      x <- as.numeric(input$tpm_cutoff) 
    }
    
    return(x)
  })
  
  
  plotHist <- reactive({
    
    expTable <- makeTable()
    x_intercept <- xint()
    
    rnge <- c(min(expTable$Expression, na.rm = T), max(expTable$Expression, na.rm = T))
    print(rnge)
    
    bincount <- rnge[2] - rnge[1]
    bincount <- case_when(bincount >= 200 ~ 100,
                          # bincount > 50 & bincount < 100 ~ 70,
                          bincount < 1 ~ 1,
                          TRUE ~ bincount)
    print(bincount)
    
    hist(expTable$Expression, 
         breaks = bincount, 
         xlab = "TPM value", 
         ylab = "Frequency", 
         main = "Distribution of TPM values")
    abline(v = x_intercept, col = "red") # Need this to be able to display the median or Q4 cutoff as well!
  })
  
  
  
  plotPie <- reactive({
    
    expTable <- makeTable() %>%
      group_by(Filter.Category) %>%
      count(Filter.Category)
    
    ref_lvl <- unique(grep("\\+", expTable$Filter.Category, value = T))
    
    print(ref_lvl)
    
    expTable <- expTable %>%
      mutate(prop = n / sum(expTable$n, na.rm = T)*100) %>%
      mutate(ypos = cumsum(prop)- 0.5*prop) %>%
      mutate(lab = paste0("N=", n, "\n(", round(prop, 2), "%)")) %>%
      ungroup() %>% # across() apparently can't find columns that are grouped, yikes! This took me forever to figure out!
      mutate(across(Filter.Category, as.factor)) %>% 
      mutate(across(Filter.Category, ~relevel(., ref = ref_lvl)))
    
    plot <-  ggplot(expTable, aes(x = "", y = prop, fill = Filter.Category)) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar("y", start = 0) +
      theme_void() + 
      # Placing text labels on a pie chart:
      # https://stackoverflow.com/questions/16184188/ggplot-facet-piechart-placing-text-in-the-middle-of-pie-chart-slices
      geom_text(aes(label = lab), position = position_stack(vjust = 0.5), color = "white", size = 6) + # Change labels to an N value & % of total!!!!
      labs(title = paste0("Proportion of ", gene(), "+ cases\n(using selected criteria)")) +
      theme(legend.position = "bottom", 
            legend.title = element_blank(),
            legend.text = element_text(size = 13),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    
    plot
    
  })
  
  
  #################################################################
  #-------------------- FINAL MODULE OUTPUTS ---------------------#
  #################################################################
  
  # https://glin.github.io/reactable/articles/examples.html#conditional-styling
  output$rankedTable <- DT::renderDataTable({
    t <- makeTable() %>% 
      filter(grepl("\\+", Filter.Category))
    
    DT::datatable(t, 
                  class = "compact nowrap hover row-border order-column", # Defines the CSS formatting of the final table, can string multiple options together
                  extensions = 'Buttons', # See https://rstudio.github.io/DT/extensions.html for more extensions & features
                  options = list(scrollY = "70vh",
                                 dom = 'Bfrtip',
                                 buttons = list(
                                   list(extend = 'excel', filename = paste0(dataset(), "_AML_", gene(), "pos_patientList_generated_", format(Sys.time(), "%m.%d.%Y"), ".xlsx"))),
                                 scrollX = TRUE,
                                 # fixedColumns = list(leftColumns = 1),
                                 searchHighlight = TRUE,
                                 pageLength = 50), 
                  escape = F) %>%
      DT::formatStyle(columns = c(1,2,4), fontSize = "100%")
  })
  
  output$histogram <- renderPlot({
    plotHist()
  })
  
  output$pieChart <- renderPlot({
    plotPie()
  })
}