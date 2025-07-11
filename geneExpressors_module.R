geneExpUI <- function(id, label = "Identifying expressors") {
  ns <- NS(id)
  
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
          overflow-y: auto;
        }
        .main-content {
          flex-grow: 1;
          padding: 15px;
        }
        .plotdwnld, .custom-file-upload .btn {
          background-color: #2096f6 !important;
          color: white !important;
          width: 100% !important;
          padding: 3px 7px;
          font-size: 1.2rem;
        }
        .custom-file-upload .form-control {
          font-size: 1.2rem;
          padding: 3px 7px;
        }
        .ui-slider-range {
          background: #2096f6 !important;
        }
        .ui-slider-handle {
          background: #2096f6 !important;
          border-color: #2096f6 !important;
        }
      "))
    ),
    fluidPage(
      theme = shinythemes::shinytheme("paper"),
      div(class = "sidebar-container",
          
          #------------------ SIDEBAR ------------------#
          div(class = "custom-sidebar",
              helpText("In the absence of protein expression data, it can be useful to identify expressors of a gene (or miRNA species) using available RNA-seq data. However, in order to determine samples that are positive or negative for the gene of interest, the application of arbitrary cutoff points to the continuous transcript expression data is required."),
              helpText("Recommended cutoff values include 5, or 10 TPM, which are higher than typical sequencing noise, or 1 TPM if the majority of cases express the gene/miRNA species near 0. 
          \nHowever, these values are arbitrary and should be adjusted based on summary statistics and/or the distribution of expression values seen in the histogram."), 
              helpText("The colored line represents the selected expression cutoff."), 
              br(),
              
              textInput(ns("tpm_cutoff"),                                                  
                        label = "Use TPM cutoff value",
                        placeholder = "Example: 5"),
              
              checkboxInput(ns("other_cutoff"), 
                            label = "Use different selection method", 
                            value = FALSE),
              
              conditionalPanel(
                condition = paste0("input['", ns("other_cutoff"), "'] == 1"),
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
                            choices = list(
                              "KMT2A/MLL rearranged" = "MLL|KMT2A-",
                              "inv(16)" = "inv\\(16\\)",
                              "t(8;21)" = "t\\(8\\;21\\)",
                              "FLT3-ITD" = "FLT3-ITD",
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
              helpText("NOTE: Please keep in mind that transcript expression does not always correlate with protein expression.")
          ),
          
          #------------------ MAIN PANEL ------------------#
          div(class = "main-content",
              tabsetPanel(
                tabPanel("Figures",
                         br(),
                         fluidRow(
                           column(12, align = "left",        
                                  plotOutput(ns("histogram"), width = "100%"),
                                  br(),
                                  plotOutput(ns("pieChart"), width = "100%"))
                         )
                ),
                tabPanel("Patient list",
                         fluidRow(
                           box(width = 12, div(style = 'overflow-x: scroll',
                                               DT::dataTableOutput(ns("rankedTable"))))
                         )
                )
              )
          )
      )
    )
  )
}

geneExp <- function(input, output, session, clinData, expData, gene, dataset) {
  
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
    if (dataset() %in% c("SWOG", "BeatAML", "TARGET", "TCGA", "StJude", "GMKF", "PCGP AML", "PCGP", "LEUCEGENE")) {
      expTable <- expData() %>%
        rownames_to_column("Gene") %>%
        filter(Gene == gene()) %>%
        pivot_longer(names_to = "Sample.ID", values_to = "Expression", -Gene) %>%
        filter(!is.na(as.numeric(Expression))) %>%
        mutate(Alterations = clinData()$Filter.Code[match(Sample.ID, clinData()$PatientID)],
               AML.Sample = clinData()$AML.Sample[match(Sample.ID, clinData()$PatientID)]) %>%
        filter(AML.Sample %in% c("AML", "BALL", "TALL"))
      
    } else if (dataset() == "CCLE") {
      expTable <- expData() %>%
        rownames_to_column("Gene") %>%
        filter(Gene == gene()) %>%
        pivot_longer(names_to = "Sample.ID", values_to = "Expression", -Gene) %>%
        filter(!is.na(Expression)) %>%
        mutate(Tissue = clinData()$Tissue[match(Sample.ID, clinData()$PatientID)],
               Malignancy = clinData()$Malignancy[match(Sample.ID, clinData()$PatientID)],
               AML.Sample = clinData()$AML.Sample[match(Sample.ID, clinData()$PatientID)],
               Alterations = clinData()$Filter.Code[match(Sample.ID, clinData()$PatientID)])
      }
    
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
        mutate(Filter.Category = case_when(Quartile == "Q4" ~ paste0(gene(), "+"), 
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
    if (dataset() %in% c("SWOG", "BeatAML", "TARGET", "TCGA", "StJude", "GMKF", "PCGP", "PCGP AML", "LEUCEGENE")) {
      expTable <- expTable %>%
        arrange(desc(Expression)) %>% 
        dplyr::select(Sample.ID, AML.Sample, matches("Protocol"), Gene, Expression, Filter.Category, Alterations)
    } else if (dataset() == "CCLE") {
      expTable <- expTable %>%
        arrange(desc(Expression)) %>%
        dplyr::select(Sample.ID, AML.Sample, Tissue, Malignancy, Gene, Expression, Filter.Category)
    }
    
    return(expTable)
    
  })
  
  
  xint <- reactive({
    
    finalTable <- makeTable()
    
    if (chooseMethod() == "quartile") {
      x <- filter(finalTable, grepl("\\+", Filter.Category)) %>%
        pull(Expression) %>%
        min(., na.rm = T)
    } else if (chooseMethod() == "median") {
      x <- med
    } else {
      # This needs to be cast as a numeric, since the textInput() widget creates a character string
      x <- as.numeric(input$tpm_cutoff) 
    }
    
    return(x)
  })
  
  
  plotHist <- reactive({
    
    expTable <- makeTable()
    x_intercept <- xint()
    
    # rnge <- c(min(expTable$Expression, na.rm = T), max(expTable$Expression, na.rm = T))
    span <- c(min(expTable$Expression, na.rm = T), max(expTable$Expression, na.rm = T))
    span <- span[2] - span[1]
    
    bincount <- case_when(span >= 200 ~ 100,
                          span < 1 ~ 1,
                          TRUE ~ span)
    
    label <- case_when(chooseMethod() == "median" ~ paste0(round(x_intercept, 2), " TPM (median)"),
                       chooseMethod() == "quartile" ~ paste0(round(x_intercept, 2), " TPM (4th quartile)"),
                       TRUE ~ paste0(round(x_intercept, 2), " TPM"))
    
    h <- hist(expTable$Expression, 
         breaks = bincount, 
         xlab = "TPM value", 
         ylab = "Frequency", 
         main = "Distribution of TPM values")
    abline(v = x_intercept, col = "red")
    text(x = x_intercept, 
         y = max(h$counts)*0.75, 
         labels = label, 
         col = "red", pos = 4)
  })
  
  
  
  plotPie <- reactive({
    
    expTable <- makeTable() %>%
      group_by(Filter.Category) %>%
      count(Filter.Category)
    
    ref_lvl <- unique(grep("\\+", expTable$Filter.Category, value = T))
    
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
      labs(title = paste0("Proportion of ", gene(), "+ cases\n(using selected criteria)\n\n")) +
      theme(legend.position = "top", 
            legend.title = element_blank(),
            legend.text = element_text(size = 14),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      scale_fill_manual(values = c("#2096f6", "#FD7370"))
    
    plot
    
  })
  
  #################################################################
  #-------------------- FINAL MODULE OUTPUTS ---------------------#
  #################################################################
  
  # https://glin.github.io/reactable/articles/examples.html#conditional-styling
  output$rankedTable <- DT::renderDataTable({
    t <- makeTable()
    
    DT::datatable(t, 
                  class = "compact nowrap hover row-border order-column", # Defines the CSS formatting of the final table, can string multiple options together
                  extensions = 'Buttons', # See https://rstudio.github.io/DT/extensions.html for more extensions & features
                  options = list(scrollY = "70vh",
                                 dom = 'Bfrtip',
                                 buttons = list(
                                   list(extend = 'excel', filename = paste0(dataset(), "_AML_", gene(), "pos_patientList_generated_", format(Sys.time(), "%m.%d.%Y")))),
                                 scrollX = TRUE,
                                 searchHighlight = TRUE,
                                 pageLength = nrow(t)), 
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
