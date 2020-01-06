# UI function for the waterfall plot module
wfPlotUI <- function(id, label = "Gene expression plot parameters"){
  
  library(DT)
  ns <- NS(id) # Setting a unique namespace for this module
  
  tagList(
            
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
                                             "Box plots" = "bx")),
                
                helpText("The grouping variable will be used to arrange patients along the x axis (for waterfall plots) 
                          or to group patients together (for box plots) based on a common clinical characteristic 
                          to help highlight expression patterns within the groups." ),
                
                helpText("NOTE: If cell lines is selected, please reference the 'Summary Stats' tab for expression data, 
                         as only one sample is available for most cell lines.") 
              ),
              
              
              ###############################################################
              #----------------------- MAIN PLOT PANEL ---------------------#
              ###############################################################
              
              mainPanel(
                
                position = "right",
                
                tabsetPanel(
                  
                  #-------------------- Waterfall plot -----------------------#
                  tabPanel("Waterfall plot", # This is the title of the tab panel, NOT the name of the plot object!
                           br(),   
                           br(),   # Linebreaks to help center the plot on the page
                           fluidRow(
                             column(10, offset = 0, align = "left",                   # This will be a reactive object that is linked to an item in the output list,
                                    plotOutput(ns("plot"), width = "100%")               # created in the "server" script
                             ), 
                             column(1, offset = 0, align = "right", # One column that takes up the entire width of the panel (all 12 columns)
                                    downloadButton(ns("plot_download"), 
                                                   label = "Download plot", class = NULL)
                             )
                           )
                  ),
                  
                  #-------------------- Summary table -----------------------#
                  tabPanel("Summary stats", 
                           br(),
                           br(),
                           fluidRow(
                             column(10, offset = 0, align = "left", 
                                    DT::dataTableOutput(ns("table")) # Table of summary stats for plot
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
}



# Server dunction for the waterfall plot module
wfPlot <- function(input, output, session, clinData, countsData, gene){
  
  library(tidyverse)
  library(DT)
  `%then%` <- shiny:::`%OR%` # See https://shiny.rstudio.com/articles/validation.html for details on the %then% operator
  
  
  #-------------------- Data preparation -----------------------#
  
    # Cleaning up the clinical data grouping columns so they're displayed nicely in the plots
    clinData <- clinData %>%
      rename(Primary.Cytogenetic.Code = `Primary Cytogenetic Code`) %>%
      rename(Cyto.Risk = `Cyto/Fusion/Molecular Risk`) %>%
      rename(CNS.Disease.Category = `CNS disease at on-study`) %>%
      rename(Event.Type = `EFS event type ID`) %>%
      mutate_at(vars(Cytogenetic.Category.2), ~gsub("\\.|\\_", " ", .)) %>%
      mutate_at(vars(Rare.Fusions), ~gsub("\\.", "\\-", .)) %>%
      mutate_at(vars(SNVs), ~gsub("\\.", " ", .)) %>%
      mutate_at(vars(Cytogenetic.Category.2, Rare.Fusions, SNVs), ~gsub("OtherAML", "Other AML", .)) %>%
      
      # Adding columns to use for re-categorizing MLL and other fusions if they occur in less than 10 patients 
      mutate(Primary.Fusion = ifelse(str_detect(pattern = "KMT2A-", string = `Primary Fusion/CNV`), "KMT2A-X", `Primary Fusion/CNV`)) %>%
      mutate(MLL.Fusion = ifelse(str_detect(pattern = "KMT2A", string = `Primary Fusion/CNV`), `Primary Fusion/CNV`, "Other AML")) %>%
      group_by(MLL.Fusion) %>%
      mutate(Fusion.Count = n()) %>%
      ungroup() %>%
      mutate_at(vars(MLL.Fusion), ~ifelse(Fusion.Count < 10, "Other MLL", .)) %>%
      mutate_at(vars(MLL.Fusion), ~forcats::fct_relevel(., "Other AML", after = Inf)) %>%
      group_by(Primary.Fusion) %>%
      mutate(Fusion.Count = n()) %>%
      ungroup() %>%
      mutate_at(vars(Primary.Fusion), ~ifelse(Fusion.Count < 10, "Other AML", .)) %>%
      
      # Adding rows for the cell lines, which are not typically included in the clinical data
      tibble::add_row(USI = c("K562.01", "K562.02", "ME1", "MO7E", "NOMO1", "Kasumi.D1", "MV4.11.D1"), Group = NA) %>%
      tibble::add_row(USI = grep("^RO|^BM", colnames(countsData), value = T), Group = "NBM") %>%
      tibble::add_row(USI = grep("34POS", colnames(countsData), value = T), Group = "CD34+ NBM")
    
  # Filtering the counts data to only retain the gene of interest & throwing
  # errors if a non-existent gene is provided
  expData <- reactive({
    validate(
      need(gene(), "Please enter a gene symbol in the text box.") %then%
        need(gene() %in% rownames(countsData), "That gene symbol does not exist in the counts data! \nDouble-check the symbol, or try an alias."))
    
    countsData[gene(), intersect(colnames(countsData), clinData$USI)]
  })
  
  # Transforming the counts into a long-format dataframe to use with ggplot
  plotData <- reactive({
    
    plotData <- expData() %>%
      gather(USI, Expression) %>%
      mutate_at(vars(Expression), ~as.numeric(.)) %>%
      left_join(., clinData, by = "USI") %>%
      mutate(Cell.Lines = ifelse(USI %in% c("K562.01", "K562.02", "ME1", "MO7E", "NOMO1", "Kasumi.D1", "MV4.11.D1"), USI, "AML and NBM samples")) %>%
      mutate_at(vars(Cell.Lines), ~forcats::fct_relevel(., "AML and NBM samples", after = Inf)) %>% # Moving the "AML and NBM" samples to 
      group_by_(input$grouping_var) %>%                                                             # the end of the factor order, so it appears greyed out
      arrange_(input$grouping_var, "Expression") # Reordering patients so that the specified groups are 
    # grouped together and ordered by increasing expression
    
    # Setting the patient order using factor levels, so they won't be rearranged 
    # alphabetically by ggplot (this step is required for ggplot waterfall plots)
    plotData$USI <- factor(plotData$USI, levels = plotData$USI)
    
    plotData
  })
  
  # Making a function that will generate the waterfall plot and 
  # can be called from multiple places in the script
  plotFun <- reactive({ 
    
    if(input$plot_type == "bx"){
      
      plotData() %>% drop_na(input$grouping_var) %>%
        ggplot(aes_string(x = input$grouping_var, y = "Expression", fill = input$grouping_var)) +
        theme_classic() +
        labs(x = "Patients", y = "Expression (TPM) \n", fill = gsub("\\.", " ", input$grouping_var)) +
        theme(axis.text.x = element_blank(),
              plot.title = element_text(size = 15, hjust = 0.5),
              axis.ticks = element_blank(),
              legend.position = "bottom") +
        ggtitle(paste0("Expression of ", gene())) +
        geom_boxplot()
      
    }else{
      
      plotData() %>% drop_na(input$grouping_var) %>%
        ggplot(aes_string(x = "USI", y = "Expression", fill = input$grouping_var)) +
        theme_classic() +
        labs(x = "Patients", y = "Expression (TPM) \n", fill = gsub("\\.", " ", input$grouping_var)) +
        theme(axis.text.x = element_blank(),
              plot.title = element_text(size = 15, hjust = 0.5),
              axis.ticks = element_blank(),
              legend.position = "bottom") +
        ggtitle(paste0("Expression of ", gene())) +
        geom_bar(stat = "identity", width = 1, position = position_dodge(width = 0.4))
    }
  })
  
  # Making a function to generate an expression summary table from the plot data
  tableFun <- reactive({
    plotData() %>%
      drop_na(input$grouping_var) %>%
      group_by_(input$grouping_var) %>%
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
    filename = function(){
      paste0("TARGET_AAML1031_", gene(), "_Waterfall_Plot_generated_", format(Sys.time(), "%m.%d.%Y"), ".png")
    }, 
    content = function(file){
      ggsave(filename = file, plot = plotFun(), width = 6, height = 4, device = "png")
    }
  )
  
  #-------------------- Summary table tab -----------------------#
  
  output$table <- DT::renderDataTable({
    DT::datatable(tableFun(), options = list(dom = "t"), autoHideNavigation = T, rownames = F)
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
