# UI function for the differentially expressed gene table module
deTableUI <- function(id, label = "Differentially expressed gene table") {
  
  library(DT)
  library(shinyWidgets)
  ns <- NS(id) # Setting a unique namespace for this module
  
  fluidPage(theme = shinytheme("lumen"),
            tags$head(tags$style(HTML('.shiny-output-error-validation { color: #93C54B; }'))), # Custom CSS to modify the app error messages
            
            sidebarLayout(
              position = "left", 
              
              sidebarPanel(
                helpText("The table on the right contains the results of differential expression (DE) analysis performed on pediatric AML subgroups against normal marrow or other AML patients."), 
                helpText("The compared groups are specified in the 'Contrast' column. The reference group used is specified in the 'Reference' column."),
                helpText("The dataset has been filtered to only display DE results for the transcript/gene of interest, aka the RNA species entered by the user in the left sidebar text box."),
                
                br(),
                
                actionButton(ns("col_key"), "Click here for column descriptions", style = 'padding:5px; font-size:80%')
              ),
            
            mainPanel(position = "right", 
                      
                      fluidRow(
                        box(width = 12, div(style = 'overflow-x: scroll',
                                            DT::dataTableOutput(ns("filteredTable")))
                        )
                      )
            )
            
  )
  )
}


# Server function for the differentially expressed gene table module
deTable <- function(input, output, session, gene, table, parent) {
  
  library(tidyverse)
  library(DT)
  library(shinyWidgets)
  
  # Making a function to generate a summary table with outcome data
  tableFun <- reactive({
    
    validate(
      need(gene(), "Please enter a gene symbol or miRNA in the text box to the left.")
    )
    
    table %>%
      mutate(across(where(is.numeric), ~round(., 2))) %>%
      filter(Gene_Symbol == gene()) %>%
      arrange(desc(logFC))
    
  })
  
  # https://glin.github.io/reactable/articles/examples.html#conditional-styling
  output$filteredTable <- DT::renderDataTable({
      DT::datatable(tableFun(), 
                    class = "compact nowrap hover row-border order-column", # Defines the CSS formatting of the final table, can string multiple options together
                    extensions = 'Buttons', # See https://rstudio.github.io/DT/extensions.html for more extensions & features
                    options = list(scrollY = "70vh",
                                   dom = 'Bfrtip',
                                   buttons = c('excel'),
                                   scrollX = TRUE,
                                   fixedColumns = list(leftColumns = 1),
                                   searchHighlight = TRUE,
                                   pageLength = 25), 
                    escape = F) %>%
      DT::formatStyle(columns = c(1,2,4), fontSize = "100%")
  })
  
  observeEvent(input$col_key, {
    showModal( # Tells Shiny to create a modal popup when the corresponding button is activated
      modalDialog(
        title = "Differential Expression Analysis Results Key",
        
        DT::renderDataTable(# Interactively renders the datatable
          DT::datatable(deColKey, # Actually produces the datatable table object
                        class = "stripe row-border", 
                        options = list(dom = "t", pageLength = nrow(deColKey)), 
                        width = "500",
                        rownames = F) %>%
            DT::formatStyle(columns = c(1,2), fontSize = "90%"),
        ),
        easyClose = TRUE)
    )
  })
  
}