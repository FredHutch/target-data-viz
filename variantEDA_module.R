# UI function for the differentially expressed gene table module
wgsTriosUI <- function(id, label = "Differentially expressed gene table") {
  
  library(DT)
  library(shinyWidgets)
  ns <- NS(id) # Setting a unique namespace for this module
  
  fluidPage(theme = shinytheme("lumen"),
            tags$head(tags$style(HTML('.shiny-output-error-validation { color: #93C54B; }'))), # Custom CSS to modify the app error messages
            
            sidebarLayout(
              position = "left", 
              
              sidebarPanel(
                helpText("placeholder"),
              ),
              
              mainPanel(position = "right", 
                        
                        fluidRow(
                          box(width = 12, div(style = 'overflow-x: scroll'))
                          )
                        )
              )
              
            )
  )
}


# Server function for the differentially expressed gene table module
wgsTriosTable <- function(input, output, session, gene, table, parent) {
  
  library(tidyverse)
  library(DT)
  library(shinyWidgets)
  `%then%` <- shiny:::`%OR%` # See https://shiny.rstudio.com/articles/validation.html for details on the %then% operator
  
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