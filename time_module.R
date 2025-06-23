# UI function for the waterfall plot module
timePlotUI <- function(id, label = "Gene expression timepoint plot parameters") {
  
  ns <- NS(id) # Setting a unique namespace for this module
  
  # Creating a list of dropdown choices for the plot type selection
  choices <- c("All Samples", "Fusions", "SNVs", "CR EOI1")
  
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
      theme = shinythemes::shinytheme(theme = "paper"),
      div(class = "sidebar-container",
          
          # Sidebar
          div(class = "custom-sidebar",
              
              radioButtons(ns("timepoint_subset"),
                           label = "Select timepoint group:",
                           choices = c("All (Dx + Rem + Rel) only" = "all",
                                       "Dx vs Rem only" = "dxrem",
                                       "Dx vs Rel only" = "dxrel"),
                           selected = "all"),
              
              br(),
              
              selectInput(ns("grouping_var"),
                          label = "Select a grouping variable",
                          choices = choices),
              
              textInput(ns("tpm_cutoff"),                                                  
                        label = "Use TPM cutoff value",
                        placeholder = "Example: 5"),
              
              checkboxGroupInput(ns("timepoint_cutoff"), 
                                 label = "Select one or both options:",
                                 choices = list("Positive at Dx" = "Dx",
                                                "Positive at Rel" = "Rel"),
                                 selected = "Dx"),
              
              br(),
              
              checkboxInput(ns("avg"), "Show Average Expression", FALSE),
              
              br(),
              
              # Plot download button
              div(style = "margin-bottom: 10px;",
                  downloadButton(ns("plot_download"),
                                 label = "Download Plot",
                                 class = "btn-primary btn-sm w-100")),
              shinyBS::bsTooltip(ns("plot_download"),
                                 title = "Click here to download the plot",
                                 placement = "right",
                                 trigger = "hover")
          ),
          
          # Main plot panel
          div(class = "main-content",
              mainPanel(
                position = "right",
                width = 12,
                tabsetPanel(
                  
                  # Plot tab
                  tabPanel("Plot",
                           br(),
                           br(),
                           fluidRow(
                             column(12, offset = 0, align = "left",
                                    plotlyOutput(ns("plot"), height = "70vh")
                             )
                           )
                  ),
                  
                  # Summary table tab
                  tabPanel("Summary stats",
                           br(),
                           br(),
                           fluidRow(
                             column(12, offset = 0, align = "left",
                                    DT::dataTableOutput(ns("table"))
                             )
                           )
                  )
                )
              )
          )
      )
    )
  )
}


# Server function for the waterfall plot module
timePlot <- function(input, output, session, clinData, expData, gene, aligner, dataset){

  plotData <- reactive({
    validate(
      need(dataset() == "TARGET" && aligner() == "star", 
           "This module is only available for TARGET data aligned with STAR"),
      need(gene(), "Please enter a gene symbol or miRNA in the text box to the left."),
      need(gene() %in% target_expData38_STAR_Dx$name, 
           paste0(gene(), " does not exist in the counts data!"))
    )
    
    gene_id <- gene()
    
    dx_gene <- target_expData38_STAR_Dx %>% filter(name == gene_id)
    rem_gene <- target_expData38_STAR_Rem %>% filter(name == gene_id)
    rel_gene <- target_expData38_STAR_Rel %>% filter(name == gene_id)
    
    dx_long <- dx_gene %>%
      dplyr::select(-id, -name) %>%
      pivot_longer(cols = everything(), names_to = "PatientID", values_to = "TPM") %>%
      mutate(Timepoint = "Diagnostic")
    
    rem_long <- rem_gene %>%
      dplyr::select(-id, -name) %>%
      pivot_longer(cols = everything(), names_to = "PatientID", values_to = "TPM") %>%
      mutate(Timepoint = "Remission")
    
    rel_long <- rel_gene %>%
      dplyr::select(-id, -name) %>%
      pivot_longer(cols = everything(), names_to = "PatientID", values_to = "TPM") %>%
      mutate(Timepoint = "Relapse")
    
    gene_df <- bind_rows(dx_long, rem_long, rel_long)
    
    if (input$timepoint_subset == "all") {
    
    manifest <- clinData() %>%
      filter(Disease.Group == "AML", ExpData == "Yes") %>%
      group_by(USI) %>%
      filter(all(c("Diagnostic", "Relapse", "Remission") %in% Timepoint)) %>%
      ungroup() %>%
      group_by(USI, Timepoint) %>%
      filter(!(Timepoint == "Remission" & n() > 1 & !grepl("EOI1", PatientID))) %>%
      ungroup() %>%
      dplyr::select(USI, PatientID, Timepoint, CR_EOI1, Primary.Fusion, WT1.Mutation, NPM1.Mutation, FLT3.ITD, CEBPA.Mutation)
    
    } else if (input$timepoint_subset == "dxrem") {
      
      manifest <- clinData() %>%
        filter(Disease.Group == "AML", ExpData == "Yes") %>%
        group_by(USI) %>%
        filter(all(c("Diagnostic", "Remission") %in% Timepoint)) %>%
        ungroup() %>%
        group_by(USI, Timepoint) %>%
        filter(!(Timepoint == "Remission" & n() > 1 & !grepl("EOI1", PatientID))) %>%
        ungroup() %>%
        dplyr::select(USI, PatientID, Timepoint, CR_EOI1, Primary.Fusion, WT1.Mutation, NPM1.Mutation, FLT3.ITD, CEBPA.Mutation)
      
      print(table(manifest$Timepoint))
      
    } else {
      
      manifest <- clinData() %>%
        filter(Disease.Group == "AML", ExpData == "Yes") %>%
        group_by(USI) %>%
        filter(all(c("Diagnostic", "Relapse") %in% Timepoint)) %>%
        ungroup() %>%
        dplyr::select(USI, PatientID, Timepoint, CR_EOI1, Primary.Fusion, WT1.Mutation, NPM1.Mutation, FLT3.ITD, CEBPA.Mutation)
      
    }
    
    if (input$grouping_var == "SNVs") {
      mutation_cols <- c("WT1.Mutation", "NPM1.Mutation", "FLT3.ITD", "CEBPA.Mutation")
      
      # Step 1: Pivot longer
      manifest_long <- manifest %>%
        pivot_longer(cols = all_of(mutation_cols), names_to = "SNV_raw", values_to = "MutationStatus") %>%
        mutate(SNVs = gsub("\\.Mutation", "", SNV_raw)) %>%
        filter(MutationStatus == "Yes")
      
      # 2. Drop existing USI (if present) to prevent conflicts
      manifest_long <- manifest_long %>%
        dplyr::select(-any_of("USI"))
      
      # Step 2: Count mutations per PatientID + Timepoint
      mutation_counts <- manifest_long %>%
        group_by(PatientID, Timepoint) %>%
        summarise(SNV_count = n(), .groups = "drop")
      
      # Step 3: Join counts to manifest_long
      manifest_long <- manifest_long %>%
        left_join(mutation_counts, by = c("PatientID", "Timepoint"))
      
      # Step 4: Join USI
      manifest_long <- manifest_long %>%
        left_join(manifest %>% dplyr::select(PatientID, USI), by = "PatientID") %>%
        mutate(USI_with_SNV = paste0(USI, "_", SNVs))
      
      # Step 5: Join gene_df
      df <- gene_df %>%
        filter(PatientID %in% manifest$PatientID) %>%
        left_join(manifest_long %>%
                    dplyr::select(PatientID, Timepoint, SNVs, USI_with_SNV), 
                  by = c("PatientID", "Timepoint")) %>%
        left_join(manifest %>% dplyr::select(PatientID, USI), by = "PatientID") %>%
        mutate(
          SNVs = ifelse(is.na(SNVs), "Other AML", SNVs),
          USI = ifelse(is.na(USI_with_SNV), USI, USI_with_SNV),
          Gene = gene_id
        ) %>%
        dplyr::select(-USI_with_SNV)
      
    } else {
      df <- gene_df %>%
        filter(PatientID %in% manifest$PatientID) %>%
        left_join(manifest, by = c("PatientID", "Timepoint")) %>%
        mutate(Gene = gene_id)
    }

    tpm_cut <- suppressWarnings(as.numeric(input$tpm_cutoff))
    
    # Map user timepoint selection codes to full Timepoint names
    tp_map <- c("Dx" = "Diagnostic", "Rel" = "Relapse")
    
    selected_timepoints <- tp_map[input$timepoint_cutoff]
    
    if (!is.na(tpm_cut) && length(selected_timepoints) > 0) {
      # For each patient, check if TPM >= cutoff at ALL selected timepoints
      qualifying_patients <- df %>%
        filter(Timepoint %in% selected_timepoints) %>%
        group_by(USI) %>%
        summarize(all_pass = all(TPM >= tpm_cut)) %>%
        filter(all_pass) %>%
        pull(USI)
      
      validate(
        need(length(qualifying_patients) > 0, "No samples meet the threshold.")
      )
      
      # Keep all timepoints for qualifying patients
      df <- df %>% filter(USI %in% qualifying_patients)
    }
    
    return(df)
  })
  
  output$plot <- renderPlotly({
    df <- plotData()
    
    # Font settings
    font_family <- "Helvetica"
    base_size <- 16
    legend_text_size <- base_size - 4
    axis_title_size <- base_size
    axis_text_size <- base_size
    plot_title_size <- base_size
    
    # Ensure Timepoint is factor with correct order
    df$Timepoint <- factor(df$Timepoint, levels = switch(input$timepoint_subset,
                                                         "all" = c("Diagnostic", "Remission", "Relapse"),
                                                         "dxrem" = c("Diagnostic", "Remission"),
                                                         c("Diagnostic", "Relapse")
    ))
    
    if (input$grouping_var == "Fusions") {
      df$color_var <- df$Primary.Fusion
      legend_title <- "Primary Fusion"
      hover_text <- paste0("USI: ", df$USI,
                           "<br>TPM: ", round(df$TPM, 2),
                           "<br>Fusion: ", df$Primary.Fusion)
    } else if (input$grouping_var == "SNVs") {
      df$color_var <- df$SNVs
      legend_title <- "SNVs"
      hover_text <- paste0("USI: ", df$USI,
                           "<br>TPM: ", round(df$TPM, 2),
                           "<br>SNV: ", df$SNVs)
    } else if (input$grouping_var == "CR EOI1") {
      df$color_var <- df$CR_EOI1
      legend_title <- "CR EOI1"
      hover_text <- paste0("USI: ", df$USI,
                           "<br>TPM: ", round(df$TPM, 2),
                           "<br>CR: ", df$CR_EOI1)
    } else {
      df$color_var <- "All"
      legend_title <- "Samples"
      hover_text <- paste0("USI: ", df$USI,
                           "<br>TPM: ", round(df$TPM, 2))
    }
    
    
    legend_title <- switch(input$grouping_var,
                           "Fusions" = "Primary Fusion",
                           "SNVs" = "SNVs",
                           "All Samples" = "Samples",
                           "CR EOI1" = "CR EOI1"
    )
    
    df$hover_text <- switch(input$grouping_var,
                            "Fusions" = paste0("USI: ", df$USI, "<br>TPM: ", round(df$TPM, 2), "<br>Fusion: ", df$Primary.Fusion),
                            "SNVs" = paste0("USI: ", df$USI, "<br>TPM: ", round(df$TPM, 2), "<br>SNV: ", df$SNVs),
                            "CR EOI1" = paste0("USI: ", df$USI, "<br>TPM: ", round(df$TPM, 2), "<br>CR: ", df$CR_EOI1),
                            "All Samples" = paste0("USI: ", df$USI, "<br>TPM: ", round(df$TPM, 2))
    )
    
    # Palette
    groups <- unique(df$color_var)
    palette_colors <- RColorBrewer::brewer.pal(min(length(groups), 8), "Set2")
    color_map <- setNames(palette_colors[seq_along(groups)], groups)
    print(color_map)
    
    # Initialize empty plot
    p <- plot_ly()
    
    # Add traces for each USI manually
    for (usi in unique(df$USI)) {
      sub_df <- df[df$USI == usi, ]
      group <- unique(sub_df$color_var)
      p <- add_trace(p,
                     data = sub_df,
                     x = ~Timepoint,
                     y = ~TPM,
                     type = 'scatter',
                     mode = 'lines+markers',
                     name = group,
                     text = ~hover_text,
                     hoverinfo = "text",
                     line = list(shape = 'linear', color = color_map[[group]]),
                     marker = list(size = 8, color = color_map[[group]]),
                     legendgroup = group,
                     showlegend = usi == df$USI[df$color_var == group][1]  # Show legend only once per group
      )
    }
    
    # Add average line if requested
    if (isTRUE(input$avg)) {
      mean_df <- aggregate(TPM ~ Timepoint, data = df, FUN = mean)
      mean_df$Timepoint <- factor(mean_df$Timepoint, levels = levels(df$Timepoint))
      mean_df <- mean_df[order(mean_df$Timepoint), ]
      p <- add_trace(p,
                     data = mean_df,
                     x = ~Timepoint,
                     y = ~TPM,
                     type = 'scatter',
                     mode = 'lines+markers',
                     name = 'Mean TPM',
                     line = list(color = 'black', width = 3, dash = 'dash'),
                     marker = list(size = 10, color = 'black'),
                     inherit = FALSE
      )
    }
    
    # Layout
    p <- p %>% layout(
      title = list(
        text = paste0("Gene Expression of ", unique(df$Gene), " Across Timepoints"),
        font = list(family = font_family, size = plot_title_size)
      ),
      xaxis = list(
        title = list(text = "Timepoint", font = list(family = font_family, size = axis_title_size)),
        tickfont = list(family = font_family, size = axis_text_size)
      ),
      yaxis = list(
        title = list(text = "TPM", font = list(family = font_family, size = axis_title_size)),
        tickfont = list(family = font_family, size = axis_text_size)
      ),
      legend = list(
        title = list(text = legend_title, font = list(family = font_family, size = legend_text_size)),
        font = list(family = font_family, size = legend_text_size),
        orientation = "v",
        x = 1.05,
        y = 0.5,
        yanchor = "middle",
        xanchor = "left"
      )
    )
    
    p
  })
  
  
  output$table <- DT::renderDataTable({
    
    DT::datatable(plotData(),
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
                  escape = F)
  })
  
  
  
  
  
  
  
  
  
}


