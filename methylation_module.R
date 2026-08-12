# UI function for the methylation module
methylPlotUI <- function(id, label = "Methylation plot parameters"){
  
  ns <- NS(id)
  
  filtered_choices <- filter(colMapping, Final_Column_Label %in% c("AML vs. normal",
                                                                   "Age category",
                                                                   "Risk classification",
                                                                   "Primary cytogenetics",
                                                                   "Other cytogenetics",
                                                                   "Fusion", 
                                                                   "CNV",
                                                                   "Rare fusions",
                                                                   "KMT2A/MLL fusions",
                                                                   "NUP98 fusions"))
  
  choices <- as.list(filtered_choices$Final_Column_Name)
  names(choices) <- filtered_choices$Final_Column_Label
  
  tagList(
    useShinyjs(),
    tags$head(
      tags$style(HTML("
          .sidebar-container { display: flex; min-height: 100vh; }
          .custom-sidebar { background-color: #f8f9fa; padding: 15px; width: 250px; flex-shrink: 0; }
          .main-content { flex-grow: 1; padding: 15px; }
          .selectize-dropdown-content .option.active { background-color: #2096f6 !important; color: white !important; }
          .selectize-dropdown-content .option.selected { background-color: #2096f6 !important; color: white !important; }
          
          /* Disable Shiny's default fade-to-translucent while recalculating */
          .shiny-plot-output.recalculating {
            opacity: 1 !important;
            transition: none !important;
          }
      "))
    ),
    fluidPage(
      theme = shinythemes::shinytheme(theme = "paper"),
      div(class = "sidebar-container",
          div(class = "custom-sidebar",
              selectInput(ns("grouping_var"),
                          label = "Select a grouping variable",
                          choices = choices),
              
              # checkboxGroupInput(ns("relation_to_island"),
              #                    label = "Include probes with Relation to Island:",
              #                    choices = c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea"),
              #                    selected = c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea")),
              
              shinyWidgets::pickerInput(
                inputId  = ns("selected_groups"),
                label    = "Select or De-select Groups to Plot",
                choices  = NULL,
                multiple = TRUE,
                options  = shinyWidgets::pickerOptions(
                  actionsBox = TRUE,
                  liveSearch = FALSE,
                  selectedTextFormat = "count > 3",
                  countSelectedText  = "{0} groups selected"
                )
              ),
              
              tags$hr(style = "margin: 15px 0;"),
              
              helpText("Shows EPIC array promoter-region (TSS200/TSS1500/1stExon) methylation
                        (Beta values) plotted against expression (TPM) for the gene entered
                        in the sidebar, faceted by block and colored by the grouping variable
                        selected above.")
          ),
          
          div(class = "main-content",
              mainPanel(
                position = "right",
                width = 10,
                tabsetPanel(
                  tabPanel("Plot",
                           br(), br(),
                           fluidRow(
                             column(12, offset = 0, align = "left",
                                    plotOutput(ns("plot"), width = "100%", height = "70vh"))
                           )
                  ),
                  tabPanel("Summary stats",
                           br(), br(),
                           fluidRow(
                             column(12, offset = 0, align = "left",
                                    h4("Per-Block Correlation (Beta vs. Expression)"),
                                    DT::dataTableOutput(ns("corTable")),
                                    br(),
                                    h4("Group Averages by Block"),
                                    DT::dataTableOutput(ns("table")))
                           )
                  ),
                )
              )
          )
      )
    )
  )
}


# Server function for the methylation module
methylPlot <- function(input, output, session, clinData, expData, gene, aligner, dataset, parent, timepoint){
  
  # ---- Step 1: find promoter-region EPIC probes for the gene, matching the
  #      user-selected Relation_to_Island categories ----
  probeAnno <- reactive({
    
    validate(
      need(gene(), "Please enter a gene symbol in the text box to the left.")
    )
    
    # First check the gene exists in the EPIC annotation at all (regardless of
    # region/island filters), so we can give a specific error if not.
    gene_present <- epic_anno %>%
      filter(grepl(paste0("(^|;)", gene(), "(;|$)"), UCSC_RefGene_Name))
    
    validate(
      need(nrow(gene_present) > 0,
           paste0(gene(), " was not found in the EPIC array annotation.\nDouble-check the symbol, or try an alias/synonym."))
    )
    
    anno_sub <- gene_present %>%
      filter(grepl("TSS200|TSS1500|1stExon", UCSC_RefGene_Group))
    
    validate(
      need(nrow(anno_sub) > 0,
           paste0("No promoter-region (TSS200/TSS1500/1stExon) probes were found for ", gene(),
                  " with the selected Relation to Island filter(s). Try adding more categories."))
    )
    
    anno_sub
  })
  
  # ---- Step 2: match those probe positions to the nearest bvals block/region ----
  blockIds <- reactive({
    anno_sub <- probeAnno()
    
    gene_gr <- GRanges(
      seqnames = anno_sub$chr,
      ranges   = IRanges(start = anno_sub$pos, width = 1)
    )
    
    bval_ids <- rownames(bvals)
    bval_gr <- GRanges(bval_ids)
    names(bval_gr) <- bval_ids
    
    nearest_hits <- distanceToNearest(gene_gr, bval_gr)
    
    ids <- unique(names(bval_gr)[subjectHits(nearest_hits)])
    
    validate(
      need(length(ids) > 0,
           paste0("No matching methylation blocks were identified for ", gene(), " in the Beta value data."))
    )
    
    ids
  })
  
  # ---- Step 3: build long-format df: Beta values x Expression x clinical data, joined by USI ----
  plotData <- reactive({
    
    validate(
      need(dataset() == "TARGET", "Methylation data is only available for the TARGET cohort."),
      need(aligner() == "star", "Methylation data is only available for the STAR-aligned TARGET data."),
      need(timepoint() == "diagnostic", "Methylation data is only available for the diagnostic timepoint."),
    )
    
    req(dataset() == "TARGET", aligner() == "star", timepoint() == "diagnostic")
    
    map <- clinData()[, c("USI", "PatientID")]
    
    ids <- blockIds()
    
    validate(
      need(gene() %in% rownames(expData()),
           paste0(gene(), " does not exist in the counts data!\nDouble-check the symbol or ID, or try an alias/synonym."))
    )
    
    beta_sub <- bvals[ids, , drop = FALSE]
    
    beta_long <- as.data.frame(beta_sub) %>%
      mutate(Block = rownames(beta_sub)) %>%
      pivot_longer(-Block, names_to = "USI", values_to = "Beta") %>%
      mutate(USI = sub("_.*", "", USI))
    
    exp_long <- expData()[gene(), , drop = FALSE] %>%
      as.data.frame() %>%
      rownames_to_column("Gene") %>%
      pivot_longer(-Gene, names_to = "PatientID", values_to = "Expression") %>%
      dplyr::select(-Gene) %>%
      mutate(Expression = as.numeric(Expression)) %>%
      inner_join(map, by = "PatientID") %>%
      dplyr::select(-PatientID)
    
    plot_df <- beta_long %>%
      inner_join(exp_long, by = "USI") %>%
      left_join(clinData(), by = "USI") %>%
      filter(!is.na(Expression), is.finite(Beta)) %>%
      mutate(Log2 = log2(Expression + 1))
    
    validate(
      need(nrow(plot_df) > 0,
           "No samples with both methylation and expression data were found for this gene.")
    )
    
    plot_df
  })
  
  # ---- Picker sync for grouping var ----
  picker_synced_to <- reactiveVal(NULL)
  
  observe({
    req(plotData())
    grps <- levels(factor(plotData()[[input$grouping_var]]))
    if (!"NBM" %in% grps) {
      grps <- c(grps, "NBM")
    }
    shinyWidgets::updatePickerInput(
      session  = session,
      inputId  = "selected_groups",
      choices  = grps,
      selected = grps
    )
    picker_synced_to(input$grouping_var)
  })
  
  filteredData <- reactive({
    req(plotData())
    req(picker_synced_to() == input$grouping_var)
    
    selected <- if (is.null(input$selected_groups) || length(input$selected_groups) == 0) {
      levels(factor(plotData()[[input$grouping_var]]))
    } else {
      input$selected_groups
    }
    
    valid_grps <- levels(factor(plotData()[[input$grouping_var]]))
    selected <- intersect(selected, valid_grps)
    req(length(selected) > 0)
    
    base <- plotData() %>%
      filter(.data[[input$grouping_var]] %in% selected)
    
    # Always include NBM samples, regardless of picker selection,
    # and force their grouping value to "NBM" so they're never NA/unlabeled
    nbm <- plotData() %>%
      filter(Disease.Group == "NBM") %>%
      mutate(!!input$grouping_var := "NBM")
    
    if (input$grouping_var == "AML.Sample") {
      base <- base %>% filter(!AML.Sample %in% c("LSNBM", "MSNBM"))
    }
    
    bind_rows(base, nbm) %>%
      distinct()
  }) %>% debounce(200)
  
  # Add/edit entries here as needed — key = column name (Final_Column_Name), 
  # value = desired level order (NOT including "NBM" — it's always appended last automatically)
  group_level_orders <- list(
    "Disease.Group" = c("AML", "NBM"),
    "Age.Category"  = c("Less than 3 years", "Between 3 and 5 years", "Between 5 and 10 years",
                        "Between 10 and 18 years", "Greater than 18 years", "Unknown"),
    "Risk"          = c("High", "Standard", "Low")
    # e.g. "MLL.Fusion" = c("KMT2A-ELL", "KMT2A-MLLT3", "KMT2A-MLLT4", "KMT2A-MLLT10"),
  )
  
  orderGroupLevels <- function(x, grouping_var) {
    x <- as.character(x)
    present <- unique(x)
    
    manual_order <- group_level_orders[[grouping_var]]
    
    if (!is.null(manual_order)) {
      ordered_levels <- intersect(manual_order, present)
      leftover <- sort(setdiff(present, c(manual_order, "NBM")))
      all_levels <- c(ordered_levels, leftover)
    } else {
      all_levels <- sort(setdiff(present, "NBM"))
    }
    
    if ("NBM" %in% present) {
      all_levels <- c(setdiff(all_levels, "NBM"), "NBM")
    }
    
    factor(x, levels = all_levels)
  }

  # ---- Per-block Beta vs Expression correlation (overall) — rounded ----
  corTable <- reactive({
    filteredData() %>%
      group_by(Block) %>%
      summarise(
        n     = n(),
        rho   = suppressWarnings(cor.test(Beta, Expression, method = "spearman")$estimate),
        p_val = suppressWarnings(cor.test(Beta, Expression, method = "spearman")$p.value),
        .groups = "drop"
      ) %>%
      mutate(
        rho   = round(rho, 2),
        p_val = round(p_val, 4)
      )
  })
  
  # ---- Per-block, per-group average Beta and TPM — Groups as rows, Blocks as columns ----
  groupAverages <- reactive({
    df <- filteredData() %>%
      group_by(Block, .data[[input$grouping_var]]) %>%
      summarise(
        n           = n(),
        `Mean Beta` = round(mean(Beta, na.rm = TRUE), 3),
        `Mean TPM`  = round(mean(Expression, na.rm = TRUE), 2),
        .groups = "drop"
      ) %>%
      dplyr::rename(Group = !!input$grouping_var)
    
    wide <- df %>%
      pivot_wider(
        id_cols     = Group,
        names_from  = Block,
        values_from = c(n, `Mean Beta`, `Mean TPM`),
        names_glue  = "{Block} {.value}"
      ) %>%
      arrange(Group)
    
    # explicitly enforce n -> Mean Beta -> Mean TPM order per Block
    blocks <- unique(df$Block)
    ordered_cols <- c("Group", unlist(lapply(blocks, function(b) paste0(b, c(" n", " Mean Beta", " Mean TPM")))))
    
    wide[, intersect(ordered_cols, colnames(wide))]
  })
  
  facetLabels <- reactive({
    corTable() %>%
      mutate(label = paste0(Block, "\nrho=", rho, ", p=", p_val))
  })
  
  # ---- Plot ----
  plotFun <- reactive({
    
    df <- filteredData() %>%
      left_join(facetLabels() %>% dplyr::select(Block, label), by = "Block")
    
    df[[input$grouping_var]] <- orderGroupLevels(df[[input$grouping_var]], input$grouping_var)
    
    grp_levels <- levels(df[[input$grouping_var]])
    req(length(grp_levels) > 0)
    
    gg_colors <- scales::hue_pal()(length(grp_levels))
    names(gg_colors) <- grp_levels
    
    centroid_df <- df %>%
      group_by(Block, label, .data[[input$grouping_var]]) %>%
      summarise(
        mean_Beta = mean(Beta, na.rm = TRUE),
        mean_log2 = mean(Log2, na.rm = TRUE),
        se_Beta   = sd(Beta, na.rm = TRUE) / sqrt(n()),
        se_log2   = sd(Log2, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      )
    
    ggplot(df, aes(x = Beta, y = Log2, color = .data[[input$grouping_var]])) +
      geom_point(size = 3, alpha = 0.75) +
      stat_ellipse(aes(group = .data[[input$grouping_var]]), type = "norm", level = 0.95, linewidth = 0.8) +
      geom_errorbar(data = centroid_df,
                    aes(x = mean_Beta, y = mean_log2,
                        ymin = mean_log2 - se_log2, ymax = mean_log2 + se_log2),
                    width = 0.02, linewidth = 0.9, inherit.aes = FALSE) +
      geom_errorbarh(data = centroid_df,
                     aes(x = mean_Beta, y = mean_log2,
                         xmin = mean_Beta - se_Beta, xmax = mean_Beta + se_Beta),
                     height = 0.2, linewidth = 0.9, inherit.aes = FALSE) +
      geom_point(data = centroid_df,
                 aes(x = mean_Beta, y = mean_log2, fill = .data[[input$grouping_var]]),
                 shape = 23, size = 5, color = "black", stroke = 1, inherit.aes = FALSE) +
      facet_wrap(~ label, scales = "free_y") +
      scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
      coord_cartesian(ylim = c(0, NA)) +
      labs(
        title = paste0(gene(), " Promoter Methylation (Beta) vs Expression, by Block"),
        subtitle = "Diamonds = group means (\u00b1 SE); ellipses = per-group 95% CI",
        x = "Beta Value",
        y = paste0("log2(", gene(), " TPM + 1)"),
        color = gsub("\\.", " ", input$grouping_var),
        fill  = gsub("\\.", " ", input$grouping_var)
      ) +
      theme_bw(base_size = bs) +
      theme(
        plot.title      = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle   = element_text(hjust = 0.5, face = "bold"),
        strip.text      = element_text(size = bs + 2, face = "bold"),
        legend.text     = element_text(size = bs + 2),
        legend.title    = element_text(size = bs + 2, face = "bold")
      ) +
      scale_color_manual(values = gg_colors) +
      scale_fill_manual(values = gg_colors)
  })
  
  output$plot <- renderPlot({
    req(filteredData())
    Sys.sleep(0.5)
    plotFun()
  })
  
  output$corTable <- DT::renderDataTable({
    DT::datatable(corTable(),
                  options = list(scrollY = "10vh", pageLength = 5, searchHighlight = TRUE, scrollX = TRUE),
                  escape = F)
  })
  output$table <- DT::renderDataTable({
    
    df <- groupAverages()
    
    metric_cols <- setdiff(colnames(df), "Group")
    blocks <- unique(sub(" (n|Mean Beta|Mean TPM)$", "", metric_cols))
    
    sketch <- htmltools::withTags(table(
      class = "display",
      thead(
        tr(
          th(rowspan = 2, "Group"),
          lapply(blocks, function(b) th(colspan = 3, b))
        ),
        tr(
          lapply(blocks, function(b) {
            tagList(
              th("n"),
              th("Mean Beta"),
              th("Mean TPM")
            )
          })
        )
      )
    ))
    
    DT::datatable(df,
                  container = sketch,
                  rownames = FALSE,
                  options = list(scrollY = "50vh", pageLength = 10, searchHighlight = TRUE, scrollX = TRUE),
                  escape = F)
  })
}
