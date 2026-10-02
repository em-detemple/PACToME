# PACT-o-ME Pathway Activation On Multiple Eukaryotes using Reactome

library(shiny)
library(DT)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(remotes)
library(shinyWidgets)

getwd() 

# ---- Load Reactome Data ----
load("reactome_data.RData")   
# provides: species_pathway_lists, species_pathway_names, target_species, species_display

# ---- Source the PAA function ----
source("PAA_function.R")

# ---- UI ----
ui <- fluidPage(
  shinyWidgets::setBackgroundColor("#89CFF0"),
  titlePanel(h1("PACToME: Pathway ACTivation Analysis on Multiple Eukaryotes with Reactome🧬")), 
  tags$head(
    tags$style(HTML("
    .nav-tabs > li > a { color: black; }
  "))
  ),
  sidebarLayout(
    sidebarPanel(
      selectInput("species", "Select Organism:",
                  choices = species_display),
      
      hr(),
      
      radioButtons("sig_method", "Significance criteria:",
                   choices = c("Default (FDR < 0.05)" = "default",
                               "Custom" = "custom"),
                   selected = "default"),
      
      conditionalPanel(
        condition = "input.sig_method == 'custom'",
        numericInput("fdr_cutoff", "FDR cutoff:", 
                     value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput("lfc_cutoff", "Absolute log2FC cutoff:", 
                     value = 0.8, min = 0, step = 0.1)
      ),
      
      hr(),
      
      fileInput("file", "Choose TSV File",
                accept = c(".tsv")),
      p("File must contain at least two columns, no header. ",
        "First column = UniProtIDs, second column = log2FC (numeric)."),
      
      # ---- Example data download button (serves the existing file) ----
      br(),
      downloadButton("downloadExample", "📥 Download Example Data (TSV)", 
                     class = "btn-info btn-sm"),
      br(), br(),
      
      hr(),
      
      strong("Processing Info:"),
      textOutput("file_status"),
      
      hr(),
      
      h5("📥 Download Results:"),
      downloadButton("downloadSig", "Significant Pathways (CSV)", 
                     class = "btn-primary btn-sm"),
      br(), br(),
      downloadButton("downloadAll", "All Results (CSV)", 
                     class = "btn-primary btn-sm"),
      br(), br(),
      downloadButton("downloadPlot", "Plot as PNG", 
                     class = "btn-success btn-sm"),
      
      # Manual download button (helper text removed)
      br(), br(),
      downloadButton("downloadManual", "Download Manual (PDF)", 
                     class = "btn-info btn-sm"),
      
      hr()
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Summary", 
                 br(),
                 h4("Run Statistics"),
                 tableOutput("summary_stats")),
        
        tabPanel("Table of Sig. Pathways",
                 br(),
                 h4("Pathways meeting significance criteria"),
                 DTOutput("sig_table")),
        
        tabPanel("Box Plot of Sig. Pathways",
                 br(),
                 h4("Box Plot of Significant Pathways"),
                 p("Pathways are ordered by median log2 fold change. ",
                   "Plot shows up to ten significant pathways sorted by lowest FDR. Box fill color reflects the median logFC (blue = repressed, orange = activated)."),
                 plotOutput("pathway_plot", height = "800px")),
        
        tabPanel("All Pathways Results",
                 br(),
                 h4("All Analyzed Pathways"),
                 DTOutput("all_results_table"))
      )
    )
  )
)

# ---- SERVER ----
server <- function(input, output, session) {
  
  selected_species <- reactive({
    input$species
  })
  
  # ---- Read input file (no header) ----
  uploaded_data <- reactive({
    req(input$file)
    ext <- tools::file_ext(input$file$name)
    if (ext == "tsv") {
      df <- read.delim(input$file$datapath, header = FALSE, stringsAsFactors = FALSE)
    } 
    
    shiny::validate(
      shiny::need(ncol(df) >= 2, "File must contain at least two columns."),
      shiny::need(nrow(df) > 0, "File is empty.")
    )
    
    # Keep first two columns, convert second to numeric
    df <- df[, 1:2, drop = FALSE]
    df[, 2] <- as.numeric(df[, 2])
    
    shiny::validate(
      shiny::need(!all(is.na(df[, 2])), "The second column contains no valid numeric values.")
    )
    
    df
  })
  
  # ---- Run PAA analysis ----
  paa_results <- reactive({
    req(uploaded_data(), selected_species())
    
    pathway_list <- species_pathway_lists[[selected_species()]]
    pathway_names <- species_pathway_names[[selected_species()]]
    
    shiny::validate(
      shiny::need(!is.null(pathway_list), 
                  "No pathway data available for the selected species.")
    )
    
    withProgress(message = 'Running Pathway Activation Analysis...', value = 0.5, {
      results <- pathway_activation_analysis(
        DE_data = uploaded_data(),
        pathway_list = pathway_list,
        pathway_names = pathway_names,
        min_genes = 4
      )
      incProgress(1, detail = "Done")
      results
    })
  })
  
  # ---- Pathway size and calculating species-specific percentages ----
  enhanced_pathway_stats <- reactive({
    req(paa_results(), selected_species())
    stats <- paa_results()$pathway_stats
    pathway_list <- species_pathway_lists[[selected_species()]]
    
    # Compute pathway size (total UniProt IDs in the pathway for this species)
    stats$Pathway_Size <- sapply(stats$Pathway_ID, function(pid) {
      length(pathway_list[[pid]])
    })
    
    # Compute percentages relative to pathway size
    stats$Pct_Input_Genes <- round((stats$Total_Genes / stats$Pathway_Size) * 100, 2)
    stats$Pct_Up <- round((stats$Genes_UP / stats$Pathway_Size) * 100, 2)
    stats$Pct_Down <- round((stats$Genes_DOWN / stats$Pathway_Size) * 100, 2)
    
    stats
  })
  
  # ---- Get user-defined significance thresholds ----
  sig_thresholds <- reactive({
    if (input$sig_method == "default") {
      list(fdr = 0.05, lfc = 0)
    } else {
      list(fdr = input$fdr_cutoff, lfc = input$lfc_cutoff)
    }
  })
  
  # ---- Filter enhanced stats based on thresholds ----
  filtered_pathway_stats <- reactive({
    req(enhanced_pathway_stats())
    stats <- enhanced_pathway_stats()
    th <- sig_thresholds()
    stats %>%
      filter(FDR < th$fdr, abs(MedianlogFC) >= th$lfc)
  })
  
  # ---- SIDEBAR STATUS ----
  output$file_status <- renderText({
    req(input$file)
    paste("File loaded:", input$file$name, 
          "| Total genes:", nrow(uploaded_data()),
          "| Species:", selected_species())
  })
  
  # ---- SUMMARY TAB (corrected) ----
  output$summary_stats <- renderTable({
    req(uploaded_data(), enhanced_pathway_stats(), filtered_pathway_stats())
    stats <- enhanced_pathway_stats()
    sig <- filtered_pathway_stats()
    
    n_sig <- nrow(sig)
    activated <- sum(sig$MedianlogFC > 0, na.rm = TRUE)
    repressed <- sum(sig$MedianlogFC < 0, na.rm = TRUE)
    
    data.frame(
      Metric = c("Total Genes in Uploaded Data",
                 "Total Pathways Analyzed", 
                 "Pathways Meeting Criteria",
                 "Activated Pathways",
                 "Repressed Pathways"),
      Value = c(
        nrow(uploaded_data()),   # fixed
        nrow(stats),             # fixed
        n_sig,                   # fixed
        activated,               # fixed
        repressed                # fixed
      )
    )
  })
  
  # ---- SIGNIFICANT PATHWAYS TABLE ----
  output$sig_table <- renderDT({
    req(filtered_pathway_stats())
    sig_df <- filtered_pathway_stats() %>%
      select(Pathway_ID, Pathway_Name, Pathway_Size, 
             Pct_Input_Genes, Pct_Up, Pct_Down,
             P_value, FDR)
    if (nrow(sig_df) == 0) {
      sig_df <- data.frame(Message = "No pathways meet the specified significance criteria.")
    }
    datatable(sig_df, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # ---- ALL RESULTS TABLE ----
  output$all_results_table <- renderDT({
    req(enhanced_pathway_stats())
    all_df <- enhanced_pathway_stats() %>%
      select(Pathway_ID, Pathway_Name, Pathway_Size, Total_Genes, 
             Pct_Input_Genes, Genes_UP, Pct_Up, Genes_DOWN, Pct_Down,
             MedianlogFC, P_value, FDR)
    datatable(all_df, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # ---- DYNAMIC GGPLOT2 PATHWAY PLOT (font size increased to 16) ----
  output$pathway_plot <- renderPlot({
    req(filtered_pathway_stats())
    sig_df <- filtered_pathway_stats()
    
    if (nrow(sig_df) == 0) {
      ggplot() + 
        annotate("text", x = 0, y = 0, label = "No pathways meet the significance criteria") +
        theme_void()
      return()
    }
    
    if (nrow(sig_df) > 10) {
      sig_df <- sig_df %>% arrange(FDR) %>% slice_head(n = 10)
    }
    
    gene_reg <- paa_results()$UniProtID_regulation
    sig_pathways <- sig_df$Pathway_ID
    
    plot_data <- lapply(sig_pathways, function(pid) {
      data.frame(
        Pathway_ID = pid,
        log2FoldChange = gene_reg$log2FoldChange[gene_reg$Pathway_ID == pid]
      )
    }) %>% bind_rows()
    
    plot_data <- plot_data %>%
      left_join(sig_df %>% select(Pathway_ID, Pathway_Name, MedianlogFC), 
                by = "Pathway_ID")
    
    pathway_order <- sig_df %>%
      arrange(MedianlogFC) %>%
      pull(Pathway_Name)
    
    plot_data$Pathway_Name <- factor(plot_data$Pathway_Name, levels = pathway_order)
    plot_data$Pathway_Name_wrapped <- str_wrap(plot_data$Pathway_Name, width = 35)
    wrapped_order <- unique(plot_data$Pathway_Name_wrapped)
    plot_data$Pathway_Name_wrapped <- factor(plot_data$Pathway_Name_wrapped, 
                                             levels = wrapped_order)
    
    max_abs_median <- max(abs(plot_data$MedianlogFC), na.rm = TRUE)
    
    p <- ggplot(plot_data, aes(x = log2FoldChange, 
                               y = Pathway_Name_wrapped, 
                               fill = MedianlogFC)) +
      geom_boxplot(outlier.shape = 21, outlier.size = 1.5, alpha = 0.8) +
      scale_fill_gradient2(
        low = "#0072B2",
        mid = "white",
        high = "#D55E00",
        midpoint = 0,
        name = "Median log2FC",
        limits = c(-max_abs_median, max_abs_median)
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", size = 0.8) +
      scale_x_continuous(expand = expansion(mult = 0.05)) +
      labs(
        x = "log2 Fold Change",
        y = NULL,
        title = "Boxplot of Significant Pathways"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        # Larger font size for pathway names (16 instead of 14)
        axis.text.y = element_text(size = 16, hjust = 1, margin = margin(r = 5)),
        axis.title.x = element_text(size = 12),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(l = 20, r = 10, t = 10, b = 10),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
      )
    
    p
    
  }, height = function() {
    req(filtered_pathway_stats())
    sig_df <- filtered_pathway_stats()
    if (nrow(sig_df) > 10) {
      sig_df <- sig_df %>% arrange(FDR) %>% slice_head(n = 10)
    }
    n <- nrow(sig_df)
    if (n == 0) return(400)
    height_px <- max(500, (5 + 0.5 * n) * 80)
    return(height_px)
  })
  
  # ---- DOWNLOAD HANDLERS ----
  output$downloadSig <- downloadHandler(
    filename = function() paste0("sig_pathways_", Sys.Date(), ".csv"),
    content = function(file) {
      req(filtered_pathway_stats())
      df <- filtered_pathway_stats() %>%
        select(Pathway_ID, Pathway_Name, Pathway_Size, 
               Pct_Input_Genes, Pct_Up, Pct_Down,
               P_value, FDR)
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  output$downloadAll <- downloadHandler(
    filename = function() paste0("all_pathways_", Sys.Date(), ".csv"),
    content = function(file) {
      req(enhanced_pathway_stats())
      df <- enhanced_pathway_stats() %>%
        select(Pathway_ID, Pathway_Name, Pathway_Size, Total_Genes, 
               Pct_Input_Genes, Genes_UP, Pct_Up, Genes_DOWN, Pct_Down,
               MedianlogFC, P_value, FDR)
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  output$downloadPlot <- downloadHandler(
    filename = function() paste0("pathway_plot_", Sys.Date(), ".png"),
    content = function(file) {
      req(filtered_pathway_stats())
      sig_df <- filtered_pathway_stats()
      
      if (nrow(sig_df) > 10) {
        sig_df <- sig_df %>% arrange(FDR) %>% slice_head(n = 10)
      }
      
      if (nrow(sig_df) == 0) {
        p <- ggplot() + 
          annotate("text", x = 0, y = 0, label = "No pathways meet the significance criteria") +
          theme_void()
        ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
        return()
      }
      
      gene_reg <- paa_results()$UniProtID_regulation
      sig_pathways <- sig_df$Pathway_ID
      
      plot_data <- lapply(sig_pathways, function(pid) {
        data.frame(
          Pathway_ID = pid,
          log2FoldChange = gene_reg$log2FoldChange[gene_reg$Pathway_ID == pid]
        )
      }) %>% bind_rows()
      
      plot_data <- plot_data %>%
        left_join(sig_df %>% select(Pathway_ID, Pathway_Name, MedianlogFC), 
                  by = "Pathway_ID")
      
      pathway_order <- sig_df %>%
        arrange(MedianlogFC) %>%
        pull(Pathway_Name)
      
      plot_data$Pathway_Name <- factor(plot_data$Pathway_Name, levels = pathway_order)
      plot_data$Pathway_Name_wrapped <- str_wrap(plot_data$Pathway_Name, width = 35)
      wrapped_order <- unique(plot_data$Pathway_Name_wrapped)
      plot_data$Pathway_Name_wrapped <- factor(plot_data$Pathway_Name_wrapped, 
                                               levels = wrapped_order)
      
      max_abs_median <- max(abs(plot_data$MedianlogFC), na.rm = TRUE)
      
      p <- ggplot(plot_data, aes(x = log2FoldChange, 
                                 y = Pathway_Name_wrapped, 
                                 fill = MedianlogFC)) +
        geom_boxplot(outlier.shape = 21, outlier.size = 1.5, alpha = 0.8) +
        scale_fill_gradient2(
          low = "#0072B2",
          mid = "white",
          high = "#D55E00",
          midpoint = 0,
          name = "Median log2FC",
          limits = c(-max_abs_median, max_abs_median)
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", size = 0.8) +
        scale_x_continuous(expand = expansion(mult = 0.05)) +
        labs(
          x = "log2 Fold Change",
          y = NULL,
          title = "Significant Pathways"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          # Larger font for downloaded plot as well (16)
          axis.text.y = element_text(size = 16, hjust = 1, margin = margin(r = 5)),
          axis.title.x = element_text(size = 13),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(l = 20, r = 10, t = 10, b = 10),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
        )
      
      n <- nrow(sig_df)
      height_in <- max(5, 5 + 0.3 * n)
      ggsave(file, plot = p, width = 10, height = height_in, dpi = 300)
    }
  )
  
  # ---- Manual download handler ----
  output$downloadManual <- downloadHandler(
    filename = function() "PACToME_UserManual.pdf",
    content = function(file) {
      if (file.exists("PACToME_UserManual.pdf")) {
        file.copy("PACToME_UserManual.pdf", file)
      } else {
        # Optionally create a temporary file with a message
        writeLines("Manual not yet available. Please contact the author.", con = file)
      }
    }
  )
  
  # ---- Example data download handler (serves the existing file) ----
  output$downloadExample <- downloadHandler(
    filename = function() "Worm_Example_Data.tsv",
    content = function(file) {
      # Copy the existing example file from the app directory.
      # If the file is not found, a fallback can be added, but we assume it exists.
      if (file.exists("Worm_Example_Data.tsv")) {
        file.copy("Worm_Example_Data.tsv", file)
      } else {
        # Optional: write an error message to the file
        writeLines("Example data file 'Worm_Example_Data.tsv' not found in the app directory.", con = file)
      }
    }
  )
  
}

shinyApp(ui, server)

