# PATHWAY ACTIVATION ANALYSIS FUNCTION
# Expects DE_data with exactly two columns:
#   Column 1: UniProt IDs (character)
#   Column 2: log2FoldChange (numeric)
# No headers are required; the function uses the first two columns as-is.

library(dplyr)
library(tibble)

pathway_activation_analysis <- function(DE_data, pathway_list, pathway_names, min_genes = 4) {
  
  # ---- Ensure DE_data has at least two columns ----
  if (ncol(DE_data) < 2) {
    stop("DE_data must contain at least two columns: IDs and log2FoldChange.")
  }
  
  # ---- Extract first two columns (no column names assumed) ----
  ids <- DE_data[, 1]
  logfc <- as.numeric(DE_data[, 2])
  
  # ---- Clean data ----
  valid <- !is.na(ids) & ids != "" & !is.na(logfc)
  ids <- ids[valid]
  logfc <- logfc[valid]
  
  if (length(ids) == 0) {
    stop("No valid data after cleaning.")
  }
  
  # ---- Pre‑calculate named vector for fast lookup ----
  log2FoldChange_vector <- setNames(logfc, ids)
  
  # ---- Initialize results containers ----
  pathway_stats <- data.frame(
    Pathway_ID = character(),
    Pathway_Name = character(),
    Total_Genes = integer(),
    Genes_UP = integer(),
    Genes_DOWN = integer(),
    Perc_Genes_UP = numeric(),
    MedianlogFC = numeric(),
    P_value = numeric(),
    UniProtIDs = character(),
    FDR = numeric(),
    stringsAsFactors = FALSE
  )
  
  UniProtID_regulation <- data.frame(
    Pathway_ID = character(),
    UniProtID = character(),
    log2FoldChange = numeric(),
    stringsAsFactors = FALSE
  )
  
  # ---- Process each pathway ----
  pathway_ids <- names(pathway_list)
  
  for (i in seq_along(pathway_list)) {
    pathway_id <- pathway_ids[i]
    pathway_UniProtIDs <- pathway_list[[pathway_id]]
    pathway_name <- pathway_names[[pathway_id]]
    
    # Find overlapping UniProtIDs (IDs that are also in the DE data)
    overlap_UniProtIDs <- intersect(pathway_UniProtIDs, ids)
    n <- length(overlap_UniProtIDs)
    
    if (n >= min_genes) {
      # Count up‑regulated UniProtIDs (log2FoldChange > 0)
      s <- length(which(log2FoldChange_vector[overlap_UniProtIDs] > 0))
      
      # Two‑sided binomial test
      binom_res <- binom.test(s, n, p = 0.5, alternative = "two.sided")
      
      # UniProtID‑level details
      pathway_data <- data.frame(
        Pathway_ID = pathway_id,
        UniProtID = overlap_UniProtIDs,
        log2FoldChange = log2FoldChange_vector[overlap_UniProtIDs],
        stringsAsFactors = FALSE
      )
      
      # Append to pathway‑level statistics
      pathway_stats <- rbind(pathway_stats, data.frame(
        Pathway_ID = pathway_id,
        Pathway_Name = pathway_name,
        Total_Genes = n,
        Genes_UP = s,
        Genes_DOWN = n - s,
        Perc_Genes_UP = round(s / n * 100, 1),
        MedianlogFC = round(median(log2FoldChange_vector[overlap_UniProtIDs], na.rm = TRUE), 2),
        P_value = binom_res$p.value,
        UniProtIDs = paste(overlap_UniProtIDs, collapse = ","),
        stringsAsFactors = FALSE
      ))
      
      # Append to UniProtID‑regulation data
      UniProtID_regulation <- rbind(UniProtID_regulation, pathway_data)
    }
  }
  
  # ---- Multiple testing correction (FDR) ----
  if (nrow(pathway_stats) > 0) {
    pathway_stats$FDR <- p.adjust(pathway_stats$P_value, method = "fdr")
  } else {
    pathway_stats$FDR <- numeric(0)
  }
  
  # ---- Return results ----
  return(list(
    pathway_stats = pathway_stats,
    UniProtID_regulation = UniProtID_regulation
  ))
}
