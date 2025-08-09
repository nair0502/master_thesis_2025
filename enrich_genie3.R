#!/usr/bin/Rscript

# Load required libraries
library(clusterProfiler)
library(enrichplot)
library(UpSetR)
library(ggplot2)
library(tidyverse)
library(readr)

# --- Fixed paths ---
term2gene_file <- "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/thesis/vaishhnavi_master-thesis/module_files/solanum_proteins.results.txt"
wgcna_module_dir <- "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/thesis/vaishhnavi_master-thesis/module_files"
louvain_dir <- "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/thesis/vaishhnavi_master-thesis/genie3_res"
output_base_dir <- "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/thesis/vaishhnavi_master-thesis/enrichment_louvain"

# --- Function to clean gene IDs ---
clean_gene_ids <- function(gene_vec) {
  gene_vec %>%
    tolower() %>%
    sub("\\..*", "", .) %>%
    unique()
}

# --- Load Mercator annotation (term2gene) ---
solanum_mercator_df <- read_tsv(term2gene_file, col_names = TRUE)
solanum_term2gene <- solanum_mercator_df %>%
  filter(!is.na(IDENTIFIER), !is.na(NAME)) %>%
  mutate(
    gene = tolower(gsub("^'+|'+$", "", IDENTIFIER)),
    gene = sub("\\..*", "", gene)
  ) %>%
  select(term = NAME, gene) %>%
  distinct()

# --- List Louvain cluster gene files ---
louvain_files <- list.files(louvain_dir, pattern = "^genes_in_louvain_cluster_.*\\.txt$", full.names = TRUE)

for (louvain_file in louvain_files) {
  # Extract cluster/module name from Louvain file name
  cluster_name <- gsub("^genes_in_louvain_cluster_|\\.txt$", "", basename(louvain_file))
  message("Processing Louvain Cluster: ", cluster_name)

  # Attempt to find corresponding WGCNA module file
  wgcna_file_path <- file.path(wgcna_module_dir, paste0(cluster_name, ".txt"))
  if (!file.exists(wgcna_file_path)) {
    warning("No matching WGCNA module file for Louvain cluster '", cluster_name, "'. Skipping.")
    next
  }

  # Create output directory for this cluster
  output_dir <- file.path(output_base_dir, paste0("wgcna_module_", cluster_name))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Read and clean gene lists
  louvain_genes <- clean_gene_ids(read_tsv(louvain_file, col_names = FALSE) %>% pull(1))
  wgcna_genes <- clean_gene_ids(read_tsv(wgcna_file_path, col_names = FALSE) %>% pull(1))

  # Find common genes between Louvain and WGCNA module
  common_genes <- intersect(louvain_genes, wgcna_genes)
  message("Number of common genes (Louvain ∩ WGCNA): ", length(common_genes))

  # Write a small summary log per cluster
  summary_log_path <- file.path(output_dir, "summary.txt")
  write_lines(
    paste("Cluster:", cluster_name,
          "\nLouvain genes:", length(louvain_genes),
          "\nWGCNA genes:", length(wgcna_genes),
          "\nCommon genes:", length(common_genes),
          "\nDate:", Sys.time()),
    summary_log_path
  )

  if (length(common_genes) == 0) {
    warning("No common genes found for cluster ", cluster_name, ". Skipping enrichment.")
    next
  }

  # Define background universe as genes in WGCNA module that are annotated
  background_genes <- intersect(wgcna_genes, solanum_term2gene$gene)

  # Filter TERM2GENE for common genes only
  term2gene_filtered <- solanum_term2gene %>% filter(gene %in% common_genes)

  # Summarize number of genes per term
  common_summary <- term2gene_filtered %>%
    group_by(term) %>%
    summarise(count = n(), Genes = paste(unique(gene), collapse = ", "))

  # Save intermediate results
  write_lines(common_genes, file.path(output_dir, "common_genes.txt"))
  write_csv(term2gene_filtered, file.path(output_dir, "matching_TERM2GENE.csv"))
  write_csv(common_summary, file.path(output_dir, "common_gene_summary.csv"))

  # ----------------------------
  # Run enrichment using enricher
  # ----------------------------
  solanum_enrich_results <- enricher(
    gene = common_genes,
    TERM2GENE = solanum_term2gene,
    universe = background_genes,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 3,
    maxGSSize = 1500,
    qvalueCutoff = 1
  )

  if (!is.null(solanum_enrich_results) && nrow(as.data.frame(solanum_enrich_results)) > 0) {
    solanum_enrich_df <- as.data.frame(solanum_enrich_results)
    solanum_enrich_df$Description <- gsub("^'+|'+$", "", solanum_enrich_df$Description)
    solanum_enrich_df$ID <- gsub("^'+|'+$", "", solanum_enrich_df$ID)
    solanum_enrich_results@result <- solanum_enrich_df

    # Save enrichment results
    write_csv(solanum_enrich_df, file.path(output_dir, "solanum_enrichment_results.csv"))
    top_enrichment_df <- solanum_enrich_df %>% arrange(p.adjust) %>% slice_head(n = 10)
    write_csv(top_enrichment_df, file.path(output_dir, "top_enrichment_results.csv"))

    num_to_plot <- min(10, nrow(solanum_enrich_df))

    # Dotplot
    plot_file_base <- file.path(output_dir, "dotplot_enrichment")
    pdf(paste0(plot_file_base, ".pdf"), width = 12, height = 12)
    print(dotplot(solanum_enrich_results, showCategory = num_to_plot, font.size = 12, label_format = 50))
    dev.off()
    png(paste0(plot_file_base, ".png"), width = 1100, height = 1000)
    print(dotplot(solanum_enrich_results, showCategory = num_to_plot, font.size = 12, label_format = 50))
    dev.off()

    # Barplot
    plot_file_base <- file.path(output_dir, "barplot_enrichment")
    pdf(paste0(plot_file_base, ".pdf"), width = 12, height = 16)
    print(barplot(solanum_enrich_results, showCategory = num_to_plot, font.size = 12))
    dev.off()
    png(paste0(plot_file_base, ".png"), width = 1100, height = 1200)
    print(barplot(solanum_enrich_results, showCategory = num_to_plot, font.size = 12))
    dev.off()
  } else {
    message("No significant enrichment results found for Louvain-WGCNA Module ", cluster_name)
  }
}
