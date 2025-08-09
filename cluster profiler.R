#!/usr/bin/Rscript

# Load required libraries
library(clusterProfiler)
library(enrichplot)
library(UpSetR)
library(ggplot2)
library(tidyverse)
library(readr)
library(igraph)
library(ggraph)

input_file <- "vaishnavi/thesis/vaishhnavi_master-thesis/module_files/solanum_proteins.results.txt"
module <- "vaishnavi/thesis/vaishhnavi_master-thesis/module_files/turquoise.txt"
output_dir_enrichment <- "vaishnavi/thesis/vaishhnavi_master-thesis/enrichments_wgcna/turquoise"


# ----------------------------
# 1. Load and Process Mercator File
# ----------------------------

# Load Mercator file correctly (tab-delimited)
solanum_mercator_df <- read_tsv(input_file, col_names = TRUE)
print(colnames(solanum_mercator_df))  # Optional: Check column names

# Build TERM2GENE by cleaning gene IDs and keeping term-gene mappings
solanum_term2gene <- solanum_mercator_df %>%
  filter(!is.na(IDENTIFIER), !is.na(NAME)) %>%
  mutate(
    gene = tolower(gsub("^'+|'+$", "", IDENTIFIER)),  # Remove surrounding single quotes
    gene = sub("\\..*", "", gene)  # Strip version suffixes from gene IDs
  ) %>%
  select(term = NAME, gene) %>%
  distinct()

# ----------------------------
# 2. Load and Process WGCNA Modules Results File
# ----------------------------

# Read WGCNA gene list and clean gene IDs
module_degs <- read_tsv(module, col_names = FALSE) %>%
  pull(1) %>%
  tolower() %>%
  sub("\\..*", "", .) %>%
  unique()


# ----------------------------
# 3. Diagnostic: Check and Print the Overlap
# ----------------------------

# Find common genes between DE list and Mercator annotation
common_genes <- intersect(module_degs, solanum_term2gene$gene)
message("Number of common genes found: ", length(common_genes))
print(head(common_genes, 20))  # Print first 20 common genes

# Filter TERM2GENE for only those genes
solanum_matching_mappings <- solanum_term2gene %>% filter(gene %in% common_genes)
message("Mapping entries for common genes:")
print(solanum_matching_mappings)

# Summarize how many genes per term
common_summary <- solanum_matching_mappings %>%
  group_by(term) %>%
  summarise(count = n(), Genes = paste(unique(gene), collapse = ", "))
message("Summary of overlapping genes by term:")
print(common_summary)

# ----------------------------
# Save Outputs
# ----------------------------

# Save list of common genes
write_lines(common_genes, file.path(output_dir_enrichment, "common_genes.txt"))

# Save the TERM2GENE entries for matching genes
write_csv(solanum_matching_mappings, file.path(output_dir_enrichment, "matching_TERM2GENE.csv"))

# Save the summary table per term
write_csv(common_summary, file.path(output_dir_enrichment, "common_gene_summary_by_term.csv"))


library(enrichplot)
library(UpSetR)

# ----------------------------
# 4. Run Enrichment Analysis with clusterProfiler
# ----------------------------

solanum_enrich_results <- enricher(
  gene = module_degs,
  TERM2GENE = solanum_term2gene,
  pvalueCutoff = 0.1,
  pAdjustMethod = "BH",
  minGSSize = 3,
  maxGSSize = 1500,
  qvalueCutoff = 0.1
)

# ----------------------------
# 5. Save Enrichment Results + Plots if results exist
# ----------------------------

if (!is.null(solanum_enrich_results) && nrow(as.data.frame(solanum_enrich_results)) > 0) {
  solanum_enrich_df <- as.data.frame(solanum_enrich_results)
  solanum_enrich_df$Description <- gsub("^'+|'+$", "", solanum_enrich_df$Description)
  solanum_enrich_df$ID <- gsub("^'+|'+$", "", solanum_enrich_df$ID)
  solanum_enrich_results@result <- solanum_enrich_df

  write.csv(solanum_enrich_df, file.path(output_dir_enrichment, "solanum_enrichment_results.csv"), row.names = FALSE)
  top_enrichment_df <- solanum_enrich_df %>% arrange(p.adjust) %>% slice_head(n = 10)
  write_csv(top_enrichment_df, file.path(output_dir_enrichment, "top_enrichment_results.csv"))

  num_to_plot <- min(10, nrow(solanum_enrich_df))

  # Dotplot
  plot_file_base <- file.path(output_dir_enrichment, "dotplot_enrichment")
  pdf(paste0(plot_file_base, ".pdf"), width = 12, height = 12)
  print(dotplot(solanum_enrich_results, showCategory = num_to_plot, font.size = 12, label_format = 50))
  dev.off()
  png(paste0(plot_file_base, ".png"), width = 1100, height = 1000)
  print(dotplot(solanum_enrich_results, showCategory = num_to_plot, font.size = 12, label_format = 50))
  dev.off()

  # Barplot
  plot_file_base <- file.path(output_dir_enrichment, "barplot_enrichment")
  pdf(paste0(plot_file_base, ".pdf"), width = 12, height = 16)
  print(barplot(solanum_enrich_results, showCategory = num_to_plot, font.size = 12))
  dev.off()
  png(paste0(plot_file_base, ".png"), width = 1100, height = 1200)
  print(barplot(solanum_enrich_results, showCategory = num_to_plot, font.size = 12))
  dev.off()

} else {
  message("No significant enrichment results found.")
}

