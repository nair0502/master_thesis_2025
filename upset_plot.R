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

# UpSet Plot
  if ("geneID" %in% colnames(solanum_enrich_df)) {
    solanum_enrich_df$GeneList <- strsplit(solanum_enrich_df$geneID, "/")
    solanum_enrich_df$ShortID <- paste0("T", seq_len(nrow(solanum_enrich_df)))
    solanum_named_lists <- setNames(solanum_enrich_df$GeneList, solanum_enrich_df$ShortID)

    write_csv(solanum_enrich_df[, c("ShortID", "ID")],
              file.path(output_dir_enrichment, "upset_term_legend.csv"))

    plot_file_base <- file.path(output_dir_enrichment, "upset_enrichment")
    pdf(paste0(plot_file_base, ".pdf"), width = 16, height = 10)
    print(upset(fromList(solanum_named_lists), order.by = "freq", text.scale = 1.5, nintersects = NA))
    dev.off()
    png(paste0(plot_file_base, ".png"), width = 1600, height = 1000)
    print(upset(fromList(solanum_named_lists), order.by = "freq", text.scale = 1.5, nintersects = NA))
    dev.off()
  }
