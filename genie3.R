#!/usr/bin/Rscript

# Load libraries
library(foreach)
library(rngtools)
library(doRNG)
library(GENIE3)
library(igraph)
library(tidyverse)
library(randomForest)
library(readr)
library(Matrix)
library(pheatmap)

# INPUT / OUTPUT FILE PATHS
vst_file <- "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/thesis/vaishhnavi_master-thesis/vst_data.csv"
TF_file <- "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/thesis/vaishhnavi_master-thesis/TF_list.txt"
output_dir <- "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/thesis/vaishhnavi_master-thesis/genie3_res"
wgcna_file <- "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/thesis/vaishhnavi_master-thesis/soft_threshold_6/Gene_to_ModuleColor.csv"
genie3_file <- "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/thesis/vaishhnavi_master-thesis/genie3_res/genie3_link_list.csv"
louvain_file <- "/dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/thesis/vaishhnavi_master-thesis/genie3_res/louvain_clusters_resolution_2_1.txt"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Load the VST data ---
vst_matrix <- read.csv(vst_file, row.names = 1, header = TRUE, check.names = FALSE)
rownames(vst_matrix) <- sub("\\.\\d+$", "", rownames(vst_matrix))  
vst_matrix <- as.matrix(vst_matrix)

# --- LOading the TF file --- 
tf_list <- read.table(TF_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tf_gene_ids <- tf_list[[1]]

# ---Check for any gene overlap between the VST data and the TF list I provided---
# So, all the TF gene_IDs should be there in VST. When printed all the TFs were there
intersecting_ids <- intersect(tf_gene_ids, rownames(vst_matrix))
cat("Intersecting TF IDs in VST matrix:", length(intersecting_ids), "\n")

# --- RUN GENIE3 and also creating the linklist ---
weightMatrix <- GENIE3(vst_matrix, regulators = tf_gene_ids, nCores = 6)
linkList <- getLinkList(weightMatrix)

# --- Save link list ---
write.csv(linkList, file.path(output_dir, "genie3_link_list.csv"), row.names = FALSE)

# --- since i have the link list i am just reading it and loading it for the plots ---
linkList <- read_csv(genie3_file)
g <- graph_from_data_frame(linkList, directed = TRUE)

# Plot network
#--- didnt carry out---
pdf(file = file.path(output_dir, "genie3_network_plot.pdf"), width = 10, height = 10)
plot(g, vertex.label = NA, vertex.size = 5, edge.arrow.size = 0.5, layout = layout_with_fr(g))
dev.off()

# Degree distribution
#--- didnt carry out---
deg <- degree(g, mode = "out")
pdf(file = file.path(output_dir, "degree_distribution.pdf"), width = 8, height = 6)
hist(deg, breaks = 50, main = "Degree Distribution (Out-degree)", xlab = "Number of Targets",
     ylab = "Number of Genes", col = "lightblue", border = "black")
dev.off()

# Edge weight distribution
#--- didnt carry out---
pdf(file = file.path(output_dir, "edge_weight_distribution.pdf"), width = 8, height = 6)
hist(linkList$weight, breaks = 50, main = "Edge Weight Distribution", xlab = "Weight",
     ylab = "Frequency", col = "lightgreen", border = "black")
abline(v = median(linkList$weight), col = "red", lty = 2)
dev.off()

# --- LOUVAIN CLUSTERING ---
linkList <- read_csv(genie3_file)
g_directed <- graph_from_data_frame(linkList, directed = TRUE)
#--- since igraph::cluster_louvain() need undirected input ---
# literature review supporting this is in my google doc
g_undirected <- as.undirected(g_directed, mode = "collapse")
# When the resolution is one the original modularity is attainted
clusters <- cluster_louvain(g_undirected)

#--- currently carrying out diff resolution to attain the same no of WGCNA modules ---
clusters <- cluster_louvain(g_undirected, resolution = 2.1)

cat("Louvain clustering produced with resolution 2.1 ", length(unique(clusters$membership)), "modules.\n")


# --- Save Louvain module assignments ---
# at the end i got 7 Louvain clusters ≠ WGCNA modules (39) --- when i used resolution 7 ---
names(clusters$membership) <- V(g_undirected)$name
louvain_membership <- data.frame(geneID = names(clusters$membership),moduleID = clusters$membership)
write.table(louvain_membership, file = file.path(output_dir, "louvain_clusters_resolution_2_1.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Module assignments saved to", file.path(output_dir, "louvain_clusters_resolution_2_1.txt"), "\n")


# --- COSINE SIMILARITY ---
#Doing it between Louvain and WGCNA modules
#Could do it between TF also

# --- Load Louvain assignments ---
louvain_df <- read.delim(file.path(output_dir, "louvain_clusters_resolution_2_1.txt"),
                         header = TRUE, stringsAsFactors = FALSE)

cat("Loaded Louvain assignments:", nrow(louvain_df), "genes\n")

# --- Load WGCNA file which has the gene_ids and the corresponding module colours ---
wgcna_raw_df <- read.csv(wgcna_file, header = TRUE, stringsAsFactors = FALSE)

cat("Loaded raw WGCNA assignments:", nrow(wgcna_raw_df), "genes\n")

# --- Clean WGCNA gene IDs (keep up to the second dot) ---
# --- got error that the louvain and the WGCNA dont match. So I am carrying out WGCNA_cleaned_df part--- 
clean_id <- function(id) sub("^(([^.]+\\.[^.]*)).*", "\\1", id)

wgcna_cleaned_df <- data.frame(
  geneID = clean_id(wgcna_raw_df$geneID),
  moduleColor = wgcna_raw_df$ModuleColor,
  stringsAsFactors = FALSE
)

# --- Clean Louvain gene IDs in the same way ---
louvain_df$geneID <- clean_id(louvain_df$geneID)

# --- Print IDs ---
cat("Sample Louvain IDs after cleaning:\n")
print(head(louvain_df$geneID, 10))

cat("Sample WGCNA IDs after cleaning:\n")
print(head(wgcna_cleaned_df$geneID, 10))

# --- Find common genes using cleaned WGCNA IDs ---
common_genes <- intersect(louvain_df$geneID, wgcna_cleaned_df$geneID)
cat("Common genes between Louvain and WGCNA:", length(common_genes), "\n")

# --- Filter both dataframes to common genes ---
louvain_df <- louvain_df[louvain_df$geneID %in% common_genes, ]
wgcna_df <- wgcna_cleaned_df[wgcna_cleaned_df$geneID %in% common_genes, ]

cat("Filtered Louvain assignments:", nrow(louvain_df), "genes\n")
cat("Filtered WGCNA assignments:", nrow(wgcna_df), "genes\n")

# --- Create membership vectors ---
louvain_membership <- setNames(louvain_df$moduleID, louvain_df$geneID)
wgcna_membership <- setNames(wgcna_df$moduleColor, wgcna_df$geneID)

# --- List unique modules ---
louvain_clusters <- sort(unique(louvain_membership))
wgcna_modules <- sort(unique(wgcna_membership))

cat("Unique Louvain modules:", length(louvain_clusters), "\n")
cat("Unique WGCNA modules:", length(wgcna_modules), "\n")

#--FOR RESOLUTION 1---
#-------- Unique Louvain modules: 7 
#-------- Unique WGCNA modules: 39 

#--FOR RESOLUTION 2.1---
#-------- Unique Louvain modules: 43 
#-------- Unique WGCNA modules: 39 


# --- Combine module assignments into a dataframe ---
combined_df <- data.frame(
  geneID = louvain_df$geneID,
  LouvainModule = louvain_membership[louvain_df$geneID],
  WGCNAModule = wgcna_membership[louvain_df$geneID],
  stringsAsFactors = FALSE
)

# --- Save combined module assignments ---
combined_file <- file.path(output_dir, "louvain_res_2_1_wgcna_combined_modules.txt")
write.table(combined_df, combined_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Combined module assignments saved to:", combined_file, "\n")

# --- Create contingency table (module overlap matrix) ---
module_table <- table(combined_df$LouvainModule, combined_df$WGCNAModule)

# --- Save module overlap table ---
overlap_file <- file.path(output_dir, "louvain_res_2_1_wgcna_module_overlap_table.txt")
write.table(module_table, overlap_file, sep = "\t", quote = FALSE)

cat("Module overlap table saved to:", overlap_file, "\n")

# --- binary matrix (presence/absence) ---
# For each gene and each module, a 1 indicates the gene is in that module, and 0 means it is not.
binary_matrix <- (module_table > 0) * 1
binary_file <- file.path(output_dir, "louvain_res2_1_wgcna_binary_matrix.txt")
write.table(binary_matrix, binary_file, sep = "\t", quote = FALSE)

cat("Binary module overlap matrix saved to:", binary_file, "\n")


# Build binary membership matrices
gene_by_louvain <- sapply(louvain_clusters, function(cl) as.integer(louvain_membership == cl))
rownames(gene_by_louvain) <- common_genes

gene_by_wgcna <- sapply(wgcna_modules, function(mod) as.integer(wgcna_membership == mod))
rownames(gene_by_wgcna) <- common_genes

# --- Writing the Cosine similarity function ---
# --- This formula measures the cosine of the angle between two vectors ---
# -------- 1 means they are identical in direction (perfect overlap). ---------
# -------- 0 means they are orthogonal (no shared genes). ---------
# -------- In between values reflect partial overlap. ---------
cosine_similarity <- function(x, y) {
  sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
}

# --- Compute similarity matrix ----
# A similarity matrix showing how much each WGCNA module overlaps with each Louvain module
# Each cell [i, j] represents the cosine similarity between WGCNA module i and Louvain module j
similarity_matrix <- matrix(NA, nrow = length(wgcna_modules), ncol = length(louvain_clusters),
                            dimnames = list(wgcna_modules, paste0("Louvain_", louvain_clusters)))

for (i in seq_along(wgcna_modules)) {
  for (j in seq_along(louvain_clusters)) {
    similarity_matrix[i, j] <- cosine_similarity(gene_by_wgcna[, i], gene_by_louvain[, j])
  }
}

# Plot heatmap
pdf(file = file.path(output_dir, "cosine_similarity_heatmap_lou_res2_1.pdf"), width = 10, height = 10)
pheatmap(similarity_matrix,
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Cosine Similarity: WGCNA vs Louvain (resolution 2.1)",
         color = colorRampPalette(c("white", "blue"))(100))
dev.off()

# Save similarity matrix
write.csv(similarity_matrix, file = file.path(output_dir, "cosine_similarity_matrix_lou_res2_1.csv"), quote = FALSE)

cat("All outputs saved to:", output_dir, "\n")
 
# -- this heatmap has the clustering ---
pdf(file = file.path(output_dir, "cosine_similarity_heatmap_clustered_lou_res2_1.pdf"), width = 10, height = 10)
pheatmap(similarity_matrix,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Cosine Similarity: WGCNA vs Louvain (resolution 2.1)",
         color = colorRampPalette(c("white", "blue"))(100))
dev.off()



## --- Arranging the heatmap based on the size ---
# Count how many genes are in each module
louvain_sizes <- sort(table(louvain_membership), decreasing = TRUE)
wgcna_sizes   <- sort(table(wgcna_membership), decreasing = TRUE)
# These give you the module names sorted by size (most genes first)
louvain_order <- names(louvain_sizes)
wgcna_order   <- names(wgcna_sizes)
# Reorder rows and columns based on size
similarity_matrix_ordered <- similarity_matrix[wgcna_order, paste0("Louvain_", louvain_order)]

#--heatmap based on size --
pdf(file = file.path(output_dir, "cosine_similarity_heatmap_size_cluster_res_2_1.pdf"), width = 10, height = 10)
pheatmap(similarity_matrix_ordered,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Cosine Similarity: WGCNA vs Louvain (Ordered by Size, resolution = 2.1)",
         color = colorRampPalette(c("white", "blue"))(100))
dev.off()

## taking the gene_ids from each louvain cluster for one particular modules who meet the cosine similarity of 0.08 and above and write a txt file for those gene_ids
#cosine similarity threshold

threshold <- 0.03

# Loop over each WGCNA module (row in similarity_matrix)
for (module in rownames(similarity_matrix)) {
  
  # Find Louvain clusters where cosine similarity ≥ threshold
  clusters_above_threshold <- colnames(similarity_matrix)[similarity_matrix[module, ] >= threshold]
  
  # Initialize empty set to collect genes
  gene_set <- character(0)
  
  # For each relevant Louvain cluster
  for (cluster_name in clusters_above_threshold) {
    
    # Extract cluster ID from name (e.g., "Louvain_5" → "5")
    cluster_id <- sub("Louvain_", "", cluster_name)
    
    # Get gene IDs in this Louvain cluster
    genes_in_cluster <- names(louvain_membership)[louvain_membership == cluster_id]
    
    # Merge gene IDs without duplicates
    gene_set <- union(gene_set, genes_in_cluster)
  }
  
  # Define output file path for this module
  output_file <- file.path(output_dir, paste0("genes_in_louvain_cluster_", module, ".txt"))
  
  # Save gene list to file (only if non-empty)
  if (length(gene_set) > 0) {
    writeLines(gene_set, output_file)
  }
}

