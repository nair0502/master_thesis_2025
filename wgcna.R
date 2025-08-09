#!/usr/bin/Rscript

# Load necessary libraries
library(readr)
library(WGCNA)
library(data.table)
library(pheatmap)
library(ggplot2)
library(BiocManager)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(ReactomePA)

output_dir <- "vaishnavi/thesis/vaishhnavi_master-thesis/soft_threshold_6"

# Set options
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Load the vst_data.csv file
file_path <- "vaishnavi/vaishnavi_masterthesis/vst_data.csv"
cat("Loading data from:", file_path, "\n")
wgcna_data <- read.csv(file = file_path, stringsAsFactors = FALSE, header = TRUE)

# Convert to data frame and set row names
wgcna_data <- as.data.frame(wgcna_data)
rownames(wgcna_data) <- wgcna_data[, 1]  # Use first column as row names
wgcna_data <- wgcna_data[, -1]           # Remove the first column

# Convert all remaining data to numeric
wgcna_data[] <- lapply(wgcna_data, function(x) as.numeric(as.character(x)))

# Transpose the data so that samples are rows and genes are columns
datExpr <- as.data.frame(t(wgcna_data))

# Convert the data into the multi-set format for WGCNA
multiExpr <- list()
multiExpr[[1]] <- list(data = datExpr)

# Check for good samples and genes
gsg <- goodSamplesGenesMS(multiExpr, verbose = 3)

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    cat("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", "), "\n")
  }
  if (sum(!gsg$goodSamples[[1]]) > 0) {
    cat("Removing samples:", paste(rownames(multiExpr[[1]]$data)[!gsg$goodSamples[[1]]], collapse = ", "), "\n")
  }
  multiExpr[[1]]$data <- multiExpr[[1]]$data[gsg$goodSamples[[1]], gsg$goodGenes]
}

# Update datExpr for single dataset use
datExpr <- multiExpr[[1]]$data

# Sample clustering
sampleTree <- hclust(dist(datExpr), method = "average")

# Save the sample dendrogram plot as a PDF
pdf(file.path(output_dir, "sample_clustering.pdf"))
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering", sub = "", xlab = "", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# Determine the soft-thresholding power
spt <- pickSoftThreshold(datExpr, networkType = "signed hybrid", verbose = 5)
cat("Selected soft thresholding power results:\n")
print(spt)

# Plot scale independence and mean connectivity, then save as PDFs
pdf(file.path(output_dir, "soft_thresholding.pdf"))
par(mfrow = c(1, 2))
plot(spt$fitIndices[, 1], spt$fitIndices[, 2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(spt$fitIndices[, 1], spt$fitIndices[, 2], col = "red")
abline(h = 0.85, col = "red")

plot(spt$fitIndices[, 1], spt$fitIndices[, 5],
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(spt$fitIndices[, 1], spt$fitIndices[, 5], labels = spt$fitIndices[, 1], col = "red")
dev.off()

# soft-thresholding power 
softPower <- 6

# Compute adjacency and TOM similarity
adjacency <- adjacency(datExpr, power = softPower, type = "signed hybrid")
TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1 - TOM

# Compute gene clustering on TOM-based dissimilarity
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")
pdf(file.path(output_dir, "dynamic_module_gene_tree.pdf"))
sizeGrWindow(12, 9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# Identify modules using dynamic tree cut
Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, 
                          deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
cat("Module counts:\n")
print(table(Modules))

ModuleColors <- labels2colors(Modules)
cat("Module colors:\n")
print(table(ModuleColors))

# Save module colors table
ModuleColors_table <- as.data.frame(table(ModuleColors))
head(ModuleColors_table)
write.csv(ModuleColors_table, file = "ModuleColors_table.csv", row.names = FALSE)

# Plot gene dendrogram with module colors
pdf(file.path(output_dir, "gene dendrogram and module colours.pdf"))
plotDendroAndColors(geneTree, ModuleColors, "Module",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
dev.off()

# Compute module eigengenes
MElist <- moduleEigengenes(datExpr, colors = ModuleColors)
MEs <- MElist$eigengenes
ME_table <- as.data.frame(MEs)
write.csv(ME_table, file = "ME_table.csv", row.names = FALSE)

# Calculate eigengene dissimilarity and cluster eigengenes
ME.dissimilarity <- 1 - cor(MElist$eigengenes, use = "complete")
METree <- hclust(as.dist(ME.dissimilarity), method = "average")
pdf(file.path(output_dir, "module_eigengene_clustering.pdf"))
par(mar = c(0, 4, 2, 0))
par(cex = 0.6)
plot(METree, main = "Clustering of module eigengenes")
dev.off()

# Module merging
merge <- mergeCloseModules(datExpr, ModuleColors, cutHeight = 0.40)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs
pdf(file.path(output_dir, "gene dendrogram and modules merged.pdf"))
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
dev.off()
write.csv(mergedMEs, file = file.path(output_dir, "ModuleColors_table.csv"), row.names = FALSE)
head(mergedMEs)

#----------------------- saving the info about the gene_IDs in every merged modules -----------------------#
# 1. Save full gene-to-module color mapping
names(mergedColors) <- colnames(datExpr)
moduleColorDF <- data.frame(ModuleColor = mergedColors)
rownames(moduleColorDF) <- names(mergedColors)
write.csv(moduleColorDF, file = file.path(output_dir, "Gene_to_ModuleColor.csv"))
cat("Saved full gene-to-module color mapping to: Gene_to_ModuleColor.csv\n")

#----------------------- Metadta introduced -----------------------#
# Module trait relationship analysis
metadata_path <- "vaishnavi/thesis/vaishhnavi_master-thesis/filteredmetadatamasters.csv"
cat("Loading metadata from:", metadata_path, "\n")
metadata <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)
head(metadata)

# Remove extra columns and set traits
allTraits <- metadata[, -c(2,3)]
Samples <- rownames(datExpr)
traitRows <- match(Samples, allTraits$run)
datTraits <- allTraits[traitRows, ]
rownames(datTraits) <- allTraits[traitRows, "run"]
head(datTraits)
#----------------------- Module–Trait Relationship with Intercation -----------------------#
datTraits$interaction <- factor(datTraits$interaction, levels = c("Biotrophic", "Symbiotic", "Wild type", "necrotrophic", "drought"))

# Create the interaction matrix (one-hot encoding)
datInteractionNumeric <- data.frame(model.matrix(~ interaction - 1, data = datTraits))
colnames(datInteractionNumeric) <- c("Biotrophic", "Symbiotic", "Wild.type", "necrotrophic", "drought")

# Compute correlation between merged module eigengenes and interaction types
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
module_trait_correlation <- cor(mergedMEs, datInteractionNumeric, use = "p")
module_trait_Pvalue <- corPvalueStudent(module_trait_correlation, nSamples)
textMatrix <- paste(signif(module_trait_correlation, 2), "\n(", signif(module_trait_Pvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(module_trait_correlation)

pdf(file.path(output_dir, "Module_Trait_Relationship_with_interaction.pdf"), 
    width = 18, height = 20, useDingbats = FALSE)
datInteractionNumeric <- data.frame(model.matrix(~ interaction - 1, data = datTraits))
MET <- orderMEs(cbind(mergedMEs, datInteractionNumeric))
par(mar = c(6, 10, 3, 1))
labeledHeatmap(Matrix = module_trait_correlation,
               xLabels = colnames(datInteractionNumeric),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1, 1),
               main = "Module-Trait Relationship with different types of interactions")
dev.off()

# Eigengene network heatmap for interaction traits
pdf(file.path(output_dir, "Eigengene_heatmap_interaction.pdf"), 
    width = 20, height = 15, useDingbats = FALSE)
datInteractionNumeric <- data.frame(model.matrix(~ interaction - 1, data = datTraits))
colnames(datInteractionNumeric) <- c("Biotrophic", "Symbiotic", "Wild.type", "necrotrophic", "drought")
MET <- orderMEs(cbind(mergedMEs, datInteractionNumeric))
par(cex = 0.8)
plotEigengeneNetworks(
  MET, 
  "", 
  marDendro = c(0, 6, 1, 2),  
  marHeatmap = c(8, 6, 1, 2),  
  cex.lab = 1.0,  
  xLabelsAngle = 45,  
  zlim = c(-1, 1)
)
dev.off()


#----------------------- Module–Trait Relationship with Pathogens -----------------------#

# Step 1: Clean the pathogen labels
datTraits$pathogen <- trimws(datTraits$pathogen)                         
datTraits <- datTraits[datTraits$pathogen != "", ]                       
datTraits$pathogen <- droplevels(factor(datTraits$pathogen))            

# Alignment
commonSamples <- intersect(rownames(mergedMEs), rownames(datTraits))
mergedMEs_pathogen <- mergedMEs[commonSamples, ]
datTraits_pathogen <- datTraits[commonSamples, ]

datPathogenNumeric <- data.frame(model.matrix(~ pathogen - 1, data = datTraits_pathogen))
rownames(datPathogenNumeric) <- rownames(datTraits_pathogen)

# remove "pathogen" prefix from X-axis labels
colnames(datPathogenNumeric) <- gsub("^pathogen", "", colnames(datPathogenNumeric))

# check
stopifnot(all(rownames(mergedMEs_pathogen) == rownames(datPathogenNumeric)))

# Compute correlation and p-values
nSamples <- nrow(mergedMEs_pathogen)
module_pathogen_correlation <- cor(mergedMEs_pathogen, datPathogenNumeric, use = "p")
module_pathogen_Pvalue <- corPvalueStudent(module_pathogen_correlation, nSamples)

#Create text matrix for heatmap
textMatrix <- paste(signif(module_pathogen_correlation, 2),
                    "\n(", signif(module_pathogen_Pvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(module_pathogen_correlation)

#Plot module–trait heatmap
pdf(file.path(output_dir, "Module_Trait_Relationship_with_Pathogen.pdf"),
    width = 18, height = 20, useDingbats = FALSE)
MET <- orderMEs(cbind(mergedMEs_pathogen, datPathogenNumeric))
par(mar = c(6, 10, 3, 1))
labeledHeatmap(Matrix = module_pathogen_correlation,
               xLabels = colnames(datPathogenNumeric),
               yLabels = names(mergedMEs_pathogen),
               ySymbols = names(mergedMEs_pathogen),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1, 1),
               main = "Module-Trait Relationship with Pathogen Types")
dev.off()

# Plot eigengene network heatmap
pdf(file.path(output_dir, "Eigengene_heatmap_with_Pathogen.pdf"),
    width = 20, height = 15, useDingbats = FALSE)
par(cex = 0.8)
plotEigengeneNetworks(
  MET,
  "",
  marDendro = c(0, 6, 1, 2),
  marHeatmap = c(8, 6, 1, 2),
  cex.lab = 1.0,
  xLabelsAngle = 45,
  zlim = c(-1, 1)
)
dev.off()



#--------------------- Module Membership vs. Gene Significance (Interaction) ---------------------#


datTraits$interaction <- factor(datTraits$interaction, levels = unique(datTraits$interaction))
interactionNumeric <- data.frame(model.matrix(~ interaction - 1, data = datTraits))

# Extract module names and calculate module membership (MM)
modNames <- substring(names(mergedMEs), 3)

# Align datExpr and interactionNumeric to only include common samples
commonSamples <- intersect(rownames(datExpr), rownames(interactionNumeric))
datExpr_aligned <- datExpr[commonSamples, ]
interactionNumeric_aligned <- interactionNumeric[commonSamples, ]

# Print dimensions
cat("Dimensions of datExpr_aligned: ", dim(datExpr_aligned), "\n")
cat("Dimensions of interactionNumeric_aligned: ", dim(interactionNumeric_aligned), "\n")

# Check for mismatch
if (nrow(datExpr_aligned) != nrow(interactionNumeric_aligned)) {
  stop("Error: The number of rows in datExpr_aligned and interactionNumeric_aligned do not match.")
}

# Calculate gene significance (GS) and p-values
geneTraitSignificance <- as.data.frame(cor(datExpr_aligned, interactionNumeric_aligned, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(datExpr_aligned)))

# Name the columns
names(geneTraitSignificance) <- paste("GS.", colnames(interactionNumeric_aligned), sep = "")
names(GSPvalue) <- paste("p.GS.", colnames(interactionNumeric_aligned), sep = "")

# Check result
print(head(geneTraitSignificance))
print(head(GSPvalue))

# Module membership correlation
geneModuleMembership <- as.data.frame(cor(datExpr_aligned, mergedMEs[commonSamples, ], use = "p"))
names(geneModuleMembership) <- paste("MM", modNames, sep = ".")

#------------------ Plotting Function ------------------#
plotModuleSignificance <- function(module, module_name) {
  column <- match(module, modNames)
  moduleGenes <- mergedColors == module
  interaction_types <- colnames(interactionNumeric_aligned)

  geneTraitSignificance[] <- lapply(geneTraitSignificance, function(x) as.numeric(as.character(x)))

  for (interaction_type in interaction_types) {
    column_name <- paste("GS.", interaction_type, sep = "")
    if (!(column_name %in% names(geneTraitSignificance))) {
      warning(paste("Warning: Column", column_name, "not found in geneTraitSignificance. Skipping..."))
      next
    }
    interaction_column <- match(column_name, names(geneTraitSignificance))
    if (sum(moduleGenes) == 0) {
      warning(paste("Warning: No genes found in module", module, ". Skipping..."))
      next
    }

    pdf(file.path(output_dir, paste0("gene_significance_for_", interaction_type, "_", module_name, ".pdf")))
    verboseScatterplot(
      abs(geneModuleMembership[moduleGenes, column]),
      abs(geneTraitSignificance[moduleGenes, interaction_column]),
      xlab = paste("Module Membership in", module_name, "module"),
      ylab = paste("Gene significance for", interaction_type),
      main = paste("Module Membership vs. Gene Significance\n(", interaction_type, ")"),
      cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_name
    )
    dev.off()
  }
}

#------------------ Dynamic Module Loop ------------------#

# Automatically detect and use unique module colors (excluding "grey" if you wish)
modules_to_plot <- unique(mergedColors)
modules_to_plot <- setdiff(unique(mergedColors), "grey")  # optional: exclude 'grey'

for (mod in modules_to_plot) {
  plotModuleSignificance(mod, mod)
}

# Save p-values for each module
for (module in modules_to_plot) {
  genes_in_module <- names(mergedColors)[mergedColors == module]
  if (length(genes_in_module) == 0) next
  pvals_in_module <- GSPvalue[genes_in_module, , drop = FALSE]
  output_file <- file.path(output_dir, paste0("PValues_", module, ".csv"))
  write.csv(pvals_in_module, file = output_file, row.names = TRUE)
  cat("Saved p-values for module:", module, "->", output_file, "\n")
}

