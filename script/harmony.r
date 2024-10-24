# Load required libraries
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(cowplot)

# # Set working directory
# setwd("c:/Users/ChangBi/Code/cellAnnotation/mclustCellAnnotation_bk/code2")

# Load count matrices
train_counts <- t(readRDS(train_path)[,-1])
test_counts <- t(readRDS(test_path)[,-1])

# Find common genes
common_genes <- intersect(rownames(train_counts), rownames(test_counts))

# Subset the count matrices to include only common genes
train_counts <- train_counts[common_genes, ]
test_counts <- test_counts[common_genes, ]

# Create Seurat objects
train_seurat <- CreateSeuratObject(counts = train_counts, project = "train")
test_seurat <- CreateSeuratObject(counts = test_counts, project = "test")

# Add source information to metadata
train_seurat$source <- "train"
test_seurat$source <- "test"

# Merge datasets
combined <- merge(train_seurat, test_seurat, add.cell.ids = c("train", "test"))

# Pre-processing
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)

npcs <- 50
combined <- RunPCA(combined, npcs = npcs)

# Run UMAP before Harmony
combined <- RunUMAP(combined, dims = 1:npcs, reduction = "pca")

# Plot before Harmony
before_plot <- DimPlot(combined, reduction = "umap", group.by = "source") +
  ggtitle("Before Harmony")

# Run Harmony
combined <- RunHarmony(combined, group.by.vars = "source")

# Run UMAP after Harmony
combined <- RunUMAP(combined, dims = 1:npcs, reduction = "harmony")

# Plot after Harmony
after_plot <- DimPlot(combined, reduction = "umap", group.by = "source") +
  ggtitle("After Harmony")

# Combine plots
combined_plot <- before_plot + after_plot
ggsave(file.path(output_path, "harmony_comparison_umap.png"), combined_plot, width = 12, height = 6)

# Extract harmonized data
harmonized_data <- Embeddings(object = combined, reduction = "harmony")

# Split the data into train and test parts
train_cells <- WhichCells(combined, expression = source == "train")
test_cells <- WhichCells(combined, expression = source == "test")

harmonized_data_train <- harmonized_data[train_cells, ]
harmonized_data_test <- harmonized_data[test_cells, ]

# Convert to data frames if they aren't already
harmonized_data_train <- as.data.frame(harmonized_data_train)
harmonized_data_test <- as.data.frame(harmonized_data_test)

# Return the harmonized data
list(train = harmonized_data_train, test = harmonized_data_test)

# # Save to CSV files
# write.csv(harmonized_data_train, file = "harmonized_features_train.csv")
# write.csv(harmonized_data_test, file = "harmonized_features_test.csv")