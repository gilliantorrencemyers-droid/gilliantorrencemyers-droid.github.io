---
title: "Ancestors log 2 fold change"
permalink: /scripts-ancestors-lg2fc/
categories: r-undergrad
---

Download it here!

- [Script_01_ancestors_log2foldchangegenerator_August2025.R](/Script_01_ancestors_log2foldchangegenerator_August2025.R)


```r
# DESeq2 Differential Expression Analysis for Yeast Strains
# Comparing GOB21 and GOB8_Gly against GOB8 (control)

# Install and load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

library(DESeq2)
library(dplyr)
library(tibble)

# Set working directory (adjust as needed)
setwd("C:/Users/gmyer/Fall2025 UR/")

# Read the count data
count_data <- read.csv("round_01_tximport_aggregated_counts.csv", 
                       row.names = 1, 
                       header = TRUE)

# Display basic info about the data
cat("Dimensions of count data:", dim(count_data), "\n")
cat("First few gene IDs:\n")
print(head(rownames(count_data)))

# Select only the samples we need for analysis
# GOB21: 4 replicates, GOB8_Gly: 3 replicates, GOB8: 4 replicates (control)
samples_to_keep <- c("GOB21_rep1", "GOB21_rep2", "GOB21_rep3", "GOB21_rep4",
                     "GOB8_Gly_rep1", "GOB8_Gly_rep2", "GOB8_Gly_rep3",
                     "GOB8_rep1", "GOB8_rep2", "GOB8_rep3", "GOB8_rep4")

# Subset the count data to include only our samples of interest
count_matrix <- count_data[, samples_to_keep]

# Convert to integer matrix (DESeq2 requires integer counts)
count_matrix <- round(as.matrix(count_matrix))

cat("Dimensions after filtering:", dim(count_matrix), "\n")
cat("Sample names:\n")
print(colnames(count_matrix))

# Create sample information dataframe
sample_info <- data.frame(
  sample = colnames(count_matrix),
  condition = c(rep("GOB21", 4),      # GOB21 replicates
                rep("GOB8_Gly", 3),   # GOB8_Gly replicates  
                rep("GOB8", 4)),      # GOB8 replicates (control)
  replicate = c(1:4, 1:3, 1:4)
)

# Set GOB8 as the reference level (control)
sample_info$condition <- factor(sample_info$condition, levels = c("GOB8", "GOB21", "GOB8_Gly"))

# Display sample information
cat("\nSample information:\n")
print(sample_info)

# Verify that sample names match between count matrix and sample info
stopifnot(colnames(count_matrix) == sample_info$sample)

# Remove genes with very low counts across all samples
# Keep genes with at least 10 counts in at least 3 samples (minimum replicate number)
keep_genes <- rowSums(count_matrix >= 10) >= 3
count_matrix_filtered <- count_matrix[keep_genes, ]

cat("\nGenes before filtering:", nrow(count_matrix), "\n")
cat("Genes after filtering:", nrow(count_matrix_filtered), "\n")

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix_filtered,
                              colData = sample_info,
                              design = ~ condition)

cat("\nDESeq2 dataset created successfully!\n")
cat("Design formula: ~ condition\n")
cat("Reference level (control):", levels(dds$condition)[1], "\n")

# Display summary of the DESeq2 object
print(dds)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Check the results names to see available comparisons
resultsNames(dds)

# Extract results for specific comparisons with alpha=0.05
# GOB21 vs GOB8 (control)
res_GOB21_vs_GOB8 <- results(dds, contrast=c("condition","GOB21","GOB8"), alpha=0.05)

# GOB8_Gly vs GOB8 (control)
res_GOB8_Gly_vs_GOB8 <- results(dds, contrast=c("condition","GOB8_Gly","GOB8"), alpha=0.05)

# View summary of results
summary(res_GOB21_vs_GOB8)
summary(res_GOB8_Gly_vs_GOB8)

# Export results to CSV files
write.csv(as.data.frame(res_GOB21_vs_GOB8), 
          file="GOB21_vs_GOB8_results.csv", row.names=TRUE)

write.csv(as.data.frame(res_GOB8_Gly_vs_GOB8), 
          file="GOB8_Gly_vs_GOB8_results.csv", row.names=TRUE)
```
