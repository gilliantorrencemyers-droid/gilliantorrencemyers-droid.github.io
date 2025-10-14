---
layout: single
title: "R Scripts & Pictures"
permalink: /r-undergrad-scripts-2PWC/
---

# R Scripts & Pictures

## R Scripts
- [Data_Analysis_2PWC_F25.R](gilliantorrencemyers-droid.github.io/Data_Analysis_2PWC_F25.R)

```r
Data_Analysis_2PWC_F25.R
#============ Data Analysis for Ancestors =================
# Created by: Gillian Myers
# Date Started: 8/23/2025
#==========================================================

#------------------------Notes-----------------------------
# Comparing GOB21 and GOB8_Gly against GOB8 (control)
# Creating Heatmap
# GOB21_vs_GOB8_results.csv
# GOB8_Gly_vs_GOB8_results.csv
#----------------------------------------------------------



#THESE ARE OLD PACKAGES FOR BIOCONDUCTOR AND R, I HAD TO MOVE THEM TO MAKE ROOM FOR UPDATES ONES
#Updated packages are in normal place: C:\Users\lmyers\AppData\Local\Temp\RtmpMdLwXW\downloaded_packages
.libPaths("C:/Users/lmyers/Desktop/Gillian/Downloaded Packages for R")

file.edit("~/.Rprofile")


# Install if not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RColorBrewer")     # Color palettes for visualization

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")           # Pretty heatmaps
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")          # Collection of data manipulation packages
}

# Load packages
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

#set wd
setwd("C:/Users/lmyers/Desktop/Gillian/Spring 2025 - Fall 2025/")

#making data frame
#GOB21_vs_GOB8_results.df <-read.csv("C:/Users/lmyers/Desktop/Gillian/Spring 2025 - Fall 2025/GOB21_vs_GOB8_results.csv", check.names = FALSE)
#head(GOB21_vs_GOB8_results.df)


#-------------Using all columns------------
#Make dataframe from the csv
GOB21_vs_GOB8_results <-read.csv("C:/Users/lmyers/Desktop/Gillian/Spring 2025 - Fall 2025/GOB21_vs_GOB8_results.csv", check.names = FALSE, row.names = 1)


#want to remove the _mrna at the end of the rows
gene_id_test <- rownames_to_column(GOB21_vs_GOB8_results, var = "Gene_ID") # Convert row names to column
gene_id_test$Gene_ID <- str_remove_all(gene_id_test$Gene_ID, "_mRNA")     # Remove "_mRNA" string
view(gene_id_test)
#Keep only that row
gene_id_test <- gene_id_test[,1] 

#Update first row names without _mrna
GOB21_vs_GOB8_results <- remove_rownames(GOB21_vs_GOB8_results)
rownames(GOB21_vs_GOB8_results) <- gene_id_test 
#-----------------------------------------------------------------------------------


#-----------isolating specific columns-----------------

#Make dataframe from the csv
GOB21_vs_GOB8_iso <-read.csv("C:/Users/lmyers/Desktop/Gillian/Spring 2025 - Fall 2025/GOB21_vs_GOB8_results.csv", check.names = FALSE, row.names = 1)


#want to remove the _mrna at the end of the rows
#labeling column GENE_ID
gene_id_test2 <- rownames_to_column(GOB21_vs_GOB8_iso, var = "Gene_ID") # Convert row names to column
#Isolating 2 columns
gene_id_test2 <- gene_id_test2[, c("Gene_ID", "log2FoldChange")]
#Verify
view(gene_id_test2)
#Removing MRNA from GENE_ID columns
gene_id_test2$Gene_ID <- str_remove_all(gene_id_test2$Gene_ID, "_mRNA")     # Remove "_mRNA" string
view(gene_id_test2)
#Keep only that row
single_row <- gene_id_test2[,1] 

#Update first row names without _mrna
gene_id_test2 <- remove_rownames(gene_id_test2)
rownames(gene_id_test2) <- single_row 

#delete first column
FCD <- gene_id_test2[, -1]

#-----------------------------------------------------------------------


#====================== works!! ========================
library(tibble)      # for rownames_to_column()
library(stringr)     # for str_remove_all()
library(pheatmap)    # for heatmap

# Step 1: Read CSV and use row 1 as row names
GOB21_vs_GOB8_iso <- read.csv(
  "C:/Users/lmyers/Desktop/Gillian/Spring 2025 - Fall 2025/GOB21_vs_GOB8_results.csv",
  check.names = FALSE,
  row.names = 1
)

# Step 2: Move row names (gene IDs) into a column
gene_id_test2 <- rownames_to_column(GOB21_vs_GOB8_iso, var = "Gene_ID")

# Step 3: Select only Gene_ID and log2FoldChange
gene_id_test2 <- gene_id_test2[, c("Gene_ID", "log2FoldChange")]

# Step 4: Remove "_mRNA" from gene IDs
gene_id_test2$Gene_ID <- str_remove_all(gene_id_test2$Gene_ID, "_mRNA")

# Step 5: Set gene names as row names and drop Gene_ID column
rownames(gene_id_test2) <- gene_id_test2$Gene_ID
FCD <- gene_id_test2[, "log2FoldChange", drop = FALSE]  # Keep as data frame

# Step 6: Plot heatmap (disable column clustering since it's only 1 column)
pheatmap(as.matrix(FCD), cluster_cols = FALSE)

#----------- P-value isolation ---------------



```

Download the R scripts by clicking the links above.
