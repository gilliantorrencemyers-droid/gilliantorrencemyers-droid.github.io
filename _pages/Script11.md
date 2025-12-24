---
title: "Pathview Schematic GLY vs GOB8 enzymes"
permalink: /scripts-PS-Gv8-enzymes/
categories: r-undergrad
---

Download it here!

- [Script_11_Pathview_GLYv8_enzymes.R](/Script_11_Pathview_GLYv8_enzymes.R)


```r
################## Pathway Schematics for G8 v 8 ####################

#===========================================================================
# Script by: Gillian Myers
# Date: 9/20/2025
# Goal: Generate PNGs for each significant pathway with gene-level labels
#===========================================================================

rm(list = ls())

# Load required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("pathview", "org.Sc.sgd.db", "KEGGREST", "AnnotationHub"), force = TRUE)

library(pathview)
library(KEGGREST)
library(org.Sc.sgd.db)
library(AnnotationDbi)
library(stringr)

# Set working directory
setwd("C:/Users/gmyer/Fall2025 UR/Pathview analysis 8v8/")

# Ensure your enrichment results object exists
head(as.data.frame(kegg_enrich_G8GLYV8))

# Convert enrichment object to a dataframe
kegg_enrich_G8GLYV8_df <- as.data.frame(kegg_enrich_G8GLYV8)

# Extract KEGG pathway IDs
pathway.id_8 <- kegg_enrich_G8GLYV8_df$ID
head(pathway.id_8)

# Extract ORF IDs (gene identifiers) and log2 fold change values
ORF_ID_8 <- as.character(kegg_enrich_G8GLYV8_df$geneID)
lgfc_8   <- gene_matrix_PWC2  # your named vector of log2FC values

# Split any multi-gene entries like "YAL001C/YBR160W"
sep_ORF_8 <- strsplit(ORF_ID_8, "/")
all_ORF_8 <- unlist(sep_ORF_8)
all_ORF_8 <- trimws(all_ORF_8)

# Build the gene.data vector
gene.data_8 <- lgfc_8
names(gene.data_8) <- all_ORF_8

# Map ORFs → Entrez IDs
mapped_8 <- AnnotationDbi::select(
  org.Sc.sgd.db,
  keys     = all_ORF_8,
  columns  = "ENTREZID",
  keytype  = "ORF"
)

mapped_8 <- na.omit(mapped_8)
mapped_8 <- mapped_8[!duplicated(mapped_8$ORF), ]

# Build final named vector
gene.data.mapped_8 <- gene.data_8[mapped_8$ORF]
names(gene.data.mapped_8) <- mapped_8$ENTREZID

head(gene.data.mapped_8)

# Choose the KEGG pathways to plot (e.g., glycolysis and TCA)
selected_pathways <- c("00010", "00020")

# Generate pathway plots (non-native, gene-label view)
for (pathway_id in selected_pathways) {
  cat("Plotting pathway:", pathway_id, "\n")
  
  tryCatch({
    pathview(
      gene.data   = gene.data.mapped_8,
      pathway.id  = pathway_id,
      species     = "sce",
      out.suffix  = paste0("G8GLYV8_GENE_", pathway_id),
      kegg.native = FALSE   # <<< gene-level schematic
    )
  }, error = function(e) {
    cat("  → Failed for:", pathway_id, "\n  Reason:", e$message, "\n")
  })
}

```
