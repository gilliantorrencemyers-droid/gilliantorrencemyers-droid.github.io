---
title: "Pathview Schematic GOB21 vs GOB8 Enzyme Map"
permalink: /scripts-PS-21v8-EM/
categories: r-undergrad
---

Download it here!
- [Script_09_Pathway_G21V8_enzyme.R](/Script_09_Pathway_G21V8_enzyme.R)

```r
################## Pathway Schematics for G21 v 8 ####################

#===========================================================================
# Script by: Gillian Myers
# Date: 10/14/2025
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
setwd("C:/Users/gmyer/Fall2025 UR/Pathview analysis 8V21/")

# Ensure your enrichment results object exists
# (must come from your previous enrichment step)
head(as.data.frame(kegg_enrich_G21V8))

# Convert to a dataframe so we can extract columns
kegg_enrich_G21V8_df <- as.data.frame(kegg_enrich_G21V8)

# Extract KEGG pathway IDs
pathway.id_21 <- kegg_enrich_G21V8_df$ID
head(pathway.id_21)

# Extract ORF IDs (gene identifiers) and log2 fold change values
ORF_ID_21 <- as.character(kegg_enrich_G21V8_df$geneID)
lgfc_21   <- gene_matrix_PWC  # this should be your named vector of log2FC values

# Split any multi-gene entries like "YAL001C/YBR160W"
sep_ORF_21 <- strsplit(ORF_ID_21, "/")
all_ORF_21 <- unlist(sep_ORF_21)
all_ORF_21 <- trimws(all_ORF_21)

# Build the gene.data vector
gene.data_21 <- lgfc_21
names(gene.data_21) <- all_ORF_21

# Map ORFs → Entrez IDs
mapped_21 <- AnnotationDbi::select(
  org.Sc.sgd.db,
  keys     = all_ORF_21,
  columns  = "ENTREZID",
  keytype  = "ORF"
)

mapped_21 <- na.omit(mapped_21)
mapped_21 <- mapped_21[!duplicated(mapped_21$ORF), ]

# Build final named vector
gene.data.mapped_21 <- gene.data_21[mapped_21$ORF]
names(gene.data.mapped_21) <- mapped_21$ENTREZID

head(gene.data.mapped_21)

# Choose the KEGG pathways to plot (e.g. glycolysis and TCA)
selected_pathways <- c("00010", "00020")

# Generate pathway plots (non-native, gene-label view)
for (pathway_id in selected_pathways) {
  cat("Plotting pathway:", pathway_id, "\n")
  
  tryCatch({
    pathview(
      gene.data   = gene.data.mapped_21,
      pathway.id  = pathway_id,
      species     = "sce",
      out.suffix  = paste0("G21V8_GENE_Revised", pathway_id),
      kegg.native = FALSE   # <<< gene-level schematic
    )
  }, error = function(e) {
    cat("  → Failed for:", pathway_id, "\n  Reason:", e$message, "\n")
  })
}

```
