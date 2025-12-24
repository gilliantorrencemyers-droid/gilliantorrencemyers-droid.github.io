---
title: "Evolved Schematics"
permalink: /scripts-evolved-schematics/
categories: r-undergrad
---

Download it here!
[Script_08_Evolved_Schematics.R](/Script_08_Evolved_Schematics.R)

```r
####################### Script 08##########################

#==========================================================
#Script by: Gillian Myers
#Date: 10/28/2025
#Objective: Generate Pathview Schematics for all evolved strains
#==========================================================

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview", force = TRUE) 
BiocManager::install("org.Sc.sgd.db", force = TRUE) 
BiocManager::install("KEGGREST", force = TRUE)
BiocManager::install("AnnotationHub")

#Load in pathview
library(pathview)
library(KEGGREST)
library(org.Sc.sgd.db)
library(AnnotationDbi)
library(AnnotationHub)
#set wd
setwd("C:/Users/gmyer/Fall2025 UR/Pathview Analysis Evolved/")


# ------------------------------
run_pathview_filtered <- function(res_df, comp_name, pathway_id) {
  cat("Processing", comp_name, "for pathway", pathway_id, "...\n")
  
  # Check required columns
  if (!all(c("log2FoldChange", "padj") %in% colnames(res_df))) {
    warning(paste("Missing log2FoldChange or padj columns in", comp_name))
    return(NULL)
  }
  
  # Filter by significance
  sig_df <- res_df[!is.na(res_df$padj) & res_df$padj <= 0.05, ]
  if (nrow(sig_df) == 0) {
    cat("  ⚠️ No significant genes (padj ≤ 0.05). Skipping.\n")
    return(NULL)
  }
  
  # Extract clean ORF names (systematic yeast IDs)
  gene_names <- gsub("_mRNA$", "", rownames(sig_df))
  gene_names <- gsub("^sce:", "", gene_names)  # remove any old prefixes just in case
  
  # Map ORF IDs → Entrez IDs using org.Sc.sgd.db
  mapped <- AnnotationDbi::select(
    org.Sc.sgd.db,
    keys     = gene_names,
    columns  = "ENTREZID",
    keytype  = "ORF"
  )
  
  # Merge LFCs with mapped Entrez IDs
  mapped <- mapped[!is.na(mapped$ENTREZID), ]
  if (nrow(mapped) == 0) {
    cat("  ⚠️ No valid Entrez ID mappings found. Skipping.\n")
    return(NULL)
  }
  
  sig_df$ORF <- gene_names
  merged <- merge(sig_df, mapped, by.x = "ORF", by.y = "ORF", all.x = FALSE)
  
  # Build numeric LFC vector named by Entrez ID
  lfc_vector <- suppressWarnings(as.numeric(merged$log2FoldChange))
  names(lfc_vector) <- as.character(merged$ENTREZID)
  
  # Filter invalid entries
  valid_idx <- !is.na(lfc_vector) & is.finite(lfc_vector)
  lfc_vector <- lfc_vector[valid_idx]
  
  cat("  Genes significant:", nrow(sig_df), "| mapped Entrez IDs:", length(lfc_vector), "\n")
  
  if (length(lfc_vector) == 0) {
    cat("  ⚠️ No valid Entrez IDs to map. Skipping Pathview.\n")
    return(NULL)
  }
  
  # Run Pathview safely
  tryCatch({
    pathview(
      gene.data   = lfc_vector,
      pathway.id  = pathway_id,
      species     = "sce",
      out.suffix  = paste0(comp_name, "_", pathway_id, "_sig"),
      kegg.native = TRUE,
      low         = list(gene = "green"),
      mid         = list(gene = "gray"),
      high        = list(gene = "red")
    )
    cat("  ✅ Pathview complete for", comp_name, "(", pathway_id, ")\n")
  }, error = function(e) {
    cat("  ❌ Pathview error:", e$message, "\n")
  })
}


# ------------------------------
# Load differential expression results
cat("Loading DESeq2 results...\n")

res_PA5_t1000_vs_GOB8    <- read.csv("PA5_t1000_vs_GOB8_results.csv", row.names = 1)
res_PA5_t1000_vs_GOB21   <- read.csv("PA5_t1000_vs_GOB21_results.csv", row.names = 1)
res_PM5_t1000_vs_GOB8    <- read.csv("PM5_t1000_vs_GOB8_results.csv", row.names = 1)
res_PO5_t1000_L_vs_GOB8  <- read.csv("PO5_t1000_L_vs_GOB8_results.csv", row.names = 1)
res_PO5_t1000_S_vs_GOB8  <- read.csv("PO5_t1000_S_vs_GOB8_results.csv", row.names = 1)

# Combine into a named list for looping
comparison_list <- list(
  "PA5_t1000_vs_GOB8"   = res_PA5_t1000_vs_GOB8,
  "PA5_t1000_vs_GOB21"  = res_PA5_t1000_vs_GOB21,
  "PM5_t1000_vs_GOB8"   = res_PM5_t1000_vs_GOB8,
  "PO5_t1000_L_vs_GOB8" = res_PO5_t1000_L_vs_GOB8,
  "PO5_t1000_S_vs_GOB8" = res_PO5_t1000_S_vs_GOB8
  )

cat("Loaded", length(comparison_list), "comparison result files.\n\n")

# ------------------------------
# KEGG pathways to plot
kegg_pathways <- c("00010", "00020")  # Glycolysis & TCA

# ------------------------------
# Run Pathview for each dataset and pathway
for (comp_name in names(comparison_list)) {
  res_df <- comparison_list[[comp_name]]
  
  for (path_id in kegg_pathways) {
    tryCatch({
      run_pathview_filtered(res_df, comp_name, path_id)
    }, error = function(e) {
      cat("❌ Error in", comp_name, "for pathway", path_id, ":", e$message, "\n")
    })
  }
}

cat("\n🎉 All Pathview analyses complete! Check your working directory for PNG outputs.\n")
```
