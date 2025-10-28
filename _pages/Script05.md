---
title: "Pathview Schematic GLY vs GOB8"
permalink: /scripts-PS-glyv8/
categories: r-undergrad
---

Download it here!

- [Script_05_Pathway_SchematicG8V8.R](/Script_05_Pathway_SchematicG8V8.R)


```r
# from this website:
##https://pathview.r-forge.r-project.org/
################## Pathway Schematics for G8 v 8 ####################


#===========================================================================
#Script by: Gillian Myers
#Date: 9/20/2025
#Goal: Get png for each significant pathways, ensure no corruption

#NOTE:
# Needs log 2 fold change from kegg enrichment not the fold enrichment, 
# have to go back into previous script and get the vector

#Issue:
#When using the log 2 fold change, the values for 8 v gly and 21 v gly are 
#extremely similar, so the images generated looks exactly the same

#Solution:
# Verify to make sure the correct data sets were used for gene.data and pathway.id
#===========================================================================
rm(list = ls())


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
setwd("C:/Users/gmyer/Fall2025 UR/Pathview analysis 8v8/")

head(as.data.frame(kegg_enrich_G8GLYV8))

#make pathway.id, [1] at the end is to only return the first ID
pathway.id <- as.data.frame(kegg_enrich_G8GLYV8)$ID
head(pathway.id)

#make gene.id, currently ORF, needs to be converted to KEGG gene ID's
kegg_enrich_G8GLYV8_df <- as.data.frame(kegg_enrich_G8GLYV8)
colnames(kegg_enrich_G8GLYV8_df)

# Get the original geneID and FoldEnrichment from your enrichment results
ORF_ID <- as.character(kegg_enrich_G8GLYV8_df$geneID)
lgfc_8 <- gene_matrix_PWC2

# Split any ORF lists like "YAL001C/YBR160W"
sep_ORF <- strsplit(ORF_ID, "/")
all_ORF <- unlist(sep_ORF)

# Remove any whitespace, just in case
all_ORF <- trimws(all_ORF)

# Make gene.data
gene.data <- lgfc_8
names(gene.data) <- paste0("sce:", all_ORF)

head(gene.data)
str(gene.data)

pathway.id <- as.data.frame(kegg_enrich_G8GLYV8)$ID
if (!grepl("^sce", pathway.id[1])) {
  pathway.id <- paste0("sce", pathway.id)
}


library(stringr)
# Remove the 'sce:' prefix from gene names before mapping
plain_ORFs <- str_remove(names(gene.data), "^sce:")

head(plain_ORFs)

library(org.Sc.sgd.db)
library(AnnotationDbi)

mapped <- AnnotationDbi::select(org.Sc.sgd.db,
                                keys = plain_ORFs,
                                columns = "ENTREZID",
                                keytype = "ORF")  # or "SYSTEMATIC"

mapped <- na.omit(mapped)
mapped <- mapped[!duplicated(mapped$ORF), ]

# Rebuild gene.data vector with Entrez IDs as names
gene.data.mapped <- gene.data
names(gene.data.mapped) <- plain_ORFs  # temporarily remove sce: for subsetting
gene.data.mapped <- gene.data.mapped[mapped$ORF]
names(gene.data.mapped) <- mapped$ENTREZID

head(gene.data.mapped)

for (i in seq_along(pathway.id)) {
  cat("Plotting pathway:", pathway.id[i], "\n")
  
  tryCatch({
    pathview(
      gene.data = gene.data.mapped,
      pathway.id = pathway.id[i],
      species = "sce",
      out.suffix = paste0("pathview_8_", i),
      kegg.native = TRUE
    )
  }, error = function(e) {
    cat("  → Failed for:", pathway.id[i], "\n  Reason:", e$message, "\n")
  })
}









#Othe method untested for this so wait to do this
selected_pathways <- c("01200")


# Loop only over selected pathways
for (pathway_id in selected_pathways) {
  cat("Plotting pathway:", pathway_id, "\n")
  
  tryCatch({
    pathview(
      gene.data = gene.data.mapped,
      pathway.id = pathway_id,
      species = "sce",
      out.suffix = paste0("G8GLYV8_", pathway_id),
      kegg.native = TRUE
    )
  }, error = function(e) {
    cat("  → Failed for:", pathway_id, "\n  Reason:", e$message, "\n")
  })
}

# 1. Start from your gene_matrix_PWC, which has ORFs as names and log2FC as values
# Example: names(gene_matrix_PWC) = ORFs like "YAL038W"

# 2. Remove any 'sce:' prefix if previously added
gene_names <- names(gene_matrix_PWC2)

# 3. Map ORFs to ENTREZID using org.Sc.sgd.db
mapped_ids <- AnnotationDbi::select(
  org.Sc.sgd.db,
  keys = gene_names,
  columns = "ENTREZID",
  keytype = "ORF"
)

# Remove NAs and duplicates
mapped_ids <- na.omit(mapped_ids)
mapped_ids <- mapped_ids[!duplicated(mapped_ids$ORF), ]

# 4. Subset gene_matrix_PWC to only mapped ORFs
gene_matrix_mapped <- gene_matrix_PWC2[mapped_ids$ORF]

# 5. Rename with ENTREZID
names(gene_matrix_mapped) <- mapped_ids$ENTREZID

# 6. This is your FINAL gene.data for pathview
gene.data.final <- gene_matrix_mapped

# 7. Now define only the desired pathways
selected_pathways <- c("00010", "00020")

# 8. Loop and plot with the correct data
for (pathway_id in selected_pathways) {
  cat("Plotting pathway:", pathway_id, "\n")
  
  tryCatch({
    pathview(
      gene.data = gene.data.final,
      pathway.id = pathway_id,
      species = "sce",
      out.suffix = paste0("GGLYV8_", pathway_id),
      kegg.native = TRUE
    )
  }, error = function(e) {
    cat("  → Failed for:", pathway_id, "\n  Reason:", e$message, "\n")
  })
}

```
