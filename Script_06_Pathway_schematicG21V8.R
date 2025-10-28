################## Pathway Schematics for G21 v 8 ####################


#===========================================================================
#Script by: Gillian Myers
#date: 10/14/2025
#Goal: Get png for each significant pathways, ensure no corruption

#NOTE:
#The use log 2 fold change instead of fold enrichment values

#Issue: graphs are identical

#Potential Solution: gene_matrix_PWC in script 01 has the gene names (no prefix),
#                     May need to use that instead of value listed below
#                     May also need to generate new script for this we r so cooked
#===========================================================================
rm(list = ls())
#Install necessary Packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview", force = TRUE) 
BiocManager::install("org.Sc.sgd.db", force = TRUE) 
BiocManager::install("KEGGREST", force = TRUE)
BiocManager::install("AnnotationHub")

#Load in packages
library(pathview)
library(KEGGREST)
library(org.Sc.sgd.db)
library(AnnotationDbi)
library(AnnotationHub)
library(stringr)

#set wd
setwd("C:/Users/gmyer/Fall2025 UR/Pathview analysis 8V21/")

#Make sure its correct
head(as.data.frame(kegg_enrich_G21V8))

#make pathway.id, [1] at the end is to only return the first ID
pathway.id_21 <- as.data.frame(kegg_enrich_G21V8)$ID
head(pathway.id_21)

#make gene.id, currently ORF, needs to be converted to KEGG gene ID's
kegg_enrich_G21V8_df <- as.data.frame(kegg_enrich_G21V8)
colnames(kegg_enrich_G21V8_df)

#from script kegg enrichment 21 v 8:
#       gene_matrix_PWC <- results_GOB21_GOB8_sig$log2FoldChange
# Get the original geneID and Fold change from your enrichment results
ORF_ID_21 <- as.character(kegg_enrich_G21V8_df$geneID)
lgfc_21 <- gene_matrix_PWC

# Split any ORF lists like "YAL001C/YBR160W"
sep_ORF_21 <- strsplit(ORF_ID_21, "/")
all_ORF_21 <- unlist(sep_ORF_21)

# Remove any whitespace, just in case
all_ORF_21 <- trimws(all_ORF_21)

# Make gene.data
gene.data_21 <- lgfc_21
names(gene.data_21) <- paste0("sce:", all_ORF_21)

#check gene.data
head(gene.data_21)
str(gene.data_21)


pathway.id_21 <- as.data.frame(kegg_enrich_G21V8)$ID
if (!grepl("^sce", pathway.id_21[1])) {
  pathway.id_21 <- paste0("sce", pathway.id_21)
}


#======Realized KeGG does not like the sce: prefix I just added to now removing it=======

# Remove the 'sce:' prefix from gene names before mapping
plain_ORFs_21 <- str_remove(names(gene.data_21), "^sce:")
head(plain_ORFs_21)

mapped_21 <- AnnotationDbi::select(org.Sc.sgd.db,
                                keys = plain_ORFs_21,
                                columns = "ENTREZID",
                                keytype = "ORF")  # or "SYSTEMATIC"

mapped_21 <- na.omit(mapped_21)
mapped_21 <- mapped_21[!duplicated(mapped_21$ORF), ]

# Rebuild gene.data vector with Entrez IDs as names
gene.data.mapped_21 <- gene.data_21
names(gene.data.mapped_21) <- plain_ORFs_21  # temporarily remove sce: for subsetting
gene.data.mapped_21 <- gene.data.mapped_21[mapped_21$ORF]
names(gene.data.mapped_21) <- mapped_21$ENTREZID

head(gene.data.mapped_21)


#Pathview to loop through alll sig genes to generate correct graphs
for (i in seq_along(pathway.id_21)) {
  cat("Plotting pathway:", pathway.id_21[i], "\n")
  
  tryCatch({
    pathview(
      gene.data = gene.data.mapped_21,
      pathway.id = pathway.id_21[i],
      species = "sce",
      out.suffix = paste0("G21V8_", i),
      kegg.native = TRUE
    )
  }, error = function(e) {
    cat("  → Failed for:", pathway.id_21[i], "\n  Reason:", e$message, "\n")
  })
}






# Define only the desired pathways, untested wait to run
selected_pathways <- c("01200")

# Loop only over selected pathways
for (pathway_id in selected_pathways) {
  cat("Plotting pathway:", pathway_id, "\n")
  
  tryCatch({
    pathview(
      gene.data = gene.data.mapped_21,
      pathway.id = pathway_id,
      species = "sce",
      out.suffix = paste0("G21V8_", pathway_id),
      kegg.native = TRUE
    )
  }, error = function(e) {
    cat("  → Failed for:", pathway_id, "\n  Reason:", e$message, "\n")
  })
}

# 2. Remove any 'sce:' prefix if previously added
gene_names <- names(gene_matrix_PWC)

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
selected_pathways <- c("01200")

# 8. Loop and plot with the correct data
for (pathway_id in selected_pathways) {
  cat("Plotting pathway:", pathway_id, "\n")
  
  tryCatch({
    pathview(
      gene.data = gene.data.final,
      pathway.id = pathway_id,
      species = "sce",
      out.suffix = paste0("G21V8_", pathway_id),
      kegg.native = TRUE
    )
  }, error = function(e) {
    cat("  → Failed for:", pathway_id, "\n  Reason:", e$message, "\n")
  })
}

