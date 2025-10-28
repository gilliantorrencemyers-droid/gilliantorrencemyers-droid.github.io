---
title: "Gly8 vs GOB8 KEGG & GO Enrichment"
permalink: /scripts-gly-8-enrichment/
---
Download it here!

- [GOB8_GLY_vs_GOB8_KEGG_Enrichment.R](/GOB8_GLY_vs_GOB8_KEGG_Enrichment.R)

```r
#############################################################################
#Gillians Additions for KEGG Enrichment Analysis
#8/26/2025

#----------------------------KEGG Enrichment Analysis---------------------
#Load in some packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", force = TRUE)           # Core differential expression analysis
BiocManager::install("org.Sc.sgd.db", force = TRUE)    # Saccharomyces cerevisiae annotation database
BiocManager::install("tibble", force = TRUE) 
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")                        # Collection of data manipulation packages
}
BiocManager::install("KEGGREST", force = TRUE)         # Interface to KEGG REST API
BiocManager::install("clusterProfiler", force=TRUE)    #KEGG



library(DESeq2)                # For differential expression analysis
library(org.Sc.sgd.db)         # S. cerevisiae gene annotation
library(tidyverse)             # For data manipulation and visualization
library(tibble)                # Gene id section
library(stringr)               # row names and stuff I think
library(KEGGREST)              #I think this is also for KEGG
library(clusterProfiler)       #For KEGG visualization


#set wd
setwd("C:/Users/gmyer/Fall2025 UR/")

#Create data frame with read inn csv
GOB8GLY_vs_GOB8_results <- read.csv("C:/Users/gmyer/Fall2025 UR/GOB8_GLY_vs_GOB8_results.csv", check.names = FALSE, row.names = 1)

#want to remove the _mrna at the end of the rows
gene_id_8_v_8 <- rownames_to_column(GOB8GLY_vs_GOB8_results, var = "Gene_ID") # Convert row names to column
gene_id_8_v_8$Gene_ID <- str_remove_all(gene_id_8_v_8$Gene_ID, "_mRNA")     # Remove "_mRNA" string
view(gene_id_8_v_8)

#Keep only that row
gene_id_8_v_8 <- gene_id_8_v_8[,1] 

#Update first row names without _mrna
GOB8GLY_vs_GOB8_results <- remove_rownames(GOB8GLY_vs_GOB8_results)
rownames(GOB8GLY_vs_GOB8_results) <- gene_id_8_v_8 

#Check again to make sure only the _mrna portion is deleted
head(GOB8GLY_vs_GOB8_results)

#Checking to see of the rownames (very far left) is labeled anything
colnames(GOB8GLY_vs_GOB8_results)

#naming and creating a separate column with all the ORFS so I can convert them
GOB8GLY_vs_GOB8_results$ORF <- rownames(GOB8GLY_vs_GOB8_results)

#Convert the ORF to ENTREZID using mapIds from sc database
GOB8GLY_vs_GOB8_results$ENTREZID <- mapIds(org.Sc.sgd.db,
                                         keys = GOB8GLY_vs_GOB8_results$ORF,
                                         column = "ENTREZID",
                                         keytype = "ORF",
                                         multiVals = "first")

# Subset significant genes and ensure they have Entrez IDs
results_GOB8GLY_GOB8_sig <- subset(GOB8GLY_vs_GOB8_results, padj < 0.05)
results_GOB8GLY_GOB8_sig_entrez <- subset(results_GOB8GLY_GOB8_sig, is.na(ENTREZID) == FALSE)

# Create a vector of log2 fold changes for pathway analysis
gene_matrix_PWC2 <- results_GOB8GLY_GOB8_sig$log2FoldChange
names(gene_matrix_PWC2) <- rownames(results_GOB8GLY_GOB8_sig)
head(gene_matrix_PWC2)  # Check format

# KEGG pathway enrichment analysis
kegg_enrich_G8GLYV8 <- enrichKEGG(gene = names(gene_matrix_PWC2),
                                organism = 'yeast',          # Specify yeast
                                pvalueCutoff = 1,         # p-value cutoff
                                qvalueCutoff = 1)         # q-value cutoff

# KEGG module enrichment (broader functional groupings)
kegg_M_enrich_G8GLYV8 <- enrichMKEGG(gene = names(gene_matrix_PWC2),
                                   organism = 'yeast',
                                   pvalueCutoff = 1,         # More lenient cutoff
                                   qvalueCutoff = 1)

# Visualize KEGG enrichment results
barplot(kegg_enrich_G8GLYV8, 
        drop = TRUE,              # Drop empty categories
        showCategory = 10,        # Show top 10 categories
        title = "KEGG Enrichment Pathways",
        font.size = 8)
dotplot(kegg_enrich_G8GLYV8)
barplot(kegg_M_enrich_G8GLYV8, 
        drop = TRUE, 
        showCategory = 10,        # Show top 50 categories
        title = "KEGG Module Enrichment Pathways",
        font.size = 7)
dotplot(kegg_M_enrich_G8GLYV8)

#---------------------------GO Analysis----------------------------
#keytypes(org.Sc.sgd.db)


#I think the symbols (Dungs original code) is just the ORF

# Prepare data for GO enrichment (use only genes with symbols)
GOB8GLY_vs_GOB8_results_go <- subset(results_GOB8GLY_GOB8_sig_entrez, is.na(ENTREZID) == FALSE)

# Create a vector of log2 fold changes with Entrez IDs as names
gene_matrix_go_2 <- GOB8GLY_vs_GOB8_results_go$log2FoldChange

# Original code by DUNG:
####names(gene_matrix_go) <- GOB21_vs_GOB8_results_go$listData[["entrez"]]

names(gene_matrix_go_2) <- GOB8GLY_vs_GOB8_results_go$ENTREZID


# Perform GO enrichment analysis for Molecular Function (MF) terms
go_enrich_G8GLY_G8 <- enrichGO(gene = names(gene_matrix_go_2),
                            OrgDb = 'org.Sc.sgd.db',     # Yeast annotation database
                            keyType = "ENTREZID",        # Input ID type
                            ont = "BP",                  # Ontology: "MF"-Molecular Function "BP"-Biological pathway "CC"-Cellular component "ALL" all three
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.10)
#readable = TRUE

# Visualize GO enrichment results
barplot(go_enrich_G8GLY_G8, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Pathways Enrichment\nG8GLY/G8",
        font.size = 8)
dotplot(go_enrich_G8GLY_G8,
        showCategory = 10,
        title = "GO Pathways Enrichment\nG8GLY/G8",)
barplot(go_enrich_G8GLY_G8, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Pathways Enrichment\nG8GLY/G8",
        font.size = 5)
```
