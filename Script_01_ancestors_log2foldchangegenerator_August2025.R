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

#############################################################################
#Gillians Additions for KEGG Enrichment Analysis
#8/24/2025

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
GOB21_vs_GOB8_results <- read.csv("C:/Users/gmyer/Fall2025 UR/GOB21_vs_GOB8_results.csv", check.names = FALSE, row.names = 1)

#want to remove the _mrna at the end of the rows
gene_id_test <- rownames_to_column(GOB21_vs_GOB8_results, var = "Gene_ID") # Convert row names to column
gene_id_test$Gene_ID <- str_remove_all(gene_id_test$Gene_ID, "_mRNA")     # Remove "_mRNA" string
view(gene_id_test)

#Keep only that row
gene_id_test <- gene_id_test[,1] 

#Update first row names without _mrna
GOB21_vs_GOB8_results <- remove_rownames(GOB21_vs_GOB8_results)
rownames(GOB21_vs_GOB8_results) <- gene_id_test 

#Check again to make sure only the _mrna portion is deleted
head(GOB21_vs_GOB8_results)

#Checking to see of the rownames (very far left) is labeled anything
colnames(GOB21_vs_GOB8_results)

#naming and creating a separate column with all the ORFS so I can convert them
GOB21_vs_GOB8_results$ORF <- rownames(GOB21_vs_GOB8_results)

#Convert the ORF to ENTREZID using mapIds from sc database
GOB21_vs_GOB8_results$ENTREZID <- mapIds(org.Sc.sgd.db,
                                          keys = GOB21_vs_GOB8_results$ORF,
                                          column = "ENTREZID",
                                          keytype = "ORF",
                                          multiVals = "first")

# Subset significant genes and ensure they have Entrez IDs
results_GOB21_GOB8_sig <- subset(GOB21_vs_GOB8_results, padj < 0.05)
results_GOB21_GOB8_sig_entrez <- subset(results_GOB21_GOB8_sig, is.na(ENTREZID) == FALSE)

# Create a vector of log2 fold changes for pathway analysis
gene_matrix_PWC <- results_GOB21_GOB8_sig$log2FoldChange
names(gene_matrix_PWC) <- rownames(results_GOB21_GOB8_sig)
head(gene_matrix_PWC)  # Check format

# KEGG pathway enrichment analysis
kegg_enrich_G21V8 <- enrichKEGG(gene = names(gene_matrix_PWC),
                          organism = 'yeast',          # Specify yeast
                          pvalueCutoff = 0.05,         # p-value cutoff
                          qvalueCutoff = 0.10)         # q-value cutoff

# KEGG module enrichment (broader functional groupings)
kegg_M_enrich_G21V8 <- enrichMKEGG(gene = names(gene_matrix_PWC),
                             organism = 'yeast',
                             pvalueCutoff = 1,         # More lenient cutoff
                             qvalueCutoff = 1)

# Visualize KEGG enrichment results
barplot(kegg_enrich_G21V8, 
        drop = TRUE,              # Drop empty categories
        showCategory = 10,        # Show top 10 categories
        title = "KEGG Enrichment Pathways",
        font.size = 8)
dotplot(kegg_enrich_G21V8)
barplot(kegg_M_enrich_G21V8, 
        drop = TRUE, 
        showCategory = 50,        # Show top 50 categories
        title = "KEGG Module Enrichment Pathways",
        font.size = 8)
dotplot(kegg_M_enrich_G21V8)

#---------------------------GO Analysis----------------------------
#keytypes(org.Sc.sgd.db)


#I think the symbols (Dungs original code) is just the ORF

# Prepare data for GO enrichment (use only genes with symbols)
GOB21_vs_GOB8_results_go <- subset(results_GOB21_GOB8_sig_entrez, is.na(ENTREZID) == FALSE)

# Create a vector of log2 fold changes with Entrez IDs as names
gene_matrix_go <- GOB21_vs_GOB8_results_go$log2FoldChange

# Original code by DUNG:
####names(gene_matrix_go) <- GOB21_vs_GOB8_results_go$listData[["entrez"]]

names(gene_matrix_go) <- GOB21_vs_GOB8_results_go$ENTREZID


# Perform GO enrichment analysis for Molecular Function (MF) terms
go_enrich_G21G8 <- enrichGO(gene = names(gene_matrix_go),
                      OrgDb = 'org.Sc.sgd.db',     # Yeast annotation database
                      keyType = "ENTREZID",        # Input ID type
                      ont = "BP",                  # Ontology: "MF"-Molecular Function "BP"-Biological pathway "CC"-Cellular component "ALL" all three
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
                      #readable = TRUE

# Visualize GO enrichment results
barplot(go_enrich_G21G8, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Pathways Enrichment\nG21/G8",
        font.size = 8)
dotplot(go_enrich_G21G8,
        showCategory = 10,
        title = "GO Pathways Enrichment\nG21/G8",)
barplot(go_enrich_G21G8, 
        drop = TRUE, 
        showCategory = 25, 
        title = "GO Pathways Enrichment\nG21/G8",
        font.size = 8)




####--------------------- SECTION 15B: GO ENRICHMENT ANALYSIS PA VS PO-----------------------####
# Gene Ontology enrichment analysis

# Prepare data for GO enrichment (use only genes with symbols)
results_sig_PA_PO_go <- subset(results_PA_PO_sig_entrez, is.na(symbol) == FALSE)

# Create a vector of log2 fold changes with Entrez IDs as names
gene_matrix_PA_PO_2 <- results_sig_PA_PO_go$log2FoldChange
names(gene_matrix_PA_PO_2) <- results_sig_PA_PO_go@listData[["entrez"]]

# Perform GO enrichment analysis for Molecular Function (MF) terms
go_enrich_PA_PO <- enrichGO(gene = names(gene_matrix_PA_PO_2),
                            OrgDb = 'org.Sc.sgd.db',     # Yeast annotation database
                            keyType = "ENTREZID",        # Input ID type
                            ont = "BP",                  # Ontology: Biological pathways
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.10)

# Visualize GO enrichment results
barplot(go_enrich_PA_PO, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Pathways Enrichment\nPA/PO ",
        font.size = 8)

barplot(go_enrich_PA_PO, 
        drop = TRUE, 
        showCategory = 25, 
        title = "GO Pathways Enrichment\nPA/PO",
        font.size = 8)
enrichplot::dotplot(go_enrich_PA_PO,
                    showCategory = 20,
                    title = "GO Pathways Enrichment\nPA/PO")