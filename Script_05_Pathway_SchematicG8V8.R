# from this website:
##https://pathview.r-forge.r-project.org/

#Script by: Gillian Myers
#Date: 9/20/2025
#

#Download bioClite, does not exits bc too old I think, use BiocManager?
#source("http://bioconductor.org/biocLite.R")
#biocLite("pathview")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview", force = TRUE) 
BiocManager::install("org.Sc.sgd.db", force = TRUE) 
BiocManager::install("KEGGREST", force = TRUE)
BiocManager::install("AnnotationHub")
#Install pathview through Rforge
#install.packages("pathview",repos="http://R-Forge.R-project.org")

#Load in pathview
library(pathview)
library(KEGGREST)
library(org.Sc.sgd.db)
library(AnnotationDbi)
library(AnnotationHub)
#set wd
setwd("C:/Users/gmyer/Fall2025 UR/")

#check previously KEGG enrich data to see if I have gene.data dn pathway.id things


head(as.data.frame(kegg_enrich_G8GLYV8))

#make pathway.id, [1] at the end is to only return the first ID
pathway.id <- as.data.frame(kegg_enrich_G8GLYV8)$ID
head(pathway.id)

#make gene.id, currently ORF, needs to be converted to KEGG gene ID's
kegg_enrich_G8GLYV8_df <- as.data.frame(kegg_enrich_G8GLYV8)
colnames(kegg_enrich_G8GLYV8_df)

#Adding "sce:" prefic so it recognizes then as KEGG names

# Get the original geneID and FoldEnrichment from your enrichment results
ORF_ID <- as.character(kegg_enrich_G8GLYV8_df$geneID)
Fold_enrichment <- kegg_enrich_G8GLYV8_df$FoldEnrichment

# Split any ORF lists like "YAL001C/YBR160W"
sep_ORF <- strsplit(ORF_ID, "/")
FE_ID <- rep(Fold_enrichment, times = sapply(sep_ORF, length))
all_ORF <- unlist(sep_ORF)

# Remove any whitespace, just in case
all_ORF <- trimws(all_ORF)

# Make gene.data
gene.data <- FE_ID
names(gene.data) <- paste0("sce:", all_ORF)

head(gene.data)
str(gene.data)


#names(gene.data) <- all_ORF



pathway.id <- as.data.frame(kegg_enrich_G8GLYV8)$ID
if (!grepl("^sce", pathway.id[1])) {
  pathway.id <- paste0("sce", pathway.id)
}

#names(gene.data) <- sub("^sce:", "", names(gene.data))
#names(gene.data) <- kegg_enriched_ids


i <- 1  # or any index
pv.out <- pathview(gene.data = gene.data, pathway.id = pathway.id[i],
                   species = "sce",
                   out.suffix = paste0("sce_pathview_analysis_", i),
                   kegg.native = TRUE)

##################SOLUTION##########################
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
      out.suffix = paste0("pathview_", i),
      kegg.native = TRUE
    )
  }, error = function(e) {
    cat("  → Failed for:", pathway.id[i], "\n  Reason:", e$message, "\n")
  })
}


