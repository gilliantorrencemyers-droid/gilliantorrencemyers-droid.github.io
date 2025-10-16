---
title: "Evolved Heatmaps"
permalink: /r-undergrad/r-undergrad-scripts-Evolved-Heatmaps/
---

## R Script

Download it here!

- [Script_04_evolved_heatmaps_August2025.R](/assets/Script_04_evolved_heatmaps_August2025.R)

################### Heatmaps for Evolved Strands ######################

#======================================================================
#Script By: Gillian Myers
#Date: 09/10/2025
#======================================================================

# Yeast Metabolic Pathway Heatmap Generator
setwd("C:/Users/gmyer/Fall2025 UR/")

# Load required libraries

library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)

# Function to create pathway heatmap
create_pathway_heatmap <- function(pathway_name, pathway_gene_mapping, res_PA5vGOB8, res_PA5vGOB21, res_GOB21vGOB8, res_PM5vGOB8, res_PO5LvGOB8, res_PO5SvGOB8, res_GOBGLYvGOB8) {
  cat("\nCreating heatmap for", pathway_name, "...\n")
  
  # Extract systematic names and add _mRNA suffix
  systematic_names_evolved <- unlist(pathway_gene_mapping)
  systematic_names_with_suffix_evolved <- paste0(systematic_names_evolved, "_mRNA")
  
  # Find genes in either dataset (FIX: union instead of intersect)
  found_PA5vGOB8<- intersect(rownames(res_PA5vGOB8), systematic_names_with_suffix_evolved)
  found_PA5vGOB21<- intersect(rownames(res_PA5vGOB21), systematic_names_with_suffix_evolved)
  found_PM5vGOB8<- intersect(rownames(res_PM5vGOB8), systematic_names_with_suffix_evolved)
  found_GOB21vGOB8<- intersect(rownames(res_GOB21vGOB8), systematic_names_with_suffix_evolved)
  
  found_PO5LvGOB8<- intersect(rownames(res_PO5LvGOB8), systematic_names_with_suffix_evolved)
  found_PO5SvGOB8<- intersect(rownames(res_PO5SvGOB8), systematic_names_with_suffix_evolved)
  found_GOBGLYvGOB8<- intersect(rownames(res_GOBGLYvGOB8), systematic_names_with_suffix_evolved)
  
  #More than 2 arguements used, cannot use union without reduced
  common_genes_evolved <- Reduce(union, list(found_PA5vGOB21, found_PA5vGOB8, found_PM5vGOB8, found_PO5LvGOB8, found_PO5SvGOB8, found_GOB21vGOB8, found_GOBGLYvGOB8))
  
  
  if(length(common_genes_evolved) == 0) {
    cat("Warning: No genes from", pathway_name, "found in either dataset\n")
    return()
  }
  
  cat("Found", length(common_genes_evolved), "genes\n")
  
  # Create matrix with log2FoldChange values
  df_evolved <- data.frame(
    PA5_t1000_vs_GOB8 = res_PA5vGOB8[common_genes_evolved, "log2FoldChange"],
    PA5_t1000_vs_GOB21 = res_PA5vGOB21[common_genes_evolved, "log2FoldChange"],
    GOB21_vs_GOB8 = res_GOB21vGOB8[common_genes_evolved, "log2FoldChange"],
    
    PM5_t1000_vs_GOB8 = res_PM5vGOB8[common_genes_evolved, "log2FoldChange"],
    PO5_t1000_L_vs_GOB8 = res_PO5LvGOB8[common_genes_evolved, "log2FoldChange"],
    PO5_t1000_S_vs_GOB8 = res_PO5SvGOB8[common_genes_evolved, "log2FoldChange"],
    GOBGLY_vs_GOB8 = res_GOBGLYvGOB8[common_genes_evolved, "log2FoldChange"],
    
    row.names = common_genes_evolved
  )
  
  # Remove rows with all NA values
  df_evolved <- df_evolved[rowSums(is.na(df_evolved)) < ncol(df_evolved), ]
  if(nrow(df_evolved) == 0) {
    cat("No valid data after removing NAs\n")
    return()
  }
  
  # Convert systematic names → common names
  gene_symbols_evolved <- sapply(rownames(df_evolved), function(gene_name) {
    sys_name_clean_evolved <- gsub("_mRNA$", "", gene_name)
    common_name_evolved <- names(pathway_gene_mapping)[pathway_gene_mapping == sys_name_clean_evolved]
    if(length(common_name_evolved) > 0) {
      return(common_name_evolved[1])
    } else {
      return(sys_name_clean_evolved)
    }
  })
  
  # Create matrix for heatmap
  l2_val_evolved <- as.matrix(df_evolved)
  rownames(l2_val_evolved) <- gene_symbols_evolved
  
  # Get p-values for masking
  p_PA5vGOB8<- res_PA5vGOB8[common_genes_evolved, "padj"]
  p_PA5vGOB21<- res_PA5vGOB21[common_genes_evolved, "padj"]
  p_GOB21vGOB8<- res_GOB21vGOB8[common_genes_evolved, "padj"]
  
  p_PM5vGOB8<- res_PM5vGOB8[common_genes_evolved, "padj"]
  p_PO5LvGOB8<- res_PO5LvGOB8[common_genes_evolved, "padj"]
  p_PO5SvGOB8<- res_PO5SvGOB8[common_genes_evolved, "padj"]
  p_GOBGLYvGOB8<- res_GOBGLYvGOB8[common_genes_evolved, "padj"]
  
  # Mask non-significant (padj > 0.05) to 0 (white)
  if("PA5_t1000_vs_GOB8" %in% colnames(l2_val_evolved)) {
    l2_val_evolved[, "PA5_t1000_vs_GOB8"][is.na(p_PA5vGOB8) | p_PA5vGOB8 > 0.05] <- 0
  }
  if("PA5_t1000_vs_GOB21" %in% colnames(l2_val_evolved)) {
    l2_val_evolved[, "PA5_t1000_vs_GOB21"][is.na(p_PA5vGOB21) | p_PA5vGOB21 > 0.05] <- 0
  }
  if("GOB21_vs_GOB8" %in% colnames(l2_val_evolved)) {
    l2_val_evolved[, "GOB21_vs_GOB8"][is.na(p_GOB21vGOB8) | p_GOB21vGOB8 > 0.05] <- 0
  }
  
  
  if("PM5_t1000_vs_GOB8" %in% colnames(l2_val_evolved)) {
    l2_val_evolved[, "PM5_t1000_vs_GOB8"][is.na(p_PM5vGOB8) | p_PM5vGOB8 > 0.05] <- 0
  }
  if("PO5_t1000_L_vs_GOB8" %in% colnames(l2_val_evolved)) {
    l2_val_evolved[, "PO5_t1000_L_vs_GOB8"][is.na(p_PO5LvGOB8) | p_PO5LvGOB8 > 0.05] <- 0
  }
 
  if("GOB8_Gly_vs_GOB8" %in% colnames(l2_val_evolved)) {
    l2_val_evolved[, "GOB8_Gly_vs_GOB8"][is.na(p_GOBGLYvGOB8) | p_GOBGLYvGOB8 > 0.05] <- 0
  }
  
  # --- Consistent color scale (white = 0) ---
  max_abs <- max(abs(l2_val_evolved), na.rm = TRUE)
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  colors <- colorRampPalette(c("darkblue", "blue", "lightblue", "white",
                               "lightcoral", "red", "darkred"))(100)
  
  # Output file
  output_file_evolved <- paste0(pathway_name, "_heatmap.png")
  
  # Create heatmap
  p_evolved <- pheatmap(l2_val_evolved,
                color = colors,
                breaks = breaks,
                scale = "none",
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                show_rownames = TRUE,
                show_colnames = TRUE,
                fontsize = 8,
                fontsize_row = 8,
                fontsize_col = 9,
                cellwidth = 35,
                cellheight = 12,
                display_numbers = TRUE,
                number_format = "%.2f",
                main = paste0(gsub("_", " ", tools::toTitleCase(pathway_name))),
                angle_col = 45,
                labels_col = c("PA5 t1000 vs GOB8", "PA5 t1000 vs GOB21", "GOB21 vs GOB8", "PM5 t1000 vs GOB8", "PO t1000 L vs GOB8","PO t1000 S vs GOB8", "GOB8 GLY vs GOB8"),
                border_color = "white")
  
  # Save file with dynamic sizing
  plot_width  <- max(3, nrow(l2_val_evolved) * 0.15 + 2.5)
  plot_height <- max(2, nrow(l2_val_evolved) * 0.25 + 1.5)
  ggsave(output_file_evolved, plot = p_evolved, width = plot_width, height = plot_height,
         dpi = 300, units = "in")
  
  cat("Saved to", output_file_evolved, "\n")
}

# ------------------------------
# Define metabolic pathways
pathway_genes_evolved <- list(
  mitochondrial_crista = list(
    "MGM1" = "YOR211C", "MIC60" = "YKR016W", "OXA1" = "YER154W", "QCR2" = "YPR191W"
  ),
  mitochondrial_DNA = list(
    "COX1" = "Q0045", "ATP8" = "Q0080", "ATP6" = "Q0085", "COB" = "Q0105",
    "OLI1" = "Q0130", "VAR1" = "Q0140", "SCEI" = "Q0160", "COX2" = "Q0250", "Q0255" = "Q0255"
  ),
  threonine_methionine_biosynthesis = list(
    "AAT1" = "YKL106W", "AAT2" = "YLR027C", "HOM2" = "YDR158W", "HOM3" = "YER052C",
    "HOM6" = "YJR139C", "MET17" = "YLR303W", "MET2" = "YNL277W", "MET6" = "YER091C",
    "MET7" = "YOR241W", "THR1" = "YHR025W", "THR4" = "YCR053W"
  ),
  sulfur_aa_biosynthesis = list(
    "CYS3" = "YAL012W", "CYS4" = "YGR155W", "HOM2" = "YDR158W", "HOM3" = "YER052C",
    "HOM6" = "YJR139C", "MET10" = "YFR030W", "MET14" = "YKL001C", "MET16" = "YPR167C",
    "MET17" = "YLR303W", "MET2" = "YNL277W", "MET3" = "YJR010W", "MET5" = "YJR137C",
    "MET6" = "YER091C", "MET7" = "YOR241W", "SAM1" = "YLR180W", "SAM2" = "YDR502C",
    "STR2" = "YJR130C", "STR3" = "YGL184C"
  ),
  branched_chain_aa_biosynthesis = list(
    "BAT1" = "YHR208W", "BAT2" = "YJR148W", "ILV1" = "YER086W", "ILV2" = "YMR108W",
    "ILV3" = "YJR016C", "ILV5" = "YLR355C", "ILV6" = "YCL009C", "LEU1" = "YGL009C",
    "LEU2" = "YCL018W", "LEU4" = "YNL104C", "LEU9" = "YOR108W"
  ),
  aromatic_aa_biosynthesis = list(
    "ARO1" = "YDR127W", "ARO2" = "YGL148W", "ARO3" = "YDR035W", "ARO4" = "YBR249C",
    "ARO7" = "YPR060C", "ARO8" = "YGL202W", "ARO9" = "YHR137W", "PHA2" = "YNL316C",
    "TRP1" = "YDR007W", "TRP2" = "YER090W", "TRP3" = "YKL211C", "TRP4" = "YDR354W",
    "TRP5" = "YGL026C", "TYR1" = "YBR166C"
  ),
  ethanol_degradation = list(
    "ACS1" = "YAL054C", "ACS2" = "YLR153C", "ADH2" = "YMR303C", "ALD2" = "YMR170C"
  ),
  glycerol_degradation = list(
    "GUT1" = "YHL032C", "GUT2" = "YIL155C"
  ),
  aerobic_respiration_ETC = list(
    "COB" = "Q0105", "COR1" = "YBL045C", "COX1" = "Q0045", "COX12" = "YLR038C",
    "COX13" = "YGL191W", "COX2" = "Q0250", "COX3" = "Q0275", "COX4" = "YGL187C",
    "COX5A" = "YNL052W", "COX6" = "YHR051W", "COX7" = "YMR256C", "COX8" = "YLR395C",
    "COX9" = "YDL067C", "CYT1" = "YOR065W", "NDI1" = "YML120C", "QCR2" = "YPR191W",
    "QCR6" = "YFR033C", "QCR7" = "YDR529C", "QCR8" = "YJL166W", "QCR9" = "YGR183C",
    "RIP1" = "YEL024W", "SDH1" = "YKL148C", "SDH2" = "YLL041C", "SDH3" = "YKL141W",
    "SDH4" = "YDR178W"
  ),
  superoxide_degradation = list(
    "CTA1" = "YDR256C", "CTT1" = "YGR088W", "SOD1" = "YJR104C", "SOD2" = "YHR008C"
  ),
  TCA_aerobic_respiration = list(
    "ACO1" = "YLR304C", "ACO2" = "YJL200C", "CIT1" = "YNR001C", "CIT3" = "YPR001W",
    "FUM1" = "YPL262W", "IDH1" = "YNL037C", "IDH2" = "YOR136W", "KGD1" = "YIL125W",
    "KGD2" = "YDR148C", "LPD1" = "YFL018C", "LSC1" = "YOR142W", "LSC2" = "YGR244C",
    "MAE1" = "YKL029C", "MDH1" = "YKL085W", "PYC1" = "YGL062W", "PYC2" = "YBR218C",
    "SDH1" = "YKL148C", "SDH2" = "YLL041C", "SDH3" = "YKL141W", "SDH4" = "YDR178W"
  ),
  TCA_glyoxylate_cycle = list(
    "ACO1" = "YLR304C", "ACO2" = "YJL200C", "CIT1" = "YNR001C", "CIT2" = "YCR005C",
    "CIT3" = "YPR001W", "DAL7" = "YIR031C", "FUM1" = "YPL262W", "ICL1" = "YER065C",
    "IDH1" = "YNL037C", "IDH2" = "YOR136W", "KGD1" = "YIL125W", "KGD2" = "YDR148C",
    "LPD1" = "YFL018C", "LSC1" = "YOR142W", "LSC2" = "YGR244C", "MAE1" = "YKL029C",
    "MDH1" = "YKL085W", "MDH2" = "YOL126C", "MDH3" = "YDL078C", "MLS1" = "YNL117W",
    "PYC1" = "YGL062W", "PYC2" = "YBR218C", "SDH1" = "YKL148C", "SDH2" = "YLL041C",
    "SDH3" = "YKL141W", "SDH4" = "YDR178W"
  ),
  glycolysis = list(
    "CDC19" = "YAL038W", "ENO1" = "YGR254W", "ENO2" = "YHR174W", "FBA1" = "YKL060C",
    "FBP1" = "YLR377C", "GPM1" = "YKL152C", "PFK1" = "YGR240C", "PFK2" = "YMR205C",
    "PGI1" = "YBR196C", "PGK1" = "YCR012W", "PYK2" = "YOR347C", "TDH1" = "YJL052W",
    "TDH2" = "YJR009C", "TDH3" = "YGR192C", "TPI1" = "YDR050C"
  )
)

# ------------------------------
# Load DE results
cat("Loading differential expression results...\n")
res_PA5_t1000_vs_GOB8    <- read.csv("PA5_t1000_vs_GOB8_results.csv", row.names = 1)
res_PA5_t1000_vs_GOB21  <- read.csv("PA5_t1000_vs_GOB21_results.csv", row.names = 1)
res_GOB21_vs_GOB8 <- read.csv("GOB21_vs_GOB8_results.csv", row.names = 1)

res_PM5_t1000_vs_GOB8 <- read.csv("PM5_t1000_vs_GOB8_results.csv", row.names = 1)
res_PO5_t1000_L_vs_GOB8 <- read.csv("PO5_t1000_L_vs_GOB8_results.csv", row.names = 1)
res_PO5_t1000_S_vs_GOB8 <- read.csv("PO5_t1000_S_vs_GOB8_results.csv", row.names = 1)
res_GOB8_Gly_vs_GOB8 <- read.csv("GOB8_Gly_vs_GOB8_results.csv", row.names = 1)

cat("Loaded PA5_t1000_vs_GOB8:", nrow(res_PA5_t1000_vs_GOB8), "genes\n")
cat("Loaded PA5_t1000_vs_GOB21:", nrow(res_PA5_t1000_vs_GOB21), "genes\n")
cat("Loaded GOB21_vs_GOB8:", nrow(res_GOB21_vs_GOB8), "genes\n")

cat("Loaded PM5_t1000_vs_GOB8:", nrow(res_PM5_t1000_vs_GOB8), "genes\n")
cat("Loaded PO5_t1000_L_vs_GOB8:", nrow(res_PO5_t1000_L_vs_GOB8), "genes\n")
cat("Loaded PO5_t1000_S_vs_GOB8:", nrow(res_PO5_t1000_S_vs_GOB8), "genes\n")
cat("Loaded GOB8_Gly_vs_GOB8:", nrow(res_GOB8_Gly_vs_GOB8), "genes\n")

# ------------------------------
# Generate heatmaps

cat("\nGenerating pathway heatmaps...\n")
for (pathway_name in names(pathway_genes_evolved)) {
  create_pathway_heatmap(pathway_name, pathway_genes_evolved[[pathway_name]],
                         res_PA5_t1000_vs_GOB8, res_PA5_t1000_vs_GOB21, res_GOB21_vs_GOB8, res_PM5_t1000_vs_GOB8, res_PO5_t1000_L_vs_GOB8, res_PO5_t1000_S_vs_GOB8, res_GOB8_Gly_vs_GOB8)
}
cat("\nAnalysis complete!\n")
