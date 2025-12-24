###################### BAR GRAPH ###################


####################################################
#Script by: Gillian Myers
#Date: 11/10/2025
####################################################


# Yeast Metabolic Pathway Log2 Fold Change Barplot Generator
setwd("C:/Users/gmyer/Fall2025 UR/")

# ------------------------------
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# ------------------------------
# Function to create barplots (significant genes only)
create_pathway_barplot <- function(pathway_name, pathway_gene_mapping,
                                   res_PA5vGOB8, res_PA5vGOB21, res_GOB21vGOB8,
                                   res_PM5vGOB8, res_PO5LvGOB8, res_PO5SvGOB8, res_GOBGLYvGOB8) {
  cat("\nCreating barplot for", pathway_name, "...\n")
  
  # Extract systematic names
  systematic_names <- unlist(pathway_gene_mapping)
  systematic_names_with_suffix <- paste0(systematic_names, "_mRNA")
  
  # Collect all matching genes across datasets
  gene_lists <- list(
    rownames(res_PA5vGOB8),
    rownames(res_PA5vGOB21),
    rownames(res_GOB21vGOB8),
    rownames(res_PM5vGOB8),
    rownames(res_PO5LvGOB8),
    rownames(res_PO5SvGOB8),
    rownames(res_GOBGLYvGOB8)
  )
  common_genes <- Reduce(union, lapply(gene_lists, function(x) intersect(x, systematic_names_with_suffix)))
  
  if (length(common_genes) == 0) {
    cat("No matching genes found for", pathway_name, "\n")
    return()
  }
  
  # Helper to extract significant (padj < 0.05) results
  get_sig_df <- function(res, treatment_name) {
    if (!all(c("log2FoldChange", "padj") %in% colnames(res))) return(NULL)
    res_subset <- res[common_genes, c("log2FoldChange", "padj")]
    res_subset <- res_subset[!is.na(res_subset$padj) & res_subset$padj < 0.05, , drop = FALSE]
    if (nrow(res_subset) == 0) return(NULL)
    data.frame(Gene = rownames(res_subset),
               Treatment = treatment_name,
               log2FoldChange = res_subset$log2FoldChange,
               stringsAsFactors = FALSE)
  }
  
  # Combine all significant results
  df_sig <- dplyr::bind_rows(
    get_sig_df(res_PA5vGOB8, "PA5_t1000_vs_GOB8"),
    get_sig_df(res_PA5vGOB21, "PA5_t1000_vs_GOB21"),
    get_sig_df(res_GOB21vGOB8, "GOB21_vs_GOB8"),
    get_sig_df(res_PM5vGOB8, "PM5_t1000_vs_GOB8"),
    get_sig_df(res_PO5LvGOB8, "PO5_t1000_L_vs_GOB8"),
    get_sig_df(res_PO5SvGOB8, "PO5_t1000_S_vs_GOB8"),
    get_sig_df(res_GOBGLYvGOB8, "GOB8_Gly_vs_GOB8")
  )
  
  if (is.null(df_sig) || nrow(df_sig) == 0) {
    cat("No significant genes (padj < 0.05) for", pathway_name, "\n")
    return()
  }
  
  # Clean up gene names (convert systematic → common)
  df_sig$Gene <- gsub("_mRNA$", "", df_sig$Gene)
  df_sig$Gene <- sapply(df_sig$Gene, function(g) {
    common <- names(pathway_gene_mapping)[pathway_gene_mapping == g]
    if (length(common) > 0) common else g
  })
  
  # ---------------------------
  # Treatment order
  # -------------------------
  # Set custom treatment order
  df_sig$Treatment <- factor(df_sig$Treatment, levels = c(
    "GOB21_vs_GOB8",
    "PA5_t1000_vs_GOB21", "PA5_t1000_vs_GOB8",
    "PM5_t1000_vs_GOB8",
    "PO5_t1000_L_vs_GOB8", "PO5_t1000_S_vs_GOB8",
    "GOB8_Gly_vs_GOB8"
  ))
  
  # -------------------------
  # Custom color palette (modern, colorblind-friendly)
  custom_colors <- c(
    "GOB21_vs_GOB8"        = "#1B9E77",  # teal
    "PA5_t1000_vs_GOB21"   = "#66A61E",  # olive green
    "PA5_t1000_vs_GOB8"    = "#A6D854",  # light green
    "PM5_t1000_vs_GOB8"    = "#E6AB02",  # goldenrod
    "PO5_t1000_L_vs_GOB8"  = "#FC8D62",  # soft salmon
    "PO5_t1000_S_vs_GOB8"  = "#E78AC3",  # bright pink
    "GOB8_Gly_vs_GOB8"     = "#8DA0CB"   # periwinkle
  )
  
  # -------------------------
  # Plot barplot
  p <- ggplot(df_sig, aes(x = Gene, y = log2FoldChange, fill = Treatment)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.85)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = paste("Significant log2 Fold Change (<0.05) for", gsub("_", " ", pathway_name)),
         y = "log2 Fold Change",
         x = "Gene") +
    scale_fill_manual(values = custom_colors) +
    coord_cartesian(ylim = c(-3, 3)) 
  
  
  # Save the plot
  output_file <- paste0(pathway_name, "_barplot_sig_CHECK.png")
  ggsave(output_file, plot = p,
         width = max(5, 0.4 * length(unique(df_sig$Gene))),
         height = 4,
         dpi = 300)
  
  cat("Saved significant barplot to", output_file, "\n")
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
res_PA5_t1000_vs_GOB21   <- read.csv("PA5_t1000_vs_GOB21_results.csv", row.names = 1)
res_GOB21_vs_GOB8        <- read.csv("GOB21_vs_GOB8_results.csv", row.names = 1)
res_PM5_t1000_vs_GOB8    <- read.csv("PM5_t1000_vs_GOB8_results.csv", row.names = 1)
res_PO5_t1000_L_vs_GOB8  <- read.csv("PO5_t1000_L_vs_GOB8_results.csv", row.names = 1)
res_PO5_t1000_S_vs_GOB8  <- read.csv("PO5_t1000_S_vs_GOB8_results.csv", row.names = 1)
res_GOB8_Gly_vs_GOB8     <- read.csv("GOB8_Gly_vs_GOB8_results.csv", row.names = 1)

cat("Loaded all datasets successfully.\n")

# ------------------------------
# Generate barplots
cat("\nGenerating pathway barplots (significant genes only)...\n")
for (pathway_name in names(pathway_genes_evolved)) {
  create_pathway_barplot(pathway_name, pathway_genes_evolved[[pathway_name]],
                         res_PA5_t1000_vs_GOB8, res_PA5_t1000_vs_GOB21, res_GOB21_vs_GOB8,
                         res_PM5_t1000_vs_GOB8, res_PO5_t1000_L_vs_GOB8, res_PO5_t1000_S_vs_GOB8, res_GOB8_Gly_vs_GOB8)
}

cat("\nAnalysis complete!\n")

