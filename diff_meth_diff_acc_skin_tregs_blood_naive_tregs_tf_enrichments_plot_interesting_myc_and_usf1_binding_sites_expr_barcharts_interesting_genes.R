# This script generates expression bar charts for genes from interesting
#   track plots arount c-Mac/USF1 binding sites.
# Author: Niklas Beumer


# Load required packages.
library(ggplot2)
library(tidyr)
library(testit)


# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the RNA sample mapping file and restrict to blood naive Tregs and 
# skin Tregs.
relevant_cell_types <- c("Skin_Treg", "Blood_naive_Treg")
rna_sample_mapping_path <- paste0(location, 
                                  "/sample_mapping_rnaseq_only_biol_rep.txt")
rna_sample_mapping <- read.table(rna_sample_mapping_path, header = T, 
                                 stringsAsFactors = F, sep = "\t")
rna_sample_mapping <- rna_sample_mapping[
  rna_sample_mapping$Cell_type_wo_blank %in% relevant_cell_types, ]
relevant_cell_types_names <- sapply(relevant_cell_types, FUN = function(x) {
  unique(rna_sample_mapping$Cell_type[
    rna_sample_mapping$Cell_type_wo_blank == x])
})

# Read in the TPM values.
tpm_file <- paste0(location, "/RNASeq/tpm/tpm_all.txt")
tpm <- read.table(tpm_file, header = T, stringsAsFactors = F)

# Specify the genes to plot bar charts for.
genes_to_plot <- c("MLPH", "NCALD", "GNA11", "SHANK3")





############################################################################
# Generate bar charts showing TPM values for the relevant genes.
############################################################################

# Identify the RNA samples corresponding to the relevant cell types.
relevant_rna_samples <- rna_sample_mapping$Sample

# Restrict the TPM values to the relevant samples.
tpm_restr <- tpm[, c("Gene_symbol", relevant_rna_samples)]

# Iterate over all relevant genes.
for (gene in genes_to_plot) {
  
  # Extract TPM values for this gene.
  tpm_gene <- tpm_restr[tpm_restr$Gene_symbol == gene, ]
  tpm_gene_melted <- pivot_longer(tpm_gene, cols = 2:ncol(tpm_gene), 
                                  names_to = "Sample", values_to = "TPM")
  tpm_gene_melted$Cell_type <- factor(
    sapply(tpm_gene_melted$Sample, FUN = function(x) {
      cell_type_wo_blank <- rna_sample_mapping$Cell_type_wo_blank[
        rna_sample_mapping$Sample == x]
      return(relevant_cell_types_names[cell_type_wo_blank])
      }), 
    levels = rev(relevant_cell_types_names)
  )
  
  # Compute log(100 * TPM + 1)
  tpm_gene_melted$TPM_log <- log1p(100 * tpm_gene_melted$TPM)
  
  # For each cell type, compute the mean log(100 * TPM + 1) and its standard 
  # deviation.
  means <- sapply(relevant_cell_types_names, FUN = function(x) {
    mean(tpm_gene_melted$TPM_log[tpm_gene_melted$Cell_type == x])
  })
  sds <- sapply(relevant_cell_types_names, FUN = function(x) {
    sd(tpm_gene_melted$TPM_log[tpm_gene_melted$Cell_type == x])
  })
  means_sds_comb <- data.frame(Cell_type = relevant_cell_types_names,
                               TPM_log = means,
                               SD = sds)
  
  # Compute the upper and lower bounds of the error bars.
  means_sds_comb$mean_minus_sd <- means_sds_comb$TPM_log - means_sds_comb$SD
  means_sds_comb$mean_plus_sd <- means_sds_comb$TPM_log + means_sds_comb$SD
  
  # Generate the bar chart.
  means_sds_comb$Cell_type <- factor(means_sds_comb$Cell_type, 
                                     levels = rev(relevant_cell_types_names))
  expr_plot <- ggplot() +
    aes(x = Cell_type, y = TPM_log, fill = Cell_type) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.05)),
                       name = "Gene expression [log(100 * TPM + 1)]") +
    scale_fill_manual(breaks = relevant_cell_types_names,
                      values = c("darkorchid1", "red"),
                      name = "Cell type") +
    geom_bar(data = means_sds_comb, stat = "identity", width = 0.7) +
    geom_errorbar(data = means_sds_comb, width = 0.3,
                  mapping = aes(ymin = mean_minus_sd, ymax = mean_plus_sd)) +
    geom_point(data = tpm_gene_melted, 
               position = position_jitter(width = 0.3, height = 0, seed = 5)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    xlab("Cell type") +
    coord_flip() +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "none")
  
  # Save the expression plot.
  expr_plot_pdf <- paste0(plot_outdir, 
                          "/skin_treg_blood_naive_Treg_expr_",
                          gene,
                          ".pdf")
  expr_plot_rds <- paste0(plot_rds_outdir, 
                          "/skin_treg_blood_naive_Treg_expr_",
                          gene,
                          ".rds")
  pdf(expr_plot_pdf, width = 5.5, height = 1.5)
  print(expr_plot)
  dev.off()
  saveRDS(expr_plot, file = expr_plot_rds)
  
}
