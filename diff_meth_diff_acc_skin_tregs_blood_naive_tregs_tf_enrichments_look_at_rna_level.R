# This script looks at RNA evidence regarding transcription factors that are
#   interesting regarding the differentiation from blood naive Treg cells to
#   skin Treg cells. This includes footprint analysis by decoupleR and Dorothea
#   and simple expression analysis.
# Author: Niklas Beumer



# Load required packages.
library(decoupleR)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(testit)
# Note: For this script to work, the package OmnipathR needs to be installed.


# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis",
                     format(Sys.time(), "%m-%d-%y"),
                     sep = "/")
if (!dir.exists(plot_outdir)) {
  dir.create(plot_outdir)
}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the DESeq2 results for the comparison of
# skin Tregs and blood naive Tregs.
diff_exp_file <- paste0(
  location,
  "/RNASeq/analysis_results/2022-01-14_diff_gene_expr_DESEq2_Skin_TregBlood_naive_Treg_results_filtered_with_significance.txt"
)
diff_exp <- read.table(diff_exp_file,
                       header = T,
                       stringsAsFactors = F)

# Read in the RNA sample mapping file and restrict to skin Tregs and blood 
# naive Tregs.
rna_sample_mapping_path <- paste0(location, 
                                  "/sample_mapping_rnaseq_only_biol_rep.txt")
rna_sample_mapping <- read.table(
  rna_sample_mapping_path,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)
cell_types_to_use <- c("Skin_Treg", "Blood_naive_Treg")
rna_sample_mapping <- rna_sample_mapping[
  rna_sample_mapping$Cell_type_wo_blank %in% cell_types_to_use, 
]

# Read in the TPM values.
tpm_file <- paste0(location, "/RNASeq/tpm/tpm_all.txt")
tpm <- read.table(tpm_file, header = T, stringsAsFactors = F)




##########################################################
# Perform transcription factor footprint analysis.
##########################################################

# Get the Dorothea network.
# By default, this extracts the data for th human system and regulons
# in the confidence levels A, B and C.
dorothea_data <- get_dorothea()

# Run the fgsea method for footprint analysis.
# To this end, use the Log2FCs as the statistic for differential expression.
deg_mat <- matrix(diff_exp$log2FoldChange, ncol = 1)
rownames(deg_mat) <- rownames(diff_exp)
footprint_res_fgsea <- run_fgsea(mat = deg_mat, network = dorothea_data)

# Perform multiple testing correction of the resulting P values.
# (Benjamini-Hochberg method)
footprint_res_fgsea_smaller <- footprint_res_fgsea[footprint_res_fgsea$statistic == "fgsea", ]
footprint_res_fgsea_smaller$p_value_adj_BH <- p.adjust(footprint_res_fgsea_smaller$p_value, method = "BH")
footprint_res_fgsea$p_value_adj_BH <- sapply(
  1:nrow(footprint_res_fgsea),
  FUN = function(x) {
    curr_tf <- footprint_res_fgsea$source[x]
    return(footprint_res_fgsea_smaller$p_value_adj_BH[footprint_res_fgsea_smaller$source == curr_tf])
  }
)

# Save the footprint analysis results.
footprint_res_fgsea_outfile <-
  paste0(
    location,
    "/treg_hierarchies/diff_meth_diff_acc_skin_tregs_blood_naive_tregs_tf_enrichments_corresp_rna_footprints.txt"
  )
write.table(
  footprint_res_fgsea,
  file = footprint_res_fgsea_outfile,
  sep = "\t",
  row.names = F,
  col.names = T,
  quote = F
)

# Select interesting transcription factors to highlight in the plot.
interesting_tfs <- c(
  "ETS1",
  "FLI1",
  "SPI1",
  "MYOD1",
  "TCF12",
  "MYC",
  "MYCN",
  "USF1",
  "JUN",
  "FOSL1",
  "FOSL2",
  "JUNB",
  "ATF3",
  "BACH1"
)

# Generate and save a volcano plot for the footprint analysis results.
footprint_res_fgsea_smaller$logp <- -log10(footprint_res_fgsea_smaller$p_value_adj_BH)
max_logp <- max(footprint_res_fgsea_smaller$logp[footprint_res_fgsea_smaller$logp < Inf])
replacement_val <- 1.5 * max_logp
footprint_res_fgsea_smaller$logp[footprint_res_fgsea_smaller$logp == Inf] <- replacement_val
footprint_res_fgsea_smaller$signif <- as.character(footprint_res_fgsea_smaller$p_value_adj_BH < 0.05)
text_data <- data.frame(
  TF = interesting_tfs,
  logp = sapply(
    interesting_tfs,
    FUN = function(x) {
      footprint_res_fgsea_smaller$logp[footprint_res_fgsea_smaller$source == x]
    }
  ),
  score = sapply(
    interesting_tfs,
    FUN = function(x) {
      footprint_res_fgsea_smaller$score[footprint_res_fgsea_smaller$source == x]
    }
  )
)
step_size <- 10 ^ floor(log10(max_logp) - 0.3)
y_breaks <- seq(0, max_logp, step_size)
volcano_plot <- ggplot(footprint_res_fgsea_smaller) +
  aes(x = score, y = logp) +
  scale_y_continuous(
    breaks = c(y_breaks, replacement_val),
    labels = c(y_breaks, "Inf"),
    expand = expansion(mult = c(0.01, 0.05)),
    name = "-log10(q)"
  ) +
  scale_colour_manual(breaks = c("TRUE", "FALSE"),
                      values = c("red", "grey50")) +
  geom_point(size = 0.5, mapping = aes(colour = signif)) +
  geom_text_repel(
    data = text_data,
    mapping = aes(label = TF),
    size = 3,
    min.segment.length = 0,
    force = 100,
    force_pull = 2
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab("Enrichment score") +
  theme_classic() +
  theme(
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.position = "none"
  )
outfile_pdf <- paste0(
  plot_outdir,
  "/treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_footprint_analysis_on_rna_level.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir,
  "/treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_footprint_analysis_on_rna_level.rds"
)
pdf(outfile_pdf, width = 4, height = 4)
print(volcano_plot)
dev.off()
saveRDS(volcano_plot, file = outfile_rds)








#################################################################
# Visualise TPM values for relevant transcription factors
# This is to show whether they display detectable expression
# in the two cell types.
#################################################################

# Identify the RNA samples corresponding to the relevant cell types.
relevant_rna_samples <- rna_sample_mapping$Sample

# Restrict the TPM values to the relevant samples.
tpm_restr <- tpm[, c("Gene_symbol", relevant_rna_samples)]

# Specify the transcription factors to plot values for.
# Fusion TFs in the Homer results are ignored.
# PU.1 = SPI1 (https://www.genecards.org/cgi-bin/carddisp.pl?gene=SPI1&keywords=PU.1; 21-Feb-2023)
# HIF2a = EPAS1 (https://www.genecards.org/cgi-bin/carddisp.pl?gene=EPAS1&keywords=HIF2A; 21-Feb-2023)
tfs_to_visualise <- c(
  "ETS1",
  "FLI1",
  "SPI1",
  "ETV2",
  "GABPA",
  "ETV1",
  "MYOD1",
  "MYF5",
  "TCF12",
  "TCF21",
  "MYOG",
  "MYC",
  "USF1",
  "MYCN",
  "EPAS1",
  "JUN",
  "FOSL2",
  "JUNB",
  "FOSL1",
  "ATF3",
  "BATF",
  "BACH2",
  "BACH1"
)

# Extract TPM values for the relevant transcription factors.
tpm_tfs <- tpm_restr[tpm_restr$Gene_symbol %in% tfs_to_visualise, ]
tpm_tfs_melted <- melt(
  tpm_tfs,
  id.vars = "Gene_symbol",
  variable.name = "Sample",
  value.name = "TPM"
)
tpm_tfs_melted$Cell_type <- factor(
  sapply(
    tpm_tfs_melted$Sample,
    FUN = function(x) {
      ifelse(grepl("Skin", x), yes = "Skin Treg", no = "Blood naive Treg")
    }
  ),
  levels = c("Blood naive Treg", "Skin Treg")
)

# Compute log(100 * TPM + 1)
tpm_tfs_melted$TPM_log <- log1p(100 * tpm_tfs_melted$TPM)

# For each cell type and each transcription factor, compute the mean 100 * log(TPM + 1)
# and its standard deviation.
means <- melt(sapply(
  c("Skin Treg", "Blood naive Treg"),
  FUN = function(x) {
    sapply(
      tfs_to_visualise,
      FUN = function(y) {
        mean(tpm_tfs_melted$TPM_log[tpm_tfs_melted$Cell_type == x &
                                      tpm_tfs_melted$Gene_symbol == y])
      }
    )
  }
))
sds <- melt(sapply(
  c("Skin Treg", "Blood naive Treg"),
  FUN = function(x) {
    sapply(
      tfs_to_visualise,
      FUN = function(y) {
        sd(tpm_tfs_melted$TPM_log[tpm_tfs_melted$Cell_type == x &
                                    tpm_tfs_melted$Gene_symbol == y])
      }
    )
  }
))
assert(all(means$Var1 == sds$Var1))
assert(all(means$Var2 == sds$Var2))
means_sds_comb <- cbind(means, sds[, 3])

# Improve column names.
colnames(means_sds_comb) <- c("Gene_symbol", "Cell_type", "TPM_log", "SD")

# Compute the upper and lower bounds of the error bars.
means_sds_comb$mean_minus_sd <- means_sds_comb$TPM_log - means_sds_comb$SD
means_sds_comb$mean_plus_sd <- means_sds_comb$TPM_log + means_sds_comb$SD

# Visualise the results and save the resulting plots.
tpm_tfs_melted$Gene_symbol <- factor(tpm_tfs_melted$Gene_symbol, levels = tfs_to_visualise)
means_sds_comb$Gene_symbol <- factor(means_sds_comb$Gene_symbol, levels = tfs_to_visualise)
means_sds_comb$Cell_type <- factor(means_sds_comb$Cell_type,
                                   levels = c("Blood naive Treg", "Skin Treg"))
expr_plot <- ggplot() +
  aes(x = Gene_symbol, y = TPM_log, fill = Cell_type) +
  scale_fill_manual(
    breaks = c("Blood naive Treg", "Skin Treg"),
    values = c("darkorchid1", "blue"),
    name = "Cell type"
  ) +
  geom_bar(
    data = means_sds_comb,
    stat = "identity",
    position = "dodge",
    width = 0.7
  ) +
  geom_errorbar(
    data = means_sds_comb,
    width = 0.3,
    position = position_dodge(width = 0.7),
    mapping = aes(ymin = mean_minus_sd, ymax = mean_plus_sd)
  ) +
  geom_point(
    data = tpm_tfs_melted,
    position = position_jitterdodge(
      jitter.width = 0.2,
      jitter.height = 0,
      dodge.width = 0.7,
      seed = 1
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  xlab("Transcription factor gene") +
  ylab("Gene expression [log(100 * TPM + 1)]") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )
expr_plot_pdf <- paste0(
  plot_outdir,
  "/treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_expression_of_relevant_tfs.pdf"
)
expr_plot_rds <- paste0(
  plot_rds_outdir,
  "/treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_expression_of_relevant_tfs.rds"
)
pdf(expr_plot_pdf,
    width = max(2.5, 0.5 * length(unique(
      tpm_tfs_melted$Gene_symbol
    ))),
    height = 4)
print(expr_plot)
dev.off()
saveRDS(expr_plot, file = expr_plot_rds)
