# This script optimises some figures for our publication.
# Author: Niklas Beumer



# Load required package(s).
library(ggplot2)
library(testit)


# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")





###### Correlation plot for DMR-peak-gene links.
plot_file <- paste0(
  location, 
  "/plot_rds_objects/treg_hierarchies_diff_meth_diff_acc_diff_gex_dmr_peak_gene_links_meth_acc_gex_corr.rds"
)
plot <- readRDS(plot_file)
plot_2 <- plot + 
  scale_colour_gradient2(
    low = "blue", 
    high = "red",
    name = "Expression\nlog2(FC)\n(Blood naive\nTregs vs.\nSkin Tregs)")
outfile_pdf <- paste0(
  plot_outdir, 
  "/plots_optim_f_paper_treg_hierarchies_diff_meth_diff_acc_diff_gex_dmr_peak_gene_links_meth_acc_gex_corr.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir, 
  "/plots_optim_f_paper_treg_hierarchies_diff_meth_diff_acc_diff_gex_dmr_peak_gene_links_meth_acc_gex_corr.rds"
)
pdf(outfile_pdf, width = 5.5, height = 4)
print(plot_2)
dev.off()
saveRDS(plot_2, file = outfile_rds)


###### Homer results plot for differential analysis between skin Tregs and 
###### blood naive Tregs.
plot_file <- paste0(
  location, 
  "/plot_rds_objects/treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_dot_plot_selected_motifs.rds"
)
plot <- readRDS(plot_file)
plot_data <- plot$data
plot_data_2 <- plot_data
plot_data_2$greater1 <- sapply(
  plot_data_2$Enrichment_score, FUN = function(x) {
    ifelse(x > 1, yes = "enriched", no = "depleted")
  }) 
plot_data_3 <- plot_data_2[grep("bZIP", plot_data_2$Motif.Name), ]
plot_3 <- plot
plot_3$data <- plot_data_3
plot_3 <- plot_3 + 
  aes(fill = greater1) +
  scale_fill_manual(name = "Enrichment",
                    breaks = c("depleted", "enriched"),
                    values = c("blue", "red"),
                    na.value = "black")
plot_data_4 <- plot_data_2[plot_data_2$Motif.Name %in% c(
  "c-Myc(bHLH)/LNCAP-cMyc-ChIP-Seq(Unpublished)/Homer",
  "c-Myc(bHLH)/mES-cMyc-ChIP-Seq(GSE11431)/Homer",
  "USF1(bHLH)/GM12878-Usf1-ChIP-Seq(GSE32465)/Homer",
  "n-Myc(bHLH)/mES-nMyc-ChIP-Seq(GSE11431)/Homer",
  "HIF2a(bHLH)/785_O-HIF2a-ChIP-Seq(GSE34871)/Homer"
), ]
plot_4 <- plot
plot_4$data <- plot_data_4
plot_4 <- plot_4 + 
  aes(fill = greater1) +
  scale_fill_manual(name = "Enrichment",
                    breaks = c("depleted", "enriched"),
                    values = c("blue", "red"),
                    na.value = "black")
plot_3_pdf <- paste0(
  plot_outdir, 
  "/plots_optim_f_paper_treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_dot_plot_selected_motifs_bzip.pdf"
)
plot_3_rds <- paste0(
  plot_rds_outdir, 
  "/plots_optim_f_paper_treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_dot_plot_selected_motifs_bzip.pdf"
)
pdf(plot_3_pdf, width = 6.5, height = 6)
print(plot_3)
dev.off()
saveRDS(plot_3, file = plot_3_rds)
plot_4_pdf <- paste0(
  plot_outdir, 
  "/plots_optim_f_paper_treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_dot_plot_selected_motifs_bhlh.pdf"
)
plot_4_rds <- paste0(
  plot_rds_outdir, 
  "/plots_optim_f_paper_treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_dot_plot_selected_motifs_bhlh.pdf"
)
pdf(plot_4_pdf, width = 6.5, height = 6)
print(plot_4)
dev.off()
saveRDS(plot_4, file = plot_4_rds)


###### TPM bar chart for transcription factors that were highlighted by diff.
###### analysis betw. skin Tregs and blood naive Tregs.
plot_file <- paste0(
  plot_rds_outdir, 
  "/treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_expression_of_relevant_tfs.rds"
)
plot <- readRDS(plot_file)
plot_data_bar <- plot$layers[[1]]$data
plot_data_errorbar <- plot$layers[[2]]$data
assert(all(sapply(1:ncol(plot_data_bar), FUN = function(x) {
  all(plot_data_bar[, x] == plot_data_errorbar[, x])
})))
plot_data_points <- plot$layers[[3]]$data
genes_to_show <- c("MYC", "USF1","MYCN", "EPAS1", "JUN", "FOSL2", "JUNB", 
                   "FOSL1", "ATF3", "BATF", "BACH2", "BACH1")
plot_data_bar_2 <- plot_data_bar[plot_data_bar$Gene_symbol %in% 
                                   genes_to_show, ]
plot_data_points_2 <- plot_data_points[plot_data_points$Gene_symbol %in% 
                                         genes_to_show, ]
plot_2 <- plot
plot_2$layers[[1]]$data <- plot_data_bar_2
plot_2$layers[[2]]$data <- plot_data_bar_2
plot_2$layers[[3]]$data <- plot_data_points_2
outfile_pdf <- paste0(
  plot_outdir, 
  "/plots_optim_f_paper_treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_expression_of_relevant_tfs.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir, 
  "/plots_optim_f_paper_treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_expression_of_relevant_tfs.rds"
)
pdf(outfile_pdf, width = 6, height = 4)
print(plot_2)
dev.off()
saveRDS(plot_2, file = outfile_rds)