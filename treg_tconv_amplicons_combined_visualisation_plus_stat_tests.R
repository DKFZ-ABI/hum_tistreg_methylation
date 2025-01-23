# This script generates combined visualisations for the amplicon results
#   from the comparison between blood naive Treg cells and blood naive Tconv 
#   cells. In addition, it performs statistical tests for these results.
# Author: Niklas Beumer



# Load required packages.
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(testit)


# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Specify the names of the amplicon experiments that are combined here.
experiment_names <- c("05_2024_pilot_run", "07_2024_further_donors_Treg_Tconv")

# Specify the location of data output.
output_loc <- paste0(location, 
                     "/amplicon_analysis/analysis_results/combination_", 
                     paste(experiment_names, collapse = "_"))

# Read in the separate data for the experiments.
data_paths <- paste0(location, "/amplicon_analysis/analysis_results/", 
                     experiment_names, "/relevant_methylation_results.rds")
meth_data <- lapply(data_paths, FUN = function(x) {
  temp <- read.csv(x)
  temp_filtered <- temp[, 2:ncol(temp)]
})
names(meth_data) <- experiment_names

# Specify the sampls that should be visualised. Note that (for the pilot run)
# I will only include Tamara10 because this is compatible with the experimental 
# strategy used for the other donors (according to Tamara).
samples_to_keep <- c("Tamara10", paste0("Donor", c("A", "B", "C", "D1", "D2")))

# Specify the exact locations of the amplicons.
amplicon_names <- c("R8225_TIGIT", "R8226_TIGIT", "R1198_IL2RA", 
                    "R10871#2_LRRN3", "R10425_SYTL3", "R115_TNFRSF8", 
                    "R12143-1_FOXP3", "FOXP3TSDR_FOXP3", "R12143-2_FOXP3", 
                    "R6632_CTLA4", "R132_TNFRSF1B", "R122_TNFRSF1B")
amplicon_locs <- GRanges(seqnames = c("3", "3", "10", "7", "6", "1", "X", "X", 
                                      "X", "2", "1", "1"), 
                         ranges = IRanges(start = c(114012711, 114014310, 
                                                    6102137, 110731169, 
                                                    159083980, 12188066, 
                                                    49119989, 49117054, 
                                                    49120347, 204736407, 
                                                    12263695, 12231169), 
                                          end = c(114013071, 114014730, 
                                                  6102543, 110731507, 
                                                  159084279, 12188381, 
                                                  49120374, 49117406, 
                                                  49120704, 204736768, 
                                                  12264119, 12231538)))
amplicon_locs$Amplicon_name <- amplicon_names






##############################################################################
# Visualise methylation levels in the combined data.
##############################################################################

# For each sample, extract relevant information.
rel_info_by_sample <- do.call(rbind, lapply(samples_to_keep, FUN = function(x) {
  relevant_table <- meth_data[[which(sapply(meth_data, FUN = function(y) {
    !(all(!(grepl(x, colnames(y)))))
  }))]]
  rel_cols <- paste0(x, ".", c("chr", "start", "Cell_type", "Methylation_beta", 
                               "Coverage", "Amplicon"))
  rel_df <- relevant_table[, rel_cols]
  colnames(rel_df) <- gsub(paste0(x, "."), "", colnames(rel_df))
  rel_df$y_axis_label <- paste0(rel_df$Cell_type, ": ", x)
  return(rel_df)
}))

# Save the collected data.
coll_data_outfile <- paste0(output_loc, "/methylation_results_combined.txt")
write.table(rel_info_by_sample, file = coll_data_outfile, sep = "\t", 
            col.names = T, row.names = F, quote = F)

# Compute log-transformed coverage values.
rel_info_by_sample$logcov <- log10(rel_info_by_sample$Coverage + 1)

# Visualise methylation data for all amplicons.
all_amplicons <- unique(rel_info_by_sample$Amplicon)
combined_meth_plots <- lapply(all_amplicons, FUN = function(x) {
  plotting_data <- rel_info_by_sample[rel_info_by_sample$Amplicon == x, ]
  plotting_data$y_axis_label <- factor(
    plotting_data$y_axis_label,
    levels = rev(sort(unique(plotting_data$y_axis_label)))
  )
  amplicon_plot <- ggplot(plotting_data) +
    aes(x = start, y = y_axis_label, fill = Methylation_beta, size = logcov) +
    scale_x_continuous(
      name = paste0("Position (Chr. ", unique(plotting_data$chr), ")")) +
    scale_fill_viridis_c(limits = c(0, 1), name = "Methylation", 
                         direction = -1) +
    scale_size_area(limits = c(0, 5), name = "log10(\nCoverage + 1)") +
    geom_hline(yintercept = unique(plotting_data$y_axis_label), 
               colour = "black") +
    geom_point(shape = 21, colour = "black", stroke = 1.5) +
    coord_cartesian(
      xlim = c(start(amplicon_locs)[amplicon_locs$Amplicon_name == x],
               end(amplicon_locs)[amplicon_locs$Amplicon_name == x]),
      ylim = c(0, length(unique(plotting_data$y_axis_label)) + 1),
      expand = F,
      default = T,
      clip = "off") +
    ylab("Cell type") +
    ggtitle(paste0("Amplicon: ", gsub("_", " (", x), ")")) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks.x = element_line(colour = "black"),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank())
  return(amplicon_plot)
})






###########################################################################
# Perform Wilcoxon signed rank tests assessing differences between
# the analysed cell types. Visualise test results.
###########################################################################

# Initialise an empty data frame that will later contain the relevant 
# values for plotting.
test_plotting_df <- data.frame(Amplicon = NA, chr = NA, start = NA, Pval = NA, 
                               Pval_adj_BH = NA)

# Initialise an empty list that will later hold all test results.
all_test_res <- list()

# Perform two-tailed Wilcoxon signed rank tests for each amplicon and each CpG 
# site therein. Within each amplicon, adjust these using the Benjamini-Hochberg
# method.
for (x in all_amplicons) {
  temp_df_for_BHcorr <- data.frame(Amplicon = NA, chr = NA, start = NA, 
                                   Pval = NA)
  res_this_amplicon <- list()
  rel_data_amplicon <- rel_info_by_sample[rel_info_by_sample$Amplicon  == x, ]
  all_positions <- unique(rel_data_amplicon$start)
  for (y in all_positions) {
    rel_data_position <- rel_data_amplicon[rel_data_amplicon$start == y, ]
    cell_types <- unique(rel_data_position$Cell_type)
    assert(length(cell_types) == 2)
    df_1 <- rel_data_position[rel_data_position$Cell_type == cell_types[1], ]
    df_2 <- rel_data_position[rel_data_position$Cell_type == cell_types[2], ]
    assert(all(gsub(cell_types[1], "", df_1$y_axis_label) == 
                 gsub(cell_types[2], "", df_2$y_axis_label)))
    wilcox_res <- wilcox.test(df_1$Methylation_beta, df_2$Methylation_beta,
                              alternative = "two.sided", paired = T)
    temp_df <- data.frame(Amplicon = x,
                          chr = unique(rel_data_position$chr), 
                          start = y, 
                          Pval = wilcox_res$p.value)
    temp_df_for_BHcorr <- rbind(temp_df_for_BHcorr, temp_df)
    res_this_amplicon[[paste0("Position_", y)]] <- wilcox_res
  }
  temp_df_for_BHcorr <- temp_df_for_BHcorr[-1, ] # Remove placeholder line.
  temp_df_for_BHcorr$Pval_adj_BH <- p.adjust(temp_df_for_BHcorr$Pval, 
                                             method = "BH")
  test_plotting_df <- rbind(test_plotting_df, temp_df_for_BHcorr)
  all_test_res[[x]] <- res_this_amplicon
}
test_plotting_df <- test_plotting_df[-1, ] # Remove placeholder line.

# Save test outcomes.
detailed_res_outfile <- paste0(output_loc, "/all_wilcox_test_res.rds")
saveRDS(all_test_res, file = detailed_res_outfile)
summary_outfile <- paste0(output_loc, "/wilcox_test_res_summary_w_BH_corr.txt")
write.table(test_plotting_df, file = summary_outfile, sep = "\t", col.names = T,
            row.names = F, quote = F)

# Visualise the Benjamini-Hochberg-corrected P values separately for each
# amplicon.
test_plotting_df$logp <- -log10(test_plotting_df$Pval_adj_BH)
test_plotting_df$significant <- factor(test_plotting_df$Pval_adj_BH < 0.05, 
                                       levels = c(T, F))
combined_signif_plots <- lapply(all_amplicons, FUN = function(x) {
  plotting_data <- test_plotting_df[test_plotting_df$Amplicon == x, ]
  amplicon_plot <- ggplot(plotting_data) +
    aes(x = start, y = "Significance", size = logp, colour = significant) +
    scale_x_continuous(
      name = paste0("Position (Chr. ", unique(plotting_data$chr), ")")) +
    scale_size_area(limits = c(0, 2), name = "-log10(q)") +
    scale_colour_manual(breaks = c(T, F), 
                        values = c("cyan", "black"), 
                        labels = c("q < 0.05", "q >= 0.05"), 
                        name = "Significance",
                        drop = F) +
    geom_hline(yintercept = "Significance", colour = "black") +
    geom_point(shape = 24, fill = "grey", stroke = 1.5) +
    coord_cartesian(
      xlim = c(start(amplicon_locs)[amplicon_locs$Amplicon_name == x],
               end(amplicon_locs)[amplicon_locs$Amplicon_name == x]),
      ylim = c(0, 2),
      expand = F,
      default = T,
      clip = "off") +
    ylab("") +
    ggtitle(paste0("Amplicon: ", gsub("_", " (", x), ")")) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks.x = element_line(colour = "black"),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank())
  return(amplicon_plot)
})







#############################################################################
# Generate plots that combine methylation data with significance data.
#############################################################################

# For each amplicon, combine a methylation plot with the respective significance 
# plot.
final_plots <- lapply(1:length(combined_meth_plots), FUN = function(x) {
  plot_grid(combined_meth_plots[[x]], combined_signif_plots[[x]],
            ncol = 1,
            axis = "lr",
            align = "v",
            rel_heights = c(4, 2))
})

# Save the combined plots.
outfile_snip <- "/amplicons_Treg_Tconv_results_all_combined_results"
amplicon_plot_pdf <- paste0(plot_outdir, outfile_snip, ".pdf")
amplicon_plot_rds <- paste0(plot_rds_outdir, outfile_snip, ".rds")
pdf(amplicon_plot_pdf, width = 8, height = 8)
print(final_plots)
dev.off()
object_f_rds <- lapply(1:length(combined_meth_plots), FUN = function(x) {
  list(meth_plot = combined_meth_plots[[x]], 
       signif_plot = combined_signif_plots[[x]])
})
names(object_f_rds) <- amplicon_names
saveRDS(object_f_rds, file = amplicon_plot_rds)
