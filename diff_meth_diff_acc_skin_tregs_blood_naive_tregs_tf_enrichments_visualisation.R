# This script visualises TF enrichment results (known motifs) for the comparison
#   between skin Tregs cells and blood naive Treg cells (both on the methylation
#   level and the chromatin accessibility level; Selected motifs).
# Author: Niklas Beumer



# load required package(s).
library(ggplot2)

# Define a location on /yyy.
location <- "yyy/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis",
                     format(Sys.time(), "%m-%d-%y"),
                     sep = "/")
if (!dir.exists(plot_outdir)) {
  dir.create(plot_outdir)
}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the results for known motifs from Homer.
homer_snips <- c(
  "diff_meth_skin_treg_blood_naive_treg_Blood_naive_Treg__hypomethylation",
  "diff_meth_skin_treg_blood_naive_treg_Skin_Treg__hypomethylation",
  "diff_acc_skin_treg_blood_naive_treg_Blood_naive_Treg__hyperaccessibility",
  "diff_acc_skin_treg_blood_naive_treg_Skin_Treg__hyperaccessibility"
)
homer_results <- lapply(
  homer_snips,
  FUN = function(x) {
    homer_file <- paste0(location,
                         "/treg_hierarchies/",
                         x,
                         "_homer_against_genome/knownResults.txt")
    read.table(
      homer_file,
      sep = "\t",
      header = T,
      stringsAsFactors = F,
      comment.char = ""
    )
  }
)
names(homer_results) <- c(
  "Blood_naive_Treg__hypomethylation",
  "Skin_Treg__hypomethylation",
  "Blood_naive_Treg__hyperaccessibility",
  "Skin_Treg__hyperaccessibility"
)










#############################################################################
# Compute enrichment scores.
# An enrichment score is the percentage of target sequences with motif
# divided by the percentage of background sequences with motif.
#############################################################################

# Compute enrichment scores and combine the results for the region classes.
#-- Iterate over homer's results.
homer_results_combined <- do.call(rbind, lapply(
  names(homer_results),
  FUN = function(x) {
    curr_df <- homer_results[[x]]
    #-- Compute enrichment scores.
    target_perc <- as.numeric(gsub("%", "", curr_df$X..of.Target.Sequences.with.Motif))
    background_perc <- as.numeric(gsub("%", "", curr_df$X..of.Background.Sequences.with.Motif))
    curr_df$Enrichment_score <- target_perc / background_perc
    #-- Scale the enrichment scores so that they range from 0 to 1.
    min_score <- min(curr_df$Enrichment_score)
    max_score <- max(curr_df$Enrichment_score)
    curr_df$Enrichment_score_scaled <- (curr_df$Enrichment_score - min_score) / (max_score - min_score)
    #-- Make column names consistent between the different results tables.
    colnames(curr_df)[6] <- "Number..of.Target.Sequences.with.Motif"
    colnames(curr_df)[8] <- "Number..of.Background.Sequences.with.Motif"
    #-- Add information on the region class.
    curr_df$Region_class <- x
    #-- Return value
    return(curr_df)
  }
))

# Save the combined results.
homer_results_combined_outfile <- paste0(
  location,
  "/treg_hierarchies/diff_meth_diff_acc_skin_treg_blood_naive_treg_homer_against_genome_combined_results.txt"
)
write.table(
  homer_results_combined,
  file = homer_results_combined_outfile,
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)







##############################################################################
# Visualise the results.
##############################################################################

# Specify the transcription factor motifs that should be shown in the plot.
tfs_to_show <- c(
  # Selection of bHLH factors
  "c-Myc(bHLH)/LNCAP-cMyc-ChIP-Seq(Unpublished)/Homer",
  "c-Myc(bHLH)/mES-cMyc-ChIP-Seq(GSE11431)/Homer",
  "USF1(bHLH)/GM12878-Usf1-ChIP-Seq(GSE32465)/Homer",
  "n-Myc(bHLH)/mES-nMyc-ChIP-Seq(GSE11431)/Homer",
  "HIF2a(bHLH)/785_O-HIF2a-ChIP-Seq(GSE34871)/Homer",
  # Selection of bZIP factors.
  "Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer",
  "Fosl2(bZIP)/3T3L1-Fosl2-ChIP-Seq(GSE56872)/Homer",
  "Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer",
  "JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer",
  "Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer",
  "Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer",
  "BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer",
  "Bach2(bZIP)/OCILy7-Bach2-ChIP-Seq(GSE44420)/Homer",
  "AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer",
  "Bach1(bZIP)/K562-Bach1-ChIP-Seq(GSE31477)/Homer"
)

# Collect values for these transcription factors.
relevant_values <- homer_results_combined[
  homer_results_combined$Motif.Name %in% tfs_to_show, 
]

# Clip Benjamini-Hochberg-adjusted P values.
relevant_values$qval_clipped <- sapply(
  relevant_values$q.value..Benjamini.,
  FUN = function(x) {
    max(x, 0.001)
  }
)

# Compute log-transformed q values.
relevant_values$logq <- -log10(relevant_values$qval_clipped) + 0.75

# Add a column indicating statistical significance.
relevant_values$significant <- relevant_values$q.value..Benjamini. < 0.05

# Generate a dot plot showing homer results.
relevant_values$Motif.Name <- factor(relevant_values$Motif.Name, levels = rev(tfs_to_show))
relevant_values$Region_class <- factor(relevant_values$Region_class, levels = names(homer_results))
dot_plot <- ggplot(relevant_values) +
  aes(
    x = Region_class,
    y = Motif.Name,
    size = logq,
    fill = Enrichment_score_scaled,
    colour = significant
  ) +
  scale_size_area(
    breaks = 0.75:3.75,
    labels = c("0", "1", "2", ">=3"),
    name = "-log10(q)"
  ) +
  scale_fill_viridis_c(
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c("min. f. class", "max. f. class"),
    name = "Enrichment",
    option = "C"
  ) +
  scale_colour_manual(
    breaks = c(T, F),
    values = c("cyan", "grey"),
    labels = c("q < 0.05", "q >= 0.05"),
    name = "Significance"
  ) +
  geom_point(shape = 21, stroke = 1) +
  xlab("Region class") +
  ylab("Motif") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      colour = "black",
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

# Save the plot.
dot_plot_pdf <- paste0(
  plot_outdir,
  "/treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_dot_plot_selected_motifs.pdf"
)
dot_plot_rds <- paste0(
  plot_rds_outdir,
  "/treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_dot_plot_selected_motifs.rds"
)
pdf(dot_plot_pdf, width = 8, height = 9.5)
print(dot_plot)
dev.off()
saveRDS(dot_plot, file = dot_plot_rds)
