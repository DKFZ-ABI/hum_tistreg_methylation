# This script analyses demethylation of LTR, SINE and DNA TEs in fat Treg cells.
# Author: Niklas Beumer



# Load required packages.
library(bsseq)
library(GenomicRanges)
library(tidyr)
library(ggplot2)
library(Seurat)
library(Signac)
library(testit)


# Specify a location on /xxx.
b330_space <- "/xxx/"
location <- paste0(b330_space, "nbeumer/hm_treg_bs_rgnsbg")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/yyy/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the BSseq object containing the smoothed methylation data and reduce 
# to the relevant cell types.
meth_data <- readRDS(paste(
  location, 
  "preprocessing_etc",
  "quality_control_after_alignment",
  "2022-03-14_bsseq_object_combined_all_samples_smoothed_w_fat.rds", sep = "/"))
cell_types <- c("Skin Treg", "Fat Treg", "Blood naive Treg")
meth_data_pdata <- pData(meth_data)
meth_data_rel <- meth_data[, meth_data_pdata$Cell_type %in% cell_types]
meth_data_rel_pdata <- pData(meth_data_rel)
rm(meth_data) # Free up memory.

# Specify the colour code for the cell types.
cell_type_col_palette <- c("blue", "green", "darkorchid1")
names(cell_type_col_palette) <- cell_types

# Read in the locations of repeat elements from RepeatMasker
rmsk_file <- paste0(location, 
                    "/external_data/2023-09-14_rmsk_hg19_from_ucsc.txt.gz")
rmsk <- read.table(rmsk_file, header = F, stringsAsFactors = F)

# Filter the repeat element data so that only those elements remain where I am
# confident that they are transposons.
conf_transp_cats <- c("LINE", "DNA", "SINE", "LTR", "RC", "Retroposon")
rmsk_filt <- rmsk[rmsk$V12 %in% conf_transp_cats, ]
rm(rmsk) # Free up memory.

# Generate a GRanges object for the annotated TE insertion sites.
# Take into account that this table originated from UCSC. Intervals are thus 
# very likely 0-based and closed.
# Modify seqlevels styles so that they match with what I used for the signature
# regions.
rmsk_filt_gr <- makeGRangesFromDataFrame(rmsk_filt,
                                         ignore.strand = T,
                                         seqnames.field = "V6",
                                         start.field = "V7",
                                         end.field = "V8",
                                         starts.in.df.are.0based = T,
                                         keep.extra.columns = T)
seqlevelsStyle(rmsk_filt_gr) <- "NCBI"

# Restrict the TEs to chromosomes for which I have methylation data.
seqlevels_to_keep <- intersect(seqlevels(rowRanges(meth_data_rel)),
                               seqlevels(rmsk_filt_gr))
rmsk_filt_2_gr <- rmsk_filt_gr[seqnames(rmsk_filt_gr) %in% seqlevels_to_keep]
# rmsk_filt_2_gr <- rmsk_filt_2_gr[1:100000] # For debugging purposes

# Specify the three TE classes that are relevant for this analysis. These three
# classes are the same for which we show a naive-skin comparison in our w
# manuscript.
relevant_tes <- c("LTR", "DNA", "SINE")

# Make a GRanges object containing locations of the TEs to show.
rmsk_filt_specialtes_gr <- rmsk_filt_2_gr[rmsk_filt_2_gr$V12 %in% relevant_tes]





###############################################################################
# Generate box plots of raw methylation stratified by sample, cell type
# and TE class.
###############################################################################

# Compute raw methylation values for all TE regions. Remove all TE occurrences
# where at least one bisulphite sample returns NA (to make the box plots 
# comparable between the samples). Afterwards, save this data.
rmsk_filt_specialtes <- as.data.frame(rmsk_filt_specialtes_gr)
avg_raw_meth <- getMeth(meth_data_rel, regions = rmsk_filt_specialtes_gr, 
                        type = "raw", what = "perRegion")
rows_to_remove <- which(apply(avg_raw_meth, 1, FUN = function(x) {
  length(which(is.na(x))) > 0
}))
avg_raw_meth_df <- cbind(rmsk_filt_specialtes[-rows_to_remove, ], 
                         avg_raw_meth[-rows_to_remove, ])
avg_raw_meth_df_melted <- pivot_longer(avg_raw_meth_df, 
                                       rownames(meth_data_rel_pdata),
                                       "Sample",
                                       values_to ="Avg_raw_meth")
avg_raw_meth_outfile <- paste0(
  location, 
  "/te_analysis/gen_loc_analysis_dna_sine_ltr_demethylation_in_fat_tregs_methylation_values.txt"
)
write.table(avg_raw_meth_df_melted, file = avg_raw_meth_outfile, 
            sep = "\t", row.names = F, col.names = F, quote = F)

# Generate a box plot showing how strongly the TE regions are methylated in 
# the different cell types (based on raw methylation), stratified by TE class.
avg_raw_meth_df_melted$Cell_type <- 
  sapply(avg_raw_meth_df_melted$Sample, FUN = function(x) {
    strsplit(strsplit(x, split = "(", fixed = T)[[1]][2], split = ";")[[1]][1]
  })
avg_raw_meth_df_melted$Cell_type <- factor(avg_raw_meth_df_melted$Cell_type, 
                                           levels = cell_types)
avg_raw_meth_df_melted$Sample <- factor(
  avg_raw_meth_df_melted$Sample, 
  levels = colnames(meth_data_rel))
avg_raw_meth_df_melted$V12 <- factor(avg_raw_meth_df_melted$V12, 
                                     levels = relevant_tes)
te_box_plot_meth <- ggplot(avg_raw_meth_df_melted) +
  aes(x = Sample, y = Avg_raw_meth, colour = Cell_type) +
  scale_y_continuous(limits = c(0, 1), 
                     name = "Raw methylation (within-region mean)") +
  scale_colour_manual(breaks = cell_types, 
                      values = cell_type_col_palette,
                      name = "Cell type") +
  geom_boxplot() +
  facet_wrap(~V12, ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0),
        panel.spacing.y = unit(0.1, "inches"))
te_box_plot_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_ltr_sine_dna_demethylation_in_fat_tregs.pdf"
)
te_box_plot_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_ltr_sine_dna_demethylation_in_fat_tregs.rds"
)
pdf(te_box_plot_pdf, width = 5, height = 6)
print(te_box_plot_meth)
dev.off()
saveRDS(te_box_plot_meth, file = te_box_plot_rds)




################################################################################
# Perform statistical tests on the differences between cell types.
# Two-tailed Student's t test on the means across the onsidered insertion sites.
################################################################################

# Define the output file that will later contain a summary of the test results.
test_res_summary_file <- paste0(
  location, 
  "/te_analysis/gen_loc_analysis_dna_sine_ltr_demethylation_in_fat_tregs_stat_test_summary.txt"
)

# Start capturing output.
sink(test_res_summary_file)

# Generate a unique ID for each relevant TE insertion site.
avg_raw_meth_df_melted$region_str <- paste(avg_raw_meth_df_melted$seqnames, 
                                           avg_raw_meth_df_melted$start, 
                                           avg_raw_meth_df_melted$end, 
                                           sep = "_")

# Extract values for the different TE classes and the different cell types.
separate_vals <- lapply(relevant_tes, FUN = function(x) {
  restr_table <- avg_raw_meth_df_melted[avg_raw_meth_df_melted$V12 == x, ]
  sublist <- lapply(cell_types, FUN = function(y) {
    restr_table[grepl(y, restr_table$Cell_type, fixed = T), ]
  })
  names(sublist) <- cell_types
  return(sublist)
})
names(separate_vals) <- relevant_tes

# Aggregate values into sample-level values by computing the mean across the 
# considered. TE insertion sites.
separate_vals_aggr <- lapply(separate_vals, FUN = function(x) {
  lapply(x, FUN = function(y) {
    samples <- unique(y$Sample)
    means <- sapply(samples, FUN = function(x) {
      mean(y$Avg_raw_meth[y$Sample == x])
    })
    return(means)
  })
})

# Iterate over the different TE classes and cell type comparisons and perform
# two-tailed Student's t tests on the means accrosss the considered TE insertion
# sites. Print the results to the output.
comparisons <- do.call(c, lapply(1:(length(cell_types) - 1), FUN = function(x) {
  lapply((x + 1):length(cell_types), FUN = function(y) {
    c(cell_types[x], cell_types[y])
  })
}))
test_res <- lapply(relevant_tes, FUN = function(x) {
  cat(paste("\n\n##", x))
  comparisons_res <- lapply(comparisons, FUN = function(y) {
    cat(paste("\n\n#", y[1], "vs.", y[2]))
    single_test_res <- t.test(separate_vals_aggr[[x]][[y[1]]],
                              separate_vals_aggr[[x]][[y[2]]],
                              var.equal = T,
                              alternative = "two.sided")
    print(single_test_res)
    cat(paste0("P = ", single_test_res$p.value, "\n"))
    return(single_test_res)
  })
  names(comparisons_res) <- lapply(comparisons, FUN = function(y) {
    paste(y, collapse = "_vs_")
  })
  pvals <- sapply(comparisons_res, FUN = function(y) {
    y$p.value
  })
  pvals_adj <- p.adjust(pvals, method = "BH")
  cat(paste("\n\n# Benjamini-Hochberg-adjusted P values:\n"))
  print(pvals_adj)
  comparisons_res$Pvals_adj_BH <- pvals_adj
  return(comparisons_res)
})
names(test_res) <- relevant_tes

# Save the test result.
test_outfile <- paste0(
  location, 
  "/te_analysis/gen_loc_analysis_dna_sine_ltr_demethylation_in_fat_tregs_stat_tests.rds"
)
saveRDS(test_res, file = test_outfile)

# Stop capturing output.
sink()
