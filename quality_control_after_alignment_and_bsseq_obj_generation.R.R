# This script performs some basic quality control of the TWGBS data after they 
#   have been aligned. It also generates a BSseq object. For the sample Green-8, 
#   only library 1 is considered because the other library was shown to display 
#   outlier behaviour.
# Author: Niklas Beumer



# Load all required packages.
library(ggplot2)
library(testit)
library(grid)
library(gridExtra)
library(reshape2)

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/home/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the location where the sequencing results can be found.
seq_results_loc <- 
  "/yyy/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid"

# Specify a location in my space.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Specify the prefix for output text and RDS files that will later contain the 
# processed data.
data_out_pref <- paste0(
  location, 
  "/preprocessing_etc/quality_control_after_alignment/", 
  format(Sys.time(), "%Y-%m-%d")
)

# Read in the sample mapping file.
sample_mapping <- read.table(paste0(location, "/sample_mapping_twgbs.txt"), 
                             header = T, sep = "\t", stringsAsFactors = F)

# Extract the list of analysed cell types.
cell_types <- unique(sample_mapping$Cell_type)

# Specify the directories where quality control statistics can be found.
qual_dirs <- paste0(seq_results_loc,
                    "/OE0436_Treg_cells_",
                    sample_mapping$Patient_ID,
                    "/",
                    sample_mapping$Sample_ID,
                    "/paired/merged-alignment")

# Specify sample names that will later appear in plots.
sample_names_for_plots <- paste0(sample_mapping$Sample_ID_readable,
                                 "\n",
                                 "(",
                                 sample_mapping$Cell_type,
                                 "; ",
                                 sample_mapping$Patient_ID_readable,
                                 ")")
sample_names_for_plots_2 <- sample_mapping$Sample_ID_readable

# Generate colour scales for the samples and cell types.
set.seed(1996)
sample_col_palette <- sample(rainbow(length(sample_names_for_plots)), 
                             length(sample_names_for_plots))
cell_type_col_palette <- c("blue", "cyan2", "orange", "darkorchid1", "red")






############################
# Insert size distribution.
############################

# Specify files containing information on the insert size distribution.
# According to the documentation 
# (https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows/wiki/3.-Results; 
# 06-05-2021), the insertsizes.txt contain insert sizes in the first column and 
# counts in the second column. The insertsize_plot.png_qcValues.txt files 
# contain the median insert size in the first row, the standard deviation 
# divided by the median in the second row and the standard deviation in the 
# third row.
insert_size_prefixes <- paste0(
  qual_dirs, 
  "/qualitycontrol/",
  ifelse(sample_mapping$Sample_ID == "green-8", yes ="lib1", no = "merged"),
  "/insertsize_distribution/", 
  sample_mapping$Sample_ID, 
  "_OE0436_Treg_cells_", 
  sample_mapping$Patient_ID
)
insert_size_files <- paste0(
  insert_size_prefixes, 
  ifelse(sample_mapping$Sample_ID == "green-8", yes ="_lib1", no = ""),
  "_insertsizes.txt"
)
insert_size_metrics_files <- paste0(
  insert_size_prefixes,
  ifelse(sample_mapping$Sample_ID == "green-8", yes ="_lib1", no = ""),
  "_insertsize_plot.png_qcValues.txt"
)

# Read in this data.
sample_repetitions <- c()
cell_type_repetitions <- c()
insert_sizes <- c()
counts <- c()
medians <- c()
for (i in 1:length(insert_size_prefixes)) {
  dist_table <- read.table(insert_size_files[i], header = F, 
                           stringsAsFactors = F)
  sample_repetitions <- c(sample_repetitions, 
                          rep(sample_names_for_plots[i], nrow(dist_table)))
  cell_type_repetitions <- c(cell_type_repetitions, 
                             rep(sample_mapping$Cell_type[i], nrow(dist_table)))
  insert_sizes <- c(insert_sizes, dist_table$V1)
  counts <- c(counts, dist_table$V2)
  metrics_table <- read.table(insert_size_metrics_files[i], header = F, 
                              stringsAsFactors = F, sep = "\t")
  medians <- c(medians, metrics_table[1, 1])
}
distributions <- data.frame(Sample = sample_repetitions,
                            Cell_type = cell_type_repetitions,
                            Insert_size = insert_sizes,
                            counts = counts)
distributions$Sample <- factor(distributions$Sample, 
                               levels = sample_names_for_plots)
distributions$Cell_type <- factor(distributions$Cell_type, levels = cell_types)
medians_df <- data.frame(Sample = sample_names_for_plots,
                         Median = paste("Median:", medians),
                         Insert_size = (0.8 * max(distributions$Insert_size)),
                         counts = (0.8 * max(distributions$counts)))
medians_df$Sample <- factor(medians_df$Sample, levels = sample_names_for_plots)

# Save the data on insert size distribution.
insert_size_outfile <- paste0(data_out_pref, 
                              "_insert_size_distribution_by_sample.txt")
write.table(distributions, file = insert_size_outfile, sep = "\t", 
            col.names = T, row.names = F, quote = F)

# Plot the insert size distributions and save this plot.
insert_size_plot <- ggplot(distributions) +
  aes(x = Insert_size, y = counts) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  geom_bar(stat = "identity", fill = "black", width = 1) +
  facet_wrap(facets = ~Sample, ncol = 3) +
  geom_text(data = medians_df, mapping = aes(label = Median)) +
  xlab("Insert size") +
  ylab("Count") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))
plot_width <- 7.5
plot_height <- 1.25 * ceiling(length(sample_names_for_plots) / 3)
pdf(paste0(plot_outdir, "/insert_size_distribution_by_sample.pdf"), 
    width = plot_width, height = plot_height)
print(insert_size_plot)
dev.off()
saveRDS(insert_size_plot, 
        file = paste0(plot_rds_outdir, 
                      "/insert_size_distribution_by_sample.rds"))

# Generate and save an alternative plot that shows the distribution as a smooth 
# curve.
insert_size_plot_2 <- ggplot(distributions) +
  aes(x = Insert_size, col = Sample, y = counts) +
  scale_colour_manual(breaks = sample_names_for_plots, 
                      values = sample_col_palette, name = "Sample") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), name = "Count") +
  geom_line() +
  xlab("Insert size") +
  facet_wrap(~Cell_type, ncol = 3) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.key.height = unit(0.3, "inches"))
plot_width <- 10
plot_height <- 1.5 * ceiling(length(cell_types) / 3)
pdf(paste0(plot_outdir, "/insert_size_distribution_by_sample_with_line.pdf"), 
    width = plot_width, height = plot_height)
print(insert_size_plot_2)
dev.off()
saveRDS(insert_size_plot_2, 
        file = paste0(plot_rds_outdir, 
                      "/insert_size_distribution_by_sample_with_line.rds"))










######################################################################
# M-Bias: Methylation levels as a function of the position in a read.
######################################################################

# Specify the directories that contain information on how methylation levels 
# change over different read positions. Unfortunately, there is no 
# ducumentation. But I suppose that the section "methylation vs. position" 
# contains the data I'm looking for. The column "CG.mC" should contain the 
# number of methylated Cs at a specific position and the column "CG.C" should 
# contain the number of un-methylated Cs.
mbias_dirs <- paste0(
  qual_dirs, 
  "/methylation/",
  ifelse(sample_mapping$Sample_ID == "green-8", yes ="lib1", no = "merged"),
  "/methylationCallingMetrics"
)

# Iterate over all samples to read in the data for CpG sites.
mbias_data <- data.frame(
  mate = NA, pos = NA, CG.mC = NA, CG.C = NA, Sample = NA
) # Placeholder line.
for (i in 1:length(mbias_dirs)) {
  
  # Iterate over all chromosomes to add up the information.
  temp_df <- data.frame(mate = c(rep("1", 150), rep("2", 150)),
                        pos = rep(1:150, 2), CG.mC = 0, CG.C = 0)
  files_in_dir <- list.files(mbias_dirs[[i]])
  for (file in files_in_dir) {
    data <- read.table(paste0(mbias_dirs[i], "/", file), header = T, 
                       stringsAsFactors = F, skip = 6, nrows = 300)
    if ("mate" %in% data$mate) {
      # Not all 150 positions present in file: Strip the unwanted lines.
      start_wrong_lines <- match("mate", data$mate)
      data <- data[-(start_wrong_lines:nrow(data)), ]
    }
    data$pos <- as.numeric(data$pos)
    for (j in 1:nrow(data)) {
      matching_line <- which(temp_df$mate == data$mate[j] & 
                               temp_df$pos == data$pos[j])
      temp_df$CG.mC[matching_line] <- temp_df$CG.mC[matching_line] + 
        as.numeric(data$CG.mC[j])
      temp_df$CG.C[matching_line] <- temp_df$CG.C[matching_line] + 
        as.numeric(data$CG.C[j])
    }
  }
  
  # Add sample information and merge the data.
  temp_df$Sample <- sample_names_for_plots[i]
  mbias_data <- rbind(mbias_data, temp_df)
}

# Remove the placeholder line, which contains NAs.
mbias_data <- mbias_data[-1, ]

# Compute the total number of CpG sites at each read position and the 
# methylation level at each position.
mbias_data$Total <- mbias_data$CG.mC + mbias_data$CG.C
mbias_data$Methyl_level <- mbias_data$CG.mC / mbias_data$Total

# Save the data on B-Bias.
mbias_outfile <- paste0(data_out_pref, "_m_bias_by_sample.txt")
write.table(mbias_data, file = mbias_outfile, sep = "\t", col.names = T, 
            row.names = F, quote = F)

# Generate an M-Bias plot.
mbias_data$Sample <- factor(mbias_data$Sample, levels = sample_names_for_plots)
mbias_plot <- ggplot(mbias_data) +
  aes(x = pos, y = Methyl_level, colour = mate) +
  scale_x_continuous(expand = c(0, 0), name = "Position in read") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05)), 
                     name = "Methylated fraction") +
  scale_colour_manual(breaks = c("1", "2"), values = c("blue", "orange"), 
                      name = "Read", labels = c("Read 1", "Read 2")) +
  geom_line(alpha = 0.6, size = 1.2) +
  geom_vline(xintercept = 10, linetype = "dashed") +
  facet_wrap(~Sample, ncol = 3) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom")

# Save the plot.
plot_width <- 7.5
plot_height <- 1.25 * ceiling(length(sample_names_for_plots) / 3)
pdf(paste0(plot_outdir, "/m_bias_plots_by_sample.pdf"), width = plot_width, 
    height = plot_height)
print(mbias_plot)
dev.off()
saveRDS(mbias_plot, 
        file = paste0(plot_rds_outdir, "/m_bias_plots_by_sample.rds"))









#########################################################################
# Coverage distribution for CpGs and BSSeq objects for further analysis.
#########################################################################

extract_coverage_distribution_from_df <- function(df, max_coverage, 
                                                  coverage_column) {
  # This function computes the coverage distribution for CpG sites from a data 
  #   frame stating the coverage for each single CpG site.
  # df: The data frame containing each CpG site in an individual row.
  # max_coverage: Numeric: The maximum coverage per CpG to expect.
  # coverage column: Character: The name of the column in df that contains the 
  #   coverage values.
  # Dependencies: None
  # Value: A data frame with the columns "cov", "count", "rel_count" and 
  #   "rel_count_cum" containing the distribution.
  # Author: Niklas Beumer
  
  # Initialise the data frame.
  temp_df <- data.frame(cov = 0:max_coverage, count = 0)
  
  # Check whether the maximum expected coverage is not exceeded.
  if (max(df[, coverage_column]) > max_coverage) {
    stop(paste("Coverages above", max_coverage, 
               "were observed. Use a higher maximum coverage."))
  }
  
  # Compute the coverage distribution.
  coverage_table <- table(df[, coverage_column])
  for (coverage_value in names(coverage_table)) {
    matching_row <- which(temp_df$cov == coverage_value)
    temp_df[matching_row, "count"] <- coverage_table[coverage_value]
  }
  
  # Compute the relative amount of CpG sites that are covered at each coverage 
  # value.
  total_cpgs <- sum(temp_df$count)
  temp_df$rel_count <- temp_df$count / total_cpgs
  
  # Compute the cumulative coverage distribution, i.e. the proportion of CpG 
  # sites covered at a particular coverage or at greater coverage.
  temp_df$rel_count_cum <- sapply(1:nrow(temp_df), FUN = function(x) {
    sum(temp_df[x:nrow(temp_df), ]$rel_count)
  })
  
  # Return the coverage distribution.
  return(temp_df)
  
}


plot_coverage_distribution_by_metavar <- function(
    dist_df, cumulative = F, metavar_colname, metavar_order, curve = F,
    metavar_2_colname = NA, metavar_2_order = NA, metavar_2_palette = NA,
    colour_name = NA, cov_to_highlight = 5, x_limit, x_name, y_name, 
    outfile_pdf, outfile_rds) {
  # This function plots the CpG coverage distribution based on a distribution 
  #   data frame.
  # dist_df: A data frame containing the coverage distribution.
  #   Must contain the columns "cov", >>metavar_colname<<, "count" and 
  #   "rel_count" or "rel_count_cum".
  # cumulative: Logical. Whether to plot the cumulative or the non-cumulative 
  #   distribution.
  # metavar_colname: Character: The name of the column according to which the 
  #   data will be split into facets.
  # metavar_order: Character; A vector specifying the order in which the 
  #   different facets will appear in the plot.
  # curve: Logical. If set to TRUE and cumulative = F, the distribution will be 
  #   shown as a curve, not as a bar chart.
  # metavar_2_colname: character; If cumulative = FALSE and curve = TRUE, this 
  #   specifies the name of a column in dist_df that is used to colour the 
  #   curves.
  # metavar_2_order: Character: If cumulative = FALSE and curve = TRUE, this 
  #   gives the order in ehich the levels in metavar_2_colname will be displayed 
  #   in the legend.
  # metavar_2_palette: Character: If cumulative = FALSE and curve = TRUE, this 
  #   vector gives the colours to use for metavar_2_colname.
  # colour_name: Character: If cumulative = FALSE and curve = TRUE, the name of 
  #   the colour legend.
  # cov_to_highlight: Numeric; If cumulative = T, the coverage value to 
  #   specifically highlight in the plot.
  # x_limit: Numeric; The maximum coverage to show on the x axis.
  # x_name: Character; The x axis label to use.
  # y_name: Character; The y axis label to use.
  # outfile_pdf: Character: Path to a PDF file where the plot will be saved.
  # outfile_rds: Character: Path to an RDS file where the plot will be saved.
  # Dependency: ggplot2
  # Value: None, the function just generates and saves the plot.
  # Author: Niklas Beumer
  
  # Collect the plotting data.
  plotting_data <- dist_df[
    , c("cov", metavar_colname, 
        ifelse(cumulative, yes = "rel_count_cum", no = "rel_count"))
  ]
  colnames(plotting_data) <- c("cov", "Metavar", "dist")
  if (!cumulative & curve) {
    plotting_data$Metavar_2 <- factor(dist_df[, metavar_2_colname], 
                                      levels = metavar_2_order)
  }
  
  # If the distribution is not cumulative and a bar chart wil be generated, 
  # process the plotting data so that all coverage values above x_limit are 
  # aggregated into a single category.
  # If the plot should show a curve, remove these coverage values completely.
  if (!cumulative) {
    rows_above_xlim <- plotting_data[plotting_data$cov > x_limit, ]
    plotting_data <- plotting_data[plotting_data$cov <= x_limit, ]
    if (!curve) {
      aggregated_line <- data.frame(cov = 1.15 * x_limit, 
                                    Metavar = unique(plotting_data$Metavar), 
                                    dist = sum(rows_above_xlim$dist))
      plotting_data <- rbind(plotting_data, aggregated_line)
    }
  }
  
  # If the distribution is not cumulative and shall be displayed as a bar chart, 
  # compute the median coverage.
  if (!cumulative & !curve) {
    single_cpg_covs <- lapply(unique(dist_df[, metavar_colname]), 
                              FUN = function(x) {
      metavar_subs <- dist_df[dist_df[, metavar_colname] == x, ]
      unlist(sapply(1:nrow(metavar_subs), FUN = function(y) {
        rep(metavar_subs[y, "cov"], metavar_subs[y, "count"])
      }))
    })
    median_covs <- sapply(single_cpg_covs, FUN = median)
    median_covs_df <- data.frame(Metavar = unique(dist_df[, metavar_colname]), 
                                 med = paste("Median:", median_covs))
    median_covs_df$Metavar <- factor(median_covs_df$Metavar, 
                                     levels = metavar_order)
  }
  
  # If the distribution is cumulative, compute the proportion of CpGs covered #
  # with a minimum coverage of cov_to_highlight.
  if (cumulative) {
    highlight <- dist_df[dist_df$cov == cov_to_highlight, 
                         c("cov", metavar_colname, "rel_count_cum")]
    colnames(highlight) <- c("cov", "Metavar", "dist")
    highlight$text_for_plot <- round(highlight$dist, 2)
    highlight$Metavar <- factor(highlight$Metavar, levels = metavar_order)
  }

  # Generate the plot.
  plotting_data$Metavar <- factor(plotting_data$Metavar, levels = metavar_order)
  cov_dist_plot <- ggplot(plotting_data) +
    aes(x = cov, y = dist) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), name = y_name) +
    facet_wrap(~Metavar, ncol = 3) +
    theme_classic()
  if (cumulative) {
    cov_dist_plot <- cov_dist_plot + 
      scale_x_continuous(limits = c(0, x_limit), 
                         expand = expansion(mult = c(0, 0.1)), name = x_name) +
      geom_line(size = 1.2) +
      geom_segment(data = highlight, 
                   mapping = aes(x = cov, xend = cov, yend = dist), 
                   y = 0, linetype = "dotted") +
      geom_segment(data = highlight, 
                   mapping = aes(xend = cov, y = dist, yend = dist), 
                   x = 0, linetype = "dotted") +
      geom_text(data = highlight, 
                mapping = aes(x = cov, y = dist, label = text_for_plot),
                nudge_x = 0.1 * x_limit) +
      theme(axis.text = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"))
  } else if (!curve) {
    cov_dist_plot <- cov_dist_plot + 
      scale_x_continuous(limits = c(-0.5, 1.3 * x_limit), 
                         expand = c(0, 0), name = x_name,
                         breaks = c(seq(0, x_limit, 10), 1.15 * x_limit),
                         labels = c(seq(0, x_limit, 10), 
                                    paste0(">", x_limit))) +
      geom_bar(stat = "identity", fill = "black", width = 1) +
      geom_vline(xintercept = x_limit + 0.5, linetype = "dotted") +
      geom_text(data = median_covs_df, 
                mapping = aes(label = med), x = 0.6 * x_limit, 
                y = 0.8 * max(plotting_data$dist)) +
      theme(axis.text = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"))
  } else {
    cov_dist_plot <- cov_dist_plot + 
      scale_x_continuous(limits = c(0, x_limit), 
                         expand = expansion(mult = c(0, 0.1)), name = x_name) +
      scale_colour_manual(breaks = metavar_2_order, 
                          values = metavar_2_palette, name = colour_name) +
      geom_line(mapping = aes(colour = Metavar_2)) +
      theme(axis.text = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            legend.key.height = unit(0.3, "inches"))
  }
  
  # Save the plot.
  plot_width <- ifelse(curve, yes = 10, no = 7.5)
  plot_height <- 1.25 * ceiling(length(levels(plotting_data$Metavar)) / 3)
  pdf(outfile_pdf, width = plot_width, height = plot_height)
  print(cov_dist_plot)
  dev.off()
  saveRDS(cov_dist_plot, file = outfile_rds)
  
}



# Specify the directories that contain the coverage values for the CPGs.
# I suppose that methylation calling has been performed by MethylCtools (at 
# least that's what Charles told me). According to the corresponding methyl 
# calling script 
# (https://github.com/hovestadt/methylCtools/blob/master/bcall.py; 20 May 2021),
# the second-to-last column of the files in these directories contain the 
# unconverted cytosines while the last column contains converted cytosines.
# The coordinates of the CpG sites are 0-based. Convert them to 1-based 
# coordinates.
coverage_dirs <- paste0(
  qual_dirs, 
  "/methylation/",
  ifelse(sample_mapping$Sample_ID == "green-8", yes ="lib1", no = "merged"),
  "/methylationCalling")
names(coverage_dirs) <- sample_names_for_plots

# Initialise the data frames that will later contain the coverage distributions.
coverage_dist_samples <- data.frame(cov = NA, count = NA, rel_count = NA, 
                                    rel_count_cum = NA, Sample = NA)
coverage_dist_cell_types <- data.frame(cov = NA, count = NA, rel_count = NA, 
                                       rel_count_cum = NA, Cell_type = NA)

# Iterate over all cell types.
num_processed_samples <- 0
for (cell_type in cell_types) { 
 
  # Find samples corresponding to this cell type.
  corresp_samples <- 
    sample_names_for_plots[sample_mapping$Cell_type == cell_type]
  
  # Initialise the data frame that will later be used to compute combined 
  # coverages in all samples of this cell type.
  combined_info_cell_type <- data.frame(V1 = NA, V2 = NA, cov = NA)
  
  # Iterate over all samples of this cell type.
  num_cpgs_already_printed <- F
  for (sample in corresp_samples) {

    # Identify the files in the directory that can be read as text files.
    # Exclude the file corresponding to the Y chromosome as all donors are 
    # female.
    text_files_in_dir <- list.files(coverage_dirs[sample])[
      grep("[^Y].CG.bed.gz$", list.files(coverage_dirs[sample]), perl = T)
    ]
  
    # Iterate over all files to combine their information.
    combined_info_sample <- data.frame(V1 = NA, V2 = NA, V6 = NA, V7 = NA, 
                                       cov = NA) # Placeholder line
    for (file in text_files_in_dir) {
      
      # Read the data in this file.
      data <- read.table(paste0(coverage_dirs[sample], "/", file), header = F, 
                         stringsAsFactors = F)
      assert(all(data$V4 == "CG"))
      
      # Ensure that all blocks of two consecutive rows correspond to the two 
      # complementary strands of the same CpG site.
      assert(nrow(data) %% 2 == 0)
      assert(all(sapply(seq(1, nrow(data), 2), FUN = function(x) { 
        data$V2[x + 1] - data$V2[x] 
      }) == 1))
      assert(all(sapply(seq(1, nrow(data), 2), FUN = function(x) { 
        data$V3[x] == "+" & data$V3[x + 1] == "-" 
      })))
       
      # Merge the two strands of each CpG site by summing up their methylated 
      # and unmethylated counts.
      assert(length(unique(data$V1)) == 1)
      data_merged <- data.frame(
        V1 = unique(data$V1),
        V2 = sapply(seq(1, nrow(data), 2), FUN = function(x) { 
          data$V2[x] 
        }) + 1, # The +1 is needed to convert the positions to 1-based 
        #         coordinates.
        V6 = sapply(seq(1, nrow(data), 2), FUN = function(x) { 
          sum(data$V6[c(x, x + 1)]) 
        }),
        V7 = sapply(seq(1, nrow(data), 2), FUN = function(x) { 
          sum(data$V7[c(x, x + 1)]) 
        })
      )
      
      # Compute the coverage at each merged CpG.
      data_merged$cov <- data_merged$V6 + data_merged$V7
      
      # Add the merged CpG sites to the combined data for this sample.
      combined_info_sample <- rbind(combined_info_sample, data_merged)

    }
    
    # Remove placeholder line.
    combined_info_sample <- combined_info_sample[-1, ]
    
    # Print a message stating the number of CpGs in the data.
    if (!num_cpgs_already_printed) {
      num_cpgs <- nrow(combined_info_sample)
      cat(paste0(num_cpgs, 
                 " CpG sites are considered in the data for '", 
                 cell_type, 
                 "'.\n"))
      num_cpgs_already_printed <- T
    }
    
    # If this is the first sample from this cell type, use the CpG sites from 
    # this sample as the CpG sites for the combined data set of this cell type.
    # Else, ensure that the CpG sites in this sample are the same as in the 
    # combined data set for the cell type and
    # add the coverage information from this sample to the combined coverage 
    # information for the cell type.
    if (nrow(combined_info_cell_type) == 1) {
      combined_info_cell_type <- rbind(
        combined_info_cell_type, 
        combined_info_sample[, c("V1", "V2", "cov")]
      )
      combined_info_cell_type <- 
        combined_info_cell_type[-1, ] # Remove placeholder line
    } else {
      assert(all(combined_info_cell_type$V1 == combined_info_sample$V1))
      assert(all(combined_info_cell_type$V2 == combined_info_sample$V2))
      combined_info_cell_type$cov <- 
        combined_info_cell_type$cov + combined_info_sample$cov
    }
    
    # Compute the coverage distribution.
    temp_df <- extract_coverage_distribution_from_df(
      combined_info_sample, max_coverage = 30000, coverage_column = "cov"
    )
    temp_df$Sample <- sample

    # Add the distribution info to the combined data for all samples.
    coverage_dist_samples <- rbind(coverage_dist_samples, temp_df)
    
    # Generate BSSeq object for this sample.
    m_matrix <- matrix(combined_info_sample$V6, ncol = 1)
    cov_matrix <- matrix(combined_info_sample$cov, ncol = 1)
    sample_names_vect <- sample
    chr_vect <- combined_info_sample$V1
    pos_vect <- combined_info_sample$V2
    bsseq_obj_sample <- bsseq::BSseq(
      M = m_matrix, Cov = cov_matrix, sampleNames = sample_names_vect,
      chr = chr_vect, pos = pos_vect, 
      pData = data.frame(Cell_type = cell_type, row.names = sample_names_vect)
    )
    
    # Append the new BSseq object to the combined BSSeq object for all sample.
    # If this doesn't exist yet, create it.
    if (exists("bsseq_obj_combined")) {
      bsseq_obj_combined <- bsseq::combine(bsseq_obj_combined, bsseq_obj_sample)
    } else {
      bsseq_obj_combined <- bsseq_obj_sample
    }
    
    # Update the number of samples that has already been processed.
    num_processed_samples <- num_processed_samples + 1
    
    # If the combined BSseq object contains two samples, save this a reduced 
    # BSseq object for later debugging purposes.
    if (num_processed_samples == 2) {
      bsseq_obj_reduced <- bsseq_obj_combined[1:3000000, 1:2]
      bsseq_reduced_outfile <- 
        paste0(data_out_pref, "_bsseq_object_combined_all_samples_reduced.rds")
      saveRDS(bsseq_obj_reduced, file = bsseq_reduced_outfile)
    }

  }
  
  # Compute the coverage distribution on the cell type level using the same 
  # approach as for the sample level.
  temp_df <- extract_coverage_distribution_from_df(
    combined_info_cell_type, max_coverage = 50000, coverage_column = "cov"
  )
  temp_df$Cell_type <- cell_type
  coverage_dist_cell_types <- rbind(coverage_dist_cell_types, temp_df)

}

# Save the BSSeq object for all samples as an RDS file.
bsseq_outfile <- paste0(data_out_pref, "_bsseq_object_combined_all_samples.rds")
saveRDS(bsseq_obj_combined, file = bsseq_outfile)

# Remove the first rows, which contains NAs.
coverage_dist_samples <- coverage_dist_samples[-1, ]
coverage_dist_cell_types <- coverage_dist_cell_types[-1, ]

# Save the data on CpG coverage distribution.
coverage_outfile_samples <- paste0(
  data_out_pref, "_cpg_coverage_distribution_wo_chr_y_by_sample.txt"
)
coverage_outfile_cell_types <- paste0(
  data_out_pref, "_cpg_coverage_distribution_wo_chr_y_by_cell_type.txt"
)
write.table(coverage_dist_samples, file = coverage_outfile_samples, sep = "\t", 
            col.names = T, row.names = F, quote = F)
write.table(coverage_dist_cell_types, file = coverage_outfile_cell_types, 
            sep = "\t", col.names = T, row.names = F, quote = F)

# Add cell type information to the sample-wise distribution.
coverage_dist_samples$Cell_type <- sapply(
  coverage_dist_samples$Sample, FUN = function(x) {
    strsplit(strsplit(x, split = "(", fixed = T)[[1]][2], split = ";")[[1]][1]
  }
)

# Plot and save the coverage distribution by sample.
plot_coverage_distribution_by_metavar(
  coverage_dist_samples, 
  metavar_colname = "Sample",
  metavar_order = sample_names_for_plots,
  x_limit = 50,
  x_name = "Coverage",
  outfile_pdf = paste0(plot_outdir, 
                       "/cpg_coverage_distribution_wo_chr_y_by_sample.pdf"),
  outfile_rds = paste0(plot_rds_outdir, 
                       "/cpg_coverage_distribution_wo_chr_y_by_sample.rds")
)

# Plot and save the coverage distribution by sample and cell type, displayed as 
# a line.
plot_coverage_distribution_by_metavar(
  coverage_dist_samples, 
  metavar_colname = "Cell_type",
  metavar_order = cell_types,
  curve = T,
  metavar_2_colname = "Sample", 
  metavar_2_order = sample_names_for_plots,
  metavar_2_palette = sample_col_palette,
  colour_name = "Sample",
  x_limit = 50,
  x_name = "Coverage",
  y_name = "Proportion of CpGs (excluding chromosome Y)",
  outfile_pdf = paste0(
    plot_outdir, "/cpg_coverage_distribution_wo_chr_y_by_sample_with_line.pdf"
  ),
  outfile_rds = paste0(
    plot_rds_outdir, 
    "/cpg_coverage_distribution_wo_chr_y_by_sample_with_line.rds"
  )
)

# Plot and save the cumulative coverage distribution by sample.
plot_coverage_distribution_by_metavar(
  coverage_dist_samples, 
  metavar_colname = "Sample",
  metavar_order = sample_names_for_plots,
  cumulative = T,
  x_limit = 50,
  x_name = "Minimum coverage",
  y_name = "Proportion of CpGs (excluding chromosome Y)",
  outfile_pdf = paste0(
    plot_outdir, "/cumulative_cpg_coverage_distribution_wo_chr_y_by_sample.pdf"
  ),
  outfile_rds = paste0(
    plot_rds_outdir, 
    "/cumulative_cpg_coverage_distribution_wo_chr_y_by_sample.rds"
  )
)

# Plot and save the coverage distribution by cell type.
plot_coverage_distribution_by_metavar(
  coverage_dist_cell_types, 
  metavar_colname = "Cell_type",
  metavar_order = cell_types,
  x_limit = 100,
  x_name = "Coverage",
  y_name = "Proportion of CpGs (excluding chromosome Y)",
  outfile_pdf = paste0(
    plot_outdir, "/cpg_coverage_distribution_wo_chr_y_by_cell_type.pdf"
  ),
  outfile_rds = paste0(
    plot_rds_outdir, "/cpg_coverage_distribution_wo_chr_y_by_cell_type.rds"
  )
)

# Plot and save the cumulative coverage distribution by cell type.
plot_coverage_distribution_by_metavar(
  coverage_dist_cell_types, 
  metavar_colname = "Cell_type",
  metavar_order = cell_types,
  cumulative = T,
  cov_to_highlight = 20,
  x_limit = 100,
  x_name = "Minimum coverage",
  y_name = "Proportion of CpGs (excluding chromosome Y)",
  outfile_pdf = paste0(
    plot_outdir, 
    "/cumulative_cpg_coverage_distribution_wo_chr_y_by_cell_type.pdf"
  ),
  outfile_rds = paste0(
    plot_rds_outdir, 
    "/cumulative_cpg_coverage_distribution_wo_chr_y_by_cell_type.rds"
  )
)
