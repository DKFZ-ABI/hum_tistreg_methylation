# This script analyses global methylation on the cell type and sample level.
# Author: Niklas Beumer



# Load required packages.
library(bsseq)
library(ggplot2)
library(parallel)

# Define a location where input is saved and where output will be saved.
in_out_location <- "/yyy/hm_treg_bs_rgnsbg"

# Define the prefix of the output data file.
data_out_pref <- paste0(
  in_out_location,
  "/preprocessing_etc/quality_control_after_alignment/",
  format(Sys.time(), "%Y-%m-%d", sep = "/")
)

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(in_out_location, "/plot_rds_objects")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis",
                     format(Sys.time(), "%m-%d-%y"),
                     sep = "/")
if (!dir.exists(plot_outdir)) {
  dir.create(plot_outdir)
}

# Read in the BSseq object containing the methylation data.
meth_data <- readRDS(
  paste(
    in_out_location,
    "preprocessing_etc/quality_control_after_alignment/2022-03-11_bsseq_object_combined_all_samples.rds",
    sep = "/"
  )
)

# Generate colour scales for the cell types.
cell_type_col_palette <- c("blue", "cyan2", "orange", "darkorchid1", "red")





#################################################
# Compute global methylation for single samples.
#################################################

remove_newlines_from_string <- function(string_name) {
  # This function removes any newline characters from a character string and 
  # replaces them by a blank space.
  # string_name: Character: The string to remove newline characters from.
  # Value: The character string with the newline character removed.
  # Dependencies: None
  # Author: Niklas Beumer
  
  return(paste(strsplit(
    string_name, split = "\n", fixed = T
  )[[1]], collapse = " "))
  
}



# Extract the coverage levels for each CpG.
coverage_levels <- getCoverage(meth_data, type = "Cov", what = "perBase")

# Extract the methylation levels.
methylation_levels <- getMeth(meth_data, type = "raw", what = "perBase")

# In order to only use accurate methylation values, only consider CpG sites with 
# a coverage of at least 2 in all samples. The constraint "all samples" is used 
# to ensure consideration of the same sites between the samples, thereby 
# ensuring full comparability between the samples.
# Also filter the coverage data in order to include the same CpG sites in 
# subsequent coverage analysis.
loci_to_keep <- which(apply(
  coverage_levels,
  1,
  FUN = function(x) {
    all(x >= 2)
  }
))
methylation_levels_filtered <- methylation_levels[loci_to_keep, ]
coverage_levels_filtered <- coverage_levels[loci_to_keep, ]

# Extract genomic positions.
# Note to myself: I was able to confirm (for some selected regions) that the order of the
# CpG sites returned by rowRanges coresponds to the oder of CpG sites returned by getMeth.
cpg_loci <- rowRanges(meth_data)
cpg_loci_filtered <- cpg_loci[loci_to_keep]

# Separate the CpG loci into autosomes and chromosome X.
autosome_loci_inds <- which(seqnames(cpg_loci_filtered) %in% as.character(1:22))
methylation_levels_filtered_aut <- methylation_levels_filtered[autosome_loci_inds, ]
coverage_levels_filtered_aut <- coverage_levels_filtered[autosome_loci_inds, ]
x_loci_inds <- which(seqnames(cpg_loci_filtered) == "X")
methylation_levels_filtered_x <- methylation_levels_filtered[x_loci_inds, ]
coverage_levels_filtered_x <- coverage_levels_filtered[x_loci_inds, ]

# Print a status message stating how many CpGs are considered for global methylation analysis.
cat(
  paste(
    nrow(methylation_levels_filtered_aut) + nrow(methylation_levels_filtered_x),
    "out of",
    nrow(methylation_levels),
    "CpG sites are considered for the computation of global methylation.\n"
  )
)

# Compute global methylation in each sample using the retained CpG.
# Global methylation is defined as the mean methylation level across these CpGs.
global_methylation_samples_aut <- apply(methylation_levels_filtered_aut, 2, FUN = mean)
global_methylation_samples_x <- apply(methylation_levels_filtered_x, 2, FUN = mean)

# Save the global methylation values for each sample.
global_meth_samp_aut_outfile <- paste0(data_out_pref, "_global_methylation_samples_autosomes.rds")
global_meth_samp_x_outfile <- paste0(data_out_pref, "_global_methylation_samples_chr_x.rds")
saveRDS(global_methylation_samples_aut, file = global_meth_samp_aut_outfile)
saveRDS(global_methylation_samples_x, file = global_meth_samp_x_outfile)

# Prepare the global methylation values to be considered in the plot.
global_meth_df_samples <- data.frame(
  Sample = rep(names(global_methylation_samples_aut), 2),
  Meth = c(
    global_methylation_samples_aut,
    global_methylation_samples_x
  ),
  class = rep(
    c("Autosomes", "Chromosome X"),
    each = length(global_methylation_samples_aut)
  )
)

# Remove newline characters from the sample names to enable nicer axis labels.
colnames(methylation_levels_filtered_aut) <- sapply(colnames(methylation_levels_filtered_aut), FUN = remove_newlines_from_string)
colnames(coverage_levels_filtered_aut) <- sapply(colnames(coverage_levels_filtered_aut), FUN = remove_newlines_from_string)
colnames(methylation_levels_filtered_x) <- sapply(colnames(methylation_levels_filtered_x), FUN = remove_newlines_from_string)
colnames(coverage_levels_filtered_x) <- sapply(colnames(coverage_levels_filtered_x), FUN = remove_newlines_from_string)
global_meth_df_samples$Sample <-  sapply(global_meth_df_samples$Sample, FUN = remove_newlines_from_string)

# Melt the methylation data so that they can be plotted.
methylation_levels_filtered_aut_melted <- reshape2::melt(
  as.data.frame(methylation_levels_filtered_aut),
  variable.name = "Sample",
  value.name = "Meth"
)
methylation_levels_filtered_aut_melted$class <- "Autosomes"
methylation_levels_filtered_x_melted <- reshape2::melt(
  as.data.frame(methylation_levels_filtered_x),
  variable.name = "Sample",
  value.name = "Meth"
)
methylation_levels_filtered_x_melted$class <- "Chromosome X"
methylation_levels_filtered_comb_melted <- rbind(methylation_levels_filtered_aut_melted,
                                                 methylation_levels_filtered_x_melted)

# Generate a box plot showing the methylation in each sample.
meth_plot_samples <- ggplot(methylation_levels_filtered_comb_melted) +
  aes(x = Sample, y = Meth, colour = class) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     limits = c(0, 1),
                     name = "Methylation") +
  scale_colour_manual(values = c("red", "blue"), name = "Chromosome\nclass") +
  scale_fill_manual(values = c("red", "blue"), name = "Chromosome\nclass") +
  geom_boxplot() +
  geom_point(
    data = global_meth_df_samples,
    mapping = aes(y = Meth, fill = class),
    shape = 24,
    position = position_dodge(width = 0.8)
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      colour = "black",
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

# Save the plot.
plot_width = 100 * ceiling(ncol(methylation_levels_filtered) / 3)
plot_height = 500
png(
  paste0(plot_outdir, "/global_methylation_by_sample_boxplot.png"),
  width = plot_width,
  height = plot_height
)
print(meth_plot_samples)
dev.off()
saveRDS(
  meth_plot_samples,
  file = paste0(
    plot_rds_outdir,
    "/global_methylation_by_sample_boxplot.rds"
  )
)







#################################################
# Compute global methylation for the cell types.
#################################################

extract_cell_types_from_sample_names <- function(sample_names) {
  # This function extracts the cell types from sample names that were generated according to my convention
  #   in the Treg-WGBS project.
  # sample_names: Character; A vector of sample names.
  # Dependencies: None.
  # Value: A vector of cell types, the order of which corresponds to the order of sample_names.
  # Author: Niklas Beumer
  
  # Extract the cell type names.
  cell_type_names <- sapply(
    sample_names,
    FUN = function(x) {
      strsplit(strsplit(x, split = "(", fixed = T)[[1]][2], split = ";")[[1]][1]
    }
  )
  
  # Return the cell type names.
  return(cell_type_names)
  
}

sample_level_to_cell_type_level <- function(global_meth_samples, cell_types) {
  # This function takes global methylation values on sample levels and returns global methylation
  #   values on cell type level.
  # global_meth_samples: Numeric; Global methylation computed in each sample. The names must
  #   be the sample names.
  # cell_types: Character; The cell types each sample corresponds to. The names must be the
  #   sample names.
  # Dependencies: None
  # Value: A matrix with the mean global methylation across samples in the first row, the
  #   lower limit of the SD error bar in the second row and the upper limit of the SD error bar
  #   in the third row.
  # Author: Niklas Beumer
  
  # Iterate over the cell types and compute mean global methylation as well as the SD.
  cell_type_values <- sapply(
    unique(cell_types),
    FUN = function(x) {
      corresp_samples <- names(cell_types[cell_types == x])
      meth_values <- global_meth_samples[corresp_samples]
      mean_value <- mean(meth_values)
      sd_value <- sd(meth_values)
      return(
        c(
          mean = mean_value,
          sd_lower = mean_value - sd_value,
          sd_upper = mean_value + sd_value
        )
      )
    }
  )
  
  return(cell_type_values)
}


# Extract cell type names from sample names.
cell_types <- extract_cell_types_from_sample_names(names(global_methylation_samples_aut))

# Compute global methylation and the corresponding SD for each sample.
# Here, global methylation is defined as the mean global methylation across all samples from a cell type.
global_methylation_cell_types_aut <- sample_level_to_cell_type_level(global_methylation_samples_aut, cell_types)
global_methylation_cell_types_x <- sample_level_to_cell_type_level(global_methylation_samples_x, cell_types)

# Save global methylation on the cell type level.
global_meth_cell_type_aut_outfile <- paste0(data_out_pref,
                                            "_global_methylation_cell_types_autosomes.rds")
global_meth_cell_type_x_outfile <- paste0(data_out_pref, "_global_methylation_cell_types_chr_x.rds")
saveRDS(global_methylation_cell_types_aut, file = global_meth_cell_type_aut_outfile)
saveRDS(global_methylation_cell_types_x, file = global_meth_cell_type_x_outfile)

# Perform Welch's t test comparing each pair of cell types.
# Afterwards, perform Benjamini-Hochberg correction.
cell_types_unique <- unique(cell_types)
p_val_df <- do.call(rbind, lapply(
  1:(length(cell_types_unique) - 1),
  FUN = function(x) {
    do.call(rbind, lapply((x + 1):length(cell_types_unique),
                          FUN = function(y) {
                            cell_type_1 <- cell_types_unique[x]
                            corresp_samples_1 <- names(cell_types)[cell_types == cell_type_1]
                            cell_type_2 <- cell_types_unique[y]
                            corresp_samples_2 <- names(cell_types)[cell_types == cell_type_2]
                            pvals <- sapply(
                              list(
                                global_methylation_samples_aut,
                                global_methylation_samples_x
                              ),
                              FUN = function(z) {
                                test_result <- t.test(z[corresp_samples_1], z[corresp_samples_2], alternative = "two.sided")
                                return(c(test_result$statistic, test_result$p.value))
                              }
                            )
                            temp_df <- data.frame(
                              Cell_type_1 = cell_type_1,
                              Cell_type_2 = cell_type_2,
                              Chromosome_category = c("Autosomes", "Chromosome X"),
                              welch_t_statistic = pvals[1, ],
                              welch_pval = pvals[2, ]
                            )
                            return(temp_df)
                          }
    ))
  }
))
p_val_df$welch_pval_adj_BH <- p.adjust(p_val_df$welch_pval, method = "BH")
p_val_outfile <- paste0(data_out_pref,
                        "_global_methylation_cell_types_stat_test_results.txt")
write.table(
  p_val_df,
  file = p_val_outfile,
  sep = "\t",
  row.names = F,
  col.names = T,
  quote = F
)

# Prepare the data on global methylation on the cell type level for plotting.
global_meth_cell_types_aut_df <- as.data.frame(t(global_methylation_cell_types_aut))
global_meth_cell_types_aut_df$Cell_type <- rownames(global_meth_cell_types_aut_df)
global_meth_cell_types_x_df <- as.data.frame(t(global_methylation_cell_types_x))
global_meth_cell_types_x_df$Cell_type <- rownames(global_meth_cell_types_x_df)
global_meth_cell_types_comb_df <- rbind(global_meth_cell_types_aut_df, global_meth_cell_types_x_df)
global_meth_cell_types_comb_df$Cell_type <- factor(global_meth_cell_types_comb_df$Cell_type,
                                                   levels = unique(cell_types))
global_meth_cell_types_comb_df$class <- rep(c("Autosomes", "Chromosome X"), each = length(unique(cell_types)))
single_values_df <- data.frame(
  Cell_type = c(cell_types[names(global_methylation_samples_aut)], cell_types[names(global_methylation_samples_x)]),
  Meth = c(
    global_methylation_samples_aut,
    global_methylation_samples_x
  ),
  class = rep(
    c("Autosomes", "Chromosome X"),
    each = length(global_methylation_samples_aut)
  )
)
single_values_df$Cell_type <- factor(single_values_df$Cell_type, levels = unique(cell_types))

# Generate a bar chart showing global methylation in the cell types including SD.
meth_plot_cell_types <- ggplot(global_meth_cell_types_comb_df) +
  aes(
    x = Cell_type,
    y = mean,
    ymin = sd_lower,
    ymax = sd_upper,
    fill = class
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     limits = c(0, 1),
                     name = "Global methylation") +
  scale_colour_manual(values = c("red", "blue"), name = "Chromosome\nclass") +
  scale_fill_manual(values = c("red", "blue"), name = "Chromosome\nclass") +
  geom_bar(stat = "identity",
           colour = "black",
           position = "dodge") +
  geom_point(
    data = single_values_df,
    mapping = aes(x = Cell_type, y = Meth, fill = class),
    inherit.aes = F,
    colour = "black",
    size = 0.75,
    position = position_jitterdodge(
      dodge.width = 0.8,
      jitter.width = 0.2,
      seed = 1
    )
  ) +
  geom_errorbar(width = 0.4, position = position_dodge(width = 0.8)) +
  xlab("Cell type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      colour = "black",
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

# Save the plot.
plot_width = ceiling(nrow(global_meth_cell_types_comb_df) / 3)
plot_height = 3.5
pdf(
  paste0(plot_outdir, "/global_methylation_by_cell_type.pdf"),
  width = plot_width,
  height = plot_height
)
print(meth_plot_cell_types)
dev.off()
saveRDS(
  meth_plot_cell_types,
  file = paste0(plot_rds_outdir, "/global_methylation_by_cell_type.rds")
)





#################################################################################################
# Visualise the methylation value distribution and coverage distribution on the cell type level.
#################################################################################################

# For each CpG site compute the mean methylation and mean coverage across the samples from a cell type.
methylation_levels_cell_types_aut <- unlist(mclapply(
  unique(cell_types),
  FUN = function(x) {
    corresp_samples <- sapply(names(cell_types)[cell_types == x], FUN = remove_newlines_from_string)
    apply(
      methylation_levels_filtered_aut[, corresp_samples],
      1,
      FUN = function(x) {
        mean(x, na.rm = T)
      }
    )
  },
  mc.cores = 7
))
coverage_levels_cell_types_aut <- unlist(mclapply(
  unique(cell_types),
  FUN = function(x) {
    corresp_samples <- sapply(names(cell_types)[cell_types == x], FUN = remove_newlines_from_string)
    apply(
      coverage_levels_filtered_aut[, corresp_samples],
      1,
      FUN = function(x) {
        mean(x, na.rm = T)
      }
    )
  },
  mc.cores = 7
))
methylation_levels_cell_types_x <- unlist(mclapply(
  unique(cell_types),
  FUN = function(x) {
    corresp_samples <- sapply(names(cell_types)[cell_types == x], FUN = remove_newlines_from_string)
    apply(
      methylation_levels_filtered_x[, corresp_samples],
      1,
      FUN = function(x) {
        mean(x, na.rm = T)
      }
    )
  },
  mc.cores = 7
))
coverage_levels_cell_types_x <- unlist(mclapply(
  unique(cell_types),
  FUN = function(x) {
    corresp_samples <- sapply(names(cell_types)[cell_types == x], FUN = remove_newlines_from_string)
    apply(
      coverage_levels_filtered_x[, corresp_samples],
      1,
      FUN = function(x) {
        mean(x, na.rm = T)
      }
    )
  },
  mc.cores = 7
))

# Prepare the data to plot the methylation value distribution in the cell types.
meth_val_distr_data_aut <- data.frame(
  Cell_type = rep(unique(cell_types), each = nrow(methylation_levels_filtered_aut)),
  level = methylation_levels_cell_types_aut,
  class = "Autosomes"
)
meth_val_distr_data_x <- data.frame(
  Cell_type = rep(unique(cell_types), each = nrow(methylation_levels_filtered_x)),
  level = methylation_levels_cell_types_x,
  class = "Chromosome X"
)
meth_val_distr_data_comb <- rbind(meth_val_distr_data_aut, meth_val_distr_data_x)
meth_val_distr_data_comb$Cell_type <- factor(meth_val_distr_data_comb$Cell_type, levels = unique(cell_types))

# Plot methylation value distribution in the different cell types.
# Note, I do not use a density geom here because I noticed that stat_density removes the
# peaks at 0 and 1 from the X chromosome data but not from the autosome data. This is
# misleading, in my opinion. This may mean that the degree of smoothing depends on
# the number of CpG sites and may thus be different for the two categories.
# Instead, I use stacked histograms.
meth_val_distr_plot <- ggplot(meth_val_distr_data_comb) +
  aes(x = level, fill = Cell_type) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), name = "CpG number") +
  scale_fill_manual(breaks = unique(cell_types),
                    values = cell_type_col_palette,
                    name = "Cell type") +
  geom_histogram(binwidth = 0.01) +
  facet_wrap( ~ class, ncol = 2, scales = "free_y") +
  xlab("Mean methylation value (across samples)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

# Save the plot.
pdf(
  paste0(
    plot_outdir,
    "/methylation_value_distribution_by_cell_type.pdf"
  ),
  width = 7,
  height = 3
)
print(meth_val_distr_plot)
dev.off()
saveRDS(
  meth_val_distr_plot,
  file = paste0(
    plot_rds_outdir,
    "/methylation_value_distribution_by_cell_type.rds"
  )
)



# Prepare the data to plot the coverage value distribution in the cell types.
cov_val_distr_data_aut <- data.frame(
  Cell_type = rep(unique(cell_types), each = nrow(coverage_levels_filtered_aut)),
  level = coverage_levels_cell_types_aut,
  class = "Autosomes"
)
cov_val_distr_data_x <- data.frame(
  Cell_type = rep(unique(cell_types), each = nrow(coverage_levels_filtered_x)),
  level = coverage_levels_cell_types_x,
  class = "Chromosome X"
)
cov_val_distr_data_comb <- rbind(cov_val_distr_data_aut, cov_val_distr_data_x)
cov_val_distr_data_comb$level[cov_val_distr_data_comb$level > 15] <- 20
cov_val_distr_data_comb$Cell_type <- factor(cov_val_distr_data_comb$Cell_type, levels = unique(cell_types))

# Plot coverage value distribution in the different cell types.
cov_val_distr_plot <- ggplot(cov_val_distr_data_comb) +
  aes(x = level, fill = Cell_type) +
  scale_x_continuous(
    name = "Mean coverage (across samples)",
    breaks = c(5, 10, 15, 20),
    labels = c("5", "10", "15", ">15")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), name = "CpG number") +
  scale_fill_manual(breaks = unique(cell_types),
                    values = cell_type_col_palette,
                    name = "Cell type") +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = 15.5, linetype = "dotted") +
  facet_wrap( ~ class, ncol = 2, scales = "free_y") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

# Save the plot.
pdf(
  paste0(
    plot_outdir,
    "/coverage_value_distribution_by_cell_type.pdf"
  ),
  width = 7,
  height = 3
)
print(cov_val_distr_plot)
dev.off()
saveRDS(
  cov_val_distr_plot,
  file = paste0(
    plot_rds_outdir,
    "/coverage_value_distribution_by_cell_type.rds"
  )
)
