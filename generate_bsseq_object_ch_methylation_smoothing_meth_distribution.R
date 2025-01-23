# This script reads in the CH methylation data for the single sequencing
#   libraries and generates a BSseq object out of them. Afterwards, it smoothes
#   this data and plots the distribution of smoothed methylation values.
# Author: Niklas Beumer



# Load required packages.
library(bsseq)
library(testit)
library(ggplot2)
library(ggridges)
library(foreach)
library(doParallel)
library(BiocParallel)

# Specify the location where the sequencing results can be found.
seq_results_loc <- "/yyy/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid"

# Specify a location in my space.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Specify the prefix for output text and RDS files that will later contain the 
# processed data.
data_out_pref <- paste0(
  location,
  "/preprocessing_etc/quality_control_after_alignment/",
  format(Sys.time(), "%Y-%m-%d")
)

# Read in the sample mapping file.
sample_mapping <- read.table(
  paste0(location, "/sample_mapping_twgbs.txt"),
  header = T,
  sep = "\t",
  stringsAsFactors = F
)

# Extract the list of analysed cell types.
cell_types <- unique(sample_mapping$Cell_type)

# Specify the sample-level directories where methylation values can be found.
meth_dirs_sample_level <- paste0(
  seq_results_loc,
  "/OE0436_Treg_cells_",
  sample_mapping$Patient_ID,
  "/",
  sample_mapping$Sample_ID,
  "/paired/merged-alignment/methylation"
)

# Specify sample names that will later appear in the bsseq object meta data.
sample_names <- paste0(
  sample_mapping$Sample_ID_readable,
  "\n",
  "(",
  sample_mapping$Cell_type,
  "; ",
  sample_mapping$Patient_ID_readable,
  ")"
)
names(meth_dirs_sample_level) <- sample_names

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis",
                     format(Sys.time(), "%m-%d-%y"),
                     sep = "/")
if (!dir.exists(plot_outdir)) {
  dir.create(plot_outdir)
}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")






####################################################################################################
# Read in the methylation data for the different libraries and generate a BSseq object out of them.
####################################################################################################

# Iterate over all sample-library combinations.
# Parallelise this process on 21 cores.
samples_and_libs <- do.call(c, lapply(
  sample_names,
  FUN = function(x) {
    list(c(x, "lib1"), c(x, "lib2"))
  }
))
registerDoParallel(21)
bsseq_objects <- foreach(comb = samples_and_libs) %dopar% {
  # Specify the directory containing the methylation values.
  # I suppose that methylation calling has been performed by MethylCtools (at least that's what Charles told me).
  # According to the corresponding methyl calling script (https://github.com/hovestadt/methylCtools/blob/master/bcall.py; 20 May 2021),
  # the second-to-last column of the files in these directories contain the unconverted cytosines while the last column contains converted cytosines.
  # The coordinates a 0-based. Convert them to 1-based coordinates.
  sample <- comb[1]
  library <- comb[2]
  meth_data_dir <- paste(meth_dirs_sample_level[sample],
                         library,
                         "methylationCalling",
                         sep = "/")
  
  # Identify the files in the directory that contain data on CH methylation and can be read as text files.
  # Restrict this list to chromosome 1 in order to stay within reasonable memory limits.
  text_file_chr_1 <- list.files(meth_data_dir)[grep("bam.1.CH.bed.gz$", list.files(meth_data_dir), perl = T)]
  # text_file_chr_1 <- list.files(meth_data_dir)[grep("bam.MT.CH.bed.gz$", list.files(meth_data_dir), perl = T)] # For debugging purposes
  
  # Read the data in this file.
  data <- read.table(
    paste0(meth_data_dir, "/", text_file_chr_1),
    header = F,
    stringsAsFactors = F
  )
  assert(all(data$V4 == "CH"))
  
  # Restrict the data to a maximum of 30,000,000 CH sites.
  data <- data[1:min(nrow(data), 3000000), ]
  
  # Compute the coverage at each CH.
  data$cov <- data$V6 + data$V7
  
  # Generate a BSSeq object for this sample-library combination.
  m_matrix <- matrix(data$V6, ncol = 1)
  cov_matrix <- matrix(data$cov, ncol = 1)
  sample_name_split <- strsplit(sample, split = "\n", fixed = T)[[1]]
  lib_names_vect <- paste0(
    sample_name_split[1],
    "; Library ",
    strsplit(library, split = "ib")[[1]][2],
    "\n",
    sample_name_split[2]
  )
  chr_vect <- data$V1
  pos_vect <- data$V2 + 1 # The +1 is needed here to convert the positions from 
  #                         0-based coordinates to 1-based coordinates.
  bsseq_obj <- BSseq(
    M = m_matrix,
    Cov = cov_matrix,
    sampleNames = lib_names_vect,
    chr = chr_vect,
    pos = pos_vect,
    pData = data.frame(
      Sample = sample,
      Cell_type = strsplit(strsplit(
        sample, split = "(", fixed = T
      )[[1]][2], split = ";")[[1]][1],
      Library = library,
      row.names = lib_names_vect
    )
  )
  
  # Return the BSseq object for the current sample-library combination so that 
  # it can be appended to the list returned from the parallel process.
  bsseq_obj
  
}
stopImplicitCluster()

# Combine the generated BSseq objects together.
meth_data_combined <- do.call(bsseq::combine, bsseq_objects)





###################################################
# Smooth the methylation data in the BSseq object.
###################################################

# Smooth methylation estimates.
meth_data_combined <- BSmooth(meth_data_combined, BPPARAM = MulticoreParam(
  workers = ncol(meth_data_combined),
  progressbar = T
))

# Save the BSseq object with the smoothed methylation estimates.
smoothed_obj_outfile <- paste0(
  data_out_pref,
  "_bsseq_object_combined_all_samples_separate_for_libraries_ch_methylation_chr1_3_mil_smoothed.rds"
)
saveRDS(meth_data_combined, file = smoothed_obj_outfile)



##########################################################################
# Plot the distribution of smoothed methylation estimates in the samples.
##########################################################################

# Extract smoothed methylation values.
meth_values_smoothed <- as.data.frame(getMeth(meth_data_combined, type = "smooth", what = "perBase"))

# Prepare the methylation values for plotting.
beta_dist_plotting_data <- reshape2::melt(meth_values_smoothed,
                                          variable.name = "Sample_Lib_Comb",
                                          value.name = "Meth")
beta_dist_plotting_data$Sample <- factor(pData(meth_data_combined)[beta_dist_plotting_data$Sample_Lib_Comb, "Sample"], levels = sample_names)
beta_dist_plotting_data$Library <- factor(pData(meth_data_combined)[beta_dist_plotting_data$Sample_Lib_Comb, "Library"], levels = c("lib2", "lib1"))

# Generate and save a beta value distribution plot that uses a density geom.
beta_dist_plot <- ggplot(beta_dist_plotting_data) +
  aes(
    x = Meth,
    y = Library,
    colour = Library,
    fill = Library
  ) +
  scale_x_continuous(limits = c(0, 1),
                     name = "Methylation value",
                     expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(expand = expansion(mult = c(0.1, 2))) +
  scale_colour_manual(breaks = c("lib1", "lib2"),
                      values = c("red", "blue")) +
  scale_fill_manual(breaks = c("lib1", "lib2"),
                    values = c("red", "blue")) +
  geom_density_ridges() +
  facet_wrap( ~ Sample, ncol = 3, scales = "fixed") +
  ylab("Density by sequencing library") +
  theme_classic() +
  theme(
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.position = "none"
  )
plot_width <- 8
plot_height <- 1.25 * ceiling(length(sample_names) / 3)
pdf(
  paste0(
    plot_outdir,
    "/beta_value_distribution_ch_methylation_chr1_3_mil_smoothed_values_by_sample_with_density.pdf"
  ),
  width = plot_width,
  height = plot_height
)
print(beta_dist_plot)
dev.off()
saveRDS(
  beta_dist_plot,
  file = paste0(
    plot_rds_outdir,
    "/beta_value_distribution_ch_methylation_chr1_3_mil_smoothed_values_by_sample_with_density.rds"
  )
)
