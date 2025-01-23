# This script performs PCA on the merged WGBS samples.
# Author: Niklas Beumer



# Load required packages.
library(bsseq)
library(ggplot2)
library(plot3D)

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
# meth_data <- readRDS(paste(in_out_location, "preprocessing_etc/quality_control_after_alignment/2022-03-11_bsseq_object_combined_all_samples_reduced.rds", sep = "/")) # For debugging purposes







####################################
# Generate 5-kb bins of the genome.
####################################

# Get lengths of the chromosomes in the reference genome used for alignment.
chr_lengths <- GenomeInfoDb::seqlengths(
  BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
)[c(as.character(1:22), "X", "MT")]

# Generate Granges specifying 5-kb bins.
# Coordinates are 1-based and intervals are closed.
bin_length <- 5000
chr_names <- unlist(sapply(
  1:length(chr_lengths),
  FUN = function(x) {
    rep(names(chr_lengths)[x], ceiling(chr_lengths[x] / bin_length))
  }
))
bin_starts <- unlist(sapply(
  chr_lengths,
  FUN = function(x) {
    seq(1, x, bin_length)
  }
))
bin_ends <- bin_starts + bin_length - 1
genome_bins <- GenomicRanges::GRanges(
  seqnames = chr_names,
  ranges = IRanges::IRanges(start = bin_starts, end = bin_ends))





############################################################
# Perform principal component analysis (PCA) on the samples.
############################################################

PCA_binned_methylation <- function(bsseq_obj, bins, outfile_rds) {
  # This function performs PCA on binned methylation data. It also generates a score plot.
  # bsseq_obj: The BSSeq object containing the methylation values to use.
  # bins: A GRanges object containing the genomic bins to use.
  # outfile_rds: Character; Path to an RDS file where the PCA output will be saved.
  # Dependencies: bsseq
  # Value: The PCA results returned by prcomp. Additionally, it saves the PCA results in the specified RDS file.
  # Author: Niklas Beumer
  
  # Compute average methylation in the different bins.
  binned_methylation <- getMeth(bsseq_obj,
                                regions = bins,
                                type = "raw",
                                what = "perRegion")
  
  # Get the mean coverage in the bins.
  binned_coverage <- getCoverage(bsseq_obj,
                                 regions = bins,
                                 type = "Cov",
                                 what = "perRegionAverage")
  
  # Only keep methylation in bins that have a coverage of at least 5 in every sample.
  bins_to_keep <- which(apply(
    binned_coverage,
    1,
    FUN = function(x) {
      all(!is.na(x)) & all(x >= 2)
    }
  ))
  binned_methylation_filtered <- binned_methylation[bins_to_keep, ]
  
  # Remove bins in which at least one sample contains a missing methylation value.
  bins_to_keep_2 <- which(apply(
    binned_methylation_filtered,
    1,
    FUN = function(x) {
      all(!is.na(x) & !is.nan(x))
    }
  ))
  binned_methylation_filtered <- binned_methylation_filtered[bins_to_keep_2, ]
  
  # Print a status message stating how many bins are retained.
  total_num_bins <- nrow(binned_methylation)
  num_bins_retained <- nrow(binned_methylation_filtered)
  cat(
    paste(
      num_bins_retained,
      "out of",
      total_num_bins,
      "bins passed filtering.\nPCA will be performed on the 2,000 most variable regions.\n"
    )
  )
  
  # Identify the 2,000 most variable bins.
  meth_vars <- apply(binned_methylation_filtered, 1, FUN = var)
  meth_most_variable <- binned_methylation_filtered[order(meth_vars, decreasing = T)[1:2000], ]
  
  # Perform PCA.
  pca_res <- prcomp(t(na.omit(meth_most_variable)), scale. = T)
  
  # Save the PCA results.
  saveRDS(pca_res, file = outfile_rds)
  
  # Return the PCA results.
  return(pca_res)
  
}

# Run PCA on different subsets of the data.
# Print a status message before every PCA analysis.
cat("Running PCA on all samples.\n")
pca_res_all <- PCA_binned_methylation(
  meth_data,
  bins = genome_bins,
  outfile_rds = paste0(data_out_pref, "_pca_results_samples.rds")
)







#################################
# Plot the samples in PCA space.
#################################

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


prepare_pca_plot <- function(pca_res,
                             colour_categories,
                             colour_palette,
                             colour_name,
                             outfile_pdf,
                             outfile_rds) {
  # This function generates a score plot out of prcomp PCA results.
  # pca_res: An object returned by prcomp.
  # colour_categories: Character; A vector containing categories to colour the samples by.
  #   The order of this vector must correspond to the order of rownames(pca_res)
  # colour_palette: Character; Vector containing the colours to use for the categories in colour_categories.
  #   Its order must correspond to the order of unique(colour_categories).
  # colour_name: Character; The name that will appear in the colour legend.
  # outfile_pdf: Character; The path to a PDF file where the plot will be saved.
  # outfile_rds: Character; The path to an RDS file where the plot will be saved.
  # Dependency: ggplot2
  # Value: None, the function just generates and saves the plot.
  # Author: Niklas Beumer
  
  # Prepare the PCA coordinates and colour categories for plotting.
  plotting_data <- as.data.frame(pca_res$x[, c(1, 2)])
  plotting_data$colour_cat <- colour_categories
  
  # Extract the proportion of explained variance.
  variance_props <- round(summary(pca_res)$importance[2, c(1, 2)], 3)
  variance_percentages <- paste0("(", 100 * variance_props, "%)")
  
  # Generate the plot.
  pca_plot <- ggplot(plotting_data) +
    aes(x = PC1, y = PC2, fill = colour_cat) +
    scale_fill_manual(
      breaks = unique(colour_categories),
      values = colour_palette,
      name = colour_name
    ) +
    geom_point(shape = 21, colour = "black") +
    xlab(paste("Principal component 1", variance_percentages[1])) +
    ylab(paste0("Principal\ncomponent 2\n", variance_percentages[2])) +
    coord_fixed() +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"))
  
  # Save the plot.
  pdf(outfile_pdf, width = 5.5, height = 3)
  print(pca_plot)
  dev.off()
  saveRDS(pca_plot, file = outfile_rds)
  
}


# Specify the colour code for the cell types.
cell_type_col_palette <- c("blue",
                           "cyan2",
                           "orange",
                           "darkorchid1",
                           "red")

# Generate score plots for the different generated PCAs.
prepare_pca_plot(
  pca_res_all,
  colour_categories = extract_cell_types_from_sample_names(rownames(pca_res_all$x)),
  colour_palette = cell_type_col_palette,
  colour_name = "Cell type",
  outfile_pdf = paste0(plot_outdir, "/PCA_plot_samples.pdf"),
  outfile_rds = paste0(plot_rds_outdir, "/PCA_plot_samples.rds")
)










############################################################################
# Generate a three-dimensional plot showing the PCA .
############################################################################

# Extract locations in principal component space.
pc1 <- pca_res$x[, 1]
pc2 <- pca_res$x[, 2]
pc3 <- pca_res$x[, 3]

# Extract the proportion of explained variance and generate axis titles.
variance_props <- round(summary(pca_res)$importance[2, 1:3], 3)
axis.titles <- paste0("PC", 1:3, " (", 100 * variance_props, "%)")

# Specify the colour categories for the dots.
colour_cats <- factor(
  extract_cell_types_from_sample_names(rownames(pca_res$x)),
  levels = unique(extract_cell_types_from_sample_names(rownames(pca_res$x)))
)

# Open the PDF file where the plot will be saved.
outfile_pdf <- paste0(plot_outdir, "/PCA_plot_samples_3d.pdf")
pdf(outfile_pdf, width = 8, height = 4)

# Generate a plot layout with two columns.
layout(matrix(c(1, 2), ncol = 2), widths = c(0.68, 0.32))

# Adjust plot margins.
par(mar = c(1, 3, 1, 1))

# Generate the three-dimensional PCA plot.
scatter3D(
  x = pc1,
  y = pc2,
  z = pc3,
  ticktype = "detailed",
  xlab = axis.titles[1],
  ylab = axis.titles[2],
  zlab = axis.titles[3],
  colvar = as.integer(colour_cats),
  col = cell_type_col_palette,
  pch = 16,
  colkey = F,
  theta = 20,
  phi = 20,
  bty = "b2",
  type = "h",
  expand = 0.3,
  xlim = c(min(pc1) - 5, max(pc1) + 5),
  ylim = c(min(pc2) - 1, max(pc2) + 1),
  zlim = c(min(pc3) - 2, max(pc3) + 2)
)

# Move to the next column of the layout.
plot.new()

# Adjust plot margins.
par(mar = c(0, 0, 0, 0))

# Generate a legend.
legend(
  x = "center",
  legend = levels(colour_cats),
  col = cell_type_col_palette,
  pch = 16,
  title = "Cell type",
  bty = "n",
  title.adj = 0
)

# Close the PDF file.
dev.off()
