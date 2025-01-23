# This script generates a methylation circos plot.
# Author: Niklas Beumer



# Load required packages.
library(bsseq)
library(circlize)
library(ComplexHeatmap)
library(viridis)

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

# Read in global methylation values.
global_meth_aut_file <- paste0(
  in_out_location,
  "/preprocessing_etc/quality_control_after_alignment/2022-03-14_global_methylation_cell_types_autosomes.rds"
)
global_meth_aut <- readRDS(global_meth_aut_file)
global_meth_x_file <- paste0(
  in_out_location,
  "/preprocessing_etc/quality_control_after_alignment/2022-03-14_global_methylation_cell_types_chr_x.rds"
)
global_meth_x <- readRDS(global_meth_x_file)



##################################################################
# Reproduce genomic bins containing 1,000,000 base pairs.
##################################################################

# Get lengths of the chromosomes in the reference genome used for alignment.
chr_lengths <- GenomeInfoDb::seqlengths(
  BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
)[c(as.character(1:22), "X", "MT")]

# Generate Granges specifying 1-Mb bins (1-based coordinates, closed intervals).
# Ensure that the last bin on each chromosome does not extend beyond the end of the chromosome.
bin_length <- 1000000
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
bin_ends <- sapply(
  1:length(bin_starts),
  FUN = function(x) {
    min(bin_starts[x] + bin_length - 1, chr_lengths[chr_names[x]])
  }
)
genome_bins <- GenomicRanges::GRanges(
  seqnames = chr_names,
  ranges = IRanges::IRanges(start = bin_starts, end = bin_ends))






##########################################################################
# Calculate mean methylation in the bins on the cell type level.
##########################################################################

# Extract average methylation values in the bins on the sample level.
binned_methylation_samples <- getMeth(meth_data,
                                      regions = genome_bins,
                                      type = "raw",
                                      what = "perRegion")

# In each bin, compute average methyation estimates across the samples of each cell type.
cell_types <- unique(pData(meth_data)$Cell_type)
binned_methylation_cell_types <- sapply(
  cell_types,
  FUN = function(x) {
    corresp_samples <- rownames(pData(meth_data))[pData(meth_data)$Cell_type == x]
    return(rowMeans(binned_methylation_samples[, corresp_samples], na.rm = T))
  }
)
colnames(binned_methylation_cell_types) <- paste0(colnames(binned_methylation_cell_types), "_Meth")

# In each bin, compute how a cell type deviates from blood naive Tconv cells.
binned_methylation_cell_types_diff <- t(apply(
  binned_methylation_cell_types,
  1,
  FUN = function(x) {
    x - x["Blood naive Tconv_Meth"]
  }
))
colnames(binned_methylation_cell_types_diff) <- paste0(colnames(binned_methylation_cell_types_diff), "_diff")

# Collect all the computed methylation values together with information on the bins.
binned_methylation_final <- data.frame(
  chr = seqnames(genome_bins),
  start = start(genome_bins),
  end = end(genome_bins)
)
binned_methylation_final <- cbind(
  binned_methylation_final,
  as.data.frame(binned_methylation_cell_types),
  as.data.frame(binned_methylation_cell_types_diff)
)

# Remove bins on chromosome MT.
binned_methylation_final <- binned_methylation_final[binned_methylation_final$chr != "MT", ]

# Save these data as a text file.
binned_meth_outfile <- paste0(data_out_pref,
                              "_binned_methylation_values_for_circos_plot.txt")
write.table(
  binned_methylation_final,
  file = binned_meth_outfile,
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)




##################################################################
# Add global methylation values.
##################################################################

# Add information regarding global methylation on autosomes and on chromosome X.
cell_types_w_glob_meth <- sapply(
  cell_types,
  FUN = function(x) {
    glob_meth_aut <- global_meth_aut["mean", x]
    glob_meth_x <- global_meth_x["mean", x]
    glob_meth_str <- paste0(round(glob_meth_aut, 2), ", ", round(glob_meth_x, 2))
    return(paste0(x, " (", glob_meth_str, ")"))
  }
)

# Include this information into the column names of the methylation matrix.
colnames(binned_methylation_final)[4:13] <- sapply(
  colnames(binned_methylation_final)[4:13],
  FUN = function(x) {
    colname_split <- strsplit(x, split = "_")[[1]]
    orig_cell_type_name <- colname_split[1]
    new_cell_type_name <- cell_types_w_glob_meth[orig_cell_type_name]
    new_colname <- paste0(new_cell_type_name,
                          "_",
                          paste(colname_split[2:length(colname_split)], collapse = "_"))
    return(new_colname)
  }
)






##################################################################
# Generate the circos plot.
##################################################################

generate_meth_circos_plot <- function(binned_meth_values,
                                      cell_types,
                                      exclude_diff_track_for,
                                      chr_to_show,
                                      genome,
                                      outfile_pdf) {
  # This function generates a circos plot for methylation values and corresponding differences
  #   in several cell types across the genome.
  # binned_meth_values: Data frame containing the methylation values in bins. Must contain the columns
  #   "chr", "start", "end" as well as a column ">cell type<_Meth" for each cell type and
  #   a column ">cell type<_Meth_diff" for each cell type.
  # cell_types: Character; The collection of cell types to use.
  # exclude_diff_track_for: Character; Cell types for which no methylation difference track
  #   should be plotted.
  # chr_to_show: Character; The names of the chromosomes to show. The chromosome names
  #   must start with "chr".
  # genome: Character; The genome to use (e.g. "hg19").
  # outfile_pdf: Character; the path to a PDF file where the plot will be saved.
  # Depdendencies: circlize, ComplexHeatmap
  # Value: None, the function just saves the plot and prints the extreme meth difference values observed
  #   in the data (To check that the automatically chosen colour scale is reasonable).
  # Author: Niklas Beumer
  
  # Get an idea of the data range.
  meth_columns <- grep("Meth$", colnames(binned_meth_values))
  meth_values <- c(as.matrix(binned_meth_values[, meth_columns]))
  max_meth <- max(meth_values, na.rm = T)
  min_meth <- min(meth_values, na.rm = T)
  meth_diff_columns <- grep("Meth_diff", colnames(binned_meth_values))
  meth_diff_values <- c(as.matrix(binned_meth_values[, meth_diff_columns]))
  max_meth_diff <- max(meth_diff_values, na.rm = T)
  min_meth_diff <- min(meth_diff_values, na.rm = T)
  
  # Generate the colour functions.
  colour_fun_meth <- colorRamp2(
    breaks = seq(min_meth, max_meth, length.out = 200),
    colors = viridis(200, direction = -1)
  )
  colour_fun_meth_diff <- colorRamp2(
    breaks = c(min_meth_diff, 0, max_meth_diff),
    colors = c("red", "white", "darkgoldenrod3")
  )
  
  # Find optimal angle degrees for the separation between chromosomes.
  angle_degrees <- c(rep(1, length(chr_to_show) - 1), 70)
  
  # Open the output device.
  pdf(outfile_pdf, width = 6, height = 6)
  
  # Specify general graphics parameters.
  # Note to myself: Using 89 degrees instead of 90 degrees helps to make the
  # cell type labels reasonably straight.
  circos.par(
    start.degree = 89,
    cell.padding = c(0, 0),
    gap.degree = angle_degrees
  )
  
  # Initialise with cytoband visualisation.
  circos.initializeWithIdeogram(
    species = genome,
    plotType = c("ideogram", "labels"),
    chromosome.index = chr_to_show
  )
  
  # Iterate over all cell types.
  for (cell_type in cell_types) {
    # Extract relevant data.
    cell_type_data <- binned_meth_values[, c(
      "chr",
      "start",
      "end",
      paste0(cell_type, "_Meth"),
      paste0(cell_type, "_Meth_diff")
    )]
    cell_type_data <- cell_type_data[!(is.na(cell_type_data[, paste0(cell_type, "_Meth")])), ]
    
    # Generate a track showing methylation values for this cell type.
    circos.track(
      sectors = cell_type_data$chr,
      track.height = 1 / (3 * length(cell_types)),
      track.margin = c(0, 0.01),
      ylim = c(0, 1),
      panel.fun = function(x, y) {
        x_start <- cell_type_data$start[cell_type_data$chr == CELL_META$sector.index]
        x_end <- cell_type_data$end[cell_type_data$chr == CELL_META$sector.index]
        colours <- colour_fun_meth(cell_type_data[cell_type_data$chr == CELL_META$sector.index, paste0(cell_type, "_Meth")])
        circos.rect(x_start, 0, x_end, 1, col = colours, border = NA)
      }
    )
    
    # Add a x axis.
    circos.yaxis(
      side = "left",
      at = 0.3,
      labels = cell_type,
      tick = F,
      sector.index = chr_to_show[1],
      labels.cex = 0.8,
      tick.length = 0
    )
    
    # Generate a track showing methylation differences for this cell type.
    # Do not do this if the user specified the current cell types under "exclude_diff_track_for".
    if (!(cell_type %in% exclude_diff_track_for)) {
      circos.track(
        sectors = cell_type_data$chr,
        ylim = c(0, 1),
        track.margin = c(0, 0),
        track.height = 1 / (5 * length(cell_types)),
        panel.fun = function(x, y) {
          x_start <- cell_type_data$start[cell_type_data$chr == CELL_META$sector.index]
          x_end <- cell_type_data$end[cell_type_data$chr == CELL_META$sector.index]
          colours <- colour_fun_meth_diff(cell_type_data[cell_type_data$chr == CELL_META$sector.index, paste0(cell_type, "_Meth_diff")])
          circos.rect(x_start, 0, x_end, 1, col = colours, border = NA)
        }
      )
    }
    
  }
  
  # Generate the colour legends.
  meth_legend <- Legend(
    col_fun = colour_fun_meth,
    title = "Mean\nmethylation",
    at = seq(0, 1, 0.2),
    legend_height = unit(0.14, "npc")
  )
  meth_diff_legend <- Legend(
    col_fun = colour_fun_meth_diff,
    title = "Mean methylation\ndifference",
    legend_height = unit(0.14, "npc")
  )
  draw(
    meth_legend,
    x = unit(0.1, "npc"),
    y = unit(0.23, "npc"),
    just = "top"
  )
  draw(
    meth_diff_legend,
    x = unit(0.9, "npc"),
    y = unit(0.23, "npc"),
    just = "top"
  )
  
  # Print the extremes of the difference scale in order to check whether the colour bar indeed
  # covers the whole range that needs to be covered.
  print("Extremes of the methylation difference colour scale:")
  print(min_meth_diff)
  print(max_meth_diff)
  
  # Close the output device.
  dev.off()
  
}


# Prepare the chromosome names so that they can be used by the circos plot function.
binned_methylation_final$chr <- paste0("chr", binned_methylation_final$chr)

# Generate the circos plot.
# Do not include a mthylation difference track for blood naive Tconvs (because this
# cell type was the base line in methylation difference computation).
generate_meth_circos_plot(
  binned_methylation_final,
  cell_types = cell_types_w_glob_meth,
  exclude_diff_track_for = grep("Blood naive Tconv", cell_types_w_glob_meth, value = T),
  chr_to_show = paste0("chr", c(1:22, "X")),
  genome = "hg19",
  outfile_pdf = paste0(plot_outdir, "/methylation_circos_plot.pdf")
)
