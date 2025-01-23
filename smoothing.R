# This script smoothes the methylation data and saves the BSseq object with the
#   smoothing results.
# Author: Niklas Beumer.



# Load required packages.
library(bsseq)
library(BiocParallel)

# Define a location where input is saved and where output will be saved.
in_out_location <- "/yyy/hm_treg_bs_rgnsbg/preprocessing_etc/quality_control_after_alignment"

# Define the prefix of the output data file.
data_out_pref <- paste(in_out_location, format(Sys.time(), "%Y-%m-%d"), sep = "/")

# Read in the BSseq objects containing the methylation data.
meth_data <- readRDS(
  paste(
    in_out_location,
    "2022-03-11_bsseq_object_combined_all_samples.rds",
    sep = "/"
  )
)
meth_data_reduced <- readRDS(
  paste(
    in_out_location,
    "2022-03-11_bsseq_object_combined_all_samples_reduced.rds",
    sep = "/"
  )
)




##################################
# Smooth data for single samples.
##################################

smooth_and_save_meth_data <- function(bsseq_obj, outfile_rds) {
  # This function smoothes methylation data stored in a BSseq object and saves the smoothed object as an RDS file.
  # bsseq_obj: The BSseq object to smooth.
  # outfile_rds: Character; The path to an RDS file where the smoothed BSseq object will be saved.
  # Dependencies: bsseq, BiocParallel
  # Value: None, the function just generates and saves the smoothed BSseq object.
  # Author: Niklas Beumer
  
  # Smooth the methylation data.
  meth_data_smoothed <- BSmooth(bsseq_obj, BPPARAM = MulticoreParam(workers = ncol(bsseq_obj), progressbar = T))
  
  # Save the smoothed BSseq object.
  saveRDS(meth_data_smoothed, file = outfile_rds)
  
}


# Smooth the methylation data and save the smoothed objects.
smooth_and_save_meth_data(
  meth_data,
  outfile_rds = paste0(
    data_out_pref,
    "_bsseq_object_combined_all_samples_smoothed.rds"
  )
)
smooth_and_save_meth_data(
  meth_data_reduced,
  outfile_rds = paste0(
    data_out_pref,
    "_bsseq_object_combined_all_samples_reduced_smoothed.rds"
  )
)





####################################################
# Smooth the full data set on the cell type level.
####################################################

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

# Collapse samples from the same cell type into joint columns.
cell_types_in_columns <- extract_cell_types_from_sample_names(colnames(meth_data))
meth_data_collapsed <- collapseBSseq(meth_data, group = cell_types_in_columns)

# Smooth the collapsed data and save the smoothed BSseq object.
smooth_and_save_meth_data(
  meth_data_collapsed,
  outfile_rds = paste0(
    data_out_pref,
    "_bsseq_object_combined_all_samples_collapsed_smoothed.rds"
  )
)
