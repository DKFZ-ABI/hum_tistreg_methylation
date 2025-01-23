# This R script looks at regions that are differentially methylated between
#   skin Treg cells and blood naive Treg cells and analyses where blood CCR8+
#   Treg cells are positioned in this comparison.
# Author: Niklas Beumer



# Load required packages.
library(GenomicRanges)
library(bsseq)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(ggplot2)


# Define a location on /yyy.
location <- "/yyy/nbeumer/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis",
                     format(Sys.time(), "%m-%d-%y"),
                     sep = "/")
if (!dir.exists(plot_outdir)) {
  dir.create(plot_outdir)
}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Specify the relevant cell types.
relevant_cell_types <- c("Blood naive Treg", "Blood CCR8+ Treg", "Skin Treg")

# Read in the list of differentially methylated regions between skin Tregs
# and blood naive Tregs.
meth_reg_file <- paste0(
  location,
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_gap_0.15_sterrratio_0.5_signature_regions_filtered.txt"
)
meth_reg <- read.table(
  meth_reg_file,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)
# meth_reg <- meth_reg[1:50000, ] # For debugging purposes.

# Read in the BSseq object containing the smoothed methylation data
# and restrict to the relevant cell types.
meth_data <- readRDS(
  paste(
    location,
    "/preprocessing_etc/quality_control_after_alignment/2022-03-14_bsseq_object_combined_all_samples_smoothed.rds",
    sep = "/"
  )
)
meth_data <- meth_data[, pData(meth_data)$Cell_type %in% relevant_cell_types]





########################################################################
# Compute cell-type-wise average raw methylation values in the relevant
# regions.
########################################################################

# Generate a GRanges object for the differentially methylated regions.
meth_reg_gr <- makeGRangesFromDataFrame(meth_reg, keep.extra.columns = T)

# Identify average raw methylation values (on the sample level).
avg_meth_sample_level <- getMeth(meth_data,
                                 regions = meth_reg_gr,
                                 type = "raw",
                                 what = "perRegion")

# Compute cell-type-level averages of these sample-level methylation values.
avg_meth_celltype_level <- sapply(
  relevant_cell_types,
  FUN = function(x) {
    celltype_subs <- avg_meth_sample_level[, grep(x, colnames(avg_meth_sample_level), fixed = T)]
    return(rowMeans(celltype_subs, na.rm = T))
  }
)



##########################################################################################
# Exclude regions for which the differential methylation tendencies reported by my multi-
# class approach do not match with what I see in my raw methylation values.
# Note: This does not exclude any regions, which also makes sense since I included
# this criterion in the methylation calling procedure. However, the code block
# is still retained for consistency with the ATAC and RNA level and to highlight
# in the output data that no regions had to be excluded.
###########################################################################################

# For each region, check whether the differential methylation tendency reported by my
# multi-class approach matches with what I see in the raw methylation values.
exclude_bools <- sapply(
  1:nrow(meth_reg),
  FUN = function(x) {
    tend <- meth_reg$automatic_annotation[x]
    skin_meth <- avg_meth_celltype_level[x, "Skin Treg"]
    blood_naive_meth <- avg_meth_celltype_level[x, "Blood naive Treg"]
    return(
      ifelse(
        tend == "Skin_Treg__hypomethylation",
        yes = skin_meth > blood_naive_meth,
        no = skin_meth < blood_naive_meth
      )
    )
  }
)

# Identify regions to exclude.
regions_to_exclude <- meth_reg[exclude_bools, 1:3]

# Save information on which regions were excluded.
output_file <- paste0(
  location,
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_region_exclusions.txt"
)
if (file.exists(output_file)) {
  file.remove(output_file)
}
capture.output(cat(
  rep(
    "*************************************************************************************************\n",
    2
  )
), file = output_file, append = T)
capture.output(
  print(
    "Regions that were excluded because the diff. meth. tendencies reported by my"
  ),
  file = output_file,
  append = T
)
capture.output(
  print(
    "multi-class approach did not match with what I saw in the raw methylation values:"
  ),
  file = output_file,
  append = T
)
capture.output(cat(
  rep(
    "*************************************************************************************************\n",
    2
  )
), file = output_file, append = T)
capture.output(cat("\n\n"), file = output_file, append = T)
if (nrow(regions_to_exclude) > 0) {
  capture.output(print(regions_to_exclude),
                 file = output_file,
                 append = T)
} else {
  capture.output(print("No regions were excluded."),
                 file = output_file,
                 append = T)
}

# Perform the exclusions.
avg_meth_sample_level <- avg_meth_sample_level[!exclude_bools, ]
avg_meth_celltype_level <- avg_meth_celltype_level[!exclude_bools, ]
meth_reg <- meth_reg[!exclude_bools, ]






##############################################################################
# Scale the methylation data so that they range from 0 to 1
# for each region.
##############################################################################

scale_matr_lines_betw_0_and_1 <- function(matr) {
  # This function scales the lines of a matrix so that values range between
  #   0 and 1.
  # matr: The matrix to scale.
  # Dependencies: none.
  # Value: The scaled matrix.
  # Author: Niklas Beumer
  
  # Iterate over the rows.
  scaled_matr <- t(apply(
    matr,
    1,
    FUN = function(x) {
      # Identify the minimum and maximum value.
      min_val <- min(x)
      max_val <- max(x)
      
      # Subtract the minimum.
      x <- x - min_val
      
      # Divide values by the difference between the maximum and the minimum.
      divisor <- max_val - min_val
      x <- x / divisor
      
      # Return the scaled values.
      return(x)
      
    }
  ))
  
  # Return the scaled matrix.
  return(scaled_matr)
  
}

# Scale the values between 0 and 1.
avg_meth_sample_level_scaled <-
  scale_matr_lines_betw_0_and_1(avg_meth_sample_level)
avg_meth_celltype_level_scaled <-
  scale_matr_lines_betw_0_and_1(avg_meth_celltype_level)






#########################################################################
# For each region, compute whether blood CCR8+ Tregs are closer to
# skin Tregs or to blood naive Tregs (based on cell-type-level
# methylation values).
#########################################################################

# Quantify distances in scaled methylation between blood CCR8+ Tregs and
# the two extreme cell types.
ccr8_naive_diffs <- avg_meth_celltype_level_scaled[, "Blood naive Treg"] -
  avg_meth_celltype_level_scaled[, "Blood CCR8+ Treg"]
ccr8_naive_dists <- abs(ccr8_naive_diffs)
ccr8_skin_diffs <- avg_meth_celltype_level_scaled[, "Skin Treg"] -
  avg_meth_celltype_level_scaled[, "Blood CCR8+ Treg"]
ccr8_skin_dists <- abs(ccr8_skin_diffs)

# For each region, determine whether blood CCR8+ Tregs are closer to
# blood naive Tregs or closer to skin Tregs.
closest_extrema <- sapply(
  1:nrow(avg_meth_celltype_level_scaled),
  FUN = function(x) {
    if (ccr8_naive_dists[x] < ccr8_skin_dists[x]) {
      return("Blood naive Treg")
    } else if (ccr8_naive_dists[x] > ccr8_skin_dists[x]) {
      return("Skin Treg")
    } else {
      return("Exactly in middle")
    }
  }
)





#########################################################################
# Collect all the computed values.
#########################################################################

# Prepare methylation values.
colnames(avg_meth_celltype_level) <-
  paste0("Raw_meth__", colnames(avg_meth_celltype_level))
raw_meth_df <- as.data.frame(avg_meth_celltype_level)
colnames(avg_meth_celltype_level_scaled) <-
  paste0("Scaled_meth__", colnames(avg_meth_celltype_level_scaled))
scaled_meth_df <- as.data.frame(avg_meth_celltype_level_scaled)

# Prepare results from distance analysis to extreme cell types.
dist_df <- data.frame(
  Scaled_meth_diff_blood_ccr8_treg_blood_naive_treg = ccr8_naive_diffs,
  Scaled_meth_dist_blood_ccr8_treg_blood_naive_treg = ccr8_naive_dists,
  Scaled_meth_diff_blood_ccr8_treg_skin_treg = ccr8_skin_diffs,
  Scaled_meth_dist_blood_ccr8_treg_skin_treg = ccr8_skin_dists
)
closest_extrema_df <- data.frame(Cell_type_closest_to_blood_ccr8_treg  = closest_extrema)

# Combine all information into a large data frame.
combined_results <- cbind(meth_reg,
                          raw_meth_df,
                          scaled_meth_df,
                          dist_df,
                          closest_extrema_df)

# Save all computed values.
combined_results_outfile <- paste0(
  location,
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_meth_discrep_regions_info.txt"
)
write.table(
  combined_results,
  file = combined_results_outfile,
  sep = "\t",
  row.names = F,
  col.names = T,
  quote = F
)





############################################################################
# Plot heat maps showing scaled methylation in the differentially
# methylated regions. Once, use cell-type-level values and once use
# sample-level values.
############################################################################

# Generate a version of the methylation information, ordered by the difference between
# blood CCR8+ Tregs and skin Tregs.
order_to_use <- order(
  combined_results$Scaled_meth_diff_blood_ccr8_treg_skin_treg *
    sapply(
      combined_results$automatic_annotation,
      FUN = function(x) {
        switch(
          x,
          Blood_naive_Treg__hypomethylation = -1,
          Skin_Treg__hypomethylation = 1
        )
      }
    )
)
combined_results_ordered <- combined_results[order_to_use, ]
avg_meth_sample_level_scaled_ordered <- avg_meth_sample_level_scaled[order_to_use, ]


###############
# Heat map on the cell type level.

# Iterate over two large categories (decreasing methylation during differentiation,
# increasing methylation during differentiation).
for (category in c("Blood_naive_Treg__hypomethylation",
                   "Skin_Treg__hypomethylation")) {
  # Extract values for this category.
  combined_results_cat <- combined_results_ordered[combined_results_ordered$automatic_annotation == category, ]
  
  # Get the matrix for heatmap plotting.
  meth_matr <- as.matrix(combined_results_cat[, grep("Scaled_meth__", colnames(combined_results_cat))])
  colnames(meth_matr) <- gsub("Scaled_meth__", "", colnames(meth_matr))
  rownames(meth_matr) <- combined_results_cat$region_ID
  
  # Generate the colour function for the heatmaps.
  meth_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                             colors = viridis(200, direction = -1))
  
  # Generate a row annotation showing to which cell type blood CCR8+ Tregs are closest
  # for each region.
  combined_results_cat$anno_values <- sapply(
    combined_results_cat$Cell_type_closest_to_blood_ccr8_treg,
    FUN = function(x) {
      ifelse(x == "Skin Treg", yes = "Closer to skin Tregs", no = "Closer to blood naive Tregs")
    }
  )
  row_anno <- rowAnnotation(
    `Blood CCR8+ Treg position` =
      combined_results_cat$anno_values,
    col = list(
      `Blood CCR8+ Treg position` = c(
        `Closer to skin Tregs` = "blue",
        `Closer to blood naive Tregs` = "darkorchid1",
        `Exactly in middle` = "grey"
      )
    )
  )
  
  # Generate the heat map.
  meth_heatmap <- Heatmap(
    meth_matr,
    show_row_names = F,
    cluster_rows = F,
    cluster_columns = F,
    col = meth_col_fun,
    right_annotation = row_anno,
    heatmap_legend_param = list(
      title = "Methylation",
      at = c(0, 1),
      labels = c("min", "max")
    ),
    column_names_max_height = unit(10, "cm")
  )
  
  
  # Save the heat map.
  meth_heatmap_pdf <- paste0(
    plot_outdir,
    "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_meth_discrep_heatmap_",
    category,
    ".pdf"
  )
  meth_heatmap_rds <- paste0(
    plot_rds_outdir,
    "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_meth_discrep_heatmap_",
    category,
    ".rds"
  )
  pdf(meth_heatmap_pdf, width = 6, height = 10)
  draw(meth_heatmap)
  dev.off()
  saveRDS(meth_heatmap, file = meth_heatmap_rds)
  
}


###############
# Heat map on the sample level.

# Iterate over two large categories (decreasing methylation during differentiation,
# increasing methylation during differentiation).
for (category in c("Blood_naive_Treg__hypomethylation",
                   "Skin_Treg__hypomethylation")) {
  # Extract values for this category.
  meth_matr <- avg_meth_sample_level_scaled_ordered[combined_results_ordered$automatic_annotation == category, ]
  rownames(meth_matr) <- combined_results$region_ID[combined_results_ordered$automatic_annotation == category]
  
  # Re-order columns according to how they should be displayed in the
  # publication.
  cell_type_order <- c("Blood naive Treg", "Blood CCR8+ Treg", "Skin Treg")
  meth_matr_ordered <- do.call(cbind, lapply(
    cell_type_order,
    FUN = function(x) {
      meth_matr[, grep(x, colnames(meth_matr), fixed = T)]
    }
  ))
  
  # Generate the colour function for the heatmaps.
  meth_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                             colors = viridis(200, direction = -1))
  
  # Generate a row annotation showing to which cell type blood CCR8+ Tregs are closest
  # for each region.
  anno_values <- sapply(
    combined_results_ordered$Cell_type_closest_to_blood_ccr8_treg[combined_results_ordered$automatic_annotation == category],
    FUN = function(x) {
      ifelse(x == "Skin Treg", yes = "Closer to skin Tregs", no = "Closer to blood naive Tregs")
    }
  )
  row_anno <- rowAnnotation(
    `Blood CCR8+ Treg position` = anno_values,
    col = list(
      `Blood CCR8+ Treg position` = c(
        `Closer to skin Tregs` = "blue",
        `Closer to blood naive Tregs` = "darkorchid1",
        `Exactly in middle` = "grey"
      )
    )
  )
  
  # Generate the heat map.
  meth_heatmap <- Heatmap(
    meth_matr_ordered,
    show_row_names = F,
    cluster_rows = F,
    cluster_columns = F,
    col = meth_col_fun,
    right_annotation = row_anno,
    heatmap_legend_param = list(
      title = "Methylation",
      at = c(0, 1),
      labels = c("min", "max")
    ),
    column_names_max_height = unit(10, "cm")
  )
  
  
  # Save the heat map.
  meth_heatmap_pdf <- paste0(
    plot_outdir,
    "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_meth_discrep_heatmap_",
    category,
    "_by_sample.pdf"
  )
  meth_heatmap_rds <- paste0(
    plot_rds_outdir,
    "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_meth_discrep_heatmap_",
    category,
    "_by_sample.rds"
  )
  pdf(meth_heatmap_pdf, width = 6, height = 10)
  draw(meth_heatmap)
  dev.off()
  saveRDS(meth_heatmap, file = meth_heatmap_rds)
  
}





#####################################################################
# Plot histograms of the differences between blood CCR8+ Tregs
# and skin Tregs.
#####################################################################

# Iterate over two large categories (decreasing methylation during diffrentiation,
# increasing methylation during differentiation).
for (category in c("Blood_naive_Treg__hypomethylation",
                   "Skin_Treg__hypomethylation")) {
  # Extract differences between skin Tregs and blood CCR8+ Tregs.
  diffs <- combined_results$Scaled_meth_diff_blood_ccr8_treg_skin_treg[combined_results$automatic_annotation == category]
  
  # Generate the histogram.
  histogram <- ggplot() +
    aes(x = diffs) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_histogram(binwidth = 0.02, fill = "black") +
    xlab("Scaled meth. skin Tregs minus scaled meth. blood CCR8+ Tregs") +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"))
  
  # Save the histogram.
  histogram_pdf <- paste0(
    plot_outdir,
    "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_meth_discrep_meth_diff_hist_",
    category,
    ".pdf"
  )
  histogram_rds <- paste0(
    plot_rds_outdir,
    "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_meth_discrep_meth_diff_hist_",
    category,
    ".rds"
  )
  pdf(histogram_pdf, width = 5, height = 4)
  print(histogram)
  dev.off()
  saveRDS(histogram, file = histogram_rds)
}





############################################################################################
# Count in how many regions blood CCR8+ Tregs are closer to skin Tregs / blood naive Tregs.
############################################################################################

# Specify the file name under which output from this section will be saved.
output_file <- paste0(
  location,
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_meth_discrep_regions_counts.txt"
)
file.remove(output_file)

# Print a status message.
capture.output(cat(
  rep(
    "*************************************************************************************************\n",
    2
  )
), file = output_file, append = T)
capture.output(
  print(
    "Numbers of regions in which blood CCR8+ Tregs are closer to skin Tregs / blood naive Tregs:"
  ),
  file = output_file,
  append = T
)
capture.output(cat(
  rep(
    "*************************************************************************************************\n",
    2
  )
), file = output_file, append = T)
capture.output(cat("\n\n"), file = output_file, append = T)

# Iterate over two large categories (decreasing methylation during diffrentiation,
# increasing methylation during differentiation).
for (category in c("Blood_naive_Treg__hypomethylation",
                   "Skin_Treg__hypomethylation")) {
  # Print what category is currently quantified.
  capture.output(
    cat("############################################\n"),
    file = output_file,
    append = T
  )
  capture.output(print(category), file = output_file, append = T)
  capture.output(cat("\n"), file = output_file, append = T)
  
  # Restrict the data to those regions that correspond to the current category.
  combined_results_cat <- combined_results[combined_results$automatic_annotation == category, ]
  
  # Compute some statistics for all regions.
  total_reg_count <- nrow(combined_results_cat)
  total_cpg_count <- sum(combined_results_cat$cpg_num)
  total_base_count <- sum(combined_results_cat$end - combined_results_cat$start + 1)
  
  # Compute statistics for the identified region categories.
  closer_to_naive_reg_count <- length(
    which(
      combined_results_cat$Cell_type_closest_to_blood_ccr8_treg == "Blood naive Treg"
    )
  )
  closer_to_naive_reg_perc <- round(100 * closer_to_naive_reg_count / total_reg_count, 2)
  closer_to_naive_cpg_count <- sum(combined_results_cat$cpg_num[combined_results_cat$Cell_type_closest_to_blood_ccr8_treg == "Blood naive Treg"])
  closer_to_naive_cpg_perc <- round(100 * closer_to_naive_cpg_count / total_cpg_count, 2)
  closer_to_naive_base_count <- sum((combined_results_cat$end - combined_results_cat$start + 1)[combined_results_cat$Cell_type_closest_to_blood_ccr8_treg == "Blood naive Treg"])
  closer_to_naive_base_perc <- round(100 * closer_to_naive_base_count / total_base_count, 2)
  closer_to_skin_reg_count <- length(which(
    combined_results_cat$Cell_type_closest_to_blood_ccr8_treg == "Skin Treg"
  ))
  closer_to_skin_reg_perc <- round(100 * closer_to_skin_reg_count / total_reg_count, 2)
  closer_to_skin_cpg_count <- sum(combined_results_cat$cpg_num[combined_results_cat$Cell_type_closest_to_blood_ccr8_treg == "Skin Treg"])
  closer_to_skin_cpg_perc <- round(100 * closer_to_skin_cpg_count / total_cpg_count, 2)
  closer_to_skin_base_count <- sum((combined_results_cat$end - combined_results_cat$start + 1)[combined_results_cat$Cell_type_closest_to_blood_ccr8_treg == "Skin Treg"])
  closer_to_skin_base_perc <- round(100 * closer_to_skin_base_count / total_base_count, 2)
  
  # Print the results.
  capture.output(print(
    paste0(
      "Closer to blood naive Tregs: ",
      closer_to_naive_reg_count,
      " regions (",
      closer_to_naive_reg_perc,
      "%), ",
      closer_to_naive_cpg_count,
      " CpGs (",
      closer_to_naive_cpg_perc,
      "%), ",
      closer_to_naive_base_count,
      " bases (",
      closer_to_naive_base_perc,
      "%)."
    )
  ), file = output_file, append = T)
  capture.output(print(
    paste0(
      "Closer to skin Tregs: ",
      closer_to_skin_reg_count,
      " regions (",
      closer_to_skin_reg_perc,
      "%), ",
      closer_to_skin_cpg_count,
      " CpGs (",
      closer_to_skin_cpg_perc,
      "%), ",
      closer_to_skin_base_count,
      " bases (",
      closer_to_skin_base_perc,
      "%)."
    )
  ), file = output_file, append = T)
  capture.output(cat("\n\n"), file = output_file, append = T)
  
}
