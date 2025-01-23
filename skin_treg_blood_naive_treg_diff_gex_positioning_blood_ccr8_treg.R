# This R script looks at genes that are differentially expressed between skin
#   Treg cells and blood naive Treg cells and analyses where blood CCR8+ Treg
#   cells are positioned in this comparison.
# Author: Niklas Beumer



# Load required packages.
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(ggplot2)
library(GenomicRanges)
library(testit)


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

# Read in the sample mapping file.
rna_sample_mapping_path <- paste0(location, 
                                  "/sample_mapping_rnaseq_only_biol_rep.txt")
rna_sample_mapping <- read.table(
  rna_sample_mapping_path,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)

# Specify the relevant cell types.
relevant_cell_types <- c("Blood naive Treg", "Blood CCR8+ Treg", "Skin Treg")

# Read in the list of differentially expressed genes between skin Tregs and
# blood naive Tregs.
diff_exp_file <- paste0(
  location,
  "/RNASeq/analysis_results/2022-01-14_diff_gene_expr_DESEq2_Skin_TregBlood_naive_Treg_results_filtered_with_significance.txt"
)
diff_exp <- read.table(diff_exp_file,
                       header = T,
                       stringsAsFactors = F)
diff_exp <- diff_exp[diff_exp$significant, ]

# Read in the TPM values.
tpm_file <- paste0(location, "/RNASeq/tpm/tpm_all.txt")
tpm <- read.table(tpm_file, header = T, stringsAsFactors = F)









########################################################################
# Compute cell-type-wise average TPM values for the relevant genes.
########################################################################

# Restrict the TPM data to the differentially expressed genes.
tpm_diff_exp <- tpm[tpm$Gene_symbol %in% rownames(diff_exp), ]

# For each cell type and each differential gene, compute the mean expression across the cells.
#-- Iterate over the three relevant cell types.
avg_tpm_celltype_level <- sapply(
  relevant_cell_types,
  FUN = function(x) {
    #-- Identify the samples from this cell type.
    corresp_samples <- rna_sample_mapping$Sample[rna_sample_mapping$Cell_type == x]
    #-- Extract TPM values for these samples.
    corresp_tpm <- tpm_diff_exp[, corresp_samples]
    rownames(corresp_tpm) <- tpm_diff_exp$Gene_symbol
    #-- Transform the TPM values by log(100 * x + 1)
    for (colname in colnames(corresp_tpm)) {
      corresp_tpm[, colname] <- log1p(100 * corresp_tpm[, colname])
    }
    #-- Compute within-cell-type means.
    means <- rowMeans(corresp_tpm)
    #-- Return the means in the order in which the corresponding genes appear in the gene expression table.
    return(means[rownames(diff_exp)])
  }
)
assert(all(rownames(avg_tpm_celltype_level) == rownames(diff_exp)))

# Also collect the sample-level TPM values.
tpm_sample_level <- do.call(cbind, lapply(
  relevant_cell_types,
  FUN = function(x) {
    #-- Identify the samples from this cell type.
    corresp_samples <- rna_sample_mapping$Sample[rna_sample_mapping$Cell_type == x]
    #-- Extract TPM values for these samples.
    corresp_tpm <- tpm_diff_exp[, corresp_samples]
    rownames(corresp_tpm) <- tpm_diff_exp$Gene_symbol
    #-- Transform the TPM values by log(100 * x + 1)
    for (colname in colnames(corresp_tpm)) {
      corresp_tpm[, colname] <- log1p(100 * corresp_tpm[, colname])
    }
    #-- Return the sample-level values in the order in which the corresponding
    #-- genes appear in the gene expression table.
    return(corresp_tpm[rownames(diff_exp), ])
  }
))
assert(all(rownames(tpm_sample_level) == rownames(diff_exp)))





##########################################################################################
# Exclude genes for which the differential expression tendencies reported by DESeq2
# do not match with what I see in my log-transformed TPM values. These cases are likely
# due to differences in normalisation between DESeq2 and my log-transformed TPMs.
###########################################################################################

# For each gene, check whether the differential expression tendency reported by DESeq2
# matches with what I see in the log-transformed TPM values.
exclude_bools <- sapply(
  1:nrow(diff_exp),
  FUN = function(x) {
    diff_tend_sign <- sign(diff_exp$log2FoldChange[x])
    skin_log_tpm <- avg_tpm_celltype_level[x, "Skin Treg"]
    blood_naive_log_tpm <- avg_tpm_celltype_level[x, "Blood naive Treg"]
    return(
      ifelse(
        diff_tend_sign == -1,
        yes = skin_log_tpm < blood_naive_log_tpm,
        no = skin_log_tpm > blood_naive_log_tpm
      )
    )
  }
)

# Identify genes to exclude.
genes_to_exclude <- rownames(diff_exp)[exclude_bools]

# Save information on which genes were excluded.
output_file <- paste0(
  location,
  "/treg_hierarchies/diff_expr_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_gene_exclusions.txt"
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
    "Genes that were excluded because the diff. expr. tendencies reported by DESeq2"
  ),
  file = output_file,
  append = T
)
capture.output(
  print("did not match with what I saw in the TPM data:"),
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
capture.output(print(genes_to_exclude), file = output_file, append = T)

# Perform the exclusions.
avg_tpm_celltype_level_filtered <- avg_tpm_celltype_level[!exclude_bools, ]
tpm_sample_level_filtered <- tpm_sample_level[!exclude_bools, ]
diff_exp_filtered <- diff_exp[!exclude_bools, ]




##############################################################################
# Scale the gene expression data so that they range from 0 to 1
# for each gene.
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
avg_tpm_celltype_level_filtered_scaled <-
  scale_matr_lines_betw_0_and_1(avg_tpm_celltype_level_filtered)
tpm_sample_level_filtered_scaled <-
  scale_matr_lines_betw_0_and_1(tpm_sample_level_filtered)






#########################################################################
# For each gene, compute whether blood CCR8+ Tregs are closer to
# skin Tregs or to blood naive Tregs.
#########################################################################

# Quantify distances in scaled expression between blood CCR8+ Tregs and
# the two extreme cell types.
ccr8_naive_diffs <- avg_tpm_celltype_level_filtered_scaled[, "Blood naive Treg"] -
  avg_tpm_celltype_level_filtered_scaled[, "Blood CCR8+ Treg"]
ccr8_naive_dists <- abs(ccr8_naive_diffs)
ccr8_skin_diffs <- avg_tpm_celltype_level_filtered_scaled[, "Skin Treg"] -
  avg_tpm_celltype_level_filtered_scaled[, "Blood CCR8+ Treg"]
ccr8_skin_dists <- abs(ccr8_skin_diffs)

# For each gene, determine whether blood CCR8+ Tregs are closer to
# blood naive Tregs or closer to skin Tregs.
closest_extrema <- sapply(
  1:nrow(avg_tpm_celltype_level_filtered_scaled),
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

# Prepare gene names.
gene_names_df <- data.frame(Gene = rownames(diff_exp_filtered))

# Prepare gene expression values.
colnames(avg_tpm_celltype_level_filtered) <-
  paste0("Log_TPM__", colnames(avg_tpm_celltype_level_filtered))
tpm_df <- as.data.frame(avg_tpm_celltype_level_filtered)
colnames(avg_tpm_celltype_level_filtered_scaled) <-
  paste0("Scaled_log_TPM__",
         colnames(avg_tpm_celltype_level_filtered_scaled))
scaled_tpm_df <- as.data.frame(avg_tpm_celltype_level_filtered_scaled)

# Prepare results from distance analysis to extreme cell types.
dist_df <- data.frame(
  Scaled_log_TPM_diff_blood_ccr8_treg_blood_naive_treg = ccr8_naive_diffs,
  Scaled_log_TPM_dist_blood_ccr8_treg_blood_naive_treg = ccr8_naive_dists,
  Scaled_log_TPM_diff_blood_ccr8_treg_skin_treg = ccr8_skin_diffs,
  Scaled_log_TPM_dist_blood_ccr8_treg_skin_treg = ccr8_skin_dists
)
closest_extrema_df <- data.frame(Cell_type_closest_to_blood_ccr8_treg  = closest_extrema)

# Combine all information into a large data frame.
combined_results <- cbind(
  gene_names_df,
  diff_exp_filtered,
  tpm_df,
  scaled_tpm_df,
  dist_df,
  closest_extrema_df
)

# Save all computed values.
combined_results_outfile <- paste0(
  location,
  "/treg_hierarchies/diff_expr_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_expr_discrep_genes_info.txt"
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
# Plot heat maps showing scaled expression of the differentially expressed
# genes. Once, use cell-type-level values and once use sample-level values.
############################################################################

# Generate a version of the expression information, ordered by the difference between
# blood CCR8+ Tregs and skin Tregs.
order_to_use <- order(
  combined_results$Scaled_log_TPM_diff_blood_ccr8_treg_skin_treg *
    sapply(
      combined_results$log2FoldChange,
      FUN = function(x) {
        ifelse(x < 0, yes =  -1, no = 1)
      }
    )
)
combined_results_ordered <- combined_results[order_to_use, ]
tpm_sample_level_filtered_scaled_ordered <- tpm_sample_level_filtered_scaled[order_to_use, ]


###############
# Heat map on the cell type level.

# Iterate over two large categories (decreasing expression during differentiation,
# increasing expression during differentiation).
for (category in c("<", ">")) {
  # Extract values for this category.
  combined_results_cat <- combined_results_ordered[eval(parse(
    text = paste0("combined_results_ordered$log2FoldChange", category, " 0")
  )), ]
  
  # Get the matrix for heatmap plotting.
  expr_matr <- as.matrix(combined_results_cat[, grep("Scaled_log_TPM__", colnames(combined_results_cat))])
  colnames(expr_matr) <- gsub("Scaled_log_TPM__", "", colnames(expr_matr))
  rownames(expr_matr) <- combined_results_cat$Gene
  
  # Generate the colour function for the heatmaps.
  expr_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                             colors = viridis(200, option = "G"))
  
  # Generate a row annotation showing to which cell type blood CCR8+ Tregs are closest
  # for each gene.
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
  expr_heatmap <- Heatmap(
    expr_matr,
    show_row_names = F,
    cluster_rows = F,
    cluster_columns = F,
    col = expr_col_fun,
    right_annotation = row_anno,
    heatmap_legend_param = list(
      title = "Expression",
      at = c(0, 1),
      labels = c("min", "max")
    ),
    column_names_max_height = unit(10, "cm"),
  )
  
  
  # Save the heat map.
  file_snip <- ifelse(category == "<", yes = "Skin_Treg__upregulation", no = "Blood_naive_Treg__upregulation")
  expr_heatmap_pdf <- paste0(
    plot_outdir,
    "/treg_hierarchies_diff_expr_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_heatmap_filtered_f_expr_discrep_",
    file_snip,
    ".pdf"
  )
  expr_heatmap_rds <- paste0(
    plot_rds_outdir,
    "/treg_hierarchies_diff_expr_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_heatmap_filtered_f_expr_discrep_",
    file_snip,
    ".rds"
  )
  pdf(expr_heatmap_pdf, width = 6, height = 10)
  draw(expr_heatmap)
  dev.off()
  saveRDS(expr_heatmap, file = expr_heatmap_rds)
  
}


###############
# Heat map on the sample level.

# Iterate over two large categories (decreasing expression during differentiation,
# increasing expression during differentiation).
for (category in c("<", ">")) {
  # Extract values for this category.
  expr_matr <- tpm_sample_level_filtered_scaled_ordered[eval(parse(
    text = paste0("combined_results_ordered$log2FoldChange", category, " 0")
  )), ]
  
  # Generate the colour function for the heatmaps.
  expr_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                             colors = viridis(200, option = "G"))
  
  # Generate a row annotation showing to which cell type blood CCR8+ Tregs are closest
  # for each gene.
  anno_values <- sapply(
    combined_results_ordered$Cell_type_closest_to_blood_ccr8_treg[eval(parse(
      text = paste0("combined_results_ordered$log2FoldChange", category, " 0")
    ))],
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
  expr_heatmap <- Heatmap(
    expr_matr,
    show_row_names = F,
    cluster_rows = F,
    cluster_columns = F,
    col = expr_col_fun,
    right_annotation = row_anno,
    heatmap_legend_param = list(
      title = "Expression",
      at = c(0, 1),
      labels = c("min", "max")
    ),
    column_names_max_height = unit(10, "cm"),
  )
  
  
  # Save the heat map.
  file_snip <- ifelse(category == "<", yes = "Skin_Treg__upregulation", no = "Blood_naive_Treg__upregulation")
  expr_heatmap_pdf <- paste0(
    plot_outdir,
    "/treg_hierarchies_diff_expr_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_heatmap_filtered_f_expr_discrep_",
    file_snip,
    "_by_sample.pdf"
  )
  expr_heatmap_rds <- paste0(
    plot_rds_outdir,
    "/treg_hierarchies_diff_expr_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_heatmap_filtered_f_expr_discrep_",
    file_snip,
    "_by_sample.rds"
  )
  pdf(expr_heatmap_pdf, width = 6, height = 10)
  draw(expr_heatmap)
  dev.off()
  saveRDS(expr_heatmap, file = expr_heatmap_rds)
  
}





#####################################################################
# Plot histograms of the differences between blood CCR8+ Tregs
# and skin Tregs.
#####################################################################

# Iterate over two large categories (decreasing expression during differentiation,
# increasing expression during differentiation).
for (category in c("<", ">")) {
  # Extract differences between skin Tregs and blood CCR8+ Tregs.
  diffs <- combined_results$Scaled_log_TPM_diff_blood_ccr8_treg_skin_treg[eval(parse(text = paste0(
    "combined_results$log2FoldChange", category, " 0"
  )))]
  
  # Generate the histogram.
  histogram <- ggplot() +
    aes(x = diffs) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_histogram(binwidth = 0.02, fill = "black") +
    xlab("Scaled logTPM skin Tregs\nminus scaled logTPM blood CCR8+ Tregs") +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"))
  
  # Save the histogram.
  file_snip <- ifelse(category == "<", yes = "Skin_Treg__upregulation", no = "Blood_naive_Treg__upregulation")
  histogram_pdf <- paste0(
    plot_outdir,
    "/treg_hierarchies_diff_expr_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_expr_diff_hist_filtered_f_expr_discrep_",
    file_snip,
    ".pdf"
  )
  histogram_rds <- paste0(
    plot_rds_outdir,
    "/treg_hierarchies_diff_expr_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_expr_diff_hist_filtered_f_expr_discrep_",
    file_snip,
    ".rds"
  )
  pdf(histogram_pdf, width = 5, height = 4)
  print(histogram)
  dev.off()
  saveRDS(histogram, file = histogram_rds)
}





############################################################################################
# Count in how many genes blood CCR8+ Tregs are closer to skin Tregs / blood naive Tregs.
############################################################################################

# Specify the file name under which output from this section will be saved.
output_file <- paste0(
  location,
  "/treg_hierarchies/diff_expr_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_expr_discrep_genes_counts.txt"
)
if (file.exists(output_file)) {
  file.remove(output_file)
}

# Print a status message.
capture.output(cat(
  rep(
    "*************************************************************************************************\n",
    2
  )
), file = output_file, append = T)
capture.output(
  print(
    "Numbers of genes in which blood CCR8+ Tregs are closer to skin Tregs / blood naive Tregs:"
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

# Iterate over two large categories (decreasing expression during differentiation,
# increasing expression during differentiation).
for (category in c("<", ">")) {
  # Print what category is currently quantified.
  capture.output(
    cat("############################################\n"),
    file = output_file,
    append = T
  )
  capture.output(print(
    ifelse(category == "<", yes = "Skin_Treg__upregulation", no = "Blood_naive_Treg__upregulation")
  ), file = output_file, append = T)
  capture.output(cat("\n"), file = output_file, append = T)
  
  # Restrict the data to those genes that correspond to the current category.
  combined_results_cat <- combined_results[eval(parse(text = paste0(
    "combined_results$log2FoldChange", category, " 0"
  ))), ]
  
  # Compute the total number of genes.
  total_gene_count <- nrow(combined_results_cat)
  
  # Compute gene numbers for the identified gene categories.
  closer_to_naive_gene_count <- length(
    which(
      combined_results_cat$Cell_type_closest_to_blood_ccr8_treg == "Blood naive Treg"
    )
  )
  closer_to_naive_gene_perc <- round(100 * closer_to_naive_gene_count / total_gene_count, 2)
  closer_to_skin_gene_count <- length(which(
    combined_results_cat$Cell_type_closest_to_blood_ccr8_treg == "Skin Treg"
  ))
  closer_to_skin_gene_perc <- round(100 * closer_to_skin_gene_count / total_gene_count, 2)
  
  # Print the results.
  capture.output(print(
    paste0(
      "Closer to blood naive Tregs: ",
      closer_to_naive_gene_count,
      " genes (",
      closer_to_naive_gene_perc,
      "%)."
    )
  ), file = output_file, append = T)
  capture.output(print(
    paste0(
      "Closer to skin Tregs: ",
      closer_to_skin_gene_count,
      " genes (",
      closer_to_skin_gene_perc,
      "%)."
    )
  ), file = output_file, append = T)
  capture.output(cat("\n\n"), file = output_file, append = T)
  
}
