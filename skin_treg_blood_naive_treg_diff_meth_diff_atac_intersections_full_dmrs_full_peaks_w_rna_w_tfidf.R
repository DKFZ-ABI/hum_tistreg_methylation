# This R script analyses hierarchies in Treg cell differentiation, focusing on 
#   DMR-peak pairs (skin Treg cells vs. blood naive Treg cells) that are 
#   assigned to differentially expressed genes (i.e. DMR-peak-gene links).
# Author: Niklas Beumer.
# Run this script with 10 cores!!!



# Load required packages.
library(bsseq)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(gtools)
library(testit)
library(ggplot2)
library(qpdf)
library(parallel)
library(cowplot)
library(rtracklayer)
library(Signac)
library(Seurat)


# Define a location on /yyy.
b330_space <- "/yyy/"
location <- paste0(b330_space, "nbeumer/hm_treg_bs_rgnsbg")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Specify the naming conventions for the relevant cell types.
relevant_cell_types <- c("Blood naive Treg", "Blood CCR8+ Treg", "Skin Treg")
relevant_cell_types_atac_naming <- c("blood_naive_treg", "blood_ccr8_treg", 
                                     "skin_treg")

# Read in the DMR-peak links.
dmr_peak_links_file <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_diff_acc_overlapping_dmr_peak_pairs_w_ccr8_positioning_w_tfidf.txt"
)
dmr_peak_links <- read.table(dmr_peak_links_file, header = T,
                             stringsAsFactors = F, sep = "\t")
colnames(dmr_peak_links) <- gsub(".", "+", colnames(dmr_peak_links), fixed = T)

# Read in the list of differentially methylated regions between skin Tregs and 
# blood naive Tregs.
meth_reg_file <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_gap_0.15_sterrratio_0.5_signature_regions_filtered.txt"
)
meth_reg <- read.table(meth_reg_file, header = T, stringsAsFactors = F, sep = "\t")

# Read in the sample mapping file.
rna_sample_mapping_path <- paste0(location, 
                                  "/sample_mapping_rnaseq_only_biol_rep.txt")
rna_sample_mapping <- read.table(rna_sample_mapping_path, header = T, 
                                 stringsAsFactors = F, sep = "\t")

# Read in the list of differentially expressed genes between skin Tregs and 
# blood naive Tregs.
diff_exp_file <- paste0(
  location, 
  "/RNASeq/analysis_results/2022-01-14_diff_gene_expr_DESEq2_Skin_TregBlood_naive_Treg_results_filtered_with_significance.txt"
)
diff_exp <- read.table(diff_exp_file, header = T, stringsAsFactors = F)
diff_exp <- diff_exp[diff_exp$significant, ]

# Read in the TPM values.
tpm_file <- paste0(location, "/RNASeq/tpm/tpm_all.txt")
tpm <- read.table(tpm_file, header = T, stringsAsFactors = F)

# Specify the prefix of the files containing genomic positions of exons, 
# transcription start sites, genes etc.
annotation_pref <- 
  "/yyy/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/GencodeV19_"

# Read in the files containing gene and exon annotation for normal chromosomes.
genes <- read.table(paste0(annotation_pref, "Genes_plain.bed.gz"), header = T, 
                    stringsAsFactors = F, comment.char = "")
exons <- read.table(paste0(annotation_pref, "Exons_plain.bed.gz"), header = T, 
                    stringsAsFactors = F, comment.char = "")

# Read in the BSseq object containing the smoothed methylation data.
# Restrict the dtaa to the relevant cell types.
meth_data_file <- paste0(
  location, 
  "/preprocessing_etc/quality_control_after_alignment/2022-03-14_bsseq_object_combined_all_samples_smoothed.rds"
)
meth_data <- readRDS(meth_data_file)
meth_data <- meth_data[, pData(meth_data)$Cell_type %in% relevant_cell_types]

# Read in the Seurat object containing scATAC-seq data for CD4+ T cells.
sc_data_cd4_file <- paste0(
  b330_space, 
  "msimon/scATAC/final/seurat_human_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar.RDS"
)
sc_data_cd4 <- readRDS(sc_data_cd4_file)

# Read in the information on Malte's scATAC-seq data for the cell types.
atac_seq_data_file <- paste0(
  location, 
  "/materials_for_scATAC_seq_tracks/Overview_Maltes_Bam_files.csv"
)
atac_seq_data_overview <- read.csv(atac_seq_data_file, stringsAsFactors = F)






##############################################################
# Process methylation data.
##############################################################

# Specify the colour palette for the cell types.
bsseq_metadata <- pData(meth_data)
cell_types_all <- unique(bsseq_metadata$Cell_type)
cell_type_col_palette <- c("blue", "orange", "darkorchid1")
names(cell_type_col_palette) <- cell_types_all

# Append this colour_palette to the meta data of the BSseq object.
bsseq_metadata$col <- cell_type_col_palette[bsseq_metadata$Cell_type]
pData(meth_data) <- bsseq_metadata





#########################################################################
# Process gene data.
#########################################################################

# Increase start and end positions by 1 nucleotide in order to convert intervals 
# to 1-based closed intervals.
genes$chromStart <- genes$chromStart + 1
exons$start <- exons$start + 1

# Generate a GRanges object out of the gene body data.
genes_gr <- makeGRangesFromDataFrame(genes,
                                     seqnames.field = "X.chrom",
                                     start.field = "chromStart",
                                     end.field = "chromEnd",
                                     keep.extra.columns = T)


# Infer the location of the transcription start sites and prepare the 
# corresponding data to appear in plots.
tss_pos <- sapply(1:nrow(genes), FUN = function(x) {
  ifelse(genes$strand[x] == "+", 
         yes = genes$chromStart[x], no = 
           genes$chromEnd[x])
})
tss_df <- data.frame(chr = genes$X.chrom, 
                     pos = tss_pos,
                     gene = genes$name,
                     strand = genes$strand)








#####################################################################
# Identify DMR-peak-gene links.
# A gene is considered linked to a DMR-peak pair if it was assigned
# to both the DMR and the peak in that pair.
# Furthermore, the gene has to be significantly differential.
#####################################################################

dmr_peak_gene_links <- do.call(rbind, lapply(1:nrow(dmr_peak_links), 
                                             FUN = function(x) {
  #-- Get all genes assigned to the corresponding DMR.
  dmr_genes <- strsplit(dmr_peak_links$Methylation__gene_assignments[x],
                        split = ", ")[[1]]
  #-- Get all genes assigned to the correponding peak.
  peak_genes <- strsplit(dmr_peak_links$Accessibility__gene_assignments[x],
                         split = ", ")[[1]]
  #-- Identify those genes that are assigned to both the DMR and the peak.
  intersect_genes <- intersect(dmr_genes, peak_genes)
  #-- Restrict to those genes that are significantly differentially expressed.
  intersect_genes_2 <- intersect(intersect_genes, rownames(diff_exp))
  #-- For each of these genes, generate a DMR-peak-gene link.
  if (length(intersect_genes_2) > 0) {
    temp_df <- do.call(rbind, lapply(intersect_genes_2, FUN = function(y) {
      temp_df_2 <- dmr_peak_links[x, ]
      temp_df_2$Expression__gene_name <- y
      return(temp_df_2)
    }))
  } else {
    temp_df <- dmr_peak_links[c(), ]
    temp_df$Expression__gene_name <- c()
  }
  return(temp_df)
}))







########################################################################
# Compute cell-type-wise average TPM values for the relevant genes.
########################################################################

# Restrict the TPM data to the genes in DMR-peak-gene links.
tpm_links <- tpm[
  tpm$Gene_symbol %in% dmr_peak_gene_links$Expression__gene_name, 
]

# For each DMR-peak-gene link, compute the mean expression across the cells.
#-- Iterate over the three relevant cell types.
avg_tpm_celltype_level <- sapply(relevant_cell_types, FUN = function(x) {
  #-- Identify the samples from this cell type.
  corresp_samples <- rna_sample_mapping$Sample[rna_sample_mapping$Cell_type == x]
  #-- Extract TPM values for these samples.
  corresp_tpm <- tpm_links[, corresp_samples]
  rownames(corresp_tpm) <- tpm_links$Gene_symbol
  #-- Transform the TPM values by log(100 * x + 1)
  for (colname in colnames(corresp_tpm)) {
    corresp_tpm[, colname] <- log1p(100 * corresp_tpm[, colname])
  }
  #-- Compute within-cell-type means.
  means <- rowMeans(corresp_tpm)
  #-- Return the means in the order in which the corresponding genes appear 
  #-- in the DMR-peak-gene links.
  return(means[dmr_peak_gene_links$Expression__gene_name])
})
assert(all(rownames(avg_tpm_celltype_level) == 
             dmr_peak_gene_links$Expression__gene_name))









##############################################################################
# Exclude genes for which the differential expression tendencies reported by 
# DESeq2 do not match with what I see in my log-transformed TPM values. 
# Note: This does not exclude any gene.
##############################################################################

# For each gene, check whether the differential expression tendency reported by 
# DESeq2 matches with what I see in the log-transformed TPM values.
exclude_bools <- sapply(1:nrow(dmr_peak_gene_links), FUN = function(x) {
  corresp_gene <- dmr_peak_gene_links$Expression__gene_name[x]
  diff_tend_sign <- sign(diff_exp[corresp_gene, "log2FoldChange"])
  skin_log_tpm <- avg_tpm_celltype_level[corresp_gene, "Skin Treg"]
  blood_naive_log_tpm <- avg_tpm_celltype_level[corresp_gene, "Blood naive Treg"]
  return(ifelse(diff_tend_sign == -1, 
                yes = skin_log_tpm < blood_naive_log_tpm, 
                no = skin_log_tpm > blood_naive_log_tpm))
})

# Identify genes to exclude.
genes_to_exclude <- rownames(diff_exp)[exclude_bools]

# Save information on which genes were excluded.
output_file <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_diff_acc_diff_expr_dmr_peak_gene_links_gene_exclusions.txt"
)
if (file.exists(output_file)) {
  file.remove(output_file)
}
capture.output(cat(rep("*************************************************************************************************\n", 2)), file = output_file, append = T)
capture.output(print("Genes that were excluded because the diff. expr. tendencies reported by DESeq2"), file = output_file, append = T)
capture.output(print("did not match with what I saw in the TPM data:"), file = output_file, append = T)
capture.output(cat(rep("*************************************************************************************************\n", 2)), file = output_file, append = T)
capture.output(cat("\n\n"), file = output_file, append = T)
capture.output(print(genes_to_exclude), file = output_file, append = T)

# Perform the exclusions.
avg_tpm_celltype_level_filtered <- avg_tpm_celltype_level[!exclude_bools, ]





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
  scaled_matr <- t(apply(matr, 1, FUN = function(x) {
    
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
    
  }))
  
  # Return the scaled matrix.
  return(scaled_matr)
  
}

# Scale the values between 0 and 1.
avg_tpm_celltype_level_filtered_scaled <- 
  scale_matr_lines_betw_0_and_1(avg_tpm_celltype_level_filtered)







##########################################################################
# Infer the differential expression tendency for each DMR-peak-gene link.
##########################################################################

# Identify differential expression tendencies.
diff_exp_tends <- sapply(dmr_peak_gene_links$Expression__gene_name, 
                         FUN = function(x) {
  log2FC <- diff_exp[x, "log2FoldChange"]
  tend <- paste0(
    ifelse(log2FC < 0, yes = "Skin_Treg", no = "Blood_naive_Treg"),
    "__hyperexpression"
  )
  return(tend)
})






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
closest_extrema <- sapply(1:nrow(avg_tpm_celltype_level_filtered_scaled), 
                          FUN = function(x) {
  if (ccr8_naive_dists[x] < ccr8_skin_dists[x]) {
    return("Blood naive Treg")
  } else if (ccr8_naive_dists[x] > ccr8_skin_dists[x]) {
    return("Skin Treg")
  } else {
    return("Exactly in middle")
  }
})





#########################################################################
# Collect all the computed values.
#########################################################################

# Prepare differential expression tendencies.
tend_df <- data.frame(Expression__signature_class = diff_exp_tends)

# Prepare gene expression values.
colnames(avg_tpm_celltype_level_filtered) <- 
  paste0("Expression__Log_TPM__", 
         gsub(" ", "_", colnames(avg_tpm_celltype_level_filtered)))
tpm_df <- as.data.frame(avg_tpm_celltype_level_filtered)
colnames(avg_tpm_celltype_level_filtered_scaled) <- 
  paste0("Expression__Scaled_log_TPM__", 
         gsub(" ", "_", colnames(avg_tpm_celltype_level_filtered_scaled)))
scaled_tpm_df <- as.data.frame(avg_tpm_celltype_level_filtered_scaled)

# Prepare results from distance analysis to extreme cell types.
dist_df <- data.frame(
  Expression__Scaled_log_TPM_diff_blood_ccr8_treg_blood_naive_treg = 
    ccr8_naive_diffs,
  Expression__Scaled_log_TPM_dist_blood_ccr8_treg_blood_naive_treg = 
    ccr8_naive_dists, 
  Expression__Scaled_log_TPM_diff_blood_ccr8_treg_skin_treg = 
    ccr8_skin_diffs,
  Expression__Scaled_log_TPM_dist_blood_ccr8_treg_skin_treg = 
    ccr8_skin_dists
)
closest_extrema_df <- data.frame(
  Expression__Blood_CCR8_Treg_pos  = paste0("Closer to ", 
                                            gsub("treg", 
                                                 "Treg", 
                                                 tolower(closest_extrema)),
                                            "s")
)

# Combine all information into a large data frame.
combined_results <- cbind(dmr_peak_gene_links,
                          tend_df,
                          tpm_df,
                          scaled_tpm_df,
                          dist_df,
                          closest_extrema_df)

# Save all computed values.
combined_results_outfile <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_diff_acc_diff_expr_dmr_peak_gene_links_w_ccr8_positioning_w_tfidf.txt"
)
write.table(combined_results, file = combined_results_outfile, 
            sep = "\t", row.names = F, col.names = T, quote = F)








############################################################################
# Plot a heat map showing scaled methylation, scaled accessibility and 
# scaled expression in the DMR-peak-gene links.
############################################################################

# Generate a version of the combined results, in which the DMR-peak-gene links
# are ordered by the difference between blood CCR8+ Tregs and skin Tregs.
order_to_use <- order(
  (combined_results$Methylation__Scaled_meth__Skin_Treg - 
    combined_results$`Methylation__Scaled_meth__Blood_CCR8+_Treg`) *
    sapply(combined_results$Methylation__automatic_annotation,
           FUN = function(x) {
             switch(x, Blood_naive_Treg__hypomethylation = -1, 
                    Skin_Treg__hypomethylation = 1)
           }
    )
)
combined_results_ordered <- combined_results[order_to_use, ]

# Iterate over the two methylation tendencies.
meth_tendencies <- unique(combined_results$Methylation__automatic_annotation)
for (tendency in meth_tendencies) {

  # Get the matrices for heatmap plotting.
  rows_to_keep <- which(
    combined_results_ordered$Methylation__automatic_annotation == tendency
  )
  meth_matr <- as.matrix(combined_results_ordered[
    rows_to_keep, grep("Scaled_meth", colnames(combined_results_ordered))
  ])
  colnames(meth_matr) <- gsub("_", 
                              " ", 
                              gsub("Methylation__Scaled_meth__", 
                                   "", 
                                   colnames(meth_matr)))
  acc_matr <- as.matrix(combined_results_ordered[
    rows_to_keep, grep("Scaled_mean_norm_acc", 
                       colnames(combined_results_ordered))
  ])
  colnames(acc_matr) <- gsub("_", 
                             " ", 
                             gsub("Accessibility__Scaled_mean_norm_acc__", 
                                  "", 
                                  colnames(acc_matr)))
  gex_matr <- as.matrix(combined_results_ordered[
    rows_to_keep, grep("Scaled_log_TPM__", colnames(combined_results_ordered))
  ])
  colnames(gex_matr) <- gsub("_", 
                             " ", 
                             gsub("Expression__Scaled_log_TPM__", 
                                  "", 
                                  colnames(gex_matr)))
  assert(nrow(meth_matr) == nrow(acc_matr))
  assert(nrow(acc_matr) == nrow(gex_matr))
  rownames(meth_matr) <- rownames(acc_matr) <- rownames(gex_matr) <- 
    1:nrow(meth_matr)
  
  # Generate the colour functions for the heatmaps.
  meth_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                             colors = viridis(200, direction = -1))
  acc_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                            colors = viridis(200, option = "B"))
  gex_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                            colors = viridis(200, option = "G"))
  
  # Generate a row annotation showing to which cell type blood CCR8+ Tregs 
  # are closest for each feature.
  row_anno_meth <- rowAnnotation(
    `Blood CCR8+ Treg position` = 
      combined_results_ordered$Methylation__Blood_CCR8_Treg_pos[rows_to_keep],
    col = list(`Blood CCR8+ Treg position` = c(
      `Closer to skin Tregs` = "blue", 
      `Closer to blood naive Tregs` = "darkorchid1",
      `Exactly in middle` = "grey"
    ))
  )
  row_anno_acc <- rowAnnotation(
    `Blood CCR8+ Treg position` = 
      combined_results_ordered$Accessibility__Blood_CCR8_Treg_pos[rows_to_keep],
    col = list(
      `Blood CCR8+ Treg position` = c(
        `Closer to skin Tregs` = "blue", 
        `Closer to blood naive Tregs` = "darkorchid1",
        `Exactly in middle` = "grey"
      ))
  )
  row_anno_gex <- rowAnnotation(
    `Blood CCR8+ Treg position` = 
      combined_results_ordered$Expression__Blood_CCR8_Treg_pos[rows_to_keep],
    col = list(
      `Blood CCR8+ Treg position` = c(
        `Closer to skin Tregs` = "blue", 
        `Closer to blood naive Tregs` = "darkorchid1",
        `Exactly in middle` = "grey"
      ))
  )
  
  # Generate single heat maps.
  meth_heatmap <- Heatmap(meth_matr, 
                          cluster_rows = F,
                          cluster_columns = F,
                          col = meth_col_fun,
                          right_annotation = row_anno_meth,
                          heatmap_legend_param = list(title = "Methylation", 
                                                      at = c(0, 1), 
                                                      labels = c("min", "max")))
  acc_heatmap <- Heatmap(acc_matr, 
                         show_row_names = F,
                         cluster_columns = F,
                         col = acc_col_fun,
                         right_annotation = row_anno_acc,
                         heatmap_legend_param = list(title = "Chrom. accessibility", 
                                                     at = c(0, 1), 
                                                     labels = c("min", "max")))
  
  gex_heatmap <- Heatmap(gex_matr, 
                         show_row_names = F,
                         cluster_columns = F,
                         col = gex_col_fun,
                         right_annotation = row_anno_gex,
                         heatmap_legend_param = list(title = "Expression", 
                                                     at = c(0, 1), 
                                                     labels = c("min", "max")))
  
  
  
  # Combine the two heat maps into a list.
  # Note: By default, this orders the rows of the second/third heat map
  # in the same manner as in the first heat map
  # (https://jokergoo.github.io/ComplexHeatmap-reference/book/a-list-of-heatmaps.html;
  # 09 May 2023)
  combined_heatmaps <- meth_heatmap + acc_heatmap + gex_heatmap
  
  # Save the combined heat map.
  combined_heatmaps_pdf <- paste0(
    plot_outdir, 
    "/diff_meth_diff_acc_diff_gex_dmr_peak_gene_links_w_ccr8_positioning_heatmap_", 
    tendency, 
    "_w_tfidf.pdf"
  )
  combined_heatmaps_rds <- paste0(
    plot_rds_outdir, 
    "/diff_meth_diff_acc_diff_gex_dmr_peak_gene_links_w_ccr8_positioning_heatmap_", 
    tendency, 
    "_w_tfidf.rds"
  )
  pdf(combined_heatmaps_pdf, width = 10, height = 10)
  draw(combined_heatmaps)
  dev.off() 
  saveRDS(combined_heatmaps, file = combined_heatmaps_rds)
  
}




############################################################################
# Plot a heat map showing scaled methylation, scaled accessibility and 
# scaled expression in the DMR-peak-gene links, ordered by genomic position.
############################################################################

# Assert that the DMR-peak-gene links are ordered by genomic position on their 
# chromosomes.
assert(all(combined_results$Methylation__region_ID ==
             mixedsort(combined_results$Methylation__region_ID)))
peaks_gr <- StringToGRanges(combined_results$Accessibility__peak_ID)
unique_chrs <- as.character(unique(seqnames(peaks_gr)))
void <- sapply(unique_chrs, FUN = function(x) {
  chr_subset <- peaks_gr[seqnames(peaks_gr) == x]
  starts <- start(chr_subset)
  assert(all(starts == sort(starts)))
})

# Iterate over the two methylation tendencies.
for (tendency in meth_tendencies) {
  
  # Get the matrices for heatmap plotting.
  rows_to_keep <- which(
    combined_results$Methylation__automatic_annotation == tendency
  )
  meth_matr <- as.matrix(combined_results[
    rows_to_keep, grep("Scaled_meth", colnames(combined_results))
  ])
  colnames(meth_matr) <- gsub("_", 
                              " ", 
                              gsub("Methylation__Scaled_meth__", 
                                   "", 
                                   colnames(meth_matr)))
  acc_matr <- as.matrix(combined_results[
    rows_to_keep, grep("Scaled_mean_norm_acc", colnames(combined_results))
  ])
  colnames(acc_matr) <- gsub("_", 
                             " ", 
                             gsub("Accessibility__Scaled_mean_norm_acc__", 
                                  "", 
                                  colnames(acc_matr)))
  gex_matr <- as.matrix(combined_results[
    rows_to_keep, grep("Scaled_log_TPM__", colnames(combined_results))
  ])
  colnames(gex_matr) <- gsub("_", 
                             " ", 
                             gsub("Expression__Scaled_log_TPM__", 
                                  "", 
                                  colnames(gex_matr)))
  assert(nrow(meth_matr) == nrow(acc_matr))
  assert(nrow(acc_matr) == nrow(gex_matr))
  rownames(meth_matr) <- rownames(acc_matr) <- rownames(gex_matr) <- 
    1:nrow(meth_matr)
  
  # Generate the colour functions for the heatmaps.
  meth_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                             colors = viridis(200, direction = -1))
  acc_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                            colors = viridis(200, option = "B"))
  gex_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                            colors = viridis(200, option = "G"))
  
  # Generate a row annotation showing to which cell type blood CCR8+ Tregs 
  # are closest for each feature.
  row_anno_meth <- rowAnnotation(
    `Blood CCR8+ Treg position` = 
      combined_results$Methylation__Blood_CCR8_Treg_pos[rows_to_keep],
    col = list(`Blood CCR8+ Treg position` = c(
      `Closer to skin Tregs` = "blue", 
      `Closer to blood naive Tregs` = "darkorchid1",
      `Exactly in middle` = "grey"
    ))
  )
  row_anno_acc <- rowAnnotation(
    `Blood CCR8+ Treg position` = 
      combined_results$Accessibility__Blood_CCR8_Treg_pos[rows_to_keep],
    col = list(
      `Blood CCR8+ Treg position` = c(
        `Closer to skin Tregs` = "blue", 
        `Closer to blood naive Tregs` = "darkorchid1",
        `Exactly in middle` = "grey"
      ))
  )
  row_anno_gex <- rowAnnotation(
    `Blood CCR8+ Treg position` = 
      combined_results$Expression__Blood_CCR8_Treg_pos[rows_to_keep],
    col = list(
      `Blood CCR8+ Treg position` = c(
        `Closer to skin Tregs` = "blue", 
        `Closer to blood naive Tregs` = "darkorchid1",
        `Exactly in middle` = "grey"
      ))
  )
  
  # Generate a row annotation showing to chromosomes.
  col_palette <- rainbow(23)
  names(col_palette) <- unique_chrs
  row_anno_chr <- rowAnnotation(
    Chromosome = 
      as.character(seqnames(peaks_gr)[rows_to_keep]),
    col = list(Chromosome = col_palette)
  )
  
  # Generate single heat maps.
  meth_heatmap <- Heatmap(meth_matr, 
                          cluster_rows = F,
                          cluster_columns = F,
                          col = meth_col_fun,
                          left_annotation = row_anno_chr,
                          right_annotation = row_anno_meth,
                          heatmap_legend_param = list(title = "Methylation", 
                                                      at = c(0, 1), 
                                                      labels = c("min", "max")))
  acc_heatmap <- Heatmap(acc_matr, 
                         show_row_names = F,
                         cluster_columns = F,
                         col = acc_col_fun,
                         right_annotation = row_anno_acc,
                         heatmap_legend_param = list(title = "Chrom. accessibility", 
                                                     at = c(0, 1), 
                                                     labels = c("min", "max")))
  
  gex_heatmap <- Heatmap(gex_matr, 
                         show_row_names = F,
                         cluster_columns = F,
                         col = gex_col_fun,
                         right_annotation = row_anno_gex,
                         heatmap_legend_param = list(title = "Expression", 
                                                     at = c(0, 1), 
                                                     labels = c("min", "max")))
  
  
  
  # Combine the two heat maps into a list.
  # Note: By default, this orders the rows of the second/third heat map
  # in the same manner as in the first heat map
  # (https://jokergoo.github.io/ComplexHeatmap-reference/book/a-list-of-heatmaps.html;
  # 09 May 2023)
  combined_heatmaps <- meth_heatmap + acc_heatmap + gex_heatmap
  
  # Save the combined heat map.
  combined_heatmaps_pdf <- paste0(
    plot_outdir, 
    "/diff_meth_diff_acc_diff_gex_dmr_peak_gene_links_w_ccr8_positioning_heatmap_", 
    tendency, 
    "_by_position_w_tfidf.pdf"
  )
  combined_heatmaps_rds <- paste0(
    plot_rds_outdir, 
    "/diff_meth_diff_acc_diff_gex_dmr_peak_gene_links_w_ccr8_positioning_heatmap_", 
    tendency, 
    "_by_position_w_tfidf.rds"
  )
  pdf(combined_heatmaps_pdf, width = 10, height = 10)
  draw(combined_heatmaps)
  dev.off() 
  saveRDS(combined_heatmaps, file = combined_heatmaps_rds)
  
}




############################################################################
# Plot a heat map showing scaled methylation, scaled accessibility and 
# scaled expression in the DMR-peak-gene links with the original order,
# this time displaying each sample/donor independently.
############################################################################

# Generate a GRanges object of the DMRs in the DMR-peak-gene links.
links_dmrs_gr <- makeGRangesFromDataFrame(
  do.call(
    rbind, 
    lapply(combined_results_ordered$Methylation__region_ID, FUN = function(x) {
      meth_reg[meth_reg$region_ID == x, ]
    })
  )
)

# Identify average raw methylation values (on the sample level).
links_avg_meth_sample_level <- getMeth(meth_data,
                                       regions = links_dmrs_gr,
                                       type = "raw",
                                       what = "perRegion")

# To enable donor-leevel heat map visualisation, subset scATAC-seq data 
# to the relevant cell types and cells from the relevant donors.
# Afterwards, perform a new TF-IDF.
blood_donors <- c("1", "2")
skin_donors <- c("4", "5")
sc_data_cd4_restr_f_donorlevel <- subset(
  sc_data_cd4, 
  subset = (treg_tconv_annot == "blood_naive_treg" & donor %in% blood_donors) |
    (treg_tconv_annot == "blood_ccr8_treg" & donor %in% blood_donors) |
    (treg_tconv_annot == "skin_treg" & donor %in% skin_donors)
)
sc_data_cd4_restr_f_donorlevel <- RunTFIDF(sc_data_cd4_restr_f_donorlevel)

# For each cell type, each donor and each region, compute the mean accessibility 
# across the cells.
links_avg_acc_donor_level <- do.call(
  cbind,
  lapply(relevant_cell_types_atac_naming, FUN = function(x) {
    sc_data_cd4_this_cell_type <- subset(sc_data_cd4_restr_f_donorlevel, 
                                         subset = treg_tconv_annot == x)
    this_cell_type_donors <- sort(unique(sc_data_cd4_this_cell_type$donor))
    norm_acc_data_this_celltype <- sapply(
      this_cell_type_donors,
      FUN = function(y) {
        sc_data_this_donor <- subset(sc_data_cd4_this_cell_type,
                                     subset = donor == y)
        norm_acc_data <- GetAssayData(
          sc_data_this_donor, 
          assay = "scATAC_raw", slot = "data"
        )[combined_results_ordered$Accessibility__peak_ID, ]
        print("Matrix extracted")
        data_aggr <- rowMeans(norm_acc_data)
        return(data_aggr)
      }
    )
    colnames(norm_acc_data_this_celltype) <- 
      paste0(x, "_donor_", this_cell_type_donors)
    return(norm_acc_data_this_celltype)
  })
)

# Collect the sample-level TPM values.
links_tpm_sample_level <- do.call(
  cbind,
  lapply(relevant_cell_types, FUN = function(x) {
    #-- Identify the samples from this cell type.
    corresp_samples <- rna_sample_mapping$Sample[
      rna_sample_mapping$Cell_type == x
    ]
    #-- Extract TPM values for these samples.
    corresp_tpm <- tpm[, corresp_samples]
    #-- Transform the TPM values by log(100 * x + 1)
    for (colname in colnames(corresp_tpm)) {
      corresp_tpm[, colname] <- log1p(100 * corresp_tpm[, colname])
    }
    #-- Return the sample-level values in the order in which the corresponding 
    #-- genes appear in the DMR-peak-gene links.
    corresp_tpm_ordered <- do.call(
      rbind, 
      lapply(combined_results_ordered$Expression__gene_name, FUN = function(y) {
        corresp_tpm[tpm$Gene_symbol == y, ]
      })
    )
    return(corresp_tpm_ordered)
  })
)
assert(nrow(links_tpm_sample_level) == nrow(links_avg_acc_donor_level))


# Scale the values between 0 and 1.
links_avg_meth_sample_level_scaled <- 
  scale_matr_lines_betw_0_and_1(links_avg_meth_sample_level)
links_avg_acc_donor_level_scaled <- 
  scale_matr_lines_betw_0_and_1(links_avg_acc_donor_level)
links_tpm_sample_level_scaled <- 
  scale_matr_lines_betw_0_and_1(links_tpm_sample_level)

# Iterate over the two methylation tendencies.
for (tendency in meth_tendencies) {
  
  # Get the matrices for heatmap plotting.
  rows_to_keep <- which(
    combined_results_ordered$Methylation__automatic_annotation == tendency
  )
  meth_matr <- links_avg_meth_sample_level_scaled[rows_to_keep, ]
  acc_matr <- links_avg_acc_donor_level_scaled[rows_to_keep, ]
  gex_matr <- links_tpm_sample_level_scaled[rows_to_keep, ]
  assert(nrow(meth_matr) == nrow(acc_matr))
  assert(nrow(acc_matr) == nrow(gex_matr))
  rownames(meth_matr) <- rownames(acc_matr) <- rownames(gex_matr) <- 
    1:nrow(meth_matr)
  
  # Re-order columns of the methylation matrix according to how they should be 
  # displayed in the publication.
  cell_type_order <- c("Blood naive Treg", "Blood CCR8+ Treg", "Skin Treg")
  meth_matr_ordered <- do.call(cbind, lapply(cell_type_order, FUN = function(x) {
    meth_matr[, grep(x, colnames(meth_matr), fixed = T)]
  }))
  
  # Generate the colour functions for the heatmaps.
  meth_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                             colors = viridis(200, direction = -1))
  acc_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                            colors = viridis(200, option = "B"))
  gex_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                            colors = viridis(200, option = "G"))
  
  # Generate a row annotation showing to which cell type blood CCR8+ Tregs 
  # are closest for each feature.
  row_anno_meth <- rowAnnotation(
    `Blood CCR8+ Treg position` = 
      combined_results_ordered$Methylation__Blood_CCR8_Treg_pos[rows_to_keep],
    col = list(`Blood CCR8+ Treg position` = c(
      `Closer to skin Tregs` = "blue", 
      `Closer to blood naive Tregs` = "darkorchid1",
      `Exactly in middle` = "grey"
    ))
  )
  row_anno_acc <- rowAnnotation(
    `Blood CCR8+ Treg position` = 
      combined_results_ordered$Accessibility__Blood_CCR8_Treg_pos[rows_to_keep],
    col = list(
      `Blood CCR8+ Treg position` = c(
        `Closer to skin Tregs` = "blue", 
        `Closer to blood naive Tregs` = "darkorchid1",
        `Exactly in middle` = "grey"
      ))
  )
  row_anno_gex <- rowAnnotation(
    `Blood CCR8+ Treg position` = 
      combined_results_ordered$Expression__Blood_CCR8_Treg_pos[rows_to_keep],
    col = list(
      `Blood CCR8+ Treg position` = c(
        `Closer to skin Tregs` = "blue", 
        `Closer to blood naive Tregs` = "darkorchid1",
        `Exactly in middle` = "grey"
      ))
  )
  
  # Generate single heat maps.
  meth_heatmap <- Heatmap(meth_matr_ordered, 
                          cluster_rows = F,
                          cluster_columns = F,
                          col = meth_col_fun,
                          right_annotation = row_anno_meth,
                          heatmap_legend_param = list(title = "Methylation", 
                                                      at = c(0, 1), 
                                                      labels = c("min", "max")),
                          width = 1)
  acc_heatmap <- Heatmap(acc_matr, 
                         show_row_names = F,
                         cluster_columns = F,
                         col = acc_col_fun,
                         right_annotation = row_anno_acc,
                         heatmap_legend_param = list(title = "Chrom. accessibility", 
                                                     at = c(0, 1), 
                                                     labels = c("min", "max")),
                         width = 1)
  
  gex_heatmap <- Heatmap(gex_matr, 
                         show_row_names = F,
                         cluster_columns = F,
                         col = gex_col_fun,
                         right_annotation = row_anno_gex,
                         heatmap_legend_param = list(title = "Expression", 
                                                     at = c(0, 1), 
                                                     labels = c("min", "max")),
                         width = 1)
  
  
  
  # Combine the two heat maps into a list.
  # Note: By default, this orders the rows of the second/third heat map
  # in the same manner as in the first heat map
  # (https://jokergoo.github.io/ComplexHeatmap-reference/book/a-list-of-heatmaps.html;
  # 09 May 2023)
  combined_heatmaps <- meth_heatmap + acc_heatmap + gex_heatmap
  
  # Save the combined heat map.
  combined_heatmaps_pdf <- paste0(
    plot_outdir, 
    "/diff_meth_diff_acc_diff_gex_dmr_peak_gene_links_w_ccr8_positioning_heatmap_", 
    tendency, 
    "_w_tfidf_by_donor.pdf"
  )
  combined_heatmaps_rds <- paste0(
    plot_rds_outdir, 
    "/diff_meth_diff_acc_diff_gex_dmr_peak_gene_links_w_ccr8_positioning_heatmap_", 
    tendency, 
    "_w_tfidf_by_donor.rds"
  )
  pdf(combined_heatmaps_pdf, width = 10, height = 10)
  draw(combined_heatmaps)
  dev.off() 
  saveRDS(combined_heatmaps, file = combined_heatmaps_rds)
  
}





########################################################################
# Visualise the proportions of RNA-based blood CCR8+ Treg positionings
# in different classes of DMR-peak-gene links.
########################################################################

# Collect all differential accessibility and expression tendencies.
diff_acc_tends <- unique(combined_results$Accessibility__signature_class)
diff_gex_tends <- unique(combined_results$Expression__signature_class)

# Generate file snippets corresponding to the tendencies.
diff_meth_tend_snips <- c("sk_hypometh", "nv_hypometh")
names(diff_meth_tend_snips) <- meth_tendencies
diff_acc_tend_snips <- c("sk_hyperacc", "nv_hyperacc")
names(diff_acc_tend_snips) <- diff_acc_tends
diff_gex_tend_snips <- c("sk_hyperexp", "nv_hyperexp")
names(diff_gex_tend_snips) <- diff_gex_tends

# Iterate over the two methylation tendencies.
for (meth_tend in meth_tendencies) {
  
  # Iterate over the two differential accessibility tendencies.
  for (acc_tend in diff_acc_tends) {
    
    # Iterate over the two differential gene expression tendencies.
    for (gex_tend in diff_gex_tends) {
      
      # Extract the data for plotting.
      plotting_data <- combined_results[
        combined_results$Methylation__automatic_annotation == meth_tend &
          combined_results$Accessibility__signature_class == acc_tend &
          combined_results$Expression__signature_class == gex_tend, 
      ]
      
      # Improve entries so that labels in the plots become clearer.
      plotting_data$Methylation__Blood_CCR8_Treg_pos <- 
        paste0("Methylation:\n", plotting_data$Methylation__Blood_CCR8_Treg_pos)
      plotting_data$Accessibility__Blood_CCR8_Treg_pos <- 
        paste0("Chr. accessibility:\n", 
               plotting_data$Accessibility__Blood_CCR8_Treg_pos)
      
      # Count the number of DMR-peak-gene link, stratified by methylation-based 
      # and accessibility-based blood CCR8+ Treg positionings.
      text_data <- do.call(rbind, lapply(
        unique(plotting_data$Methylation__Blood_CCR8_Treg_pos), 
        FUN = function(x) {
          do.call(rbind, lapply(
            unique(plotting_data$Accessibility__Blood_CCR8_Treg_pos),
            FUN = function(y) {
              temp_df <- data.frame(
                Methylation__Blood_CCR8_Treg_pos = x,
                Accessibility__Blood_CCR8_Treg_pos = y,
                count = paste0(length(which(
                  plotting_data$Methylation__Blood_CCR8_Treg_pos == x & 
                    plotting_data$Accessibility__Blood_CCR8_Treg_pos == y
                )), " links"),
                y = 0.5
              )
            }
          ))
        }
      ))
      
      # Generate a pie chart showing RNA-based blood CCR8+ positionings,
      # stratified by methylation-based and accessibility-based blood CCR8+
      # Treg positionings.
      pie_chart <- ggplot() +
        scale_fill_manual(breaks = c("Closer to blood naive Tregs", 
                                     "Closer to skin Tregs"),
                          values = c("darkorchid1", "blue"),
                          name = "Expression-based blood\nCCR8+Treg position") +
        geom_bar(data = plotting_data,
                 mapping = aes(x = "A", fill = Expression__Blood_CCR8_Treg_pos),
                 position = "fill", 
                 colour = "black") +
        geom_label(data = text_data, 
                   mapping = aes(x = "A", y = y, label = count),
                   alpha = 0.75) +
        coord_polar(theta = "y", direction = -1) +
        facet_grid(cols = vars(Methylation__Blood_CCR8_Treg_pos),
                   rows = vars(Accessibility__Blood_CCR8_Treg_pos)) +
        ggtitle(paste(meth_tend, acc_tend, gex_tend, sep = ",\n")) +
        theme_classic() +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank())
      
      # Save the pie chart.
      file_snip <- paste0(
        "/diff_meth_diff_acc_diff_gex_dmr_peak_gene_links_w_ccr8_positioning_piecharts_", 
        diff_meth_tend_snips[meth_tend], 
        "_", 
        diff_acc_tend_snips[acc_tend], 
        "_", 
        diff_gex_tend_snips[gex_tend]
      )
      pie_outfile_pdf <- paste0(plot_outdir, file_snip, "_w_tfidf.pdf")
      pie_outfile_rds <- paste0(plot_rds_outdir, file_snip, "_w_tfidf.rds")
      pdf(pie_outfile_pdf, width = 6, height = 5)
      print(pie_chart)
      dev.off()
      saveRDS(pie_chart, file = pie_outfile_rds)
      
    }
  }
}







############################################################################
# Generate track plots fo the DMR-peak-gene links.
############################################################################

get_cpg_pos_in_region <- function(bsseq_obj, region) {
  # This function extracts the CpG sites overlapping with a specified region.
  # bsseq_obj: The BSseq object to extract CpG sites from.
  # region: A GRanges object specifying the region to look at.
  # Dependencies: bsseq, GenomicRanges (need not be loaded but must be installed)
  # Value: A GRanges object with the CpG sites.
  # Author: Niklas Beumer
  
  # Find overlaps between the region and all CpGs in the BSseq object.
  region_overlaps <- GenomicRanges::findOverlaps(query = region, subject = bsseq_obj@rowRanges)
  
  # Index the CpG positions in the BSseq object by the indices of the subject hits.
  cpgs_in_region <- bsseq_obj@rowRanges[to(region_overlaps), ]
  
  # Return the extracted CpG sites.
  return(cpgs_in_region)
  
}

read_bigwig_file_cell_type <- function(paths_table, cell_type_col, bw_path_col, paths_prefix = "", cell_type) {
  # This function reads in a normalised bigwig file for a cell type.
  # paths_table: Data frame containing paths to bigwig files.
  # cell_type_col: Character; The name of the column in paths_table that contains cell type information.
  # bw_path_col: Character; The name of the column in paths_table that contains bigwig file paths.
  # paths_prefix: Character; A prefix to put in front of every path. If specified, must end with a "/".
  # cell_type: Character; The cell type to read the bigwig file for.
  # Dependencies: rtracklayer, testit
  # Value: A GRanges object containing the data from the bigwig file. 
  # Author: Niklas Beumer
  
  # Specify the path to the bigwig file.
  bigwig_path_prelim <- unique(paths_table[paths_table[, cell_type_col] == cell_type, bw_path_col])
  bigwig_path <- paste0(paths_prefix, bigwig_path_prelim)
  if (length(grep(".bw$", bigwig_path, perl = T)) == 0) { # If the path only leads to a directory but not
    # to the file itself, look for the file name.
    bigwig_file_name <- list.files(bigwig_path, pattern = ".bw")
    bigwig_path <- paste0(bigwig_path, "/", bigwig_file_name)
  }
  
  # Read in the bigwig file.
  bw_data <- import(bigwig_path)
  
  # Check that bins whose length is not a multiple of 50 occur only once per chromosome.
  # This behaviour would be expected because the bin at the very end of a chromosome will not be
  # an exact multiple of 50.
  assert(all(table(as.character(seqnames(bw_data))[(which(width(ranges(bw_data)) %% 50 != 0))]) == 1))
  
  # Synchronise the chromosome names with the naming convention used in my other work.
  seqlevelsStyle(bw_data) <- "NCBI"
  
  # Return the GRanges object containing the data from the bigwig files.
  return(bw_data)
  
}

plot_meth_and_other_signal_many_regions_custom_titles_custom_highlight_many_cell_types <- 
  function(bsseq_obj, regions_gr, extend = rep(5000, length(regions_gr)), genes, exons, 
           tss_pos, atac_bigwig_path_table, atac_bigwig_prefix, rna_bigwig_path_table, 
           add_track_bws = NA, add_track_names = NA, meth_highlight = NA, atac_highlight = NA, 
           rna_highlight = NA, gene_highlight = NA, add_track_highlight = NA,
           plot_titles, outfile_pdf, ncores) {
    # This function plots smoothed methalytion values, scATAC-seq signals and RNA-seq signals 
    #   in several regions using custom plot titles (for more than two cell types.
    # bsseq_obj: The BSseq object to retrieve methylation data from. In the pData, this object must contain
    #   the columns "Cell_type" and "col". All cell types present in the BSseq object will be
    #   plotted.
    # regions_gr: A GRanges object containing the regions to plot.
    # extend: Numeric: The number of bases by which the regions wil be extended on either side.
    #   (The plot will range from the first to the last CpG inside this extended region).
    # genes: A GRanges object containing information on human genes. Gene names have the 
    #   meta data column name "name".
    # exons: A data.frame containing data on human exons. Must contain the columns
    #   "gene_name", "start" and "end".
    # tss_pos: A data.frame containing data on human transcription start sites. Must contain the columns 
    #   "pos", "gene" and "strand".
    # atac_bigwig_path_table: A data.frame containing information on my bigwig files with Malte's scATAC-seq data
    #   (The table follows my own convention, using the columns "Cell_type", "Bam_file_path" and "Norm_bigwig_dir").
    # atac_bigwig_prefix: Character; A prefix to put in front of all bigwig file paths mentioned in bigwig_path_table.
    #   Must end with "/".
    # rna_bigwig_path_table: A data.frame containing paths to bigwig files with RNA signals.
    # add_track_bws: Character; Vector of paths to BigWig files whose contents will be plotted as additional tracks.
    #   If no additional tracks should be plotted, set this to NA.
    # add_track_names: Character; A vector of names that should appear for the additional tracks specified
    #   in "add_track_bws". If no additional tracks should be plotted, set this to NA. If "add_track_bws" is
    #   set, this argument must contain a value.
    # meth_highlight: GRanges object containing the regions to highlight in the methylation track.
    #   May contain regions outside of the actual plotting window. Set this to NA if no highlighting is desired.
    # atac_highlight: GRanges object containing the regions to highlight in the ATAC track.
    #   May contain regions outside of the actual plotting window. Set this to NA if no highlighting is desired.
    # rna_highlight: GRanges object containing the regions to highlight in the RNA track.
    #   May contain regions outside of the actual plotting window. Set this to NA if no highlighting is desired.
    # gene_highlight: GRanges object containing the regions to highlight in the gene track.
    #   May contain regions outside of the actual plotting window. Set this to NA if no highlighting is desired.
    # add_track_highlight: GRanges object containing the regions to highlight in the additional tracks.
    #   May contain regions outside of the actual plotting window. Set this to NA if no highlighting is desired
    #   or no additional tracks should be plotted.
    # plot_titles: Character; The titles to use for the single plots. These should be in an order corresponding
    #   to the order of regions_gr.
    # outfile_pdf: Character; Path to a PDF file where the plot will be saved.
    # ncores: Numeric; The number of cores to use.
    # Dependencies: bsseq, ggplot2, cowplot, GenomicRanges, testit, rtracklayer, qpdf, parallel
    # Value: None, the function just generates and saves the plots.
    # Author: Niklas Beumer
    
    # Generate a GRanges object that contains the regions of interest extended by the specified ranges.
    regions_extend_gr <- regions_gr + extend
    
    # Extract the cell types, sample names and colours corresponding to the 
    # samples in the BSseq object.
    cell_types_for_plot <- pData(bsseq_obj)$Cell_type
    cell_types_for_plot_unique <- unique(cell_types_for_plot)
    colours_for_plot <- unique(pData(bsseq_obj)$col)
    sample_names_for_plot <- rownames(pData(bsseq_obj))
    
    # Extract smoothed methylation values for the extended regions.
    smoothed_meth_in_regions <- getMeth(bsseq_obj, regions = regions_extend_gr, 
                                        type = "smooth", what = "perBase")
    
    # Read in Malte's scATAC-seq data.
    cell_type_atac_data <- lapply(cell_types_for_plot_unique, FUN = function(x) {
      read_bigwig_file_cell_type(paths_table = atac_bigwig_path_table,
                                 cell_type_col = "Cell_type",
                                 bw_path_col = "Norm_bigwig_dir",
                                 paths_prefix = atac_bigwig_prefix,
                                 cell_type = x)
    })
    names(cell_type_atac_data) <- c(cell_types_for_plot_unique)
    
    # Read in the RNA track data for the two cell types.
    # Also specify colours to use in the RNA plot. 
    rna_data_list <- lapply(cell_types_for_plot_unique, FUN = function(x) {
      read_bigwig_file_cell_type(paths_table = rna_sample_mapping,
                                 cell_type_col = "Cell_type",
                                 bw_path_col = "Bigwig_path",
                                 cell_type = x)
    })
    names(rna_data_list) <- cell_types_for_plot_unique
    rna_colours <- colours_for_plot
    
    # If additional tracks should be plotted, read in the corresponding BigWig files.
    # If a bin is only a single nucleotide long, increase its end position by 1 so
    # that it can be displayed in the plot.
    if (!(all(is.na(add_track_bws)))) {
      additional_tracks <- lapply(add_track_bws, FUN = function(x) {
        print(paste0("Reading ", x))
        bw_data <- import(x)
        seqlevelsStyle(bw_data) <- "NCBI"
        length_1_inds <- which(width(bw_data) == 1)
        end(bw_data)[length_1_inds] <- end(bw_data)[length_1_inds] + 1
        return(bw_data)
      })
    }
    
    
    # Iterate over all specified regions and save a plot for each of them.
    void <- mclapply(1:length(regions_extend_gr), FUN = function(l) {
      
      
      #-------------- Preliminary steps -------------------
      
      # Extract the CpG positions inside the extended region from the BSseq object.
      region_to_plot <- regions_extend_gr[l]
      cpgs_in_extended_region <- get_cpg_pos_in_region(bsseq_obj, region_to_plot)
      
      # Extract the chromosome on which the region is located.
      chromosome_for_plot <- as.character(seqnames(region_to_plot))
      
      
      #-------------- Methylation data -------------------
      
      # Prepare the methylation data for plotting.
      meth_plotting_data <- data.frame(Pos = rep(start(cpgs_in_extended_region), ncol(smoothed_meth_in_regions[[l]])),
                                       Meth = c(smoothed_meth_in_regions[[l]]),
                                       Sample = rep(sample_names_for_plot, each = nrow(smoothed_meth_in_regions[[l]])),
                                       Cell_type = rep(cell_types_for_plot, each = nrow(smoothed_meth_in_regions[[l]])))
      meth_plotting_data$Cell_type <- factor(meth_plotting_data$Cell_type, levels = names(rna_data_list))
      
      # Identify all highlight regions that overlap with the region to plot.
      if (all(!(is.na(meth_highlight)))) {
        overlaps <- findOverlaps(query = meth_highlight,
                                 subject = region_to_plot,
                                 ignore.strand = T)
        regions_to_highlight <- as.data.frame(meth_highlight[from(overlaps)])
        regions_to_highlight$start <- sapply(regions_to_highlight$start, FUN = function(start) {
          max(c(start, min(meth_plotting_data$Pos)))
        })
        regions_to_highlight$end <- sapply(regions_to_highlight$end, FUN = function(end) {
          min(c(end, max(meth_plotting_data$Pos)))
        })
        regions_to_highlight <- regions_to_highlight[regions_to_highlight$end > min(meth_plotting_data$Pos) &
                                                       regions_to_highlight$start < max(meth_plotting_data$Pos), ]
      }
      
      # Generate the plot showing smoothed methylation values.
      meth_plot <- ggplot(meth_plotting_data) +
        scale_x_continuous(name = paste0("Position (Chr. ", chromosome_for_plot, ")"),
                           limits = c(min(meth_plotting_data$Pos), max(meth_plotting_data$Pos)),
                           expand = c(0, 0)) +
        scale_y_continuous(name = "Smoothed\nmethylation",
                           limits = c(0, 1),
                           expand = expansion(mult = c(0.15, 0))) +
        aes(x = Pos, y = Meth, colour = Cell_type, group = Sample) +
        scale_colour_manual(breaks = names(rna_data_list), 
                            values = rna_colours,
                            name = "Cell type",
                            drop = F,
                            guide = guide_legend(override.aes = list(size = 6), # Make line legend appear as rectangles
                                                 nrow = ifelse(length(rna_data_list) > 2, yes = 2, no = 1))) 
      if (all(!(is.na(meth_highlight)))) {
        meth_plot <- meth_plot + 
          geom_rect(data = regions_to_highlight,
                    mapping = aes(xmin = start, xmax = end),
                    ymin = -0.15, ymax = 1, fill = "yellow", alpha = 0.5,
                    inherit.aes = F)
      }
      meth_plot <- meth_plot + 
        geom_line() +
        geom_rug(y = 0, sides = "b", colour = "black", length = unit(0.1, "npc")) +
        ggtitle(plot_titles[l]) +
        theme_classic() +
        theme(axis.text.y = element_text(colour = "black"),
              axis.ticks.y = element_line(colour = "black"),
              legend.position = "top",
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0.5))
      
      
      #-------------- scATAC-seq data -------------------
      
      # Identify all highlight regions that overlap with the region to plot.
      if (all(!(is.na(atac_highlight)))) {
        overlaps <- findOverlaps(query = atac_highlight,
                                 subject = region_to_plot,
                                 ignore.strand = T)
        regions_to_highlight <- as.data.frame(atac_highlight[from(overlaps)])
        regions_to_highlight$start <- sapply(regions_to_highlight$start, FUN = function(start) {
          max(c(start, min(meth_plotting_data$Pos)))
        })
        regions_to_highlight$end <- sapply(regions_to_highlight$end, FUN = function(end) {
          min(c(end, max(meth_plotting_data$Pos)))
        })
        regions_to_highlight <- regions_to_highlight[regions_to_highlight$end > min(meth_plotting_data$Pos) &
                                                       regions_to_highlight$start < max(meth_plotting_data$Pos), ]
      }
      
      # Restrict the data to those bins that overlap with CpGs in the extended plotting region.
      atac_subset_region <- GRanges(seqnames = chromosome_for_plot,
                                    ranges = IRanges(min(meth_plotting_data$Pos), max(meth_plotting_data$Pos)))
      atac_data_relevant <- lapply(cell_type_atac_data, FUN = function(x) {
        x[from(findOverlaps(x, atac_subset_region))]
      })
      
      # Prepare the data to plot scATAC-seq tracks.
      atac_plotting_data <- do.call(rbind, lapply(1:length(atac_data_relevant), FUN = function(x) {
        temp_df <- as.data.frame(atac_data_relevant[[x]])
        temp_df$Cell_type <- names(atac_data_relevant)[x]
        return(temp_df)
      }))
      atac_plotting_data$Cell_type <- factor(atac_plotting_data$Cell_type,
                                             levels = cell_types_for_plot_unique)
      atac_plotting_data$start <- sapply(atac_plotting_data$start, FUN = function(x) {
        max(x, min(meth_plotting_data$Pos))
      })
      atac_plotting_data$end <- sapply(atac_plotting_data$end, FUN = function(x) {
        min(x, max(meth_plotting_data$Pos))
      })
      
      # Generate the ATAC-seq signal plot.
      atac_plot <- ggplot(atac_plotting_data) +
        aes(xmin = start, xmax = end, ymax = score, fill = Cell_type) +
        scale_x_continuous(name = paste0("Position (Chr. ", chromosome_for_plot, ")"),
                           limits = c(min(meth_plotting_data$Pos), max(meth_plotting_data$Pos)),
                           expand = c(0, 0)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                           name = "Binned scATAC-seq\nread count [RPKM]") +
        scale_fill_manual(breaks = levels(atac_plotting_data$Cell_type),
                          values = cell_type_col_palette,
                          name = "Cell type")
      if (all(!(is.na(atac_highlight)))) {
        atac_plot <- atac_plot +
          geom_rect(data = regions_to_highlight,
                    mapping = aes(xmin = start, xmax = end),
                    ymin = 0, ymax = 1.1 * max(atac_plotting_data$score), fill = "yellow", alpha = 0.5,
                    inherit.aes = F)
      }
      atac_plot <- atac_plot +
        geom_rect(ymin = 0) +
        facet_wrap(~Cell_type, ncol = 1) +
        theme_classic() +
        theme(axis.text.y = element_text(colour = "black"),
              axis.ticks.y = element_line(colour = "black"),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              strip.text = element_blank(),
              strip.background = element_blank(),
              legend.position = "none")
      
      
      #-------------- RNA track -------------------
      
      # Identify all highlight regions that overlap with the region to plot.
      if (all(!(is.na(rna_highlight)))) {
        overlaps <- findOverlaps(query = rna_highlight,
                                 subject = region_to_plot,
                                 ignore.strand = T)
        regions_to_highlight <- as.data.frame(rna_highlight[from(overlaps)])
        regions_to_highlight$start <- sapply(regions_to_highlight$start, FUN = function(start) {
          max(c(start, min(meth_plotting_data$Pos)))
        })
        regions_to_highlight$end <- sapply(regions_to_highlight$end, FUN = function(end) {
          min(c(end, max(meth_plotting_data$Pos)))
        })
        regions_to_highlight <- regions_to_highlight[regions_to_highlight$end > min(meth_plotting_data$Pos) &
                                                       regions_to_highlight$start < max(meth_plotting_data$Pos), ]
      }
      
      # Prepare the data to plot RNA tracks.
      rna_df_list <- lapply(1:length(rna_data_list), FUN = function(m) {
        region_data <- rna_data_list[[m]][from(findOverlaps(rna_data_list[[m]], atac_subset_region))]
        df_prelim <- as.data.frame(region_data)
        df_prelim$Cell_type <- names(rna_data_list)[m]
        return(df_prelim)
      })
      rna_plotting_data <- do.call(rbind, rna_df_list)
      rna_plotting_data$Cell_type <- factor(rna_plotting_data$Cell_type,
                                            levels = names(rna_data_list))
      rna_plotting_data$start <- sapply(rna_plotting_data$start, FUN = function(x) {
        max(x, min(meth_plotting_data$Pos))
      })
      rna_plotting_data$end <- sapply(rna_plotting_data$end, FUN = function(x) {
        min(x, max(meth_plotting_data$Pos))
      })
      
      # Generate the RNA track plot.
      rna_plot <- ggplot(rna_plotting_data) +
        aes(xmin = start, xmax = end, ymax = score, fill = Cell_type) +
        scale_x_continuous(name = paste0("Position (Chr. ", chromosome_for_plot, ")"),
                           limits = c(min(meth_plotting_data$Pos), max(meth_plotting_data$Pos)),
                           expand = c(0, 0)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                           name = "Binned RNA-seq\nread count [RPKM]") +
        scale_fill_manual(breaks = levels(rna_plotting_data$Cell_type),
                          values = rna_colours,
                          name = "Cell type")
      if (all(!(is.na(rna_highlight)))) {
        rna_plot <- rna_plot +
          geom_rect(data = regions_to_highlight,
                    mapping = aes(xmin = start, xmax = end),
                    ymin = 0, ymax = 1.1 * max(rna_plotting_data$score), fill = "yellow", alpha = 0.5,
                    inherit.aes = F)
      }
      rna_plot <- rna_plot +
        geom_rect(ymin = 0) +
        facet_wrap(~Cell_type, ncol = 1) +
        theme_classic() +
        theme(axis.text.y = element_text(colour = "black"),
              axis.ticks.y = element_line(colour = "black"),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              strip.text = element_blank(),
              strip.background = element_blank(),
              legend.position = "none")
      
      
      #-------------- Custom tracks (if requested) ----------------
      
      if (!(all(is.na(add_track_bws)))) {
        
        if (all(!(is.na(add_track_highlight)))) {
          overlaps <- findOverlaps(query = add_track_highlight,
                                   subject = region_to_plot,
                                   ignore.strand = T)
          regions_to_highlight <- as.data.frame(add_track_highlight[from(overlaps)])
          regions_to_highlight$start <- sapply(regions_to_highlight$start, FUN = function(start) {
            max(c(start, min(meth_plotting_data$Pos)))
          })
          regions_to_highlight$end <- sapply(regions_to_highlight$end, FUN = function(end) {
            min(c(end, max(meth_plotting_data$Pos)))
          })
          regions_to_highlight <- regions_to_highlight[regions_to_highlight$end > min(meth_plotting_data$Pos) &
                                                         regions_to_highlight$start < max(meth_plotting_data$Pos), ]
        }
        
        # Restrict the data to those bins that overlap with CpGs in the extended plotting region.
        add_track_subset_region <- GRanges(seqnames = chromosome_for_plot,
                                           ranges = IRanges(min(meth_plotting_data$Pos), max(meth_plotting_data$Pos)))
        add_track_data <- lapply(additional_tracks, FUN = function(x) {
          x_restr <- x[from(findOverlaps(x, add_track_subset_region))]
          if (length(x_restr) == 0) {
            x_restr <- add_track_subset_region
            x_restr$score <- 0
          }
          return(x_restr)
        })
        
        # Prepare the data to plot additional tracks.
        # If no bins are present for the plotting region, use a placeholder with a score equal to 0.
        add_track_plotting_data <- as.data.frame(add_track_data[[1]])
        add_track_plotting_data$Name <- add_track_names[1]
        if (length(additional_tracks) > 1) {
          for (new_index in 2:length(additional_tracks)) {
            add_track_plotting_data_temp <- as.data.frame(add_track_data[[new_index]])
            add_track_plotting_data_temp$Name <- add_track_names[new_index]
            add_track_plotting_data <- rbind(add_track_plotting_data, add_track_plotting_data_temp)
          }
        }
        add_track_plotting_data$Name <- factor(add_track_plotting_data$Name,
                                               levels = add_track_names)
        add_track_plotting_data$start <- sapply(add_track_plotting_data$start, FUN = function(x) {
          max(x, min(meth_plotting_data$Pos))
        })
        add_track_plotting_data$end <- sapply(add_track_plotting_data$end, FUN = function(x) {
          min(x, max(meth_plotting_data$Pos))
        })
        
        # Generate the  plot.
        add_track_plot <- ggplot(add_track_plotting_data) +
          aes(xmin = start, xmax = end, ymax = score) +
          scale_x_continuous(name = paste0("Position (Chr. ", chromosome_for_plot, ")"),
                             limits = c(min(meth_plotting_data$Pos), max(meth_plotting_data$Pos)),
                             expand = c(0, 0)) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                             name = "Signal",
                             limits = c(0, max(add_track_plotting_data$score)))
        if (all(!(is.na(add_track_highlight)))) {
          add_track_plot <- add_track_plot +
            geom_rect(data = regions_to_highlight,
                      mapping = aes(xmin = start, xmax = end),
                      ymin = 0, ymax = 1.1 * max(add_track_plotting_data$score), fill = "yellow", alpha = 0.5,
                      inherit.aes = F)
        }
        add_track_plot <- add_track_plot +
          geom_rect(ymin = 0, fill = "black") +
          facet_wrap(~Name, ncol = 1, strip.position = "right") +
          theme_classic() +
          theme(axis.text.y = element_text(colour = "black"),
                axis.ticks.y = element_line(colour = "black"),
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.line.x = element_blank(),
                strip.text = element_text(size = 6))
        
      }
      
      
      #-------------- Gene track -------------------
      
      # Identify all highlight regions that overlap with the region to plot.
      if (all(!(is.na(gene_highlight)))) {
        overlaps <- findOverlaps(query = gene_highlight,
                                 subject = region_to_plot,
                                 ignore.strand = T)
        regions_to_highlight <- as.data.frame(gene_highlight[from(overlaps)])
        regions_to_highlight$start <- sapply(regions_to_highlight$start, FUN = function(start) {
          max(c(start, min(meth_plotting_data$Pos)))
        })
        regions_to_highlight$end <- sapply(regions_to_highlight$end, FUN = function(end) {
          min(c(end, max(meth_plotting_data$Pos)))
        })
        regions_to_highlight <- regions_to_highlight[regions_to_highlight$end > min(meth_plotting_data$Pos) &
                                                       regions_to_highlight$start < max(meth_plotting_data$Pos), ]
      }
      
      # Identify genes that overlap with the extended region of interest.
      genes_in_extend_region_gr <- genes_gr[to(findOverlaps(regions_extend_gr[l], genes, type = "any"))]
      genes_in_extend_region <- genes_in_extend_region_gr$name
      
      # Prepare data on the genes and exons in the extended region of interest for plotting.
      if (length(genes_in_extend_region_gr) > 0) {
        y_pos_genes <- sapply(as.character(strand(genes_in_extend_region_gr)), FUN = function(x) {
          ifelse(x == "+", yes = 0.5, no = -0.5)
        })
        
        region_genes_df <- data.frame(x = sapply(1:length(genes_in_extend_region_gr), FUN = function(y) {
          max(start(genes_in_extend_region_gr[y]), min(meth_plotting_data$Pos))
        }),
        xend = sapply(1:length(genes_in_extend_region_gr), FUN = function(y) {
          min(end(genes_in_extend_region_gr[y]), max(meth_plotting_data$Pos))
        }),
        y = y_pos_genes,
        yend = y_pos_genes) 
        
        exons_in_extend_region <- exons[exons$gene_name %in% genes_in_extend_region_gr$name &
                                          exons$end > min(meth_plotting_data$Pos) &
                                          exons$start < max(meth_plotting_data$Pos), ]
        
        if (nrow(exons_in_extend_region) > 0) {
          region_exons_df <- data.frame(xmin = sapply(exons_in_extend_region$start, FUN = function(y) {
            max(y, min(meth_plotting_data$Pos))
          }),
          xmax = sapply(exons_in_extend_region$end, FUN = function(y) {
            min(y, max(meth_plotting_data$Pos))
          }),
          ymin = sapply(exons_in_extend_region$strand, FUN = function(x) {
            ifelse(x == "+", yes = 0.3, no = -0.7)
          }),
          ymax = sapply(exons_in_extend_region$strand, FUN = function(x) {
            ifelse(x == "+", yes = 0.7, no = -0.3)
          }))
        }
        gene_names_df <- data.frame(label = genes_in_extend_region_gr$name,
                                    x = sapply(1:nrow(region_genes_df), FUN = function(j) {
                                      mean(c(region_genes_df$x[j], region_genes_df$xend[j]))
                                    }),
                                    y = sapply(as.character(strand(genes_in_extend_region_gr)), FUN = function(x) {
                                      ifelse(x == "+", yes = 0.15, no = -0.15)
                                    }))
        tss_in_extend_region <- tss_pos[tss_pos$gene %in% genes_in_extend_region_gr$name &
                                          tss_pos$pos >= min(meth_plotting_data$Pos) &
                                          tss_pos$pos <= max(meth_plotting_data$Pos), ]
        if (nrow(tss_in_extend_region) > 0) {
          region_tss_plotting_data_1 <- data.frame(x = tss_in_extend_region$pos,
                                                   xend = tss_in_extend_region$pos,
                                                   y = sapply(1:nrow(tss_in_extend_region), FUN = function(x) {
                                                     ifelse(tss_in_extend_region$strand[x] == "+", yes = 0.5, no = -0.5)
                                                   }),
                                                   yend = sapply(1:nrow(tss_in_extend_region), FUN = function(x) {
                                                     ifelse(tss_in_extend_region$strand[x] == "+", yes = 0.85, no = -0.85)
                                                   }))
          region_tss_plotting_data_2 <- data.frame(x = tss_in_extend_region$pos,
                                                   xend = sapply(1:nrow(tss_in_extend_region), FUN = function(x) {
                                                     ifelse(tss_in_extend_region$strand[x] == "+", 
                                                            yes = tss_in_extend_region$pos[x] + 0.075 * diff(range(meth_plotting_data$Pos)),
                                                            no = tss_in_extend_region$pos[x] - 0.075 * diff(range(meth_plotting_data$Pos))
                                                     )
                                                   }),
                                                   y = sapply(1:nrow(tss_in_extend_region), FUN = function(x) {
                                                     ifelse(tss_in_extend_region$strand[x] == "+", yes = 0.85, no = -0.85)
                                                   }),
                                                   yend = sapply(1:nrow(tss_in_extend_region), FUN = function(x) {
                                                     ifelse(tss_in_extend_region$strand[x] == "+", yes = 0.85, no = -0.85)
                                                   }))
        }
      }
      
      # Generate a gene annotaton track.
      gene_plot <- ggplot() +
        scale_y_continuous(limits = c(-1, 1)) +
        scale_x_continuous(name = paste0("Position (Chr. ", chromosome_for_plot, ")"),
                           limits = c(min(meth_plotting_data$Pos), max(meth_plotting_data$Pos)),
                           expand = c(0, 0))
      if (all(!(is.na(gene_highlight)))) {
        gene_plot <- gene_plot +
          geom_rect(data = regions_to_highlight,
                    mapping = aes(xmin = start, xmax = end),
                    ymin = -1.2, ymax = 1.2, fill = "yellow", alpha = 0.5)
      }
      if (length(genes_in_extend_region_gr) > 0) {
        gene_plot <- gene_plot + 
          geom_segment(data = region_genes_df,
                       mapping = aes(x = x, y = y, xend = xend, yend = yend),
                       colour = "black") +
          geom_text(data = gene_names_df,
                    mapping = aes(x = x, y = y, label = label),
                    size = 3)
        if (nrow(exons_in_extend_region) > 0) {
          gene_plot <- gene_plot +
            geom_rect(data = region_exons_df,
                      mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                      fill = "grey50")
        }
        if (nrow(tss_in_extend_region) > 0) {
          gene_plot <- gene_plot +
            geom_segment(data = region_tss_plotting_data_1,
                         mapping = aes(x = x, y = y, xend = xend, yend = yend),
                         size = 1, lineend = "square") +
            geom_segment(data = region_tss_plotting_data_2,
                         mapping = aes(x = x, y = y, xend = xend, yend = yend),
                         arrow = arrow(type = "closed", angle = 20, length = unit(0.1, "inches")),
                         size = 1)
        }
      }
      gene_plot <- gene_plot +
        theme_classic() +
        theme(axis.text.x = element_text(colour = "black"),
              axis.ticks.x = element_line(colour = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank())
      
      
      #-------------- Combination of plots --------------------
      
      # Combine the plots into one and return the combined plot.
      if (!(all(is.na(add_track_bws)))) {
        combined_plot <- plot_grid(meth_plot, atac_plot, rna_plot, add_track_plot, gene_plot, 
                                   ncol = 1,
                                   axis = "lr",
                                   align = "v",
                                   rel_heights = c(2.5, 1.5, 1.5, 3, 1.5))
      } else {
        combined_plot <- plot_grid(meth_plot, atac_plot, rna_plot, gene_plot, 
                                   ncol = 1,
                                   axis = "l",
                                   align = "v",
                                   rel_heights = c(2.5, 1.75, 1.75, 1.5))
      }
      
      # Save the combined plot at a preliminary location.
      pdf(paste0(outfile_pdf, "_", l, ".pdf"), width = 5.5, height = ifelse(!(all(is.na(add_track_bws))), yes = 10, no = 8))
      print(combined_plot)
      dev.off()
      
      # Return an empty string so as to not inflate memory.
      return("")
      
    }, mc.cores = ncores)
    
    # Combine the single PDF files into one. Do so in chunks of 10,000 plots.
    # Remove temporary files.
    if (length(regions_gr) > 1) {
      counter <- 0
      for (l in seq(1, length(regions_gr), 10000)) {
        counter <- counter + 1
        files_to_merge <- paste0(outfile_pdf, "_", l:min((l + 9999), length(regions_gr)), ".pdf")
        chunk_outfile_name <- paste0(outfile_pdf, "_intermediate_", counter, ".pdf")
        pdf_combine(files_to_merge, output = chunk_outfile_name)
        file.remove(files_to_merge)
      }
      comb_files_to_merge <- paste0(outfile_pdf, "_intermediate_", 1:counter, ".pdf")
      pdf_combine(comb_files_to_merge, output = outfile_pdf)
      file.remove(comb_files_to_merge)
    } else {
      file.rename(from = paste0(outfile_pdf, "_1.pdf"),
                  to = outfile_pdf)
    }
    
  } 





# Generate a GRanges object containing the overlapping regions between 
# the DMR and the peak in each DMR-peak-gene link.
meth_gr <- makeGRangesFromDataFrame(do.call(
  rbind, lapply(1:nrow(combined_results), FUN = function(x) {
    corresp_reg_id <- combined_results$Methylation__region_ID[x]
    return(meth_reg[meth_reg$region_ID == corresp_reg_id, 1:3])
  })
))
acc_gr <- StringToGRanges(combined_results$Accessibility__peak_ID,
                          starts.in.df.are.0based = T)
seqlevelsStyle(acc_gr) <- "NCBI"
overl_gr <- pintersect(meth_gr, acc_gr)
assert(all(width(overl_gr) >= 1))

# Generate track plots for the DMR-peak-gene links. In the methylation track,
# highlight all DMRs in DMR-peak-gene links. In the ATAC track, highlight all
# peaks in DMR-peak-gene links. In the gene track, highlight the overlap 
# between the corresponding DMR and the corresponding peak.
plot_meth_and_other_signal_many_regions_custom_titles_custom_highlight_many_cell_types(
  bsseq_obj = meth_data,
  regions_gr = overl_gr,
  extend = 3000,
  genes = genes_gr,
  exons = exons,
  tss_pos = tss_df,
  atac_bigwig_path_table = atac_seq_data_overview,
  atac_bigwig_prefix = b330_space,
  rna_bigwig_path_table = rna_sample_mapping,
  add_track_bws = NA,
  add_track_names = NA,
  meth_highlight = meth_gr,
  atac_highlight = acc_gr,
  rna_highlight = NA,
  gene_highlight = overl_gr,
  add_track_highlight = NA,
  plot_titles = paste(combined_results$Methylation__region_ID, 
                      combined_results$Accessibility__peak_ID,
                      combined_results$Expression__gene_name,
                      sep = ";\n"),
  outfile_pdf = paste0(
    location, 
    "/treg_hierarchies/diff_meth_diff_acc_diff_expr_dmr_peak_gene_links_trackplots_w_tfidf.pdf"
  ),
  ncores = 10
)