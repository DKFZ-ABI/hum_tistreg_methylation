# This script generates a gerenal expression heat map (across all TEs and cell
#   types) based on the counts returned by TEspeX.
# Author: Niklas Beumer



# Load required package(s).
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(testit)


# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the RNA sample mapping file.
sample_mapping_file <- paste0(location, 
                              "/sample_mapping_rnaseq_only_biol_rep.txt")
sample_mapping <- read.table(sample_mapping_file, header = T, 
                             stringsAsFactors = F, sep = "\t")

# Read in the consensus TE sequences.
cons_seq_file <- paste0(
  location, 
  "/external_data/Giri_repbase_data_2023-11-03/Human_TEs_identified_by_browsing.fasta"
)
cons_seq <- readLines(cons_seq_file)

# Read in the file containing expression counts returned by TEspeX.
counts_file <- paste0(
  location, 
  "/te_analysis/cons_seq_tespex_results_bulk_rna_all/outfile.txt"
)
counts <- read.table(counts_file, header = T, stringsAsFactors = F)

# Define the relevant cell types.
relevant_cell_types <- c("Skin Treg", "Skin Tconv", "Blood CCR8+ Treg", 
                         "Blood naive Treg", "Blood naive Tconv")

# Specify a colour palette for the cell types.
cell_type_col_palette <- c("blue", "cyan2", "orange", "darkorchid1", "red")
names(cell_type_col_palette) <- relevant_cell_types

# Specify the TEs that should be highlighted in the heat maps.
tes_to_highlight <- c("HERVIP10F", "LTR45B", "L1PA12", "MER4D")





############################################################################
# Compute TPM values for the TEs.
############################################################################

# For each TE, compute the sequence length.
sequence_lengths <- data.frame(
  TE = counts$TE,
  length = sapply(counts$TE, FUN = function(x) {
    start_line_ind <- grep(paste0(">", x, "\t"), cons_seq)
    assert(length(start_line_ind) == 1)
    end_line_ind <- start_line_ind
    found <- F
    while(!found) {
      end_line_ind <- end_line_ind + 1
      if (end_line_ind > length(cons_seq)) {
        found <- T
      } else {
        found <- grepl(">", cons_seq[end_line_ind])
      }
    }
    sequence <- paste0(cons_seq[(start_line_ind + 1):(end_line_ind - 1)],
                       collapse = "")
    return(nchar(sequence))
  })
)

# Normalise the count values to transcripts per million (TPM). The base line
# here is the total number of TE "transcripts".
counts_matr <- as.matrix(counts[, 2:ncol(counts)])
rownames(counts_matr) <- counts$TE
assert(all(rownames(counts_matr) == sequence_lengths$TE))
counts_matr_norm_length <- apply(counts_matr, 2, FUN = function(x) {
  x / sequence_lengths$length
})
total_colsums <- colSums(counts_matr_norm_length)
assert(all(names(total_colsums) == colnames(counts_matr_norm_length)))
tpm_matr <- t(apply(counts_matr_norm_length, 1, FUN = function(x) {
  10^6 * x / total_colsums
}))

# Save the TPM table.
tpm_outfile <- paste0(
  location, 
  "/te_analysis/cons_seq_tespex_results_bulk_rna_all_tpm_values.txt"
)
write.table(tpm_matr, file = tpm_outfile, sep = "\t", col.names = T, 
            row.names = T, quote = F)







#############################################################################
# Generate heat maps showing the expression values of each TE in on the
# sample level or cell type level.
#############################################################################

# Transform TPM values by log(100*TPM+1).
tpm_log <- log1p(100 * tpm_matr)

# Compute cell-type-wise average values of normal TPM values.
# Take into account that one of the blood CCR8+ Treg samples had a typo in the
# Fastq file name.
tpm_celltype <- t(sapply(relevant_cell_types, FUN = function(x) {
  relevant_samples <- sample_mapping$Sample[
    sample_mapping$Cell_type == x]
  relevant_samples <- gsub("CCR8_Treg_Yellow5", "CCRR_Treg_Yellow5",
                           relevant_samples)
  relevant_colnames <- sapply(relevant_samples, FUN = function(y) {
    grep(y, colnames(tpm_matr), value = T)
  })
  means <- rowMeans(tpm_matr[, relevant_colnames])
  return(means)
}))

# Compute cell-type-wise average values of log-transformed TPM values.
# Take into account that one of the blood CCR8+ Treg samples had a typo in the
# Fastq file name.
tpm_log_celltype <- t(sapply(relevant_cell_types, FUN = function(x) {
  relevant_samples <- sample_mapping$Sample[
    sample_mapping$Cell_type == x]
  relevant_samples <- gsub("CCR8_Treg_Yellow5", "CCRR_Treg_Yellow5",
                           relevant_samples)
  relevant_colnames <- sapply(relevant_samples, FUN = function(y) {
    grep(y, colnames(tpm_log), value = T)
  })
  means <- rowMeans(tpm_log[, relevant_colnames])
  return(means)
}))



############# Sample-level: Raw counts.

# Prepare the matrix for heatmap generation.
relevant_samples <- unlist(lapply(relevant_cell_types, FUN = function(x) {
  sample_mapping$Sample[sample_mapping$Cell_type %in% x]
}))
relevant_samples <- gsub("CCR8_Treg_Yellow5", "CCRR_Treg_Yellow5", 
                         relevant_samples)
cols_to_consider <- sapply(relevant_samples, FUN = function(x) {
  grep(x, colnames(counts_matr), value = T)
})
relevant_matr <- t(counts_matr[, cols_to_consider])

# Get a clustering for the columns.
relevant_matr_dist <- dist(t(relevant_matr))
relevant_matr_hclust <- hclust(relevant_matr_dist)

# Replace 0 values by NA. This is to facilitate printing of the heat map.
relevant_matr[relevant_matr == 0] <- NA

# Generate and save a heat map showing expression values for all TEs in all cell 
# types (on the sample level).
row_anno_df <- data.frame(
  Cell_type = sapply(rownames(relevant_matr), FUN = function(x) {
    x_split <- strsplit(x, split = ".fastq")[[1]][1]
    x_split <- gsub("CCRR_Treg_Yellow5", "CCR8_Treg_Yellow5", x_split)
    corresp_cell_type <- sample_mapping$Cell_type[
      sample_mapping$Sample == x_split]
    return(corresp_cell_type)
  }))
row_anno_df$Cell_type <- factor(row_anno_df$Cell_type, 
                                levels = relevant_cell_types)
row_anno_obj <- rowAnnotation(
  df = row_anno_df, 
  annotation_label = "Cell type",
  col = list(Cell_type = cell_type_col_palette),
  annotation_legend_param = list(title_gp = gpar(fontsize = 17, 
                                                 fontface = "bold"), 
                                 labels_gp = gpar(fontsize = 15),
                                 title = "Cell type"))
mark_anno_obj <- HeatmapAnnotation(
  placeholder = anno_mark(at = sapply(tes_to_highlight, FUN = function(x) {
    which(colnames(relevant_matr) == x)
  }), 
  labels = tes_to_highlight,
  side = "bottom",
  labels_gp = gpar(fontsize = 10))
)
expr_heatmap_sample <- Heatmap(
  relevant_matr,
  col = colorRamp2(seq(min(relevant_matr, na.rm = T), 
                       max(relevant_matr, na.rm = T), 
                       length.out = 200), 
                   mako(200)),
  cluster_rows = F,
  cluster_columns = relevant_matr_hclust,
  show_row_names = F,
  show_column_names = F,
  left_annotation = row_anno_obj,
  bottom_annotation = mark_anno_obj,
  heatmap_legend_param = list(
    title = "Expression\n[raw counts]\n", 
    labels_gp = gpar(fontsize = 15), 
    title_gp = gpar(fontsize = 17, fontface = "bold"))
)
heatmap_outfile_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_cons_seq_expr_quant_tespex_te_levels_heatmap_rawcounts_sample_level.pdf"
)
heatmap_outfile_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_cons_seq_expr_quant_tespex_te_levels_heatmap_rawcounts_sample_level.rds"
)
pdf(heatmap_outfile_pdf, width = 15, height = 5)
draw(expr_heatmap_sample)
dev.off()
saveRDS(expr_heatmap_sample, file = heatmap_outfile_rds)




############# Sample-level: Normal TPM values.

# Prepare the matrix for heatmap generation.
relevant_samples <- unlist(lapply(relevant_cell_types, FUN = function(x) {
  sample_mapping$Sample[sample_mapping$Cell_type %in% x]
}))
relevant_samples <- gsub("CCR8_Treg_Yellow5", "CCRR_Treg_Yellow5", 
                         relevant_samples)
cols_to_consider <- sapply(relevant_samples, FUN = function(x) {
  grep(x, colnames(tpm_matr), value = T)
})
relevant_matr <- t(tpm_matr[, cols_to_consider])

# Get a clustering for the columns.
relevant_matr_dist <- dist(t(relevant_matr))
relevant_matr_hclust <- hclust(relevant_matr_dist)

# Replace 0 values by NA. This is to facilitate printing of the heat map.
relevant_matr[relevant_matr == 0] <- NA

# Generate and save a heat map showing expression values for all TEs in all cell 
# types (on the sample level).
row_anno_df <- data.frame(
  Cell_type = sapply(rownames(relevant_matr), FUN = function(x) {
    x_split <- strsplit(x, split = ".fastq")[[1]][1]
    x_split <- gsub("CCRR_Treg_Yellow5", "CCR8_Treg_Yellow5", x_split)
    corresp_cell_type <- sample_mapping$Cell_type[
      sample_mapping$Sample == x_split]
    return(corresp_cell_type)
  }))
row_anno_df$Cell_type <- factor(row_anno_df$Cell_type, 
                                levels = relevant_cell_types)
row_anno_obj <- rowAnnotation(
  df = row_anno_df, 
  annotation_label = "Cell type",
  col = list(Cell_type = cell_type_col_palette),
  annotation_legend_param = list(title_gp = gpar(fontsize = 17, 
                                                 fontface = "bold"), 
                                 labels_gp = gpar(fontsize = 15),
                                 title = "Cell type"))
mark_anno_obj <- HeatmapAnnotation(
  placeholder = anno_mark(at = sapply(tes_to_highlight, FUN = function(x) {
    which(colnames(relevant_matr) == x)
  }), 
  labels = tes_to_highlight,
  side = "bottom",
  labels_gp = gpar(fontsize = 10))
)
expr_heatmap_sample <- Heatmap(
  relevant_matr,
  col = colorRamp2(seq(min(relevant_matr, na.rm = T), 
                       max(relevant_matr, na.rm = T), 
                       length.out = 200), 
                   mako(200)),
  cluster_rows = F,
  cluster_columns = relevant_matr_hclust,
  show_row_names = F,
  show_column_names = F,
  left_annotation = row_anno_obj,
  bottom_annotation = mark_anno_obj,
  heatmap_legend_param = list(
    title = "Expression\n[TPM]\n", 
    labels_gp = gpar(fontsize = 15), 
    title_gp = gpar(fontsize = 17, fontface = "bold"))
)
heatmap_outfile_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_cons_seq_expr_quant_tespex_te_levels_heatmap_tpm_sample_level.pdf"
)
heatmap_outfile_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_cons_seq_expr_quant_tespex_te_levels_heatmap_tpm_sample_level.rds"
)
pdf(heatmap_outfile_pdf, width = 15, height = 5)
draw(expr_heatmap_sample)
dev.off()
saveRDS(expr_heatmap_sample, file = heatmap_outfile_rds)




############# Sample-level: Log-transformed TPM values.

# Prepare the matrix for heatmap generation.
relevant_samples <- unlist(lapply(relevant_cell_types, FUN = function(x) {
  sample_mapping$Sample[sample_mapping$Cell_type %in% x]
}))
relevant_samples <- gsub("CCR8_Treg_Yellow5", "CCRR_Treg_Yellow5", 
                         relevant_samples)
cols_to_consider <- sapply(relevant_samples, FUN = function(x) {
  grep(x, colnames(tpm_log), value = T)
})
relevant_matr <- t(tpm_log[, cols_to_consider])

# Get a clustering for the columns.
relevant_matr_dist <- dist(t(relevant_matr))
relevant_matr_hclust <- hclust(relevant_matr_dist)

# Replace 0 values by NA. This is to facilitate printing of the heat map.
relevant_matr[relevant_matr == 0] <- NA

# Generate and save a heat map showing expression values for all TEs in all cell 
# types (on the sample level).
row_anno_df <- data.frame(
  Cell_type = sapply(rownames(relevant_matr), FUN = function(x) {
    x_split <- strsplit(x, split = ".fastq")[[1]][1]
    x_split <- gsub("CCRR_Treg_Yellow5", "CCR8_Treg_Yellow5", x_split)
    corresp_cell_type <- sample_mapping$Cell_type[
      sample_mapping$Sample == x_split]
    return(corresp_cell_type)
  }))
row_anno_df$Cell_type <- factor(row_anno_df$Cell_type, 
                                levels = relevant_cell_types)
row_anno_obj <- rowAnnotation(
  df = row_anno_df, 
  annotation_label = "Cell type",
  col = list(Cell_type = cell_type_col_palette),
  annotation_legend_param = list(title_gp = gpar(fontsize = 17, 
                                                 fontface = "bold"), 
                                 labels_gp = gpar(fontsize = 15),
                                 title = "Cell type"))
mark_anno_obj <- HeatmapAnnotation(
  placeholder = anno_mark(at = sapply(tes_to_highlight, FUN = function(x) {
    which(colnames(relevant_matr) == x)
  }), 
  labels = tes_to_highlight,
  side = "bottom",
  labels_gp = gpar(fontsize = 10))
)
expr_heatmap_sample <- Heatmap(
  relevant_matr,
  col = colorRamp2(seq(min(relevant_matr, na.rm = T), 
                       max(relevant_matr, na.rm = T), 
                       length.out = 200), 
                   mako(200)),
  cluster_rows = F,
  cluster_columns = relevant_matr_hclust,
  show_row_names = F,
  show_column_names = F,
  left_annotation = row_anno_obj,
  bottom_annotation = mark_anno_obj,
  heatmap_legend_param = list(
    title = "Expression\n[log(100*TPM+1)]\n", 
    labels_gp = gpar(fontsize = 15), 
    title_gp = gpar(fontsize = 17, fontface = "bold"))
)
heatmap_outfile_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_cons_seq_expr_quant_tespex_te_levels_heatmap_logtpm_sample_level.pdf"
)
heatmap_outfile_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_cons_seq_expr_quant_tespex_te_levels_heatmap_logtpm_sample_level.rds"
)
pdf(heatmap_outfile_pdf, width = 15, height = 5)
draw(expr_heatmap_sample)
dev.off()
saveRDS(expr_heatmap_sample, file = heatmap_outfile_rds)




############# Cell-type-level: Normal TPM values.

# Get a clustering for the columns.
relevant_matr <- tpm_celltype
relevant_matr_dist <- dist(t(relevant_matr))
relevant_matr_hclust <- hclust(relevant_matr_dist)

# Replace 0 values by NA. This is to facilitate printing of the heat map.
relevant_matr[relevant_matr == 0] <- NA

# Generate and save a heat map showing expression values for all TEs in all cell 
# types (on the cell type level).
mark_anno_obj <- HeatmapAnnotation(
  placeholder = anno_mark(at = sapply(tes_to_highlight, FUN = function(x) {
    which(colnames(relevant_matr) == x)
  }), 
  labels = tes_to_highlight,
  side = "bottom",
  labels_gp = gpar(fontsize = 10))
)
expr_heatmap_celltype <- Heatmap(
  relevant_matr,
  col = colorRamp2(seq(min(relevant_matr, na.rm = T), 
                       max(relevant_matr, na.rm = T), 
                       length.out = 200), 
                   mako(200)),
  cluster_rows = F,
  cluster_columns = relevant_matr_hclust,
  show_row_names = T,
  show_column_names = F,
  bottom_annotation = mark_anno_obj,
  heatmap_legend_param = list(
    title = "Mean\nexpression\n[TPM]\n", 
    labels_gp = gpar(fontsize = 15), 
    title_gp = gpar(fontsize = 17, fontface = "bold"))
)
heatmap_outfile_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_cons_seq_expr_quant_tespex_te_levels_heatmap_tpm_celltype_level.pdf"
)
heatmap_outfile_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_cons_seq_expr_quant_tespex_te_levels_heatmap_tpm_celltype_level.rds"
)
pdf(heatmap_outfile_pdf, width = 15, height = 3)
draw(expr_heatmap_celltype)
dev.off()
saveRDS(expr_heatmap_celltype, file = heatmap_outfile_rds)




############# Cell-type-level: Log-transformed TPM values.

# Get a clustering for the columns.
relevant_matr <- tpm_log_celltype
relevant_matr_dist <- dist(t(relevant_matr))
relevant_matr_hclust <- hclust(relevant_matr_dist)

# Replace 0 values by NA. This is to facilitate printing of the heat map.
relevant_matr[relevant_matr == 0] <- NA

# Generate and save a heat map showing expression values for all TEs in all cell 
# types (on the cell type level).
mark_anno_obj <- HeatmapAnnotation(
  placeholder = anno_mark(at = sapply(tes_to_highlight, FUN = function(x) {
    which(colnames(relevant_matr) == x)
  }), 
  labels = tes_to_highlight,
  side = "bottom",
  labels_gp = gpar(fontsize = 10))
)
expr_heatmap_celltype <- Heatmap(
  relevant_matr,
  col = colorRamp2(seq(min(relevant_matr, na.rm = T), 
                       max(relevant_matr, na.rm = T), 
                       length.out = 200), 
                   mako(200)),
  cluster_rows = F,
  cluster_columns = relevant_matr_hclust,
  show_row_names = T,
  show_column_names = F,
  bottom_annotation = mark_anno_obj,
  heatmap_legend_param = list(
    title = "Mean\nexpression\n[log(100*TPM+1)]\n", 
    labels_gp = gpar(fontsize = 15), 
    title_gp = gpar(fontsize = 17, fontface = "bold"))
)
heatmap_outfile_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_cons_seq_expr_quant_tespex_te_levels_heatmap_logtpm_celltype_level.pdf"
)
heatmap_outfile_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_cons_seq_expr_quant_tespex_te_levels_heatmap_logtpm_celltype_level.rds"
)
pdf(heatmap_outfile_pdf, width = 15, height = 3)
draw(expr_heatmap_celltype)
dev.off()
saveRDS(expr_heatmap_celltype, file = heatmap_outfile_rds)

