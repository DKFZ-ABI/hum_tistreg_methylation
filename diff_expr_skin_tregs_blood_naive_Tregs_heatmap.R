# This script generates a heat map for those genes that are differentially 
#   expressed between skin Treg cells and blood CD45RA+ Treg cells.
# Author: Niklas Beumer



# Load required packages.
library(ComplexHeatmap)
library(circlize)
library(viridis)


# Specify a location on /xxx.
location <- "/xxx/nbeumer/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/yyy/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Specify the cell types.
cell_types_my_spelling <- c("Skin_Treg", "Blood_naive_Treg")
cell_types_rna_spelling <- c("Skin_Treg", "Blood_naive_Treg")

# Read in the RNA sample mapping file and restrict to the relevant cell types.
rna_sample_mapping_path <- paste0(location, 
                                  "/sample_mapping_rnaseq_only_biol_rep.txt")
rna_sample_mapping <- read.table(rna_sample_mapping_path, header = T, 
                                 stringsAsFactors = F, sep = "\t")
rna_sample_mapping <- rna_sample_mapping[
  rna_sample_mapping$Cell_type_wo_blank %in% cell_types_rna_spelling, 
]

# Read in the table with differential gene expression results between
# skin Tregs and blood naive Tregs.
rna_sig_file <- paste0(
  location, 
  "/RNASeq/analysis_results/2022-01-14_diff_gene_expr_DESEq2_Skin_TregBlood_naive_Treg_results_filtered_with_significance.txt"
)
rna_sig <- read.table(rna_sig_file, header = T, stringsAsFactors = F)

# Restrict the differential gene expression data to those genes that displayed
# statistical significance.
rna_sig <- rna_sig[rna_sig$significant, ]

# Identify the differential expression tendency.
rna_sig$Signature_category <- sapply(rna_sig$log2FoldChange, FUN = function(x) {
  paste0(ifelse(x > 0, yes = "Blood_naive_Treg", no = "Skin_Treg"), 
         "__hyperexpression")
})

# Read in the TPM values and restrict them to the relevant cell types and genes.
tpm_file <- paste0(location, "/RNASeq/tpm/tpm_all.txt")
tpm <- read.table(tpm_file, header = T, stringsAsFactors = F)
relevant_samples <- unlist(lapply(cell_types_rna_spelling, FUN = function(x) {
  rna_sample_mapping$Sample[rna_sample_mapping$Cell_type_wo_blank == x]
}))
tpm_relevant <- tpm[
  tpm$Gene_symbol %in% rownames(rna_sig), c("Gene_symbol", relevant_samples)
]
rownames(tpm_relevant) <- tpm_relevant$Gene_symbol
tpm_relevant <- tpm_relevant[, -1] # Remove gene symbol column.

# Transform the TPM values by log(100 * x + 1).
tpm_relevant <- log1p(100 * tpm_relevant)

# Re-order columns in the heat map to make the order of blood naive Treg
# cell samples match with other heat maps shown in the manuscript.
tpm_relevant <- tpm_relevant[, c(1, 2, 3, 4, 6, 8, 9, 5, 7)]

# Specify the colour code for the cell types.
cell_type_col_palette <- c("blue", "darkorchid1")
names(cell_type_col_palette) <- gsub("_", " ", cell_types_my_spelling)






###############################################################################
# Generate a heat map showing TPM values in all signature classes
# with at least 10 genes.
###############################################################################

# Quantify how many genes are present for each differential tendency.
sig_class_table <- table(rna_sig$Signature_category)

# Bring the differential tendencies into a custom order.
custom_order <- c("Skin_Treg__hyperexpression", 
                  "Blood_naive_Treg__hyperexpression")
rna_sig_ordered <- do.call(rbind, lapply(custom_order, FUN = function(x) {
  rna_sig[rna_sig$Signature_category == x, ]
}))

# Generate a colour palette for the differential tendencies.
set.seed(2023)
sig_classes_col_palette <- sample(rainbow(length(custom_order)), 
                                  length(custom_order))
names(sig_classes_col_palette) <- custom_order


# Generate a column annotation.
cell_types_by_sample <- sapply(colnames(tpm_relevant), FUN = function(x) {
  cell_types_my_spelling[
    cell_types_rna_spelling == rna_sample_mapping$Cell_type_wo_blank[
      rna_sample_mapping$Sample == x
    ]
    ]
})
col_anno_df <- data.frame(
  Cell_type = factor(gsub("_", " ", cell_types_by_sample), 
                     levels = gsub("_", " ", cell_types_my_spelling)),
  row.names = colnames(tpm_relevant)
)
col_anno_obj <- HeatmapAnnotation(
  df = col_anno_df, 
  col = list(Cell_type = cell_type_col_palette),
  annotation_label = "Cell type",
  annotation_legend_param = list(title = "Cell type", 
                                 title_gp = gpar(fontsize = 17, 
                                                 fontface = "bold"), 
                                 labels_gp = gpar(fontsize = 15), 
                                 grid_height = unit(8, "mm"),
                                 grid_width = unit(8, "mm"))
)

# Perform row-scaling of the TPM matrix.
tpm_matr_rowscaled <- t(scale(t(as.matrix(tpm_relevant))))

# Iterate over all signatures to generate accessibility heat maps.
single_heat_maps <- lapply(custom_order, FUN = function(x) {
  row_anno_df <- data.frame(
    Signature = rep(x, length(which(rna_sig_ordered$Signature_category == x)))
  )
  row_anno_df$Signature <- factor(row_anno_df$Signature, levels = custom_order)
  row_anno_obj <- rowAnnotation(
    df = row_anno_df, 
    col = list(Signature = sig_classes_col_palette),
    show_annotation_name = ifelse(x == custom_order[length(custom_order)], 
                                  yes = T, no = F),
    annotation_legend_param = list(title_gp = gpar(fontsize = 17, 
                                                   fontface = "bold"), 
                                   labels_gp = gpar(fontsize = 15),
                                   grid_height = unit(8, "mm"),
                                   grid_width = unit(8, "mm"))
    
  )
  maximum_gene_num <- max(sig_class_table)
  gene_num <- length(which(rna_sig_ordered$Signature_category == x))
  gene_num_anno_obj <- rowAnnotation(
    Placeholder = anno_block(labels = paste0(gene_num, "\nGenes"), 
                             show_name = F,
                             labels_rot = 0,
                             labels_gp = gpar(fontsize = 9)),
    show_annotation_name = x == custom_order[length(custom_order)],
    annotation_name_gp = gpar(fontsize = 12))
  matr_to_show <- tpm_matr_rowscaled[rownames(rna_sig_ordered)[
    rna_sig_ordered$Signature_category == x], 
  ]
  this_signature_heatmap <- Heatmap(
    matr_to_show,
    col = colorRamp2(seq(min(tpm_matr_rowscaled), 
                         max(tpm_matr_rowscaled), 
                         length.out = 200), 
                     mako(200)),
    cluster_rows = F,
    cluster_columns = F,
    show_row_names = F,
    show_column_names = F,
    heatmap_legend_param = list(
      title = "Gene\nexpression\n[log(100*TPM+1),\nrow-scaled]", 
      labels_gp = gpar(fontsize = 15), 
      title_gp = gpar(fontsize = 17, fontface = "bold")),
    left_annotation = row_anno_obj,
    right_annotation = gene_num_anno_obj,
    show_heatmap_legend = ifelse(x == custom_order[1], yes = T, no = F),
    height = 1
  )
  return(this_signature_heatmap)
})


# Merge the single heat maps into one and save this merged heat map.
indices <- 1:length(custom_order)
heatmap_drawing_string <- paste0(
  "draw(col_anno_obj %v% ",
  paste0("single_heat_maps[[", indices, "]]", collapse = " %v% "),
  ", merge_legends = T)"
)
signature_heatmap_outfile_pdf <- paste0(
  plot_outdir, 
  "/diff_expr_skin_treg_blood_naive_treg_heatmap.pdf"
)
signature_heatmap_outfile_rds <- paste0(
  plot_rds_outdir, 
  "/diff_expr_skin_treg_blood_naive_treg_heatmap.rds"
)
pdf(signature_heatmap_outfile_pdf, width = 9, height = 10)
eval(parse(text = heatmap_drawing_string))
dev.off()
saveRDS(single_heat_maps, file = signature_heatmap_outfile_rds)


