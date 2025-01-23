# This script generates a combined heat map of methylation, accessibility and 
#   expression for DMR-peak-gene links in the comparison between skin
#   Treg cells and blood naive Treg cells, stratified by the direction in which 
#   methylation and accessibility change.
# Author: Niklas Beumer



# Load required packages.
library(Seurat)
library(Signac)
library(bsseq)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(parallel)
library(testit)


# Define a location on /yyy.
b330_space <- "/yyy/"
location <- paste0(b330_space, "yyy/hm_treg_bs_rgnsbg")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Specify the naming conventions for the relevant cell types.
relevant_cell_types <- c("Skin Treg", "Blood naive Treg")
relevant_cell_types_atac_naming <- c("skin_treg", "blood_naive_treg")

# Read in the DMR-peak-gene links.
dmr_peak_links_file <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_diff_acc_diff_gex_dmr_peak_gene_links_meth_acc_gex_corr.txt"
)
dmr_peak_links <- read.table(dmr_peak_links_file, header = T,
                             stringsAsFactors = F, sep = "\t")

# Read in the RNA sample mapping file.
rna_sample_mapping_path <- paste0(location, 
                                  "/sample_mapping_rnaseq_only_biol_rep.txt")
rna_sample_mapping <- read.table(rna_sample_mapping_path, header = T, 
                                 stringsAsFactors = F, sep = "\t")

# Read in the TPM values.
tpm_file <- paste0(location, "/RNASeq/tpm/tpm_all.txt")
tpm <- read.table(tpm_file, header = T, stringsAsFactors = F)

# Read in the BSseq object containing the smoothed methylation data.
# Restrict the data to the relevant cell types.
meth_data_file <- paste0(
  location, 
  "/preprocessing_etc/quality_control_after_alignment/2022-03-14_bsseq_object_combined_all_samples_smoothed.rds"
)
meth_data <- readRDS(meth_data_file)
meth_data <- meth_data[, pData(meth_data)$Cell_type %in% relevant_cell_types]

# Read in the Seurat object containing scATAC-seq data of CD4+ cells.
sc_data_cd4_file <- paste0(
  b330_space, 
  "yyy/scATAC/final/seurat_human_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar.RDS"
)
sc_data_cd4 <- readRDS(sc_data_cd4_file)




##############################################################
# Process methylation data.
##############################################################

# Specify the colour palette for the cell types.
bsseq_metadata <- pData(meth_data)
cell_types_all <- unique(bsseq_metadata$Cell_type)
cell_type_col_palette <- c("blue", "darkorchid1")
names(cell_type_col_palette) <- cell_types_all

# Append this colour_palette to the meta data of the BSseq object.
bsseq_metadata$col <- cell_type_col_palette[bsseq_metadata$Cell_type]
pData(meth_data) <- bsseq_metadata






###############################################################
# Extract relevant data.
###############################################################

# Define the four categories of DMR-peak-gene links. They are defined by the 
# four quadrants of the correpation plot between differential methylation and 
# differential accessibility.
tendencies_meth <- paste0(c("Skin_Treg", "Blood_naive_Treg"),
                          "__hypomethylation")
tendencies_atac <- paste0(c("Skin_Treg", "Blood_naive_Treg"),
                          "__hyperaccessibility")

# Within the tendencies, order the DMR-peak-gene links alphabetically
# by gene name.
links_ordered <- do.call(c, lapply(tendencies_meth, FUN = function(x) {
  lapply(tendencies_atac, FUN = function(y) {
    temp_df <- 
      dmr_peak_links[dmr_peak_links$Methylation__automatic_annotation == x &
                       dmr_peak_links$Accessibility__signature_class == y, ]
    temp_df_ordered <- temp_df[order(temp_df$Expression__gene_name), ]
    return(temp_df_ordered)
  })
}))
table_order <- unlist(do.call(c, lapply(tendencies_meth, FUN = function(x) {
  lapply(tendencies_atac, FUN = function(y) {
    paste0(x, "_", y)
  })
})))




########### Methylation level.

# Generate a GRanges object for the differentially methylated regions.
meth_reg_gr <- lapply(links_ordered, FUN = function(x) {
  makeGRangesFromDataFrame(
    x[, 1:3], 
    seqnames.field = "Methylation__seqnames",
    start.field = "Methylation__start",
    end.field = "Methylation__end")
})

# Extract average raw methylation values for respective regions.
avg_raw_meth <- lapply(meth_reg_gr, FUN = function(x) {
  temp <- getMeth(meth_data, regions = x, what = "perRegion", type = "raw")
  rownames(temp) <- 1:nrow(temp)
  return(temp)
})



########### Accessibility level.

# Restrict scATAC-seq data to the relevant cell types.
sc_data_cd4_restr <- subset(
  sc_data_cd4, 
  subset = treg_tconv_annot %in% c(relevant_cell_types_atac_naming)
)

# Re-run TF-IDF.
sc_data_cd4_restr <- RunTFIDF(sc_data_cd4_restr)

# For each relevant peak and each cell type, compute mean normalised 
# accessibility scores.
mean_norm_acc <- lapply(links_ordered, FUN = function(z) {
  temp <- sapply(relevant_cell_types_atac_naming, FUN = function(x) {
    sc_data_cd4_this_cell_type <- subset(sc_data_cd4_restr, 
                                         subset = treg_tconv_annot == x)
    norm_acc_data <- GetAssayData(sc_data_cd4_this_cell_type, 
                                  assay = "scATAC_raw", slot = "data")
    real_peak_ids <- sapply(z$Accessibility__peak_ID, FUN = function(y) {
      strsplit(y, split = ".", fixed = T)[[1]][1]
    })
    unlist(mclapply(real_peak_ids, FUN = function(y) {
      mean(norm_acc_data[y, ])
    }, mc.cores = 10))
  })
  rownames(temp) <- 1:nrow(temp)
  return(temp)
})





########### Expression level.

# Restrict TPM values them to the relevant cell types and genes.
# Transform the values by log(100 * TPM + 1) and scale them row-wise.
relevant_samples <- unlist(lapply(relevant_cell_types, FUN = function(x) {
  rna_sample_mapping$Sample[rna_sample_mapping$Cell_type == x]
}))
tpm_relevant <- lapply(links_ordered, FUN = function(x) {
  temp_df <- do.call(rbind, lapply(x$Expression__gene_name, FUN = function(y) {
    tpm[tpm$Gene_symbol == y, relevant_samples]
  }))
  temp_df <- log1p(100 * temp_df)
  rownames(temp_df) <- 1:nrow(temp_df)
  return(t(scale(t(as.matrix(temp_df)))))
})







############################################################################
# Generate the common heat map including all DMR-peak-gene links.
#############################################################################

for (i in 1:length(links_ordered)) {
  
  # Perform some sanity checks.
  assert(nrow(avg_raw_meth[[i]]) == nrow(mean_norm_acc[[i]]))
  assert(nrow(mean_norm_acc[[i]]) == nrow(tpm_relevant[[i]]))
  
  # Generate the colour functions for the heatmaps.
  meth_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                             colors = viridis(200, direction = -1))
  acc_col_fun <- colorRamp2(breaks = seq(min(mean_norm_acc[[i]]), 
                                         max(mean_norm_acc[[i]]), 
                                         length.out = 200),
                            colors = viridis(200, option = "B"))
  gex_col_fun <- colorRamp2(breaks = seq(min(tpm_relevant[[i]]), 
                                         max(tpm_relevant[[i]]), 
                                         length.out = 200),
                            colors = viridis(200, option = "G"))
  
  # Generate column annotations for the heat maps.
  meth_cell_types_by_sample <- pData(meth_data)$Cell_type
  meth_col_anno_df <- data.frame(
    Cell_type = factor(meth_cell_types_by_sample, levels = relevant_cell_types),
    row.names = rownames(pData(meth_data)))
  meth_col_anno_obj <- HeatmapAnnotation(
    df = meth_col_anno_df, 
    col = list(Cell_type = cell_type_col_palette),
    annotation_label = "Cell type",
    annotation_legend_param = list(title = "Cell type"),
    show_annotation_name = F)
  acc_col_anno_df <- data.frame(
    Cell_type = relevant_cell_types,
    row.names = relevant_cell_types_atac_naming)
  acc_col_anno_obj <- HeatmapAnnotation(
    df = acc_col_anno_df, 
    col = list(Cell_type = cell_type_col_palette),
    annotation_label = "Cell type",
    annotation_legend_param = list(title = "Cell type"),
    show_annotation_name = F)
  rna_cell_types_by_sample <- sapply(
    colnames(tpm_relevant[[i]]), 
    FUN = function(x) {
      rna_sample_mapping$Cell_type[rna_sample_mapping$Sample == x]
    })
  rna_col_anno_df <- data.frame(
    Cell_type = factor(rna_cell_types_by_sample, 
                       levels = relevant_cell_types),
    row.names = colnames(tpm_relevant[[i]]))
  rna_col_anno_obj <- HeatmapAnnotation(
    df = rna_col_anno_df,
    col = list(Cell_type = cell_type_col_palette),
    annotation_label = "Cell type",
    annotation_legend_param = list(title = "Cell type"))
  
  # Generate single heat maps.
  meth_heatmap <- Heatmap(avg_raw_meth[[i]], 
                          cluster_rows = F,
                          cluster_columns = F,
                          col = meth_col_fun,
                          show_row_names = F,
                          show_column_names = F,
                          heatmap_legend_param = list(title = "Methylation"),
                          top_annotation = meth_col_anno_obj,
                          width = unit(3, "cm"))
  acc_heatmap <- Heatmap(mean_norm_acc[[i]], 
                         cluster_rows = F,
                         cluster_columns = F,
                         col = acc_col_fun,
                         show_row_names = F,
                         show_column_names = F,
                         heatmap_legend_param = list(
                           title = "Mean\nchr.\naccessibility"),
                         top_annotation = acc_col_anno_obj,
                         width = unit(3, "cm"))
  gex_heatmap <- Heatmap(tpm_relevant[[i]], 
                         cluster_rows = F,
                         cluster_columns = F,
                         col = gex_col_fun,
                         show_row_names = T,
                         show_column_names = F,
                         heatmap_legend_param = list(
                           title = "Gene\nexpression\n[log(100*TPM+1),\nrow-scaled]"),
                         width = unit(3, "cm"),
                         top_annotation = rna_col_anno_obj,
                         row_labels = links_ordered[[i]]$Expression__gene_name,
                         row_names_gp = gpar(fontsize = 6))
  
  
  
  # Combine the three heat maps into a list.
  # Note: By default, this orders the rows of the second/third heat map
  # in the same manner as in the first heat map
  # (https://jokergoo.github.io/ComplexHeatmap-reference/book/a-list-of-heatmaps.html;
  # 09 May 2023)
  combined_heatmaps <- meth_heatmap + acc_heatmap + gex_heatmap
  
  # Save the combined heat map.
  combined_heatmaps_pdf <- paste0(
    plot_outdir, 
    "/diff_meth_diff_acc_diff_gex_skin_naive_dmr_peak_gene_links_heatmap_", 
    table_order[i], 
    ".pdf"
  )
  combined_heatmaps_rds <- paste0(
    plot_rds_outdir, 
    "/diff_meth_diff_acc_diff_gex_skin_naive_dmr_peak_gene_links_heatmap_", 
    table_order[i], 
    ".rds"
  )
  pdf(combined_heatmaps_pdf, width = 10, height = 10)
  draw(combined_heatmaps)
  dev.off() 
  saveRDS(combined_heatmaps, file = combined_heatmaps_rds)
  
}
