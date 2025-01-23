# This script generates a heat map showing methylation (skin Treg cells and
#   blood naive Treg cells around predicted binding sites for selected bZIP 
#   transcription factors.
# Author: Niklas Beumer



# Load required packages.
library(Seurat)
library(Signac)
library(bsseq)
library(ComplexHeatmap)
library(circlize)
library(viridis)


# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the Seurat object for scATAC-seq data of CD4+ T cells, which also
# contains predicted binding sites for the selected bZIP factors.
sc_data_cd4_file <- paste0(
  location, 
  "/treg_hierarchies/seurat_obj_cd4_skin_treg_blood_naive_treg_updated_w_motifs_and_footprints_and_chromvar-f_homer_relevant_TFs.rds"
)
sc_data_cd4 <- readRDS(sc_data_cd4_file)

# Read in the BSseq object containing methylation data.
meth_data_file <- paste0(
  location, 
  "/preprocessing_etc/quality_control_after_alignment/2022-03-14_bsseq_object_combined_all_samples_smoothed.rds"
)
meth_data <- readRDS(meth_data_file)

# Specify the relevant transcription factor motifs. These are the bZIP motifs
# that are already shown for the differential analysis between skin Tregs and
# blood naive Tregs.
motifs_to_show <- c("Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer",
                    "Fosl2(bZIP)/3T3L1-Fosl2-ChIP-Seq(GSE56872)/Homer",
                    "Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer",
                    "JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer",
                    "Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer",
                    "Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer",
                    "BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer",
                    "Bach2(bZIP)/OCILy7-Bach2-ChIP-Seq(GSE44420)/Homer",
                    "AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer",
                    "Bach1(bZIP)/K562-Bach1-ChIP-Seq(GSE31477)/Homer")





#########################################################################
# Process methylation data.
#########################################################################

# Specify the colour palette for the cell types.
bsseq_metadata <- pData(meth_data)
cell_types_all <- unique(bsseq_metadata$Cell_type)
cell_type_col_palette <- c("blue", "cyan2", "orange", "darkorchid1", "red")
names(cell_type_col_palette) <- cell_types_all

# Append this colour_palette to the meta data of the BSseq object.
bsseq_metadata$col <- cell_type_col_palette[bsseq_metadata$Cell_type]
pData(meth_data) <- bsseq_metadata

# Restrict the data to those samples that contain skin Tregs and blood naive 
# Tregs.
cell_types <- c("Skin Treg", "Blood naive Treg")
cell_types_by_sample <- pData(meth_data)$Cell_type
samples_to_keep <- which(cell_types_by_sample %in% cell_types)
meth_data_relevant <- meth_data[, samples_to_keep]






#############################################################################
# Compute average raw methylation (sample-level) in windows around predicted
# transcription factor binding sites.
#############################################################################

# Specify how many base pairs each motif will be extended to both sides.
motif_extension <- 200

# Extract positions of predicted transcription factor binding sites.
motifs_obj <- Motifs(sc_data_cd4)
motifs_pos <- motifs_obj@positions[motifs_to_show]

# Harmonise seqlevels style with what is used in the BSSeq object.
motifs_pos_ncbi <- lapply(motifs_pos, FUN = function(x) {
  seqlevelsStyle(x) <- "NCBI"
  return(x)
})

# Extend the motif positions by the specified amount of bases.
motifs_pos_extended <- lapply(motifs_pos_ncbi, FUN = function(x) {
  x_ext <- x + motif_extension
  return(x_ext)
})

# Extract average raw methylation in extended regions around predicted 
# motif sites.
avg_raw_meth <- lapply(motifs_pos_extended, FUN = function(x) {
  meth_matr <- getMeth(meth_data_relevant, regions = x, type = "raw", 
                       what = "perRegion")
  return(meth_matr)
})

# Save these average raw methylation values together with the locations.
avg_raw_meth_outfile <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_skin_vs_naive_pred_rel_bzip_sites_ext_",
  motif_extension, 
  "_bp_locs_and_avg_raw_meth.rds")
list_to_save <- list(motif_locs_ext = motifs_pos_extended,
                     avg_raw_meth = avg_raw_meth)
saveRDS(list_to_save, file = avg_raw_meth_outfile)


# Exclude all regions that only return NA values.
avg_raw_meth_naomit <- lapply(avg_raw_meth, FUN = function(x) {
  lines_to_exclude <- which(apply(x, 1, FUN = function(y) {
  all(is.na(y))
  }))
  x_naomit <- x[-lines_to_exclude, ]
  return(x_naomit)
})







###########################################################################
# Generate a heat map showing these average methylation values.
###########################################################################

# Extract the colour palette for the cell types.
col_palette_relevant <- cell_type_col_palette[cell_types]

# Generate a colour palette for the transcription factor motifs.
motifs_col_palette <- rainbow(length(motifs_to_show))
names(motifs_col_palette) <- motifs_to_show

# Generate a column annotation.
cell_types_by_sample <- pData(meth_data_relevant)$Cell_type
col_anno_df <- data.frame(Cell_type = factor(cell_types_by_sample, 
                                             levels = cell_types),
                          row.names = rownames(pData(meth_data_relevant)))
col_anno_obj <- HeatmapAnnotation(
  df = col_anno_df, 
  col = list(Cell_type = col_palette_relevant),
  annotation_label = "Cell type",
  annotation_legend_param = list(title = "Cell type", 
                                 title_gp = gpar(fontsize = 17, 
                                                 fontface = "bold"), 
                                 labels_gp = gpar(fontsize = 15), 
                                 grid_height = unit(8, "mm"),
                                 grid_width = unit(8, "mm")))

# Iterate over the motifs.
single_heat_maps <- lapply(motifs_to_show, FUN = function(x) {
  #-- Generate the heat map for this motif. Regions are ordered by their mean
  #-- raw methylation across the samples.
  row_anno_df <- data.frame(Motif = rep(x, nrow(avg_raw_meth_naomit[[x]])))
  row_anno_df$Motif <- factor(row_anno_df$Motif, levels = motifs_to_show)
  row_anno_obj <- rowAnnotation(
    df = row_anno_df, 
    col = list(Motif = motifs_col_palette),
    show_annotation_name = ifelse(x == motifs_to_show[length(motifs_to_show)], 
                                  yes = T, no = F),
    annotation_legend_param = list(title_gp = gpar(fontsize = 17, 
                                                   fontface = "bold"), 
                                   labels_gp = gpar(fontsize = 15),
                                   grid_height = unit(8, "mm"),
                                   grid_width = unit(8, "mm")))
  ordering <- order(rowMeans(avg_raw_meth_naomit[[x]]))
  this_motif_heatmap <- Heatmap(avg_raw_meth_naomit[[x]][ordering, ],
                                col = colorRamp2(seq(1, 0, length.out = 200), 
                                                 viridis(200)),
                                cluster_rows = F,
                                cluster_columns = F,
                                show_row_names = F,
                                show_column_names = F,
                                heatmap_legend_param = list(
                                  title = "Methylation", 
                                  labels_gp = gpar(fontsize = 15), 
                                  title_gp = gpar(fontsize = 17, 
                                                  fontface = "bold")),
                                left_annotation = row_anno_obj,
                                show_heatmap_legend = ifelse(
                                  x == motifs_to_show[1], yes = T, no = F
                                ),
                                height = 1)
  return(list(this_motif_heatmap, row_anno_df))
})

# Merge the single heat maps into one and save this merged heat map.
indices <- 1:length(motifs_to_show)
heatmap_drawing_string <- paste0(
  "draw(col_anno_obj %v% ",
  paste0("single_heat_maps[[", indices, "]][[1]]", collapse = " %v% "),
  ", merge_legends = T)")
heatmap_outfile_pdf <- paste0(
  plot_outdir, 
  "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_pred_rel_bzip_sites_heatmap.pdf")
heatmap_outfile_rds <- paste0(
  plot_rds_outdir, 
  "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_pred_rel_bzip_sites_heatmap.rds")
rowanno_outfile_rds <- paste0(
  plot_rds_outdir, 
  "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_pred_rel_bzip_sites_heatmap_rowanno_dfs.rds")
colanno_outfile_rds <- paste0(
  plot_rds_outdir, 
  "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_pred_rel_bzip_sites_heatmap_colanno_df.rds")
pdf(heatmap_outfile_pdf, width = 11, height = 10)
eval(parse(text = heatmap_drawing_string))
dev.off()
saveRDS(single_heat_maps, file = heatmap_outfile_rds)
row_anno_dfs <- lapply(single_heat_maps, FUN = function(x) {x[[2]]})
names(row_anno_dfs) <- motifs_to_show
saveRDS(row_anno_dfs, file = rowanno_outfile_rds)
saveRDS(col_anno_df, file = colanno_outfile_rds)








###########################################################################
# Perform statistical tests assessing differences between skin Tregs and
# blood naive Tregs. To this end, average methylation signal in every
# region across biological replicates of each cell type and then perform
# a Wilcoxon signed-rank test.
###########################################################################

# Iterate over the motifs.
all_test_res <- lapply(motifs_to_show, FUN = function(x) {
  #-- Extract methylation values for this motif.
  rel_meth_vals <- avg_raw_meth_naomit[[x]]
  #-- Average the methylation values across the biological replicates of each 
  #-- cell type.
  rel_meth_vals_aggr <- sapply(cell_types, FUN = function(y) {
    rowMeans(rel_meth_vals[, grep(y, colnames(rel_meth_vals))],
             na.rm = T)
  })
  #-- Perform a two-tailed Wilcoxon signed-rank test between the two cell types.
  test_res <- wilcox.test(rel_meth_vals_aggr[, 1], rel_meth_vals_aggr[, 2],
                          alternative = "two.sided", paired = T)
  #-- Return the results of the test.
  return(test_res)
})
names(all_test_res) <- motifs_to_show

# Save the detailed test results.
test_res_outfile <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_skin_vs_naive_pred_rel_bzip_sites_ext_",
  motif_extension, 
  "_wilcox_test_res.rds")
saveRDS(all_test_res, file = test_res_outfile)

# Save a summary of the P values.
test_summary_outfile <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_skin_vs_naive_pred_rel_bzip_sites_ext_",
  motif_extension, 
  "_wilcox_test_summary.txt"
)
sink(test_summary_outfile)
void <- lapply(motifs_to_show, FUN = function(x) {
  print(paste0(x, ": P = ", all_test_res[[x]]$p.value))
})
sink()
