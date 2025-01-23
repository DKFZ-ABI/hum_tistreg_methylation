# This script generates a heat map showing differentially methylated regions
#   between skin Treg cells and blood naive Treg cells.
# Author: Niklas Beumer



# Load required packages.
library(GenomicRanges)
library(bsseq)
library(ComplexHeatmap)
library(circlize)
library(viridis)


# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

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
relevant_cell_types <- c("Skin Treg", "Blood naive Treg")

# Read in the list of differentially methylated regions between skin Tregs and 
# blood naive Tregs.
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
# meth_reg <- meth_reg[1:5000, ] # For debugging purposes.

# Generate a GRanges object for the differentially methylated regions.
meth_reg_gr <- makeGRangesFromDataFrame(meth_reg, keep.extra.columns = T)

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











###########################################################################
# Generate the heat map.
###########################################################################

# Specify a colour palette for the cell types.
cell_type_col_palette <- c("blue", "darkorchid1")
names(cell_type_col_palette) <- relevant_cell_types

# Generate a colour palette for the automatic annotation categories.
categories <- unique(meth_reg$automatic_annotation)
categories_col_palette <- rainbow(length(categories))
names(categories_col_palette) <- categories

# Generate a column annotation.
cell_types_by_sample <- pData(meth_data)$Cell_type
col_anno_df <- data.frame(
  Cell_type = factor(cell_types_by_sample, levels = relevant_cell_types),
  row.names = rownames(pData(meth_data))
)
col_anno_obj <- HeatmapAnnotation(
  df = col_anno_df,
  col = list(Cell_type = cell_type_col_palette),
  annotation_label = "Cell type",
  annotation_legend_param = list(
    title = "Cell type",
    title_gp = gpar(fontsize = 17, fontface = "bold"),
    labels_gp = gpar(fontsize = 15),
    grid_height = unit(8, "mm"),
    grid_width = unit(8, "mm")
  )
)

# Iterate over two automatic annotation categories.
single_heat_maps <- lapply(
  categories,
  FUN = function(x) {
    #-- Get a GRanges object for this signature.
    this_signature_gr <- meth_reg_gr[meth_reg_gr$automatic_annotation == x]
    #-- Extract average raw methylation values for this signature.
    avg_raw_meth_this_signature <- getMeth(meth_data,
                                           regions = this_signature_gr,
                                           what = "perRegion",
                                           type = "raw")
    #-- Generate the heat map for this signature.
    row_anno_df <- data.frame(Signature = rep(x, length(this_signature_gr)))
    row_anno_df$Signature <- factor(row_anno_df$Signature, levels = categories)
    row_anno_obj <- rowAnnotation(
      df = row_anno_df,
      col = list(Signature = categories_col_palette),
      show_annotation_name = ifelse(x == categories[length(categories)], yes = T, no = F),
      annotation_legend_param = list(
        title_gp = gpar(fontsize = 17, fontface = "bold"),
        labels_gp = gpar(fontsize = 15),
        grid_height = unit(8, "mm"),
        grid_width = unit(8, "mm")
      )
    )
    #-- Generate a text-based row annotation depicting the number of regions in this signature.
    region_num <- length(this_signature_gr)
    region_num_anno_obj <- rowAnnotation(
      Placeholder = anno_block(
        labels = paste0(region_num, "\nRegions"),
        show_name = F,
        labels_rot = 0,
        labels_gp = gpar(fontsize = 9)
      ),
      show_annotation_name = x == categories[length(categories)],
      annotation_name_gp = gpar(fontsize = 12)
    )
    this_signature_heatmap <- Heatmap(
      avg_raw_meth_this_signature,
      col = colorRamp2(seq(1, 0, length.out = 200), viridis(200)),
      cluster_rows = F,
      cluster_columns = F,
      show_row_names = F,
      show_column_names = F,
      heatmap_legend_param = list(
        title = "Methylation",
        labels_gp = gpar(fontsize = 15),
        title_gp = gpar(fontsize = 17, fontface = "bold")
      ),
      left_annotation = row_anno_obj,
      right_annotation = region_num_anno_obj,
      show_heatmap_legend = ifelse(x == categories[1], yes = T, no = F),
      height = 1
    )
    return(list(this_signature_heatmap, row_anno_df))
  }
)

# Merge the single heat maps into one and save this merged heat map.
indices <- 1:length(categories)
heatmap_drawing_string <- paste0(
  "draw(col_anno_obj %v% ",
  paste0("single_heat_maps[[", indices, "]][[1]]", collapse = " %v% "),
  ", merge_legends = T)"
)
signature_heatmap_outfile_pdf <- paste0(
  plot_outdir,
  "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_heatmap.pdf"
)
signature_heatmap_outfile_rds <- paste0(
  plot_rds_outdir,
  "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_heatmap.rds"
)
signature_heatmap_rowanno_outfile_rds <- paste0(
  plot_rds_outdir,
  "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_heatmap_rowanno_dfs.rds"
)
signature_heatmap_colanno_outfile_rds <- paste0(
  plot_rds_outdir,
  "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_heatmap_colanno_df.rds"
)
pdf(signature_heatmap_outfile_pdf,
    width = 9,
    height = 10)
eval(parse(text = heatmap_drawing_string))
dev.off()
saveRDS(single_heat_maps, file = signature_heatmap_outfile_rds)
row_anno_dfs <- lapply(
  single_heat_maps,
  FUN = function(x) {
    x[[2]]
  }
)
names(row_anno_dfs) <- categories
saveRDS(row_anno_dfs, file = signature_heatmap_rowanno_outfile_rds)
saveRDS(col_anno_df, file = signature_heatmap_colanno_outfile_rds)
