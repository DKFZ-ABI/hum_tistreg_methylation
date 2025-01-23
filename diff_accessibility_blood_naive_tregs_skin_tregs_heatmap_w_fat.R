# This script generates an accessibility heat map for peaks that are
#   differentially accessible between skin Treg cells and blood naive Treg
#   cells, also visualising fat Treg cells.
# Author: Niklas Beumer
# Run this script with 10 cores!!!



# Load required packages.
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(Seurat)
library(Signac)
library(parallel)


# Specify a location on /xxx.
b330_space <- "/xxx/"
location <- paste0(b330_space, "nbeumer/hm_treg_bs_rgnsbg")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/yyy/hm_treg_bs_rgnsbg/analysis",
                     format(Sys.time(), "%m-%d-%y"),
                     sep = "/")
if (!dir.exists(plot_outdir)) {
  dir.create(plot_outdir)
}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the Seurat object containing scATAC-seq data of CD4+ cells.
sc_data_cd4_file <- paste0(
  b330_space,
  "msimon/scATAC/final/seurat_human_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar.RDS"
)
sc_data_cd4 <- readRDS(sc_data_cd4_file)

# Read in the table with peaks that are differentially accessible between skin
# Tregs and blood naive Tregs.
atac_sig_file <- paste0(
  location,
  "/differential_accessibility/diff_acc_analysis_Skin_TregBlood_naive_Treg_results_with_significance.txt"
)
atac_sig <- read.table(
  atac_sig_file,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)

# Restrict the peaks to those that displayed statistical significance.
atac_sig <- atac_sig[atac_sig$significant, ]

# Generate a column specifying the differential tendency.
atac_sig$Signature_category <- sapply(
  atac_sig$avg_log2FC,
  FUN = function(x) {
    ifelse(x < 0, yes = "Skin_Treg__hyperaccessibility", no = "Blood_naive_Treg__hyperaccessibility")
  }
)

# Bring the signature peaks into a custom order.
custom_order <- c("Skin_Treg__hyperaccessibility",
                  "Blood_naive_Treg__hyperaccessibility")
atac_sig_ordered <- do.call(rbind, lapply(
  custom_order,
  FUN = function(x) {
    atac_sig[atac_sig$Signature_category == x, ]
  }
))

# Generate a colour palette for the differential tendencies.
categories_col_palette <- rainbow(length(custom_order))
names(categories_col_palette) <- custom_order

# Specify the cell types.
cell_types_my_spelling <- c("Skin_Treg", "Fat_Treg", "Blood_naive_Treg")
cell_types_atac_spelling <- c("skin_treg", "fat_treg", "blood_naive_treg")
names(cell_types_my_spelling) <- cell_types_atac_spelling

# Specify the colour code for the cell types.
cell_type_col_palette <- c("blue", "green", "darkorchid1")
names(cell_type_col_palette) <- gsub("_", " ", cell_types_my_spelling)

# Specify the identity of donors that provided blood  and skin samples.
blood_donors <- c("1", "2")
skin_donors <- c("4", "5")
fat_donors <- c("3", "4", "5")




###############################################################################
# Generate a heat map showing accessibility scores for both diffferential
# tendencies (skin Treg hyperaccessibility and blood naive Treg
# hyperaccessibility.
###############################################################################

# Restrict scATAC-seq data to the relevant cell types.
sc_data_cd4_restr <- subset(sc_data_cd4, subset = treg_tconv_annot %in% c(cell_types_atac_spelling))

# Re-run TF-IDF.
sc_data_cd4_restr <- RunTFIDF(sc_data_cd4_restr)

# For each relevant peak and each cell type, compute mean normalised accessibility scores.
mean_norm_acc <- sapply(
  cell_types_atac_spelling,
  FUN = function(x) {
    sc_data_cd4_this_cell_type <- subset(sc_data_cd4_restr, subset = treg_tconv_annot == x)
    norm_acc_data <- GetAssayData(sc_data_cd4_this_cell_type,
                                  assay = "scATAC_raw",
                                  slot = "data")
    unlist(mclapply(
      rownames(atac_sig_ordered),
      FUN = function(y) {
        mean(norm_acc_data[y, ])
      },
      mc.cores = 10
    ))
  }
)
rownames(mean_norm_acc) <- rownames(atac_sig_ordered)

# Generate a column annotation.
col_anno_df <- data.frame(Cell_type = factor(
  gsub("_", " ", cell_types_my_spelling),
  levels = gsub("_", " ", cell_types_my_spelling)
),
row.names = cell_types_atac_spelling)
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

# Iterate over all signatures to generate accessibility heat maps.
single_heat_maps <- lapply(
  custom_order,
  FUN = function(x) {
    row_anno_df <- data.frame(Signature = rep(x, length(
      which(atac_sig_ordered$Signature_category == x)
    )))
    row_anno_df$Signature <- factor(row_anno_df$Signature, levels = custom_order)
    row_anno_obj <- rowAnnotation(
      df = row_anno_df,
      col = list(Signature = categories_col_palette),
      show_annotation_name = ifelse(x == custom_order[length(custom_order)], yes = T, no = F),
      annotation_legend_param = list(
        title_gp = gpar(fontsize = 17, fontface = "bold"),
        labels_gp = gpar(fontsize = 15),
        grid_height = unit(8, "mm"),
        grid_width = unit(8, "mm")
      )
    )
    region_num <- length(which(atac_sig_ordered$Signature_category == x))
    region_num_anno_obj <- rowAnnotation(
      Placeholder = anno_block(
        labels = paste0(region_num, "\nPeaks"),
        show_name = F,
        labels_rot = 0,
        labels_gp = gpar(fontsize = 9)
      ),
      show_annotation_name = x == custom_order[length(custom_order)],
      annotation_name_gp = gpar(fontsize = 12)
    )
    matr_to_show <- mean_norm_acc[rownames(atac_sig_ordered)[atac_sig_ordered$Signature_category == x], ]
    this_signature_heatmap <- Heatmap(
      matr_to_show,
      col = colorRamp2(seq(
        min(mean_norm_acc), max(mean_norm_acc), length.out = 200
      ), inferno(200)),
      cluster_rows = F,
      cluster_columns = F,
      show_row_names = F,
      show_column_names = F,
      heatmap_legend_param = list(
        title = "Mean\nchr.\naccessibility",
        labels_gp = gpar(fontsize = 15),
        title_gp = gpar(fontsize = 17, fontface = "bold")
      ),
      left_annotation = row_anno_obj,
      right_annotation = region_num_anno_obj,
      show_heatmap_legend = ifelse(x == custom_order[1], yes = T, no = F),
      height = 1
    )
    return(this_signature_heatmap)
  }
)


# Merge the single heat maps into one and save this merged heat map.
indices <- 1:length(custom_order)
heatmap_drawing_string <- paste0(
  "draw(col_anno_obj %v% ",
  paste0("single_heat_maps[[", indices, "]]", collapse = " %v% "),
  ", merge_legends = T)"
)
signature_heatmap_outfile_pdf <- paste0(plot_outdir,
                                        "/diff_acc_skin_treg_blood_naive_treg_heatmap_w_fat_treg.pdf")
signature_heatmap_outfile_rds <- paste0(plot_rds_outdir,
                                        "/diff_acc_skin_treg_blood_naive_treg_heatmap_w_fat_treg.rds")
pdf(signature_heatmap_outfile_pdf,
    width = 9,
    height = 10)
eval(parse(text = heatmap_drawing_string))
dev.off()
saveRDS(single_heat_maps, file = signature_heatmap_outfile_rds)






##########################################################################
# Generate another heat map that shows aggregated accessibility data
# separately for the donors that donated the relevant tissue(s).
##########################################################################

# Restrict scATAC-seq data to the relevant cell types and to cells from donors
# that donated the relevant tissue(s).
sc_data_cd4_restr <- subset(
  sc_data_cd4,
  subset = (treg_tconv_annot == "blood_naive_treg" &
              donor %in% blood_donors) |
    (treg_tconv_annot == "skin_treg" & donor %in% skin_donors) |
    (treg_tconv_annot == "fat_treg" & donor %in% fat_donors)
)

# Re-run TF-IDF.
sc_data_cd4_restr <- RunTFIDF(sc_data_cd4_restr)

# For each relevant peak and each cell type, compute mean normalised
# accessibility scores separately for each donor.
mean_norm_acc <- do.call(cbind, lapply(
  cell_types_atac_spelling,
  FUN = function(x) {
    sc_data_cd4_this_cell_type <- subset(sc_data_cd4_restr, subset = treg_tconv_annot == x)
    this_cell_type_donors <- sort(unique(sc_data_cd4_this_cell_type$donor))
    norm_acc_data_this_celltype <- sapply(
      this_cell_type_donors,
      FUN = function(y) {
        sc_data_this_donor <- subset(sc_data_cd4_this_cell_type, subset = donor == y)
        norm_acc_data <- GetAssayData(sc_data_this_donor, assay = "scATAC_raw", slot = "data")
        print("Matrix extracted")
        data_aggr <- unlist(mclapply(
          rownames(atac_sig_ordered),
          FUN = function(z) {
            mean(norm_acc_data[z, ])
          },
          mc.cores = 10
        ))
        return(data_aggr)
      }
    )
    colnames(norm_acc_data_this_celltype) <-
      paste0(x, "_donor_", this_cell_type_donors)
    return(norm_acc_data_this_celltype)
  }
))
rownames(mean_norm_acc) <- rownames(atac_sig_ordered)

# Generate a column annotation.
col_anno_df <- data.frame(Cell_type = factor(
  sapply(
    colnames(mean_norm_acc),
    FUN = function(x) {
      cell_type <- strsplit(x, split = "_donor")[[1]][1]
      return(gsub("_", " ", cell_types_my_spelling[cell_type]))
    }
  ),
  levels = gsub("_", " ", cell_types_my_spelling)
),
row.names = colnames(mean_norm_acc))
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

# Iterate over all signatures to generate accessibility heat maps.
single_heat_maps <- lapply(
  custom_order,
  FUN = function(x) {
    row_anno_df <- data.frame(Signature = rep(x, length(
      which(atac_sig_ordered$Signature_category == x)
    )))
    row_anno_df$Signature <- factor(row_anno_df$Signature, levels = custom_order)
    row_anno_obj <- rowAnnotation(
      df = row_anno_df,
      col = list(Signature = categories_col_palette),
      show_annotation_name = ifelse(x == custom_order[length(custom_order)], yes = T, no = F),
      annotation_legend_param = list(
        title_gp = gpar(fontsize = 17, fontface = "bold"),
        labels_gp = gpar(fontsize = 15),
        grid_height = unit(8, "mm"),
        grid_width = unit(8, "mm")
      )
    )
    region_num <- length(which(atac_sig_ordered$Signature_category == x))
    region_num_anno_obj <- rowAnnotation(
      Placeholder = anno_block(
        labels = paste0(region_num, "\nPeaks"),
        show_name = F,
        labels_rot = 0,
        labels_gp = gpar(fontsize = 9)
      ),
      show_annotation_name = x == custom_order[length(custom_order)],
      annotation_name_gp = gpar(fontsize = 12)
    )
    matr_to_show <- mean_norm_acc[rownames(atac_sig_ordered)[atac_sig_ordered$Signature_category == x], ]
    this_signature_heatmap <- Heatmap(
      matr_to_show,
      col = colorRamp2(seq(
        min(mean_norm_acc), max(mean_norm_acc), length.out = 200
      ), inferno(200)),
      cluster_rows = F,
      cluster_columns = F,
      show_row_names = F,
      show_column_names = F,
      heatmap_legend_param = list(
        title = "Mean\nchr.\naccessibility",
        labels_gp = gpar(fontsize = 15),
        title_gp = gpar(fontsize = 17, fontface = "bold")
      ),
      left_annotation = row_anno_obj,
      right_annotation = region_num_anno_obj,
      show_heatmap_legend = ifelse(x == custom_order[1], yes = T, no = F),
      height = 1
    )
    return(this_signature_heatmap)
  }
)


# Merge the single heat maps into one and save this merged heat map.
indices <- 1:length(custom_order)
heatmap_drawing_string <- paste0(
  "draw(col_anno_obj %v% ",
  paste0("single_heat_maps[[", indices, "]]", collapse = " %v% "),
  ", merge_legends = T)"
)
signature_heatmap_outfile_pdf <- paste0(
  plot_outdir,
  "/diff_acc_skin_treg_blood_naive_treg_heatmap_w_fat_treg_by_donor.pdf"
)
signature_heatmap_outfile_rds <- paste0(
  plot_rds_outdir,
  "/diff_acc_skin_treg_blood_naive_treg_heatmap_w_fat_treg_by_donor.rds"
)
pdf(signature_heatmap_outfile_pdf,
    width = 9,
    height = 10)
eval(parse(text = heatmap_drawing_string))
dev.off()
saveRDS(single_heat_maps, file = signature_heatmap_outfile_rds)
