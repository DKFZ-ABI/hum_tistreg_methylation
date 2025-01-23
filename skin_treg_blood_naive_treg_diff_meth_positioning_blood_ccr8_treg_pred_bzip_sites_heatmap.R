# This script generates a heat map showing the positioning of blood CCR8+ Tregs
#   (methylation-wise) with respect to predicted binding sites for selected bZIP
#   transcription factors.
# Author: Niklas Beumer



# Load required packages.
library(Seurat)
library(Signac)
library(bsseq)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(GenomicRanges)
library(testit)

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

# Read in the Seurat object for scATAC-seq data of human CD4+ T cells with
# motif information for relevant bZIP motifs.
sc_data_cd4_file <- paste0(
  location,
  "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_seur_obj_w_vzip_footprints.rds"
)
sc_data_cd4 <- readRDS(sc_data_cd4_file)

# Read in the BSseq object containing methylation data.
meth_data_file <- paste0(
  location,
  "/preprocessing_etc/quality_control_after_alignment/2022-03-14_bsseq_object_combined_all_samples_smoothed.rds"
)
meth_data <- readRDS(meth_data_file)

# Read in the table of differentially methylated regions between skin Tregs
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

# Specify the names of the bZIP factors that should be visualised.
motifs_to_show <- c(
  "BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer",
  "Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer",
  "AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer",
  "Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer",
  "JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer",
  "Fosl2(bZIP)/3T3L1-Fosl2-ChIP-Seq(GSE56872)/Homer",
  "Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer",
  "Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer",
  "Bach2(bZIP)/OCILy7-Bach2-ChIP-Seq(GSE44420)/Homer",
  "HLF(bZIP)/HSC-HLF.Flag-ChIP-Seq(GSE69817)/Homer",
  "Atf1(bZIP)/K562-ATF1-ChIP-Seq(GSE31477)/Homer",
  "NFIL3(bZIP)/HepG2-NFIL3-ChIP-Seq(Encode)/Homer",
  "bZIP:IRF(bZIP,IRF)/Th17-BatF-ChIP-Seq(GSE39756)/Homer",
  "MafA(bZIP)/Islet-MafA-ChIP-Seq(GSE30298)/Homer"
)






#########################################################################
# Process methylation data.
#########################################################################

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



# Specify the colour palette for the cell types.
bsseq_metadata <- pData(meth_data)
cell_types_all <- unique(bsseq_metadata$Cell_type)
cell_type_col_palette <- c("blue",
                           "cyan2",
                           "orange",
                           "darkorchid1",
                           "red")
names(cell_type_col_palette) <- cell_types_all

# Append this colour_palette to the meta data of the BSseq object.
bsseq_metadata$col <- cell_type_col_palette[bsseq_metadata$Cell_type]
pData(meth_data) <- bsseq_metadata

# Restrict the data to those samples that contain skin Tregs, blood CCR8+ Tregs
# or blood naive Tregs.
cell_types <- c("Blood naive Treg", "Blood CCR8+ Treg", "Skin Treg")
cell_types_by_sample <- pData(meth_data)$Cell_type
samples_to_keep <- which(cell_types_by_sample %in% cell_types)
meth_data_relevant <- meth_data[, samples_to_keep]

# Specify how many base pairs each motif will be extended to both sides.
motif_extension <- 200

# Extract positions of predicted transcription factor binding sites.
motifs_obj <- Motifs(sc_data_cd4)
motifs_pos <- motifs_obj@positions[motifs_to_show]

# Harmonise seqlevels style with what is used in the BSSeq object.
motifs_pos_ncbi <- lapply(
  motifs_pos,
  FUN = function(x) {
    seqlevelsStyle(x) <- "NCBI"
    return(x)
  }
)

# Restrict to motifs that overlap with DMRs between skin Tregs and
# blood naive Tregs.
meth_reg_gr <- makeGRangesFromDataFrame(meth_reg)
motifs_pos_ncbi_overl <- lapply(
  motifs_pos_ncbi,
  FUN = function(x) {
    overl <- findOverlaps(x, meth_reg_gr)
    x_restr <- x[unique(from(overl)), ]
    return(x_restr)
  }
)

# Extend the motif positions by the specified amount of bases.
motifs_pos_extended <- lapply(
  motifs_pos_ncbi_overl,
  FUN = function(x) {
    x_ext <- x + motif_extension
    return(x_ext)
  }
)

# Identify average raw methylation values (on the sample level) in regions
# around predicted motif sites.
avg_meth_sample_level <- lapply(
  motifs_pos_extended,
  FUN = function(x) {
    getMeth(
      meth_data_relevant,
      regions = x,
      type = "raw",
      what = "perRegion"
    )
  }
)

# Compute cell-type-level averages of these sample-level methylation values.
avg_meth_celltype_level <- lapply(
  avg_meth_sample_level,
  FUN = function(y) {
    sapply(
      cell_types,
      FUN = function(x) {
        celltype_subs <- y[, grep(x, colnames(y), fixed = T)]
        return(rowMeans(celltype_subs, na.rm = T))
      }
    )
  }
)

# Scale the values between 0 and 1.
avg_meth_celltype_level_scaled <- lapply(avg_meth_celltype_level, 
                                         FUN = scale_matr_lines_betw_0_and_1)

# Quantify distances in scaled methylation between blood CCR8+ Tregs and
# the two extreme cell types.
diffs_dists <- lapply(
  avg_meth_celltype_level_scaled,
  FUN = function(x) {
    ccr8_naive_diffs <- x[, "Blood naive Treg"] - x[, "Blood CCR8+ Treg"]
    ccr8_naive_dists <- abs(ccr8_naive_diffs)
    ccr8_skin_diffs <- x[, "Skin Treg"] - x[, "Blood CCR8+ Treg"]
    ccr8_skin_dists <- abs(ccr8_skin_diffs)
    temp_matr <- matrix(c(
      ccr8_naive_diffs,
      ccr8_naive_dists,
      ccr8_skin_diffs,
      ccr8_skin_dists
    ),
    ncol = 4)
    return(temp_matr)
  }
)


# For each region, determine whether blood CCR8+ Tregs are closer to
# blood naive Tregs or closer to skin Tregs.
closest_extrema <- lapply(
  diffs_dists,
  FUN = function(y) {
    sapply(
      1:nrow(y),
      FUN = function(x) {
        if (is.na(y[x, 2]) | is.na(y[x, 4])) {
          return(NA)
        } else if (y[x, 2] < y[x, 4]) {
          return("Blood naive Treg")
        } else if (y[x, 2] > y[x, 4]) {
          return("Skin Treg")
        } else {
          return("Exactly in middle")
        }
      }
    )
  }
)

# Save these average raw methylation values together with the locations.
avg_raw_meth_outfile <- paste0(
  location,
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_pred_rel_bzip_sites_overl_w_dmrs_ext_",
  motif_extension,
  "_bp_locs_and_avg_raw_meth_and_other_info.rds"
)
list_to_save <- list(
  motif_locs_ext = motifs_pos_extended,
  avg_raw_meth_samples = avg_meth_sample_level,
  avg_raw_meth_celltypes = avg_meth_celltype_level,
  avg_raw_meth_celltypes_scaled = avg_meth_celltype_level_scaled,
  closest_extrema = closest_extrema
)
saveRDS(list_to_save, file = avg_raw_meth_outfile)


# Exclude all regions that return NA values regarding the closest extremum.
avg_meth_celltype_level_scaled_naomit <- lapply(
  1:length(avg_meth_celltype_level_scaled),
  FUN = function(x) {
    lines_to_exclude <- which(is.na(closest_extrema[[x]]))
    if (length(lines_to_exclude) > 0) {
      naomit <- avg_meth_celltype_level_scaled[[x]][-lines_to_exclude, ]
    } else {
      naomit <- avg_meth_celltype_level_scaled[[x]]
    }
    return(naomit)
  }
)
names(avg_meth_celltype_level_scaled_naomit) <-
  names(avg_meth_celltype_level_scaled)
closest_extrema_naomit <- lapply(
  closest_extrema,
  FUN = function(x) {
    lines_to_exclude <- which(is.na(x))
    if (length(lines_to_exclude) > 0) {
      naomit <- x[-lines_to_exclude]
    } else {
      naomit <- x
    }
    return(naomit)
  }
)
diffs_dists_naomit <- lapply(
  1:length(diffs_dists),
  FUN = function(x) {
    lines_to_exclude <- which(is.na(closest_extrema[[x]]))
    if (length(lines_to_exclude) > 0) {
      naomit <- diffs_dists[[x]][-lines_to_exclude, ]
    } else {
      naomit <- diffs_dists[[x]]
    }
  }
)
names(diffs_dists_naomit) <- names(diffs_dists)








############################################################################
# Plot a heat map showing scaled methylation
# around the predicted motif sites.
############################################################################

# Extract the colour palette for the cell types.
col_palette_relevant <- cell_type_col_palette[cell_types]

# Generate a colour palette for the transcription factor motifs.
motifs_col_palette <- rainbow(length(motifs_to_show))
names(motifs_col_palette) <- motifs_to_show

# Iterate over the motifs.
single_heat_maps <- lapply(
  motifs_to_show,
  FUN = function(x) {
    #-- Order the values according to their blood CCR8+ Treg positioning.
    naive_hypo_inds <- which(
      avg_meth_celltype_level_scaled_naomit[[x]][, "Skin Treg"] >
        avg_meth_celltype_level_scaled_naomit[[x]][, "Blood naive Treg"]
    )
    naive_hypo_order <- naive_hypo_inds[order(diffs_dists_naomit[[x]][naive_hypo_inds, 3])]
    equal_inds <- which(
      avg_meth_celltype_level_scaled_naomit[[x]][, "Skin Treg"] ==
        avg_meth_celltype_level_scaled_naomit[[x]][, "Blood naive Treg"]
    )
    skin_hypo_inds <- which(
      avg_meth_celltype_level_scaled_naomit[[x]][, "Skin Treg"] <
        avg_meth_celltype_level_scaled_naomit[[x]][, "Blood naive Treg"]
    )
    skin_hypo_order <- skin_hypo_inds[order(diffs_dists_naomit[[x]][skin_hypo_inds, 3])]
    order <- c(skin_hypo_order, equal_inds, naive_hypo_order)
    assert(all(
      sort(order) == 1:nrow(avg_meth_celltype_level_scaled_naomit[[x]])
    ))
    vals_ordered <- avg_meth_celltype_level_scaled_naomit[[x]][order, ]
    closest_extrema_ordered <- closest_extrema_naomit[[x]][order]
    #-- Generate the heat map for this motif.
    row_anno_motif_df <- data.frame(Motif = rep(x, nrow(
      avg_meth_celltype_level_scaled_naomit[[x]]
    )))
    row_anno_motif_df$Motif <- factor(row_anno_motif_df$Motif, levels = motifs_to_show)
    row_anno_motif_obj <- rowAnnotation(
      df = row_anno_motif_df,
      col = list(Motif = motifs_col_palette),
      show_annotation_name = ifelse(x == motifs_to_show[length(motifs_to_show)], yes = T, no = F),
      annotation_legend_param = list(
        title_gp = gpar(fontsize = 17, fontface = "bold"),
        labels_gp = gpar(fontsize = 15),
        grid_height = unit(8, "mm"),
        grid_width = unit(8, "mm")
      )
    )
    row_anno_closest_obj <- rowAnnotation(
      `Blood CCR8+ Treg position` = sapply(
        closest_extrema_ordered,
        FUN = function(x) {
          switch(x,
                 `Blood naive Treg` = "Closer to blood naive Tregs",
                 `Skin Treg` = "Closer to skin Tregs",
                 "Exactly in middle")
        }
      ),
      col = list(
        `Blood CCR8+ Treg position` = c(
          `Closer to skin Tregs` = "blue",
          `Closer to blood naive Tregs` = "darkorchid1",
          `Exactly in middle` = "grey"
        )
      ),
      show_annotation_name = ifelse(x == motifs_to_show[length(motifs_to_show)], yes = T, no = F),
      annotation_legend_param = list(
        title_gp = gpar(fontsize = 17, fontface = "bold"),
        labels_gp = gpar(fontsize = 15),
        grid_height = unit(8, "mm"),
        grid_width = unit(8, "mm")
      )
    )
    this_motif_heatmap <- Heatmap(
      vals_ordered,
      col = colorRamp2(seq(1, 0, length.out = 200), viridis(200)),
      cluster_rows = F,
      cluster_columns = F,
      show_row_names = F,
      show_column_names = T,
      heatmap_legend_param = list(
        title = "Methylation",
        labels_gp = gpar(fontsize = 15),
        title_gp = gpar(fontsize = 17, fontface = "bold"),
        at = c(0, 1),
        labels = c("min", "max")
      ),
      left_annotation = row_anno_motif_obj,
      right_annotation = row_anno_closest_obj,
      show_heatmap_legend = ifelse(x == motifs_to_show[1], yes = T, no = F),
      height = 1,
      column_names_max_height = unit(10, "cm")
    )
    return(list(
      this_motif_heatmap,
      row_anno_motif_df,
      row_anno_closest_obj
    ))
  }
)
# Merge the single heat maps into one and save this merged heat map.
indices <- 1:length(motifs_to_show)
heatmap_drawing_string <- paste0(
  "draw(",
  paste0("single_heat_maps[[", indices, "]][[1]]", collapse = " %v% "),
  ", merge_legends = T)"
)
heatmap_outfile_pdf <- paste0(
  plot_outdir,
  "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_pred_rel_bzip_sites_overl_w_dmrs_heatmap.pdf"
)
heatmap_outfile_rds <- paste0(
  plot_rds_outdir,
  "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_pred_rel_bzip_sites_overl_w_dmrs_heatmap.rds"
)
rowanno_outfile_rds <- paste0(
  plot_rds_outdir,
  "/treg_hierarchies_diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_pred_rel_bzip_sites_overl_w_dmrs_heatmap_rowanno_dfs.rds"
)
pdf(heatmap_outfile_pdf, width = 11, height = 10)
eval(parse(text = heatmap_drawing_string))
dev.off()
saveRDS(single_heat_maps, file = heatmap_outfile_rds)
row_anno_dfs <- lapply(
  single_heat_maps,
  FUN = function(x) {
    x[2:3]
  }
)
names(row_anno_dfs) <- motifs_to_show
saveRDS(row_anno_dfs, file = rowanno_outfile_rds)