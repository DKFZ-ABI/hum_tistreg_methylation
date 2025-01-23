# This script generates heat maps with raw methylation values for the track
#   plots showing DMR-peak-gene links from the comparison between blood naive
#   Treg cells and skin Treg cells.
# Author: Niklas Beumer



# Load required packages.
library(GenomicRanges)
library(bsseq)
library(testit)
library(Signac)
library(ComplexHeatmap)
library(circlize)
library(viridis)


# Define a location on /xxx.
b330_space <- "/xxx/"
location <- paste0(b330_space, "nbeumer/hm_treg_bs_rgnsbg")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/yyyy/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the list of DMR-peak-gene links.
links_file <<- paste0(
  location, 
  "/treg_hierarchies/diff_meth_diff_acc_diff_gex_dmr_peak_gene_links_meth_acc_gex_corr.txt"
)
links <- read.table(links_file, header = T, stringsAsFactors = F, sep = "\t")

# Specify the indices of links for whose track plots we need heat maps.
links_to_show_inds <- c(807, 992, 614, 11, 383)

# Generate a GRanges object containing the overlapping regions between 
# the DMR and the peak in each relevant DMR-peak-gene link.
meth_gr <- makeGRangesFromDataFrame(links[links_to_show_inds, c(1:3, 7)], 
                                    keep.extra.columns = T)
acc_gr <- StringToGRanges(links[links_to_show_inds, "Accessibility__peak_ID"],
                          starts.in.df.are.0based = T)
seqlevelsStyle(acc_gr) <- "NCBI"
overl_gr <- pintersect(meth_gr, acc_gr)
assert(all(width(overl_gr) >= 1))
assert(length(overl_gr) == length(meth_gr))

# Extend the regions in the same way that they were extended for track plot 
# generation.
overl_gr <- overl_gr + 3000

# Specify the relevant cell types.
relevant_cell_types <- c("Blood naive Treg", "Skin Treg")

# Read in the BSseq object containing the smoothed methylation data
# and restrict to the relevant cell types.
meth_data <- readRDS(paste(
  location, 
  "/preprocessing_etc/quality_control_after_alignment/2022-03-14_bsseq_object_combined_all_samples_smoothed.rds", 
  sep = "/"
))
meth_data <- meth_data[, pData(meth_data)$Cell_type %in% relevant_cell_types]



################################################################################
# Generate heat maps showing raw. methylatin for each CpG that was covered in 
# the track plots.
################################################################################

# Specify a colour palette for the cell types.
cell_type_col_palette <- c("darkorchid1", "blue")
names(cell_type_col_palette) <- relevant_cell_types

# Generate a column annotation.
cell_types_by_sample <- pData(meth_data)$Cell_type
row_anno_df <- data.frame(
  Cell_type = factor(cell_types_by_sample, levels = relevant_cell_types),
  row.names = rownames(pData(meth_data))
)
row_anno_obj <- rowAnnotation(
  df = row_anno_df, 
  col = list(Cell_type = cell_type_col_palette),
  annotation_label = "Cell type",
  annotation_legend_param = list(title = "Cell type", 
                                 title_gp = gpar(fontsize = 17, 
                                                 fontface = "bold"), 
                                 labels_gp = gpar(fontsize = 15), 
                                 grid_height = unit(8, "mm"),
                                 grid_width = unit(8, "mm"))
)

# Iterate over regions.
single_heatmaps <- lapply(1:length(overl_gr), FUN = function(x) {
  #-- Extract raw methylation values for each CpG in the corresponding region.
  raw_meth_this_region <- getMeth(meth_data, 
                                  regions = overl_gr[x], 
                                  what = "perBase", 
                                  type = "raw")[[1]]
  #-- Generate the heat map for this region.
  this_region_heatmap <- Heatmap(
    t(raw_meth_this_region),
    col = colorRamp2(seq(1, 0, length.out = 200), viridis(200)),
    cluster_rows = F,
    cluster_columns = F,
    show_row_names = F,
    show_column_names = F,
    heatmap_legend_param = list(title = "Methylation", 
                                labels_gp = gpar(fontsize = 15), 
                                title_gp = gpar(fontsize = 17, 
                                                fontface = "bold")),
    left_annotation = row_anno_obj,
    show_heatmap_legend = T,
    column_title = ifelse(
      x <= length(links_to_show_inds),
      yes = paste(links$Methylation__region_ID[links_to_show_inds[x]], 
                  links$Accessibility__peak_ID[links_to_show_inds[x]],
                  links$Expression__gene_name[links_to_show_inds[x]],
                  sep = ",\n"),
      no = "Custom\nFOXP3 region"
    )
  )
  return(this_region_heatmap)
})

# Save the heat maps.
heatmap_pdf <- paste0(
  plot_outdir, 
  "/dmr_peak_gene_links_blood_naive_treg_skin_treg_meth_heatmmaps_for_trackplots.pdf"
)
heatmap_rds <- paste0(
  plot_rds_outdir, 
  "/dmr_peak_gene_links_blood_naive_treg_skin_treg_meth_heatmmaps_for_trackplots.rds"
)
pdf(heatmap_pdf, width = 10, height = 6)
void <- lapply(single_heatmaps, FUN = draw)
dev.off()
saveRDS(single_heatmaps, file = heatmap_rds)