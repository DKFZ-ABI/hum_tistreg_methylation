# This script generates validation box plots for the hypothesis that LTR 
#   transposons are more relevant on the methylation level than on the chromatin 
#   accessibility level during skin Treg differentiation.
# Author: Niklas Beumer



# Load required packages.
library(bsseq)
library(GenomicRanges)
library(tidyr)
library(ggplot2)
library(Seurat)
library(Signac)
library(testit)


# Specify a location on /xxx.
b330_space <- "/xxx/"
location <- paste0(b330_space, "nbeumer/hm_treg_bs_rgnsbg")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/yyy/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the BSseq object containing the smoothed methylation data and reduce 
# to skin Tregs and blood naive Tregs.
meth_data <- readRDS(paste(
  location, 
  "/preprocessing_etc",
  "quality_control_after_alignment",
  "2022-03-14_bsseq_object_combined_all_samples_smoothed.rds", sep = "/"))
cell_types <- c("Skin Treg", "Blood naive Treg")
meth_data_pdata <- pData(meth_data)
meth_data_rel <- meth_data[, meth_data_pdata$Cell_type %in% cell_types]
meth_data_rel_pdata <- pData(meth_data_rel)

# Read in the Seurat object containing scATAC-seq data for CD4+ T cells
sc_data_cd4_file <- paste0(
  b330_space, 
  "msimon/scATAC/final/seurat_human_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar.RDS"
)
sc_data_cd4 <- readRDS(sc_data_cd4_file)
cell_types_atac_naming <- c("skin_treg", "blood_naive_treg")

# Specify the identity of donors that provided blood and skin samples in 
# scATACC-seq.
blood_donors <- c("1", "2")
skin_donors <- c("4", "5")

# Restrict scATAC-seq data to the relevant cell types and to cells from donors
# that donated the relevant tissue(s).
sc_data_cd4_rel <- subset(
  sc_data_cd4, 
  subset = (treg_tconv_annot == "blood_naive_treg" & donor %in% blood_donors) |
    (treg_tconv_annot == "skin_treg" & donor %in% skin_donors)
)
rm(sc_data_cd4) # Free up memory.
sc_data_cd4_rel <- RunTFIDF(sc_data_cd4_rel)
# sc_data_cd4_rel <- sc_data_cd4_rel[, 1:5000] # For debugging purposes.

# Specify the colour code for the cell types.
cell_type_col_palette <- c("blue", "darkorchid1")
names(cell_type_col_palette) <- cell_types

# Read in the locations of repeat elements from RepeatMasker
rmsk_file <- paste0(location, 
                    "/external_data/2023-09-14_rmsk_hg19_from_ucsc.txt.gz")
rmsk <- read.table(rmsk_file, header = F, stringsAsFactors = F)

# Filter the repeat element data so that only LTR transposons remain.
rmsk_filt <- rmsk[rmsk$V12 == "LTR", ]

# Generate a GRanges object for the annotated TE insertion sites.
# Take into account that this table originated from UCSC. Intervals are thus 
# very likely 0-based and closed.
# Modify seqlevels styles so that they match with what I used for the signature
# regions.
rmsk_filt_gr <- makeGRangesFromDataFrame(rmsk_filt,
                                         ignore.strand = T,
                                         seqnames.field = "V6",
                                         start.field = "V7",
                                         end.field = "V8",
                                         starts.in.df.are.0based = T,
                                         keep.extra.columns = T)
seqlevelsStyle(rmsk_filt_gr) <- "NCBI"

# Restrict the TEs to chromosomes for which I have methylation data.
seqlevels_to_keep <- intersect(seqlevels(rowRanges(meth_data_rel)),
                               seqlevels(rmsk_filt_gr))
rmsk_filt_2_gr <- rmsk_filt_gr[seqnames(rmsk_filt_gr) %in% seqlevels_to_keep]
# rmsk_filt_2_gr <- rmsk_filt_2_gr[1:100000] # For debugging purposes






###############################################################################
# Generate a box plot of raw methylation stratified by sample and cell type.
###############################################################################

# Compute raw methylation values for all TE regions. Remove all TE occurrences
# where at least one bisulphite sample returns NA (to make the box plots 
# comparable between the samples). Afterwards, save this data.
avg_raw_meth <- getMeth(meth_data_rel, regions = rmsk_filt_2_gr, 
                        type = "raw", what = "perRegion")
rows_to_remove <- which(apply(avg_raw_meth, 1, FUN = function(x) {
  length(which(is.na(x))) > 0
}))
rmsk_filt_2 <- as.data.frame(rmsk_filt_2_gr)
avg_raw_meth_df <- cbind(rmsk_filt_2[-rows_to_remove, ], 
                         avg_raw_meth[-rows_to_remove, ])
avg_raw_meth_df_melted <- pivot_longer(avg_raw_meth_df, 
                                       rownames(meth_data_rel_pdata),
                                       "Sample",
                                       values_to ="Avg_raw_meth")
avg_raw_meth_outfile <- paste0(
  location, 
  "/te_analysis/gen_loc_epigen_validation_ltr_discrep_betw_meth_and_atac_boxplots_data_meth.txt"
)
write.table(avg_raw_meth_df_melted, file = avg_raw_meth_outfile, 
            sep = "\t", row.names = F, col.names = F, quote = F)

# Generate a box plot showing how strongly the TE regions are methylated in 
# the different cell types (based on raw methylation).
avg_raw_meth_df_melted$Cell_type <- 
  sapply(avg_raw_meth_df_melted$Sample, FUN = function(x) {
    strsplit(strsplit(x, split = "(", fixed = T)[[1]][2], split = ";")[[1]][1]
  })
avg_raw_meth_df_melted$Cell_type <- factor(avg_raw_meth_df_melted$Cell_type, 
                                           levels = cell_types)
avg_raw_meth_df_melted$Sample <- factor(
  avg_raw_meth_df_melted$Sample, 
  levels = colnames(meth_data_rel))
te_box_plot_meth <- ggplot(avg_raw_meth_df_melted) +
  aes(x = Sample, y = Avg_raw_meth, colour = Cell_type) +
  scale_y_continuous(limits = c(0, 1), 
                     name = "Raw methylation\n(within-region mean)") +
  scale_colour_manual(breaks = cell_types, 
                      values = cell_type_col_palette,
                      name = "Cell type") +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0),
        panel.spacing.y = unit(0.1, "inches"))
te_box_plot_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_epigen_validation_ltr_discrep_boxplot_meth.pdf"
)
te_box_plot_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_epigen_validation_ltr_discrep_boxplot_meth.rds"
)
pdf(te_box_plot_pdf, width = 5.5, height = 3.5)
print(te_box_plot_meth)
dev.off()
saveRDS(te_box_plot_meth, file = te_box_plot_rds)






################################################################################
# Generate a box plot of normalised acsessibility stratified by donor and cell 
# type. The normalised accessibility values are based on peaks that overlap with 
# the respective TEs.
################################################################################

# Make a GRanges object out of the peak locations.
peaks_gr <- StringToGRanges(rownames(sc_data_cd4_rel), 
                            starts.in.df.are.0based = T)
seqlevelsStyle(peaks_gr) <- "NCBI"

# For each relevant TE location, find peaks that overlap with the TE insertion
# site.
te_locs_peaks_overl <- findOverlaps(rmsk_filt_2_gr, peaks_gr, ignore.strand = T)
corresp_peak_inds <- unique(to(te_locs_peaks_overl))
peaks_by_tes <- data.frame(TE = rep("LTR", length(corresp_peak_inds)),
                           Peak = rownames(sc_data_cd4_rel)[corresp_peak_inds])

# For each peak overlapping with a relevant TE, compute the average normalised 
# accessibility (per donor and cell type).
norm_acc_matr <- GetAssayData(sc_data_cd4_rel, assay = "scATAC_raw", 
                              slot = "data")
assert(all(colnames(norm_acc_matr) == colnames(sc_data_cd4_rel)))
peaks_by_tes_melted <- do.call(
  rbind, 
  lapply(cell_types_atac_naming, FUN = function(x) {
    corresp_donors <- sort(unique(
      sc_data_cd4_rel$donor[sc_data_cd4_rel$treg_tconv_annot == x]
    ))
    mean_acc_this_celltype <- do.call(
      rbind, 
      lapply(corresp_donors, FUN = function(y) {
        corresp_cells <- colnames(sc_data_cd4_rel) [
          sc_data_cd4_rel$treg_tconv_annot == x & sc_data_cd4_rel$donor == y
        ]
        norm_acc_matr_subs <- norm_acc_matr[peaks_by_tes$Peak, 
                                            corresp_cells]
        acc_vals <- rowMeans(norm_acc_matr_subs)
        temp_df <- data.frame(Cell_type = x,
                              Donor = y,
                              Celltype_donor_comb = paste0(x, "_donor", y),
                              Mean_norm_acc = acc_vals)
        return(cbind(peaks_by_tes, temp_df))
      }
    ))
    return(mean_acc_this_celltype)
  }
))

# Save the computed mean accessibility values.
peaks_by_tes_outfile <- paste0(
  location, 
  "/te_analysis/gen_loc_epigen_validation_ltr_discrep_betw_meth_and_atac_boxplots_data_atac.txt"
)
write.table(peaks_by_tes_melted, file = peaks_by_tes_outfile, 
            sep = "\t", row.names = F, col.names = F, quote = F)


# Generate a box plot showing how accessible the peaks overlapping with TE 
# regions are in the different cell types.
peaks_by_tes_melted$Celltype_donor_comb <- factor(
  peaks_by_tes_melted$Celltype_donor_comb,
  levels = unique(peaks_by_tes_melted$Celltype_donor_comb)
)
te_box_plot_atac <- ggplot(peaks_by_tes_melted) +
  aes(x = Celltype_donor_comb, y = Mean_norm_acc, colour = Cell_type) +
  # scale_x_discrete(breaks = cell_types_atac_naming,
  #                  labels = cell_types,
  #                  name = "Cell type") +
  scale_y_continuous(name = "Mean norm. acc.") +
  scale_colour_manual(breaks = cell_types_atac_naming,
                      values = unname(cell_type_col_palette),
                      labels = cell_types,
                      name = "Cell type") +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0),
        panel.spacing.y = unit(0.1, "inches"))
te_box_plot_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_epigen_validation_ltr_discrep_boxplot_atac.pdf"
)
te_box_plot_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_epigen_validation_ltr_discrep_boxplot_atac.rds"
)
pdf(te_box_plot_pdf, width = 5, height = 4.5)
print(te_box_plot_atac)
dev.off()
saveRDS(te_box_plot_atac, file = te_box_plot_rds)




