# This script generates validation box plots for the epigenetics of TEs that 
#   were shown to be relevant on the RNA level for the comparison between
#   skin Treg cells and blood naive Treg cells and also enriched in DMRs between 
#   these cell types.
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
# to the relevant cell types.
meth_data <- readRDS(paste(
  location, 
  "preprocessing_etc",
  "quality_control_after_alignment",
  "2022-03-14_bsseq_object_combined_all_samples_smoothed.rds", sep = "/"))
cell_types <- c("Skin Treg", "Blood naive Treg")
cell_types_2 <- c("Skin Treg", "Blood CCR8+ Treg", "Blood naive Treg")
meth_data_pdata <- pData(meth_data)
meth_data_rel <- meth_data[, meth_data_pdata$Cell_type %in% cell_types]
meth_data_rel_pdata <- pData(meth_data_rel)
meth_data_rel_2 <- meth_data[, meth_data_pdata$Cell_type %in% cell_types_2]
meth_data_rel_2_pdata <- pData(meth_data_rel_2)
rm(meth_data) # Free up memory.

# Read in the Seurat object containing scATAC-seq data for CD4+ T cells, reduce 
# to the relevant cell types and perform a new TF-IDF.
sc_data_cd4_file <- paste0(
  b330_space, 
  "msimon/scATAC/final/seurat_human_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar.RDS"
)
sc_data_cd4 <- readRDS(sc_data_cd4_file)
cell_types_atac_naming <- c("skin_treg", "blood_naive_treg")
cell_types_atac_naming_2 <- c("skin_treg", "blood_ccr8_treg", 
                              "blood_naive_treg")

# Specify the identity of donors that provided blood and skin samples in 
# scATAC-seq.
blood_donors <- c("1", "2")
skin_donors <- c("4", "5")

# Restrict scATAC-seq data to the relevant cell types and to cells from donors
# that donated the relevant tissue(s).
sc_data_cd4_rel <- subset(
  sc_data_cd4, 
  subset = (treg_tconv_annot == "blood_naive_treg" & donor %in% blood_donors) |
    (treg_tconv_annot == "skin_treg" & donor %in% skin_donors)
)
sc_data_cd4_rel_2 <- subset(
  sc_data_cd4, 
  subset = (treg_tconv_annot == "blood_naive_treg" & donor %in% blood_donors) |
    (treg_tconv_annot == "skin_treg" & donor %in% skin_donors) |
    (treg_tconv_annot == "blood_ccr8_treg" & donor %in% blood_donors)
)
rm(sc_data_cd4) # Free up memory.
sc_data_cd4_rel <- RunTFIDF(sc_data_cd4_rel)
sc_data_cd4_rel_2 <- RunTFIDF(sc_data_cd4_rel_2)
# sc_data_cd4_rel <- sc_data_cd4_rel[, 1:5000] # For debugging purposes.
# sc_data_cd4_rel_2 <- sc_data_cd4_rel_2[, 1:5000] # For debugging purposes.

# Specify the colour code for the cell types.
cell_type_col_palette <- c("blue", "darkorchid1")
names(cell_type_col_palette) <- cell_types
cell_type_col_palette_2 <- c("blue", "orange", "darkorchid1")
names(cell_type_col_palette_2) <- cell_types_2


# Read in the locations of repeat elements from RepeatMasker
rmsk_file <- paste0(location, 
                    "/external_data/2023-09-14_rmsk_hg19_from_ucsc.txt.gz")
rmsk <- read.table(rmsk_file, header = F, stringsAsFactors = F)

# Filter the repeat element data so that only those elements remain where I am
# confident that they are transposons.
conf_transp_cats <- c("LINE", "DNA", "SINE", "LTR", "RC", "Retroposon")
rmsk_filt <- rmsk[rmsk$V12 %in% conf_transp_cats, ]
rm(rmsk) # Free up memory.


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

# Specify the two TEs that are relevant for this analysis. These two TEs are
# overexpressed in skin Tregs with respect to blood naive Tregs and are enriched
# in the skin_treg__hypomethylation regions of this comparison. It also happens
# that they are not significnatly enriched in the skin_treg__hyperaccessibility
# regions. 
relevant_tes <- c("HERVIP10F-int", "LTR45B")

# Make a GRanges object containing locations of the TEs to show.
rmsk_filt_specialtes_gr <- rmsk_filt_2_gr[rmsk_filt_2_gr$V11 %in% relevant_tes]




###############################################################################
# Generate box plots of raw methylation stratified by sample, cell type
# and element name.
###############################################################################

# Compute raw methylation values for all TE regions. Remove all TE occurrences
# where at least one bisulphite sample returns NA (to make the box plots 
# comparable between the samples). Afterwards, save this data.
rmsk_filt_specialtes <- as.data.frame(rmsk_filt_specialtes_gr)
avg_raw_meth <- getMeth(meth_data_rel, regions = rmsk_filt_specialtes_gr, 
                        type = "raw", what = "perRegion")
rows_to_remove <- which(apply(avg_raw_meth, 1, FUN = function(x) {
  length(which(is.na(x))) > 0
}))
avg_raw_meth_df <- cbind(rmsk_filt_specialtes[-rows_to_remove, ], 
                         avg_raw_meth[-rows_to_remove, ])
avg_raw_meth_df_melted <- pivot_longer(avg_raw_meth_df, 
                                       rownames(meth_data_rel_pdata),
                                       "Sample",
                                       values_to ="Avg_raw_meth")
avg_raw_meth_outfile <- paste0(
  location, 
  "/te_analysis/gen_loc_epigen_validation_tes_w_relevance_on_rna_level_boxplots_data_by_special_elements_meth.txt"
)
write.table(avg_raw_meth_df_melted, file = avg_raw_meth_outfile, 
            sep = "\t", row.names = F, col.names = F, quote = F)

# Repeat the previous step considering blood CCR8+ Tregs as third cell type.
avg_raw_meth_2 <- getMeth(meth_data_rel_2, regions = rmsk_filt_specialtes_gr, 
                          type = "raw", what = "perRegion")
rows_to_remove <- which(apply(avg_raw_meth_2, 1, FUN = function(x) {
  length(which(is.na(x))) > 0
}))
avg_raw_meth_2_df <- cbind(rmsk_filt_specialtes[-rows_to_remove, ], 
                           avg_raw_meth_2[-rows_to_remove, ])
avg_raw_meth_2_df_melted <- pivot_longer(avg_raw_meth_2_df, 
                                         rownames(meth_data_rel_2_pdata),
                                         "Sample",
                                         values_to ="Avg_raw_meth")
avg_raw_meth_outfile <- paste0(
  location, 
  "/te_analysis/gen_loc_epigen_validation_tes_w_relevance_on_rna_level_boxplots_data_by_special_elements_meth_w_ccr8.txt"
)
write.table(avg_raw_meth_2_df_melted, file = avg_raw_meth_outfile, 
            sep = "\t", row.names = F, col.names = F, quote = F)

# Generate a box plot showing how strongly the TE regions are methylated in 
# the different cell types (based on raw methylation), stratified by TE name.
avg_raw_meth_df_melted$Cell_type <- 
  sapply(avg_raw_meth_df_melted$Sample, FUN = function(x) {
    strsplit(strsplit(x, split = "(", fixed = T)[[1]][2], split = ";")[[1]][1]
  })
avg_raw_meth_df_melted$Cell_type <- factor(avg_raw_meth_df_melted$Cell_type, 
                                           levels = cell_types)
avg_raw_meth_df_melted$Sample <- factor(
  avg_raw_meth_df_melted$Sample, 
  levels = colnames(meth_data_rel))
avg_raw_meth_df_melted$V11 <- factor(avg_raw_meth_df_melted$V11, 
                                     levels = relevant_tes)
te_box_plot_meth <- ggplot(avg_raw_meth_df_melted) +
  aes(x = Sample, y = Avg_raw_meth, colour = Cell_type) +
  scale_y_continuous(limits = c(0, 1), 
                     name = "Raw methylation (within-region mean)") +
  scale_colour_manual(breaks = cell_types, 
                      values = cell_type_col_palette,
                      name = "Cell type") +
  geom_boxplot() +
  facet_wrap(~V11, ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0),
        panel.spacing.y = unit(0.1, "inches"))
te_box_plot_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_epigen_validation_tes_w_relevance_on_rnalevel_boxplot_meth.pdf"
)
te_box_plot_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_epigen_validation_tes_w_relevance_on_rnalevel_boxplot_meth.rds"
)
pdf(te_box_plot_pdf, width = 5, height = 4.5)
print(te_box_plot_meth)
dev.off()
saveRDS(te_box_plot_meth, file = te_box_plot_rds)

# Generate the same box plot including blood CCR8+ Tregs.
avg_raw_meth_2_df_melted$Cell_type <- 
  sapply(avg_raw_meth_2_df_melted$Sample, FUN = function(x) {
    strsplit(strsplit(x, split = "(", fixed = T)[[1]][2], split = ";")[[1]][1]
  })
avg_raw_meth_2_df_melted$Cell_type <- factor(avg_raw_meth_2_df_melted$Cell_type, 
                                             levels = cell_types_2)
avg_raw_meth_2_df_melted$Sample <- factor(
  avg_raw_meth_2_df_melted$Sample, 
  levels = colnames(meth_data_rel_2))
avg_raw_meth_2_df_melted$V11 <- factor(avg_raw_meth_2_df_melted$V11, 
                                       levels = relevant_tes)
te_box_plot_meth_2 <- ggplot(avg_raw_meth_2_df_melted) +
  aes(x = Sample, y = Avg_raw_meth, colour = Cell_type) +
  scale_y_continuous(limits = c(0, 1), 
                     name = "Raw methylation (within-region mean)") +
  scale_colour_manual(breaks = cell_types_2, 
                      values = cell_type_col_palette_2,
                      name = "Cell type") +
  geom_boxplot() +
  facet_wrap(~V11, ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0),
        panel.spacing.y = unit(0.1, "inches"))
te_box_plot_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_epigen_validation_tes_w_relevance_on_rnalevel_boxplot_meth_w_ccr8.pdf"
)
te_box_plot_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_epigen_validation_tes_w_relevance_on_rnalevel_boxplot_meth_w_ccr8.rds"
)
pdf(te_box_plot_pdf, width = 6.5, height = 4.5)
print(te_box_plot_meth_2)
dev.off()
saveRDS(te_box_plot_meth_2, file = te_box_plot_rds)






################################################################################
# Generate box plots of normalised acessibility stratified by sample, cell type
# and element name. The normalised accessibility values are based on peaks that 
# overlap with the respective TEs.
################################################################################

# Make a GRanges object out of the peak locations.
peaks_gr <- StringToGRanges(rownames(sc_data_cd4_rel), 
                            starts.in.df.are.0based = T)
seqlevelsStyle(peaks_gr) <- "NCBI"

# For each relevant TE location, find peaks that overlap with the TE insertion
# site.
peaks_by_tes <- do.call(rbind, lapply(relevant_tes, FUN = function(x) {
  te_locs <- rmsk_filt_specialtes_gr[rmsk_filt_specialtes_gr$V11 == x]
  te_locs_peaks_overl <- findOverlaps(te_locs, peaks_gr, ignore.strand = T)
  corresp_peak_inds <- unique(to(te_locs_peaks_overl))
  temp_df <- data.frame(TE = rep(x, length(corresp_peak_inds)),
                        Peak = rownames(sc_data_cd4_rel)[corresp_peak_inds])
  return(temp_df)
}))

# For each peak overlapping with a relevant TE, compute the average normalised 
# accessibility (per cell type). Save these values.
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

peaks_by_tes_outfile <- paste0(
  location, 
  "/te_analysis/gen_loc_epigen_validation_tes_w_relevance_on_rna_level_boxplots_data_by_special_elements_atac.txt"
)
write.table(peaks_by_tes_melted, file = peaks_by_tes_outfile, 
            sep = "\t", row.names = F, col.names = F, quote = F)

# Repeat the previous step while also considering blood CCR8+ Tregs.
norm_acc_matr_2 <- GetAssayData(sc_data_cd4_rel_2, assay = "scATAC_raw", 
                                slot = "data")
assert(all(colnames(norm_acc_matr_2) == colnames(sc_data_cd4_rel_2)))
peaks_by_tes_melted_2 <- do.call(
  rbind, 
  lapply(cell_types_atac_naming_2, FUN = function(x) {
    corresp_donors <- sort(unique(
      sc_data_cd4_rel_2$donor[sc_data_cd4_rel_2$treg_tconv_annot == x]
    ))
    mean_acc_this_celltype <- do.call(
      rbind, 
      lapply(corresp_donors, FUN = function(y) {
        corresp_cells <- colnames(sc_data_cd4_rel_2) [
          sc_data_cd4_rel_2$treg_tconv_annot == x & sc_data_cd4_rel_2$donor == y
        ]
        norm_acc_matr_subs <- norm_acc_matr_2[peaks_by_tes$Peak, 
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

peaks_by_tes_outfile <- paste0(
  location, 
  "/te_analysis/gen_loc_epigen_validation_tes_w_relevance_on_rna_level_boxplots_data_by_special_elements_atac_w_ccr8.txt"
)
write.table(peaks_by_tes_melted_2, file = peaks_by_tes_outfile, 
            sep = "\t", row.names = F, col.names = F, quote = F)


# Generate a box plot showing how accessible the peaks overlapping with TE 
# regions are in the different cell types, stratified by TE name.
peaks_by_tes_melted$Celltype_donor_comb <- factor(
  peaks_by_tes_melted$Celltype_donor_comb, 
  levels = unique(peaks_by_tes_melted$Celltype_donor_comb)
)
peaks_by_tes_melted$TE <- factor(peaks_by_tes_melted$TE, levels = relevant_tes)
te_box_plot_atac <- ggplot(peaks_by_tes_melted) +
  aes(x = Celltype_donor_comb, y = Mean_norm_acc, colour = Cell_type) +
  scale_y_continuous(name = "Mean norm. acc.") +
  scale_colour_manual(breaks = cell_types_atac_naming,
                      values = unname(cell_type_col_palette),
                      labels = cell_types,
                      name = "Cell type") +
  geom_boxplot() +
  facet_wrap(~TE, ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0),
        panel.spacing.y = unit(0.1, "inches"))
te_box_plot_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_epigen_validation_tes_w_relevance_on_rnalevel_boxplot_atac.pdf"
)
te_box_plot_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_epigen_validation_tes_w_relevance_on_rnalevel_boxplot_atac.rds"
)
pdf(te_box_plot_pdf, width = 5, height = 5.5)
print(te_box_plot_atac)
dev.off()
saveRDS(te_box_plot_atac, file = te_box_plot_rds)

# Generate a box plot that also shows blood CCR8+ Tregs.
peaks_by_tes_melted_2$Celltype_donor_comb <- factor(
  peaks_by_tes_melted_2$Celltype_donor_comb,
  levels = unique(peaks_by_tes_melted_2$Celltype_donor_comb)
)
peaks_by_tes_melted_2$TE <- factor(peaks_by_tes_melted_2$TE, 
                                   levels = relevant_tes)
te_box_plot_atac_2 <- ggplot(peaks_by_tes_melted_2) +
  aes(x = Celltype_donor_comb, y = Mean_norm_acc, colour = Cell_type) +
  scale_y_continuous(name = "Mean norm. acc.") +
  scale_colour_manual(breaks = cell_types_atac_naming_2,
                      values = unname(cell_type_col_palette_2),
                      labels = cell_types_2,
                      name = "Cell type") +
  geom_boxplot() +
  facet_wrap(~TE, ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0),
        panel.spacing.y = unit(0.1, "inches"))
te_box_plot_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_epigen_validation_tes_w_relevance_on_rnalevel_boxplot_atac_w_ccr8.pdf"
)
te_box_plot_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_epigen_validation_tes_w_relevance_on_rnalevel_boxplot_atac_w_ccr8.rds"
)
pdf(te_box_plot_pdf, width = 5, height = 5.5)
print(te_box_plot_atac_2)
dev.off()
saveRDS(te_box_plot_atac_2, file = te_box_plot_rds)