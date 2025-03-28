# This script generates violin plots and related diagrams for methylation and 
#   chromatin accessibility of HERVIP10F-int and LTR45B insertion 
#   sites for skin Treg cells, blood CCR8+ Treg cells and blood CD45RA+ Treg 
#   cells.
# Author: Niklas Beumer
# Run this script with 10 cores!!!



# Load required packages.
library(bsseq)
library(GenomicRanges)
library(tidyr)
library(ggplot2)
library(Seurat)
library(Signac)
library(testit)
library(parallel)


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
cell_types <- c("Skin Treg", "Blood CCR8+ Treg", "Blood naive Treg")
meth_data_pdata <- pData(meth_data)
meth_data_rel <- meth_data[, meth_data_pdata$Cell_type %in% cell_types]
meth_data_rel_pdata <- pData(meth_data_rel)
rm(meth_data) # Free up memory.

# Read in the Seurat object containing scATAC-seq data for CD4+ T cells, reduce 
# to the relevant cell types and perform a new TF-IDF.
sc_data_cd4_file <- paste0(
  b330_space, 
  "msimon/scATAC/final/seurat_human_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar.RDS"
)
sc_data_cd4 <- readRDS(sc_data_cd4_file)
cell_types_atac_naming <- c("skin_treg", "blood_ccr8_treg", "blood_naive_treg")

# Specify the identity of donors that provided blood and skin samples in 
# scATAC-seq.
blood_donors <- c("1", "2")
skin_donors <- c("4", "5")

# Restrict scATAC-seq data to the relevant cell types and to cells from donors
# that donated the relevant tissue(s).
sc_data_cd4_rel <- subset(
  sc_data_cd4, 
  subset = (treg_tconv_annot == "blood_naive_treg" & donor %in% blood_donors) |
    (treg_tconv_annot == "skin_treg" & donor %in% skin_donors) |
    (treg_tconv_annot == "blood_ccr8_treg" & donor %in% blood_donors)
)
rm(sc_data_cd4) # Free up memory.
sc_data_cd4_rel <- RunTFIDF(sc_data_cd4_rel)
# sc_data_cd4_rel <- sc_data_cd4_rel[, 1:5000] # For debugging purposes.

# Specify the colour code for the cell types.
cell_type_col_palette <- c("blue", "orange", "darkorchid1")
names(cell_type_col_palette) <- cell_types

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

# Specify the two TEs that are relevant for this analysis.
relevant_tes <- c("HERVIP10F-int", "LTR45B")

# Make a GRanges object containing locations of the TEs to show.
rmsk_filt_specialtes_gr <- rmsk_filt_2_gr[rmsk_filt_2_gr$V11 %in% relevant_tes]






###############################################################################
# Generate violin plots of raw methylation stratified by sample, cell type
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
  "/te_analysis/gen_loc_hervip10f_ltr45b_vio_plot_skin_bloodccr8_bloodnaive_data_meth.txt"
)
write.table(avg_raw_meth_df_melted, file = avg_raw_meth_outfile, 
            sep = "\t", row.names = F, col.names = F, quote = F)

# Generate a violin plot showing how strongly the TE regions are methylated in 
# the different cell types (based on raw methylation), stratified by TE element.
# Highlight median values inside the violins.
avg_raw_meth_df_melted$Cell_type <- 
  sapply(avg_raw_meth_df_melted$Sample, FUN = function(x) {
    strsplit(strsplit(x, split = "(", fixed = T)[[1]][2], split = ";")[[1]][1]
  })
avg_raw_meth_df_melted$Cell_type <- factor(avg_raw_meth_df_melted$Cell_type, 
                                           levels = cell_types)
avg_raw_meth_df_melted$Sample <- factor(
  avg_raw_meth_df_melted$Sample, 
  levels = colnames(meth_data_rel))
te_violin_plot_meth <- ggplot(avg_raw_meth_df_melted) +
  aes(x = Sample, y = Avg_raw_meth, fill = Cell_type) +
  scale_y_continuous(limits = c(0, 1), 
                     name = "Raw methylation (within-region mean)",
                     expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(breaks = cell_types, 
                    values = cell_type_col_palette,
                    name = "Cell type") +
  geom_violin(draw_quantiles = 0.5, scale = "width") +
  facet_wrap(~ V11, ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0),
        panel.spacing.y = unit(0.1, "inches"))
te_violin_plot_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_hervip10f_ltr45b_vio_plot_skin_bloodccr8_bloodnaive_violinplot_meth.pdf"
)
te_violin_plot_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_hervip10f_ltr45b_vio_plot_skin_bloodccr8_bloodnaive_violinplot_meth.rds"
)
pdf(te_violin_plot_pdf, width = 5, height = 4.5)
print(te_violin_plot_meth)
dev.off()
saveRDS(te_violin_plot_meth, file = te_violin_plot_rds)



########## Related diagram 1: 
# Dot plot showing, for every sample, the median of raw methylation acrosss the 
# considered TE insertion sites.
dot_plot_data <- do.call(rbind, lapply(cell_types, FUN = function(x) {
  corresp_samp <- rownames(meth_data_rel_pdata)[
    meth_data_rel_pdata$Cell_type == x
  ]
  this_celltype_df <- do.call(
    rbind, 
    lapply(relevant_tes, FUN = function(y) {
      this_class_df <- do.call(rbind, lapply(corresp_samp, FUN = function(z) {
        meth_vals <- avg_raw_meth_df[avg_raw_meth_df$V11 == y, z]
        median_meth <- median(meth_vals)
        temp_df <- data.frame(Cell_type = x,
                              TE_element = y,
                              Sample = z,
                              Median_meth = median_meth)
      return(temp_df)
    }))
    return(this_class_df)
  }))
  return(this_celltype_df)
}))
dot_plot_data$Cell_type <- factor(dot_plot_data$Cell_type, levels = cell_types)
te_dot_plot_meth <- ggplot(dot_plot_data) +
  aes(x = Cell_type, y = Median_meth, colour = Cell_type) +
  scale_y_continuous(name = "Median raw methylation across TE sites",
                     limits = c(0, 1),
                     expand = expansion(mult = c(0.05, 0.1))) +
  scale_colour_manual(breaks = cell_types, 
                      values = cell_type_col_palette,
                      name = "Cell type") +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 20)) +
  facet_wrap(~ TE_element, ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0),
        panel.spacing.y = unit(0.1, "inches"))
dot_plot_meth_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_hervip10f_ltr45b_vio_plot_skin_bloodccr8_bloodnaive_violinplot_meth_related_dotplot.pdf"
)
dot_plot_meth_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_hervip10f_ltr45b_vio_plot_skin_bloodccr8_bloodnaive_violinplot_meth_related_dotplot.rds"
)
pdf(dot_plot_meth_pdf, width = 3, height = 4.5)
print(te_dot_plot_meth)
dev.off()
saveRDS(te_dot_plot_meth, file = dot_plot_meth_rds)
  


########## Related diagram 2: 
# For each considered TE insertion site and each combination of two cell types, 
# the difference between mean methylation in each cell type is computed. The 
# distribution of those differences is shown by a ridge plot.
cell_type_comps <- do.call(c, lapply(1:(length(cell_types) - 1), 
                                     FUN = function(x) {
  lapply((x + 1):length(cell_types), FUN = function(y) {
    return(c(cell_types[x], cell_types[y]))
  })
}))
ridge_plot_data <- do.call(rbind, lapply(cell_type_comps, FUN = function(comp) {
  cell_type_1_samp <- rownames(meth_data_rel_pdata)[
    meth_data_rel_pdata$Cell_type == comp[1]
  ]
  cell_type_2_samp <- rownames(meth_data_rel_pdata)[
    meth_data_rel_pdata$Cell_type == comp[2]
  ]
  this_comp_data <- do.call(
    rbind, 
    lapply(relevant_tes, FUN = function(x) {
      corresp_df <- avg_raw_meth_df[avg_raw_meth_df$V11 == x, ]
      meth_diffs <- do.call(
        rbind, 
        mclapply(1:nrow(corresp_df), FUN = function(y) {
          cell_type_1_mean <- mean(as.numeric(corresp_df[y, cell_type_1_samp]))
          cell_type_2_mean <- mean(as.numeric(corresp_df[y, cell_type_2_samp]))
          diff_to_show <- cell_type_1_mean - cell_type_2_mean
          temp_df <- data.frame(TE_element = x, 
                                meth_diff = diff_to_show,
                                Comparison = paste(comp, collapse = " vs. "))
          temp_df <- cbind(temp_df, corresp_df[y, ])
          return(temp_df)
        }, mc.cores = 10)
      )
      return(meth_diffs)
    })
  )
  return(this_comp_data)
}))
te_counts <- do.call(rbind, lapply(cell_type_comps, FUN = function(comp) {
  this_comp_df <- do.call(
    rbind, 
    lapply(relevant_tes, FUN = function(x) {
      num_pos <- length(which(
        ridge_plot_data$meth_diff[
          ridge_plot_data$TE_element == x & 
            ridge_plot_data$Comparison == paste(comp, collapse = " vs. ")
        ] > 0
      ))
      num_neg <- length(which(
         ridge_plot_data$meth_diff[
           ridge_plot_data$TE_element == x & 
             ridge_plot_data$Comparison == paste(comp, collapse = " vs. ")
          ] < 0
      ))
      temp_df <- data.frame(Comparison = paste(comp, collapse = " vs. "),
                            TE_element = x,
                            meth_diff = c(min(ridge_plot_data$meth_diff) + 0.3,
                                          max(ridge_plot_data$meth_diff) - 0.3),
                            text = c(num_neg, num_pos))
      return(temp_df)
    })
  )
  return(this_comp_df)
}))
te_counts$y_pos <- rep(
  seq(1, 3, length.out = length(cell_type_comps)),
  each = 2 * length(relevant_tes)
)
te_ridge_plot_meth <- ggplot(ridge_plot_data) +
  aes(x = meth_diff, colour = Comparison) +
  scale_x_continuous(
    name = "Diff. in mean methylation across samples"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    name = "Density of TE insertion sites [Arbitrary units]"
  ) +
  scale_colour_manual(values = c("grey50", "turquoise1", "lightpink1")) +
  geom_density() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = te_counts, mapping = aes(label = text, y = y_pos)) +
  facet_wrap(~ TE_element, ncol = 1) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))
ridge_plot_meth_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_hervip10f_ltr45b_vio_plot_skin_bloodccr8_bloodnaive_violinplot_meth_related_ridgeplot.pdf"
)
ridge_plot_meth_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_hervip10f_ltr45b_vio_plot_skin_bloodccr8_bloodnaive_violinplot_meth_related_ridgeplot.rds"
)
pdf(ridge_plot_meth_pdf, width = 7, height = 4.5)
print(te_ridge_plot_meth)
dev.off()
saveRDS(te_ridge_plot_meth, file = ridge_plot_meth_rds)





################################################################################
# Generate violin plots of normalised acessibility stratified by sample, cell 
# type and element name. The normalised accessibility values are based on peaks 
# that overlap with the respective TEs.
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
  "/te_analysis/gen_loc_hervip10f_ltr45b_vio_plot_skin_bloodccr8_bloodnaive_data_atac.txt"
)
write.table(peaks_by_tes_melted, file = peaks_by_tes_outfile, 
            sep = "\t", row.names = F, col.names = F, quote = F)


# Generate a violin plot showing how accessible the peaks overlapping with TE 
# regions are in the different cell types and donors, stratified by TE class.
# Highlight median values inside the violins.
peaks_by_tes_melted$Celltype_donor_comb <- factor(
  peaks_by_tes_melted$Celltype_donor_comb,
  levels = unique(peaks_by_tes_melted$Celltype_donor_comb)
)
te_violin_plot_atac <- ggplot(peaks_by_tes_melted) +
  aes(x = Celltype_donor_comb, y = Mean_norm_acc, fill = Cell_type) +
  scale_y_continuous(name = "Mean norm. acc.",
                     expand = expansion(mult = c(0.01, 0.1))) +
  scale_fill_manual(breaks = cell_types_atac_naming,
                    values = unname(cell_type_col_palette),
                    labels = cell_types,
                    name = "Cell type") +
  geom_violin(draw_quantiles = 0.5) +
  facet_wrap(~TE, ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0),
        panel.spacing.y = unit(0.1, "inches"))
te_violin_plot_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_hervip10f_ltr45b_vio_plot_skin_bloodccr8_bloodnaive_violinplot_atac.pdf"
)
te_violin_plot_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_hervip10f_ltr45b_vio_plot_skin_bloodccr8_bloodnaive_violinplot_atac.rds"
)
pdf(te_violin_plot_pdf, width = 5, height = 4.5)
print(te_violin_plot_atac)
dev.off()
saveRDS(te_violin_plot_atac, file = te_violin_plot_rds)





########## Related diagram 1: 
# Dot plot showing, for every sample, the median of normalised accessibility 
#   acrosss the considered TE insertion sites.
dot_plot_data <- do.call(rbind, lapply(
  cell_types_atac_naming, 
  FUN = function(x) {
    corresp_samp <- unique(peaks_by_tes_melted$Celltype_donor_comb[
      peaks_by_tes_melted$Cell_type == x
    ])
    this_celltype_df <- do.call(rbind, lapply(
      unique(peaks_by_tes_melted$TE), 
      FUN = function(y) {
        this_class_df <- do.call(rbind, lapply(corresp_samp, FUN = function(z) {
          acc_vals <- peaks_by_tes_melted$Mean_norm_acc[
            peaks_by_tes_melted$TE == y & 
              peaks_by_tes_melted$Celltype_donor_comb == z
          ]
          median_acc <- median(acc_vals)
          temp_df <- data.frame(Cell_type = x,
                                TE_element = y,
                                Sample = z,
                                Median_acc = median_acc)
          return(temp_df)
        }))
        return(this_class_df)
      }
    ))
    return(this_celltype_df)
  }
))
dot_plot_data$Cell_type <- factor(dot_plot_data$Cell_type, 
                                  levels = cell_types_atac_naming)
te_dot_plot_atac <- ggplot(dot_plot_data) +
  aes(x = Cell_type, y = Median_acc, colour = Cell_type) +
  scale_y_continuous(name = "Median norm. accessibility across TE sites",
                     limits = c(min(te_violin_plot_atac$data$Mean_norm_acc),
                                max(te_violin_plot_atac$data$Mean_norm_acc))) +
  scale_colour_manual(breaks = cell_types_atac_naming, 
                      values = unname(cell_type_col_palette),
                      name = "Cell type") +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 20)) +
  facet_wrap(~ TE_element, ncol = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0),
        panel.spacing.y = unit(0.1, "inches"))
dot_plot_atac_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_hervip10f_ltr45b_vio_plot_skin_bloodccr8_bloodnaive_violinplot_atac_related_dotplot.pdf"
)
dot_plot_atac_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_hervip10f_ltr45b_vio_plot_skin_bloodccr8_bloodnaive_violinplot_atac_related_dotplot.rds"
)
pdf(dot_plot_atac_pdf, width = 3, height = 4.5)
print(te_dot_plot_atac)
dev.off()
saveRDS(te_dot_plot_atac, file = dot_plot_atac_rds)



########## Related diagram 2: 
# For each considered TE insertion site and each combination of two cell types, 
# the difference between mean methylation in each cell type is computed. The 
# distribution of those differences is shown by a ridge plot.
cell_type_comps <- do.call(
  c, 
  lapply(1:(length(cell_types_atac_naming) - 1), FUN = function(x) {
    lapply((x + 1):length(cell_types_atac_naming), FUN = function(y) {
      return(c(cell_types_atac_naming[x], cell_types_atac_naming[y]))
    })
  })
)
ridge_plot_data <- do.call(rbind, lapply(cell_type_comps, FUN = function(comp) {
  cell_type_1_samp <- unique(peaks_by_tes_melted$Celltype_donor_comb[
    peaks_by_tes_melted$Cell_type == comp[1]
  ])
  cell_type_2_samp <- unique(peaks_by_tes_melted$Celltype_donor_comb[
    peaks_by_tes_melted$Cell_type == comp[2]
  ])
  this_comp_data <- do.call(
    rbind, 
    lapply(unique(peaks_by_tes_melted$TE), FUN = function(x) {
      corresp_df <- peaks_by_tes_melted[
        peaks_by_tes_melted$TE == x & peaks_by_tes_melted$Cell_type %in% comp, 
      ]
      assert(length(unique(corresp_df$Peak)) == 
               nrow(corresp_df) / 
               (length(c(cell_type_1_samp, cell_type_2_samp))))
      corresp_peaks <- unique(corresp_df$Peak)
      acc_diffs <- do.call(
        rbind, 
        mclapply(corresp_peaks, FUN = function(y) {
          cell_type_1_mean <- mean(as.numeric(corresp_df$Mean_norm_acc[
            corresp_df$Peak == y & 
              corresp_df$Celltype_donor_comb %in% cell_type_1_samp
          ]))
          cell_type_2_mean <- mean(as.numeric(corresp_df$Mean_norm_acc[
            corresp_df$Peak == y & 
              corresp_df$Celltype_donor_comb %in% cell_type_2_samp
          ]))
          diff_to_show <- cell_type_1_mean - cell_type_2_mean
          temp_df <- data.frame(TE_element = x, 
                                acc_diff = diff_to_show,
                                Comparison = paste(comp, collapse = " vs. "),
                                Peak = y)
          return(temp_df)
        }, mc.cores = 10)
      )
      return(acc_diffs)
    })
  )
  return(this_comp_data)
}))
te_counts <- do.call(rbind, lapply(cell_type_comps, FUN = function(comp) {
  this_comp_df <- do.call(
    rbind, 
    lapply(unique(peaks_by_tes_melted$TE), FUN = function(x) {
      num_pos <- length(which(
        ridge_plot_data$acc_diff[
          ridge_plot_data$TE_element == x & 
            ridge_plot_data$Comparison == paste(comp, collapse = " vs. ")
        ] > 0
      ))
      num_neg <- length(which(
        ridge_plot_data$acc_diff[
          ridge_plot_data$TE_element == x & 
            ridge_plot_data$Comparison == paste(comp, collapse = " vs. ")
        ] < 0
      ))
      temp_df <- data.frame(Comparison = paste(comp, collapse = " vs. "),
                            TE_element = x,
                            acc_diff = c(min(ridge_plot_data$acc_diff) + 0.25,
                                         max(ridge_plot_data$acc_diff) - 0.25),
                            text = c(num_neg, num_pos))
      return(temp_df)
    })
  )
  return(this_comp_df)
}))
te_counts$y_pos <- rep(
  seq(1, 3, length.out = length(cell_type_comps)),
  each = 2 * length(unique(peaks_by_tes_melted$TE))
)
te_ridge_plot_atac <- ggplot(ridge_plot_data) +
  aes(x = acc_diff, colour = Comparison) +
  scale_x_continuous(
    name = "Diff. in mean accessibility across samples"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    name = "Density of TE-overlapping peaks [Arbitrary units]"
  ) +
  scale_colour_manual(values = c("grey50", "turquoise1", "lightpink1")) +
  geom_density() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(data = te_counts, mapping = aes(label = text, y = y_pos)) +
  facet_wrap(~ TE_element, ncol = 1) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))
ridge_plot_atac_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_hervip10f_ltr45b_vio_plot_skin_bloodccr8_bloodnaive_violinplot_atac_related_ridgeplot.pdf"
)
ridge_plot_atac_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_hervip10f_ltr45b_vio_plot_skin_bloodccr8_bloodnaive_violinplot_atac_related_ridgeplot.rds"
)
pdf(ridge_plot_atac_pdf, width = 7, height = 4.5)
print(te_ridge_plot_atac)
dev.off()
saveRDS(te_ridge_plot_atac, file = ridge_plot_atac_rds)
