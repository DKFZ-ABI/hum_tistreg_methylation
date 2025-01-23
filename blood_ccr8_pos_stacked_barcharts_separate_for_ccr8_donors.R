# This script generates stacked bar charts for blood CCR8+ Treg cell 
#   positionings, stratified by blood CCR8+ Treg cell donor.
# Author: Niklas Beumer



# Load required packages.
library(GenomicRanges)
library(bsseq)
library(Seurat)
library(Signac)
library(ggplot2)


# Define a location on /xxx.
b330_space <- "/xxx/"
location <- paste0(b330_space, "nbeumer/hm_treg_bs_rgnsbg")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/yyy/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Specify the relevant cell types.
relevant_cell_types <- c("Blood naive Treg", "Blood CCR8+ Treg", "Skin Treg")
relevant_cell_types_atac_naming <- c("blood_naive_treg", "blood_ccr8_treg", 
                                     "skin_treg")

# Read in the list of differentially methylated regions between skin Tregs and 
# blood naive Tregs.
meth_reg_file <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_gap_0.15_sterrratio_0.5_signature_regions_filtered.txt"
)
meth_reg <- read.table(meth_reg_file, header = T, stringsAsFactors = F, 
                       sep = "\t")
# meth_reg <- meth_reg[1:50000, ] # For debugging purposes.

# Read in the BSseq object containing the smoothed methylation data and restrict 
# to the relevant cell types.
meth_data <- readRDS(paste(
  location, 
  "/preprocessing_etc/quality_control_after_alignment/2022-03-14_bsseq_object_combined_all_samples_smoothed.rds", 
  sep = "/"
))
meth_data <- meth_data[, pData(meth_data)$Cell_type %in% relevant_cell_types]

# Read in the list of differential peaks between skin Tregs and blood naive 
# Tregs.
atac_reg_file <- paste0(
  location, 
  "/differential_accessibility/diff_acc_analysis_Skin_TregBlood_naive_Treg_results_with_significance.txt"
)
atac_reg <- read.table(atac_reg_file, header = T, stringsAsFactors = F)
atac_reg <- atac_reg[atac_reg$significant, ]

# Read in the Seurat object containing scATAC-seq data for CD4+ T cells.
sc_data_cd4_file <- paste0(
  b330_space, 
  "msimon/scATAC/final/seurat_human_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar.RDS"
)
sc_data_cd4 <- readRDS(sc_data_cd4_file)

# Subset scATAC-seq data to the relevant cell types and cells from the relevant 
# donors.
# Afterwards, perform a new TF-IDF.
blood_donors <- c("1", "2")
skin_donors <- c("4", "5")
sc_data_cd4_restr <- subset(
  sc_data_cd4, 
  subset = (treg_tconv_annot == "blood_naive_treg" & donor %in% blood_donors) |
    (treg_tconv_annot == "blood_ccr8_treg" & donor %in% blood_donors) |
    (treg_tconv_annot == "skin_treg" & donor %in% skin_donors)
)
sc_data_cd4_restr <- RunTFIDF(sc_data_cd4_restr)

# Read in the sample mapping file for RNAseq data (only biological replicates).
rna_sample_mapping_path <- paste0(location, 
                                  "/sample_mapping_rnaseq_only_biol_rep.txt")
rna_sample_mapping <- read.table(rna_sample_mapping_path, header = T, 
                                 stringsAsFactors = F, sep = "\t")

# Read in the list of differentially expressed genes between skin Tregs and 
# blood naive Tregs.
diff_exp_file <- paste0(
  location, 
  "/RNASeq/analysis_results/2022-01-14_diff_gene_expr_DESEq2_Skin_TregBlood_naive_Treg_results_filtered_with_significance.txt"
)
diff_exp <- read.table(diff_exp_file, header = T, stringsAsFactors = F)
diff_exp <- diff_exp[diff_exp$significant, ]

# Read in the TPM values.
tpm_file <- paste0(location, "/RNASeq/tpm/tpm_all.txt")
tpm <- read.table(tpm_file, header = T, stringsAsFactors = F)





#############################################################################
# For each blood CCR8+ Treg donor and each feature, compute positioning
# of blood CCR8+ Tregs with respect to blood naive Tregs and skin Tregs.
#############################################################################

###############################
# 1. Methylation level.

# Generate a GRanges object for the differentially methylated regions.
meth_reg_gr <- makeGRangesFromDataFrame(meth_reg, keep.extra.columns = T)

# Identify average raw methylation values (on the sample level).
avg_meth_sample_level <- getMeth(meth_data,
                                 regions = meth_reg_gr,
                                 type = "raw",
                                 what = "perRegion")

# Compute cell-type-level averages of these sample-level methylation values
# (only needed for skin Tregs and bloood naive Tregs).
avg_meth_celltype_level <- sapply(
  c("Blood naive Treg", "Skin Treg"), 
  FUN = function(x) {
    celltype_subs <- avg_meth_sample_level[
      , grep(x, colnames(avg_meth_sample_level), fixed = T)
    ]
    return(rowMeans(celltype_subs, na.rm = T))
})

# To make all blood CCR8+ Treg samples fully comparable in the final plot,
# only keep regions where all blood CCR8+ Treg samples have non-missing values.
blood_ccr8_samples <- colnames(meth_data)[
  pData(meth_data)$Cell_type == "Blood CCR8+ Treg"
]
keep_bools <- apply(
  avg_meth_sample_level[, blood_ccr8_samples], 
  1, 
  FUN = function(x) {
    all(!(is.na(x)) & !(is.nan(x)))
  }
)
avg_meth_sample_level <- avg_meth_sample_level[keep_bools, ]
avg_meth_celltype_level <- avg_meth_celltype_level[keep_bools, ]
meth_reg <- meth_reg[keep_bools, ]

# For blood CCR8+ Treg sample and each region, determine whether blood CCR8+ 
# Tregs are closer to blood naive Tregs or closer to skin Tregs.
blood_ccr8_pos_meth <- do.call(rbind, lapply(
  blood_ccr8_samples, 
  FUN = function(sample) {
    blood_ccr8_vals <- avg_meth_sample_level[, sample]
    ccr8_naive_dist <- abs(
      blood_ccr8_vals - avg_meth_celltype_level[, "Blood naive Treg"]
    )
    ccr8_skin_dist <- abs(
      blood_ccr8_vals - avg_meth_celltype_level[, "Skin Treg"]
    )
    ccr8_closer_to <- sapply(1:length(blood_ccr8_vals), FUN = function(x) {
      if (ccr8_naive_dist[x] < ccr8_skin_dist[x]) {
        return("Blood naive Treg")
      } else if (ccr8_naive_dist[x] > ccr8_skin_dist[x]) {
        return("Skin Treg")
      } else {
        return("Exactly in middle")
      }
    })
    this_sample_df <- data.frame(
      Omics_level = "Methylation",
      Donor_or_sample = gsub("\n", " ", sample),
      Level_donor_comb = paste0("Meth.: ", gsub("\n", " ", sample)),
      Feature = meth_reg$region_ID,
      Sig_cat = gsub("omethylation", "", meth_reg$automatic_annotation),
      Blood_ccr8_val = blood_ccr8_vals,
      Blood_naive_val = avg_meth_celltype_level[, "Blood naive Treg"],
      Skin_val = avg_meth_celltype_level[, "Skin Treg"],
      Blood_ccr8_blood_naive_dist = ccr8_naive_dist,
      Blood_ccr8_skin_dist = ccr8_skin_dist,
      Blood_ccr8_closer_to = ccr8_closer_to
    )
    return(this_sample_df)
  }
))



###############################
# 2. Chr. accessibility level.

# Compute cell-type-level average accessibility values for each region (only
# required for blood naive Tregs and skin Tregs).
avg_acc_celltype_level <- sapply(
  c("blood_naive_treg", "skin_treg"), 
  FUN = function(x) {
    acc_values <- GetAssayData(subset(sc_data_cd4_restr, 
                                      subset = treg_tconv_annot == x), 
                               assay = "scATAC_raw", 
                               slot = "data")[rownames(atac_reg), ]
  acc_means <- rowMeans(acc_values)
  return(acc_means)
})
colnames(avg_acc_celltype_level) <- c("Blood naive Treg", "Skin Treg")

# For each blood CCR8+ Treg donor, compute dono-wise average accessibility
# values for each region.
avg_acc_donor_level <- sapply(blood_donors, FUN = function(x) {
  acc_values <- GetAssayData(
    subset(sc_data_cd4_restr, 
           subset = treg_tconv_annot == "blood_ccr8_treg" & donor == x), 
    assay = "scATAC_raw", 
    slot = "data"
  )[rownames(atac_reg), ]
  acc_means <- rowMeans(acc_values)
  return(acc_means)
})
colnames(avg_acc_donor_level) <- paste0("Blood_CCR8_Tregs_Donor", blood_donors)

# For each blood CCR8+ Treg donor and each region, determine whether blood CCR8+ 
# Tregs are closer to blood naive Tregs or closer to skin Tregs.
blood_ccr8_pos_acc <- do.call(rbind, lapply(
  colnames(avg_acc_donor_level), 
  FUN = function(sample) {
    blood_ccr8_vals <- avg_acc_donor_level[, sample]
    ccr8_naive_dist <- abs(
      blood_ccr8_vals - avg_acc_celltype_level[, "Blood naive Treg"]
    )
    ccr8_skin_dist <- abs(
      blood_ccr8_vals - avg_acc_celltype_level[, "Skin Treg"]
    )
    ccr8_closer_to <- sapply(1:length(blood_ccr8_vals), FUN = function(x) {
      if (ccr8_naive_dist[x] < ccr8_skin_dist[x]) {
        return("Blood naive Treg")
      } else if (ccr8_naive_dist[x] > ccr8_skin_dist[x]) {
        return("Skin Treg")
      } else {
        return("Exactly in middle")
      }
    })
    this_sample_df <- data.frame(
      Omics_level = "Accessibility",
      Donor_or_sample = sample,
      Level_donor_comb = paste0("Acc.: ", sample),
      Feature = rownames(atac_reg),
      Sig_cat = sapply(atac_reg$avg_log2FC, FUN = function(y) {
        ifelse(y < 0, yes = "Skin_Treg__hyp", no = "Blood_naive_Treg__hyp")
      }),
      Blood_ccr8_val = blood_ccr8_vals,
      Blood_naive_val = avg_acc_celltype_level[, "Blood naive Treg"],
      Skin_val = avg_acc_celltype_level[, "Skin Treg"],
      Blood_ccr8_blood_naive_dist = ccr8_naive_dist,
      Blood_ccr8_skin_dist = ccr8_skin_dist,
      Blood_ccr8_closer_to = ccr8_closer_to
    )
    return(this_sample_df)
  }
))



###############################
# 3. Expression level.

# Restrict the TPM data to the differentially expressed genes.
tpm_diff_exp <- tpm[tpm$Gene_symbol %in% rownames(diff_exp), ]

# Compute cell-type-wise average TPM values for the relevant genes.
avg_tpm_celltype_level <- sapply(relevant_cell_types, FUN = function(x) {
  corresp_samples <- rna_sample_mapping$Sample[
    rna_sample_mapping$Cell_type == x
  ]
  corresp_tpm <- tpm_diff_exp[, corresp_samples]
  rownames(corresp_tpm) <- tpm_diff_exp$Gene_symbol
  for (colname in colnames(corresp_tpm)) {
    corresp_tpm[, colname] <- log1p(100 * corresp_tpm[, colname])
  }
  means <- rowMeans(corresp_tpm)
  return(means[rownames(diff_exp)])
})

# Also collect the sample-level TPM values.
tpm_sample_level <- do.call(
  cbind, 
  lapply(relevant_cell_types, FUN = function(x) {
    corresp_samples <- rna_sample_mapping$Sample[rna_sample_mapping$Cell_type == x]
    corresp_tpm <- tpm_diff_exp[, corresp_samples]
    rownames(corresp_tpm) <- tpm_diff_exp$Gene_symbol
    for (colname in colnames(corresp_tpm)) {
      corresp_tpm[, colname] <- log1p(100 * corresp_tpm[, colname])
    }
    return(corresp_tpm[rownames(diff_exp), ])
  })
)

# Perform gene exclusions as in the script that originally assessed blood CCR8+
# Treg positionings.
exclude_bools <- sapply(1:nrow(diff_exp), FUN = function(x) {
  diff_tend_sign <- sign(diff_exp$log2FoldChange[x])
  skin_log_tpm <- avg_tpm_celltype_level[x, "Skin Treg"]
  blood_naive_log_tpm <- avg_tpm_celltype_level[x, "Blood naive Treg"]
  return(ifelse(diff_tend_sign == -1, 
                yes = skin_log_tpm < blood_naive_log_tpm, 
                no = skin_log_tpm > blood_naive_log_tpm))
})
genes_to_exclude <- rownames(diff_exp)[exclude_bools]
avg_tpm_celltype_level_filtered <- avg_tpm_celltype_level[!exclude_bools, ]
tpm_sample_level_filtered <- tpm_sample_level[!exclude_bools, ]
diff_exp_filtered <- diff_exp[!exclude_bools, ]

# For each blood CCR8+ Treg sample and each region, determine whether blood 
# CCR8+ Tregs are closer to blood naive Tregs or closer to skin Tregs.
blood_ccr8_samples <- rna_sample_mapping$Sample[
  rna_sample_mapping$Cell_type == "Blood CCR8+ Treg"
]
blood_ccr8_pos_expr <- do.call(rbind, lapply(
  blood_ccr8_samples, 
  FUN = function(sample) {
    blood_ccr8_vals <- tpm_sample_level_filtered[, sample]
    ccr8_naive_dist <- abs(
      blood_ccr8_vals - avg_tpm_celltype_level_filtered[, "Blood naive Treg"]
    )
    ccr8_skin_dist <- abs(
      blood_ccr8_vals - avg_tpm_celltype_level_filtered[, "Skin Treg"]
    )
    ccr8_closer_to <- sapply(1:length(blood_ccr8_vals), FUN = function(x) {
      if (ccr8_naive_dist[x] < ccr8_skin_dist[x]) {
        return("Blood naive Treg")
      } else if (ccr8_naive_dist[x] > ccr8_skin_dist[x]) {
        return("Skin Treg")
      } else {
        return("Exactly in middle")
      }
    })
    this_sample_df <- data.frame(
      Omics_level = "Expression",
      Donor_or_sample = sample,
      Level_donor_comb = paste0("Expr.: ", sample),
      Feature = rownames(diff_exp_filtered),
      Sig_cat = sapply(diff_exp_filtered$log2FoldChange, FUN = function(y) {
        ifelse(y < 0, yes = "Skin_Treg__hyp", no = "Blood_naive_Treg__hyp")
      }),
      Blood_ccr8_val = blood_ccr8_vals,
      Blood_naive_val = avg_tpm_celltype_level_filtered[, "Blood naive Treg"],
      Skin_val = avg_tpm_celltype_level_filtered[, "Skin Treg"],
      Blood_ccr8_blood_naive_dist = ccr8_naive_dist,
      Blood_ccr8_skin_dist = ccr8_skin_dist,
      Blood_ccr8_closer_to = ccr8_closer_to
    )
    return(this_sample_df)
  }
))



###############################
# 4. Combination of the three omics levels.

# Combine the results on blood CCR8+ Treg positioning from the three omics 
# levels.
blood_ccr8_pos_combined <- rbind(blood_ccr8_pos_meth,
                                 blood_ccr8_pos_acc,
                                 blood_ccr8_pos_expr)

# Save the collected data.
blood_ccr8_pos_outfile <- paste0(
  location, 
  "/treg_hierarchies/blood_ccr8_treg_positioning_by_blood_ccr8_donor_all_omics_levels.txt"
)
write.table(blood_ccr8_pos_combined, file = blood_ccr8_pos_outfile, sep = "\t",
            col.names = T, row.names = F, quote = F)







#############################################################################
# Quantify proportions of blood CCR8+ Treg positionings on the three omics
# levels and in the different blood CCR8+ Treg samples/donors.
#############################################################################

# Quantify absolute numbers and proportions of features that place blood CCR8+
# Tregs closer to skin Tregs or closer to blood naive Tregs.
omics_cats <- unique(blood_ccr8_pos_combined$Omics_level)
positioning_cats <- c("Skin Treg", "Exactly in middle", "Blood naive Treg")
tendencies_cats <- c("Blood_naive_Treg__hyp", "Skin_Treg__hyp")
blood_ccr8_pos_props <- do.call(
  rbind, 
  lapply(tendencies_cats, FUN = function(x) {
    this_tend_df <- do.call(rbind, lapply(omics_cats, FUN = function(y) {
      corresp_samples <- unique(blood_ccr8_pos_combined$Donor_or_sample[
        blood_ccr8_pos_combined$Omics_level == y &
          blood_ccr8_pos_combined$Sig_cat == x
      ])
      this_omics_df <- do.call(
        rbind, 
        lapply(corresp_samples, FUN = function(z) {
          this_sample_total_feat <- nrow(blood_ccr8_pos_combined[
            blood_ccr8_pos_combined$Omics_level == y &
              blood_ccr8_pos_combined$Sig_cat == x &
              blood_ccr8_pos_combined$Donor_or_sample == z
            , ])
          this_sample_df <- do.call(rbind, lapply(
            positioning_cats, 
            FUN  = function(i) {
              total_count <- nrow(blood_ccr8_pos_combined[
                blood_ccr8_pos_combined$Omics_level == y &
                  blood_ccr8_pos_combined$Sig_cat == x &
                  blood_ccr8_pos_combined$Donor_or_sample == z &
                  blood_ccr8_pos_combined$Blood_ccr8_closer_to == i
                , ])
              relative_count <- total_count / this_sample_total_feat
              counts_df <- data.frame(
                Omics_level = y,
                Sample_or_Donor = z,
                Sig_cat = x,
                Sig_cat_w_sample_or_donor = paste0(y, ": ", x, ": ", z),
                Blood_ccr8_closer_to = i,
                Absolute_feature_count = total_count,
                Relative_feature_proportion = relative_count
              )
              return(counts_df)
            }))
          return(this_sample_df)
        })
      )
      return(this_omics_df)
    }))
    return(this_tend_df)
  })
)






###########################################################################
# Visualise the results.
###########################################################################

# Generate stacked bar charts showing the distribution of blood CCR8+ Treg
# positionings for the different samples/donors in the three omics levels.
blood_ccr8_pos_props$Sig_cat_w_sample_or_donor <- factor(
  blood_ccr8_pos_props$Sig_cat_w_sample_or_donor,
  levels = unique(blood_ccr8_pos_props$Sig_cat_w_sample_or_donor)
)
positioning_prop_plot <- ggplot(blood_ccr8_pos_props) +
  aes(x = Sig_cat_w_sample_or_donor, y = Relative_feature_proportion,
      fill = Blood_ccr8_closer_to) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("darkorchid1", "grey", "blue")) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, 
                                   vjust = 0.5),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))

# Save the plot.
positioning_prop_pdf <- paste0(plot_outdir, 
                               "/blood_ccr8_pos_props_by_donor.pdf")
positioning_prop_rds <- paste0(plot_rds_outdir, 
                               "/blood_ccr8_pos_props_by_donor.rds")
pdf(positioning_prop_pdf, width = 5, height = 7)
print(positioning_prop_plot)
dev.off()
saveRDS(positioning_prop_plot, file = positioning_prop_rds)

