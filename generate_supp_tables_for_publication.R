# This script generates the supplementary tables for the manuscript.
# Author: Niklas Beumer



# Load required package(s).
library(writexl)
library(GenomicRanges)
library(Signac)
library(parallel)
library(testit)


# Define a location on /xxx.
location <- "/xxx/nbeumer/hm_treg_bs_rgnsbg"

# Specify a location where output will be saved.
output_loc <- paste0(location, 
                     "/publication_2023_data_files/supplementary_tables")



########
# Supp. Table 1: Methylation-level cell type signatures.

# Read in the corresponding data file(s).
data_file <- paste0(
  location, 
  "/differential_methylation/multi_class_signature_extraction_gap_0.15_largestsecondlargest_1.5_sterrratio_0.5_signature_regions_filtered.txt"
)
data <- read.table(data_file, header = T, stringsAsFactors = F, sep = "\t")

# Optimise the way in which data is displayed (e.g. column names, remove
# unnecessary columns, ...).
colnames(data)[1] <- "chr"
colnames(data)[5] <- "signature_category"
data <- data[, !(colnames(data) %in% c(
  "region_ID",
  "nearest_gene", 
  "distance_to_nearest_gene", 
  "at_least_3_cpgs", 
  "retained", 
  "cov_at_least_2_all_samples", 
  "raw_meth_identifies_same_lower_meth_group", 
  "correctly_classified_samples_log_reg_at_least_14"
))]

# Save the supplementary table.
table_outfile <- paste0(
  output_loc, 
  "/Supp_Table_1_methylation_cell_type_signatures.xlsx"
)
write_xlsx(data, path = table_outfile)




########
# Supp. Table 2: Comparison between blood naive Tregs and blood naive Tconvs.
# DMRs, differential peaks and differential genes, plus DMR-peak-gene links.

# Read in the corresponding data file(s).
dmrs_file <- paste0(
  location, 
  "/differential_methylation/diff_meth_blood_naive_tconv_blood_naive_treg_gap_0.15_sterrratio_0.5_signature_regions_filtered.txt"
)
dmrs <- read.table(dmrs_file, header = T, stringsAsFactors = F, sep = "\t")
diffpeaks_file <- paste0(
  location, 
  "/differential_accessibility/diff_acc_analysis_Blood_naive_TregBlood_naive_Tconv_results_with_significance.txt"
)
diffpeaks <- read.table(diffpeaks_file, header = T, stringsAsFactors = F)
diffgenes_file <- paste0(
  location, 
  "/RNASeq/analysis_results/2022-01-14_diff_gene_expr_DESEq2_Blood_naive_TregBlood_naive_Tconv_results_filtered_with_significance.txt"
)
diffgenes <- read.table(diffgenes_file, header = T, stringsAsFactors = F)
links_file <- paste0(
  location, 
  "/differential_methylation/diff_meth_diff_acc_diff_gex_naivetconv_naivetreg_dmr_peak_gene_links_meth_acc_gex_corr.txt"
)
links <- read.table(links_file, header = T, stringsAsFactors = F, sep = "\t")

# Optimise the way in which data is displayed (e.g. column names, remove
# unnecessary columns, ...).
colnames(dmrs)[1] <- "chr"
colnames(dmrs)[5] <- "differentiality_direction"
dmrs <- dmrs[, !(colnames(dmrs) %in% c(
  "region_ID",
  "nearest_gene", 
  "distance_to_nearest_gene", 
  "at_least_3_cpgs", 
  "retained", 
  "cov_at_least_2_all_samples", 
  "raw_meth_identifies_same_lower_meth_group", 
  "correctly_classified_samples_log_reg_at_least_6"
))]
diffpeaks$peak_interval <- rownames(diffpeaks)
diffpeaks$more_accessible_in <- sapply(diffpeaks$avg_log2FC, FUN = function(x) {
  if(x < 0) {
    return("Blood naive Treg")
  } else if (x > 0) {
    return("Blood naive Tconv")
  } else {
    return("none")
  }
})
diffgenes$gene_name <- rownames(diffgenes)
diffgenes$more_expressed_in <- 
  sapply(diffgenes$log2FoldChange, FUN = function(x) {
    if(x < 0) {
      return("Blood naive Treg")
    } else if (x > 0) {
      return("Blood naive Tconv")
    } else {
      return("none")
    }
  }
)
colnames(links)[1] <- "Methylation__chr"
colnames(links)[7] <- "Methylation__differentiality_direction"
colnames(links)[20] <- "Methylation__Avg_raw_meth_diff_baseline_tconv"
colnames(links)[28] <- "Accessibility__peak_interval"
colnames(links)[30] <- "Accessibility__gene_assignment_overlap_types"
colnames(links)[31] <- "Accessibility__avg_log2FC_baseline_tconv"
colnames(links)[35] <- "Expression__log2FoldChange_baseline_tconv"
links <- links[, !(colnames(links) %in% c(
  "Methylation__strand",
  "Methylation__region_ID", 
  "Methylation__nearest_gene", 
  "Methylation__distance_to_nearest_gene", 
  "Methylation__at_least_3_cpgs", 
  "Methylation__retained", 
  "Methylation__cov_at_least_2_all_samples", 
  "Methylation__raw_meth_identifies_same_lower_meth_group",
  "Methylation__correctly_classified_samples_log_reg_at_least_6"
))]

# Save the supplementary table.
table_outfile <- paste0(
  output_loc, 
  "/Supp_Table_2_diff_analysis_treg_tconv_w_links.xlsx"
)
list_to_save <- list(DMR = dmrs,
                     Diff_peaks = diffpeaks,
                     Diff_genes = diffgenes,
                     DMR_peak_gene_links = links)
write_xlsx(list_to_save, path = table_outfile)




########
# Supp. Table 3: Comparison between blood naive Tregs and skin Tregs.
# DMRs, differential peaks and differential genes, plus DMR-peak-gene links.

# Read in the corresponding data file(s).
dmrs_file <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_gap_0.15_sterrratio_0.5_signature_regions_filtered.txt"
)
dmrs <- read.table(dmrs_file, header = T, stringsAsFactors = F, sep = "\t")
diffpeaks_file <- paste0(
  location, 
  "/differential_accessibility/diff_acc_analysis_Skin_TregBlood_naive_Treg_results_with_significance.txt"
)
diffpeaks <- read.table(diffpeaks_file, header = T, stringsAsFactors = F)
diffgenes_file <- paste0(
  location, 
  "/RNASeq/analysis_results/2022-01-14_diff_gene_expr_DESEq2_Skin_TregBlood_naive_Treg_results_filtered_with_significance.txt"
)
diffgenes <- read.table(diffgenes_file, header = T, stringsAsFactors = F)
links_file <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_diff_acc_diff_gex_dmr_peak_gene_links_meth_acc_gex_corr.txt"
)
links <- read.table(links_file, header = T, stringsAsFactors = F, sep = "\t")

# Optimise the way in which data is displayed (e.g. column names, remove
# unnecessary columns, ...).
colnames(dmrs)[1] <- "chr"
colnames(dmrs)[5] <- "differentiality_direction"
dmrs <- dmrs[, !(colnames(dmrs) %in% c(
  "region_ID",
  "nearest_gene", 
  "distance_to_nearest_gene", 
  "at_least_3_cpgs", 
  "retained", 
  "cov_at_least_2_all_samples", 
  "raw_meth_identifies_same_lower_meth_group", 
  "correctly_classified_samples_log_reg_at_least_6"
))]
diffpeaks$peak_interval <- rownames(diffpeaks)
diffpeaks$more_accessible_in <- sapply(diffpeaks$avg_log2FC, FUN = function(x) {
  if(x < 0) {
    return("Skin Treg")
  } else if (x > 0) {
    return("Blood naive Treg")
  } else {
    return("none")
  }
})
diffgenes$gene_name <- rownames(diffgenes)
diffgenes$more_expressed_in <- 
  sapply(diffgenes$log2FoldChange, FUN = function(x) {
    if(x < 0) {
      return("Skin Treg")
    } else if (x > 0) {
      return("Blood naive Treg")
    } else {
      return("none")
    }
  }
  )
colnames(links)[1] <- "Methylation__chr"
colnames(links)[7] <- "Methylation__differentiality_direction"
colnames(links)[20] <- "Methylation__Avg_raw_meth_diff_baseline_blood_naive"
colnames(links)[28] <- "Accessibility__peak_interval"
colnames(links)[30] <- "Accessibility__gene_assignment_overlap_types" 
colnames(links)[31] <- "Accessibility__avg_log2FC_baseline_blood_naive"
colnames(links)[35] <- "Expression__log2FoldChange_baseline_blood_naive"
links <- links[, !(colnames(links) %in% c(
  "Methylation__strand",
  "Methylation__region_ID", 
  "Methylation__nearest_gene", 
  "Methylation__distance_to_nearest_gene", 
  "Methylation__at_least_3_cpgs", 
  "Methylation__retained", 
  "Methylation__cov_at_least_2_all_samples", 
  "Methylation__raw_meth_identifies_same_lower_meth_group",
  "Methylation__correctly_classified_samples_log_reg_at_least_6"
))]

# Save the supplementary table.
table_outfile <- paste0(
  output_loc, 
  "/Supp_Table_3_diff_analysis_skintreg_bloodnaivetreg_w_links.xlsx"
)
list_to_save <- list(DMR = dmrs,
                     Diff_peaks = diffpeaks,
                     Diff_genes = diffgenes,
                     DMR_peak_gene_links = links)
write_xlsx(list_to_save, path = table_outfile)




########
# Supp. Table 5: All homer results for the comparison between skin Tregs and
# blood naive Tregs.

# Read in the corresponding data files.
homer_snips <- c(
  "diff_meth_skin_treg_blood_naive_treg_Blood_naive_Treg__hypomethylation",
  "diff_meth_skin_treg_blood_naive_treg_Skin_Treg__hypomethylation",
  "diff_acc_skin_treg_blood_naive_treg_Blood_naive_Treg__hyperaccessibility",
  "diff_acc_skin_treg_blood_naive_treg_Skin_Treg__hyperaccessibility"
)
homer_results <- lapply(homer_snips, FUN = function(x) {
  homer_file <- paste0(location, "/treg_hierarchies/", x, 
                       "_homer_against_genome/knownResults.txt")
  read.table(homer_file, sep = "\t", header = T, stringsAsFactors = F, 
             comment.char = "")
})
names(homer_results) <- c("Blood_naive_Treg__hypometh", 
                          "Skin_Treg__hypometh",
                          "Blood_naive_Treg__hyperacc", 
                          "Skin_Treg__hyperacc")

# Optimise column names.
homer_results <- lapply(homer_results, FUN = function(x) {
  colnames(x) <- sapply(colnames(x), FUN = function(y) {
    colname_temp <- gsub(".", 
                         "_", 
                         gsub("X.", 
                              "Prop", 
                              gsub("X..of.Background.Sequences.with.Motif.of",
                                   "Num_Background_Sequences_with_Motif_of", 
                                   gsub("X..of.Target.Sequences.with.Motif.of",
                                        "Num_Target_Sequences_with_Motif_of",
                                        y)),
                              fixed = T), 
                         fixed = T)
    charnum <- nchar(colname_temp)
    if (substr(colname_temp, charnum, charnum) == "_") {
      colname_temp <- substr(colname_temp, 1, charnum - 1)
    }
    return(colname_temp)
  })
  return(x)
})

# Compute enrichment scores.
homer_results <- lapply(homer_results, FUN = function(x) {
  target_perc <- as.numeric(
    gsub("%", "", x$Prop_of_Target_Sequences_with_Motif)
  )
  background_perc <- as.numeric(
    gsub("%", "", x$Prop_of_Background_Sequences_with_Motif)
  )
  x$Enrichment_score <- target_perc / background_perc
  return(x)
})

# Save the tables.
table_outfile <- paste0(
  output_loc, 
  "/Supp_Table_5_motif_enrichment_results_skin_naive.xlsx"
)
write_xlsx(homer_results, path = table_outfile)





########
# Supp. Table 6: Blood CCR8+ Treg positionings on all three levels.

# Read in the corresponding data file(s).
methdata_file <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_meth_discrep_regions_info.txt"
)
methdata <- read.table(methdata_file, header = T, stringsAsFactors = F, 
                       sep = "\t")
accdata_file <- paste0(
  location, 
  "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_acc_discrep_regions_info.txt"
)
accdata <- read.table(accdata_file, header = T, stringsAsFactors = F, 
                      sep = "\t")
gexdata_file <- paste0(
  location, 
  "/treg_hierarchies/diff_expr_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_expr_discrep_genes_info.txt"
)
gexdata <- read.table(gexdata_file, header = T, stringsAsFactors = F, 
                      sep = "\t")

# Optimise the way in which data is displayed (e.g. column names, remove
# unnecessary columns, ...).
colnames(methdata)[1] <- "chr"
colnames(methdata)[5] <- "differentiality_direction"
methdata <- methdata[, !(colnames(methdata) %in% c(
  "region_ID",
  "nearest_gene", 
  "distance_to_nearest_gene", 
  "at_least_3_cpgs", 
  "retained", 
  "cov_at_least_2_all_samples", 
  "raw_meth_identifies_same_lower_meth_group", 
  "correctly_classified_samples_log_reg_at_least_6"
))]
accdata$more_accessible_in <- sapply(accdata$avg_log2FC, FUN = function(x) {
  if(x < 0) {
    return("Skin Treg")
  } else if (x > 0) {
    return("Blood naive Treg")
  } else {
    return("none")
  }
})
accdata <- accdata[, c(1:7, 22, 8:21)]
gexdata$more_expressed_in <- 
  sapply(gexdata$log2FoldChange, FUN = function(x) {
    if(x < 0) {
      return("Skin Treg")
    } else if (x > 0) {
      return("Blood naive Treg")
    } else {
      return("none")
    }
  }
  )
gexdata <- gexdata[, c(1:3, 20, 4:19)]

# Save the supplementary table.
table_outfile <- paste0(
  output_loc, 
  "/Supp_Table_6_blood_CCR8_Treg_positioning_all_levels.xlsx"
)
list_to_save <- list(Methylation = methdata,
                     Accessibility = accdata,
                     Expression = gexdata)
write_xlsx(list_to_save, path = table_outfile)




########
# Supp. Table 7: Table of all TE subfamilies assigned to classes together with 
# the numbers of insertion sites.

# Read in the locations of repeat elements from RepeatMasker
rmsk_file <- paste0(location, 
                    "/external_data/2023-09-14_rmsk_hg19_from_ucsc.txt.gz")
rmsk <- read.table(rmsk_file, header = F, stringsAsFactors = F)

# Filter the repeat element data so that only those elements remain where I am
# confident that they are transposons.
conf_transp_cats <- c("LINE", "DNA", "SINE", "LTR", "RC", "Retroposon")
rmsk_filt <- rmsk[rmsk$V12 %in% conf_transp_cats, ]
# rmsk_filt <- rmsk_filt[1:100000, ] # For debugging purposes.

# Identify all TE subfamilies.
all_subfs <- unique(rmsk_filt$V11)

# Iterate over all TE subfamilies; collect the corresponding class names and 
# the number of insertion sites.
subfamilies <- rmsk_filt$V11
classes <- rmsk_filt$V12
all_data_collected <- do.call(rbind, lapply(all_subfs, FUN = function(x) {
  this_subf_inds <- which(subfamilies == x)
  ins_count <- length(this_subf_inds)
  class <- classes[this_subf_inds]
  class_unique <- unique(class)
  assert(length(class_unique) == 1)
  temp_df <- data.frame(Subfamily = x,
                        Class = class_unique,
                        Insertion_site_count = ins_count)
  return(temp_df)
}))


# Save the supplementary table.
table_outfile <- paste0(
  output_loc, 
  "/Supp_Table_7_TE_subfamilies_w_classes_and_insertion_counts.xlsx"
)
write_xlsx(all_data_collected, path = table_outfile)




########
# Supp. Table 8: All TE insertion sites with indication of whether they overlap
# with a DMR or differential peak.

# Generate a GRanges object for the annotated TE insertion sites.
# Take into account that this table originated from UCSC. Intervals are thus 
# very likely 0-based and closed.
# Generate a version of this GRanges object that uses the same seqlevels style
# that I also used for DMRs.
rmsk_filt_gr <- makeGRangesFromDataFrame(rmsk_filt,
                                         ignore.strand = T,
                                         seqnames.field = "V6",
                                         start.field = "V7",
                                         end.field = "V8",
                                         starts.in.df.are.0based = T,
                                         keep.extra.columns = T)
rmsk_filt_ncbi_gr <- rmsk_filt_gr
seqlevelsStyle(rmsk_filt_ncbi_gr) <- "NCBI"

# Prepare a preliminary table showing TE insertion sites in a readable format.
final_table <- data.frame(
  Chr = rmsk_filt$V6,
  Start = rmsk_filt$V7,
  End = rmsk_filt$V8,
  TE_class = rmsk_filt$V12,
  TE_subfamily = rmsk_filt$V11,
  Overlapping_skintreg_bloodnaive_DMRs = "none",
  Overlapping_skintreg_bloodnaivetreg_differentialpeaks = "none"
)

# Read in the locations of DMRs between skin Tregs and blood naive Tregs.
dmrs_file <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_gap_0.15_sterrratio_0.5_signature_regions_filtered.txt"
)
dmrs <- read.table(dmrs_file, header = T, stringsAsFactors = F, sep = "\t")

# Generate a GRanges object for the DMRs.
dmrs_gr <- makeGRangesFromDataFrame(dmrs, keep.extra.columns = T)

# Read in the table with peaks that are differentially accessible between skin 
# Tregs and blood naive Tregs.
atac_sig_file <- paste0(
  location, 
  "/differential_accessibility/diff_acc_analysis_Skin_TregBlood_naive_Treg_results_with_significance.txt"
)
atac_sig <- read.table(atac_sig_file, header = T, stringsAsFactors = F, 
                       sep = "\t")

# Restrict the peaks to those that displayed statistical significance.
atac_sig <- atac_sig[atac_sig$significant, ]

# Generate a column specifying the differential tendency.
atac_sig$Signature_category <- sapply(atac_sig$avg_log2FC, FUN = function(x) {
  ifelse(x < 0, yes = "Skin_Treg__hyperaccessibility", 
         no = "Blood_naive_Treg__hyperaccessibility")
})

# Generate a GRanges object for the differential peaks.
sig_gr <- StringToGRanges(rownames(atac_sig),
                          starts.in.df.are.0based=FALSE)
sig_gr$Peak <- rownames(atac_sig)
sig_gr$Signature_category <- atac_sig$Signature_category
seqlevelsStyle(sig_gr) <- "NCBI"

# For each TE insertion site, add information on what DMRs and what differential
# peak it overlaps with.
te_dmr_overl <- findOverlaps(rmsk_filt_ncbi_gr, dmrs_gr, ignore.strand = T)
te_dmr_overl_qhits <- from(te_dmr_overl)
te_dmr_overl_shits <- to(te_dmr_overl)
dmr_chrs <- dmrs$seqnames
dmr_starts <- dmrs$start
dmr_ends <- dmrs$end
dmrs_difftends <- dmrs$automatic_annotation
final_table$Overlapping_skintreg_bloodnaive_DMRs <- 
  unlist(mclapply(1:nrow(final_table), FUN = function(x) {
    if (x %in% te_dmr_overl_qhits) {
      corresp_dmr_inds <- te_dmr_overl_shits[which(te_dmr_overl_qhits == x)]
      corresp_dmr_str <- paste(sapply(corresp_dmr_inds, FUN = function(y) {
        paste0(dmr_chrs[y], "_", dmr_starts[y], "_", dmr_ends[y], 
               " (", dmrs_difftends[y], ")")
      }), collapse = "; ")
      return(corresp_dmr_str)
    } else {
      return("none")
    }
  }, mc.cores = 10))
te_diffpeak_overl <- findOverlaps(rmsk_filt_ncbi_gr, sig_gr, ignore.strand = T)
te_diffpeak_overl_qhits <- from(te_diffpeak_overl)
te_diffpeak_overl_shits <- to(te_diffpeak_overl)
peak_intervals <- rownames(atac_sig)
peak_difftends <- atac_sig$Signature_category
final_table$Overlapping_skintreg_bloodnaivetreg_differentialpeaks <- 
  unlist(mclapply(1:nrow(final_table), FUN = function(x) {
    if (x %in% te_diffpeak_overl_qhits) {
      corresp_diffpeak_inds <- te_diffpeak_overl_shits[
        which(te_diffpeak_overl_qhits == x)
      ]
      corresp_diffpeak_str <- paste(
        sapply(corresp_diffpeak_inds, FUN = function(y) {
          paste0(peak_intervals[y], " (", peak_difftends[y], ")")
        }), collapse = "; ")
      return(corresp_diffpeak_str)
    } else {
      return("none")
    }
  }, mc.cores = 10))


# Save the supplementary table.
table_outfile <- paste0(
  output_loc, 
  "/Supp_Table_8_TE_sites_and_DMR_diffpeak_overlaps.csv"
)
write.csv(final_table, file = table_outfile, row.names = F, quote = F)




########
# Supp. Table 9: Expression levels of TEs in all analysed cell types 
# (sample level).

# Read in the corresponding data file(s).
tpm_file <- paste0(
  location, 
  "/te_analysis/cons_seq_tespex_results_bulk_rna_all_tpm_values.txt"
)
tpm <- read.table(tpm_file, header = T, stringsAsFactors = F)

# Optimise the way in which data is displayed (e.g. column names, remove
# unnecessary columns, ...).
colnames(tpm) <- sapply(colnames(tpm), FUN = function(x) {
  x_new <- gsub("p024_RNAseq_", "", x)
  x_new <- gsub("CCR8_Treg", "Blood_CCR8_Treg", x_new)
  x_new <- gsub("CCRR_Treg", "Blood_CCR8_Treg", x_new)
  x_new <- gsub("Naive_Treg", "Blood_Naive_Treg", x_new)
  x_new <- gsub("Naive_Tconv", "Blood_Naive_Tconv", x_new)
  x_new <- strsplit(x_new, split = "_S")[[1]][1]
})
colnames(tpm) <- paste0("TPM_", colnames(tpm))
tpm$TE_subfamily <- rownames(tpm)
tpm <- tpm[, c(22, 1:21)]

# Save the supplementary table.
table_outfile <- paste0(
  output_loc, 
  "/Supp_Table_9_TE_expression_TPM_by_sample.xlsx"
)
write_xlsx(tpm, path = table_outfile)




########
# Supp. Table 10: HERVIP10F-int and LTR45B insertion sites overlapping with 
# skin Treg hypomethylation DMRs.

# Read in the corresponding data file(s).
hervip_file <- paste0(
  location, 
  "/te_analysis/cons_seq_insertion_sites_HERVIP10F-int_overl_w_skin_treg_hypometh_dmrs.txt"
)
hervip <- read.table(hervip_file, header = T, stringsAsFactors = F, sep = "\t")
ltr45b_file <- paste0(
  location, 
  "/te_analysis/cons_seq_insertion_sites_LTR45B_overl_w_skin_treg_hypometh_dmrs.txt"
)
ltr45b <- read.table(ltr45b_file, header = T, stringsAsFactors = F, sep = "\t")

# Optimise the way in which data is displayed (e.g. column names, remove
# unnecessary columns, ...).
hervip <- hervip[, c(1:3, 6:8)]
ltr45b <- ltr45b[, c(1:3, 6:8)]

# Save the supplementary table.
list_to_save <- list(`HERVIP10F-int` = hervip,
                     LRT45B = ltr45b)
table_outfile <- paste0(
  output_loc, 
  "/Supp_Table_10_HERVIP10Fint_LTR45B_insertion_sites_overlapping_w_skin_treg_hypometh_DMRs.xlsx"
)
write_xlsx(list_to_save, path = table_outfile)

