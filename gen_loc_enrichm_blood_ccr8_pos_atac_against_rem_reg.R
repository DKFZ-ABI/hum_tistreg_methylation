# This script analyses whether TEs are enriched in skin Treg hyperaccessibility 
#   peaks between skin Treg cells and blood naive Treg cells stratified by the 
#   positioning of blood CCR8+ Treg cells (with remaining regions as 
#   background).
# Author: Niklas Beumer



# Load required packages.
library(parallel)
library(ggplot2)
library(GenomicRanges)
library(ggrepel)
library(Signac)


# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the locations of DMRs between skin Tregs and blood naive Tregs with 
# the positioning of blood CCR8+ Tregs.
peaks_file <- paste0(
  location, 
  "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_acc_discrep_regions_info.txt"
)
peaks <- read.table(peaks_file, header = T, stringsAsFactors = F, sep = "\t")

# Restrict the peaks to the skin Treg hyperaccessibility class.
peaks <- peaks[peaks$avg_log2FC < 0, ]

# Generate a GRanges object for the DMRs.
peaks_gr <- StringToGRanges(peaks$Peak, starts.in.df.are.0based = T)
peaks_gr$Cell_type_closest_to_blood_ccr8_treg <-
  peaks$Cell_type_closest_to_blood_ccr8_treg
seqlevelsStyle(peaks_gr) <- "NCBI"

# Generate a meta data column that nicely displays positioning of blood CCR8+
# Tregs.
peaks_gr$ccr8pos <- paste0("Closer_to_", 
                           gsub(" ", "_", 
                                peaks_gr$Cell_type_closest_to_blood_ccr8_treg))

# Read in the locations of repeat elements from RepeatMasker
rmsk_file <- paste0(location, 
                    "/external_data/2023-09-14_rmsk_hg19_from_ucsc.txt.gz")
rmsk <- read.table(rmsk_file, header = F, stringsAsFactors = F)

# Filter the repeat element data so that only those elements remain where I am
# confident that they are transposons.
conf_transp_cats <- c("LINE", "DNA", "SINE", "LTR", "RC", "Retroposon")
rmsk_filt <- rmsk[rmsk$V12 %in% conf_transp_cats, ]

# Generate a GRanges object for the annotated TE insertion sites.
# Take into account that this table originated from UCSC. Intervals are thus 
# very likely 0-based and closed.
# Modify seqlevels styles so that they match with what I I used for the DMRs.
rmsk_filt_gr <- makeGRangesFromDataFrame(rmsk_filt,
                                         ignore.strand = T,
                                         seqnames.field = "V6",
                                         start.field = "V7",
                                         end.field = "V8",
                                         starts.in.df.are.0based = T,
                                         keep.extra.columns = T)
seqlevelsStyle(rmsk_filt_gr) <- "NCBI"

# Add a meta data column specifying that all elements are TEs.
rmsk_filt_gr$all_tes <- "All_TEs"







#############################################################################
# Investigate enrichments of TEs in the DMRs.
#############################################################################

# Print a status message.
print("Computing enrichments")

# Specify the meta data columns by which TEs will be grouped.
# For my analysis. I will group TEs into four different levels:
# 1. All TEs together (regardless of their class etc.)
# 2. Repeat classes (column "V12")
# 3. Repeat family (Column "V13")
# 4. Repeat name (Column "V11")
grouping_cols <- c("all_tes", "V12", "V13", "V11")
granularity_names <- c("all", "class", "family", "name")

# Separate the peaks into the two region classes that should be compared.
peaks_closer_skin_gr <- peaks_gr[peaks_gr$ccr8pos == 
                                   "Closer_to_Skin_Treg"]
peaks_closer_naive_gr <- peaks_gr[peaks_gr$ccr8pos == 
                                    "Closer_to_Blood_naive_Treg"]

# Iterate over the region classes.
enrichment_results_complete <- lapply(grouping_cols, FUN = function(j) {
    print(j)
  
  # Extract the column containing the relevant values for granularity.
  granularity_vals <- rmsk_filt_gr@elementMetadata[, j]
  
  # Identify overlaps between query intervals and subject intervals.
  peaks_closer_skin_overl <- findOverlaps(peaks_closer_skin_gr, rmsk_filt_gr, 
                                          ignore.strand = T)
  peaks_closer_naive_overl <- findOverlaps(peaks_closer_naive_gr, rmsk_filt_gr, 
                                           ignore.strand = T)
  
  peaks_closer_skin_overl_qhits <- from(peaks_closer_skin_overl)
  peaks_closer_skin_overl_subjhits <- to(peaks_closer_skin_overl)
  peaks_closer_naive_overl_qhits <- from(peaks_closer_naive_overl)
  peaks_closer_naive_overl_subjhits <- to(peaks_closer_naive_overl)
  
  # Iterate over all groups of subject intervals according to the chosen 
  # granularity.
  unique_groups <- unique(granularity_vals)
  # unique_groups <-
  #   unique_groups[1:min(50, length(unique_groups))] # For debugging purposes
  enrichment_dat <- do.call(rbind, mclapply(unique_groups, FUN = function(y) {
    
    # Identify corresponding subject loci.
    corresp_subj_ind <- which(granularity_vals == y)
    
    # Generate the confusion matrix for Fisher's exact test.
    num_cl_skin_overl <- length(unique(peaks_closer_skin_overl_qhits[
      peaks_closer_skin_overl_subjhits %in% corresp_subj_ind]))
    num_cl_skin_nonoverl <- length(peaks_closer_skin_gr) - num_cl_skin_overl
    num_cl_naive_overl <- length(unique(peaks_closer_naive_overl_qhits[
      peaks_closer_naive_overl_subjhits %in% corresp_subj_ind]))
    num_cl_naive_nonoverl <- length(peaks_closer_naive_gr) - num_cl_naive_overl
    matrix_for_fisher <- matrix(c(num_cl_skin_overl, num_cl_skin_nonoverl, 
                                  num_cl_naive_overl, num_cl_naive_nonoverl),
                                nrow = 2, byrow = T)
    
    # Perform a two-tailed Fisher's exact test.
    fisher_res <- fisher.test(matrix_for_fisher)
    
    # Return a data frame containing the results.
    temp_df <- data.frame(
      ID = y, 
      Num_closer_to_skin_overlapping_w_te = num_cl_skin_overl,
      Num_closer_to_skin_not_overlapping_w_te = num_cl_skin_nonoverl,
      Num_closer_to_bloodnaive_overlapping_w_te = num_cl_naive_overl,
      Num_closer_to_bloodnaive_not_overlapping_w_te = num_cl_naive_nonoverl,
      OR = fisher_res$estimate,
      Fisher_Pval = fisher_res$p.value
    )
    
    # Return the collected data.
    return(temp_df)
    
  }, mc.cores = 10))
  
  # Correct P values for multiple hypothesis testing using the Benjamini-
  # Hochberg method.
  enrichment_dat$Pval_adj_BH <- p.adjust(enrichment_dat$Fisher_Pval, 
                                         method = "BH")
  
  # Return the computed enrichment data.
  return(enrichment_dat)
  
})
names(enrichment_results_complete) <- granularity_names

# Save the computed enrichment data.
enrichment_results_complete_outfile <- paste0(
  location, 
  "/te_analysis/gen_loc_enrichm_blood_ccr8_pos_atac_skin_treg_hyperacc_against_rem_reg.rds"
)
saveRDS(enrichment_results_complete, file = enrichment_results_complete_outfile)









#############################################################################
# Visualise enrichment results.
#############################################################################

# Iterate over the different levels of enrichment results.
plot_list <- lapply(granularity_names, FUN = function(granularity) {

  # Extract corresponding enrichment data.
  plotting_data <- enrichment_results_complete[[granularity]]
  
  # Compute logarithms of odds ratio and adjusted P value.
  plotting_data$logenr <- log2(plotting_data$OR)
  plotting_data$logp <- (-1) * log10(plotting_data$Pval_adj_BH)
  
  # Clip P values.
  sanity_check <- length(which(plotting_data$logp < Inf)) > 0
  if (sanity_check) {
    max_noninf_p <- max(plotting_data$logp[plotting_data$logp < Inf])
    pval_replacement <- max(1, 1.25 * max_noninf_p)
  } else {
    max_noninf_p <- 0.75
    pval_replacement <- 1
  }
  plotting_data$logp[plotting_data$logp == Inf] <- pval_replacement
  
  # Clip enrichment statistics.
  sanity_check <- length(which(plotting_data$logenr < 0 & 
                                 !(is.na(plotting_data$logenr)) &
                                 plotting_data$logenr > -Inf)) > 0
  if (sanity_check) {
    min_negat_enr <- min(
      plotting_data$logenr[plotting_data$logenr < 0 & 
                             !(is.na(plotting_data$logenr)) &
                             plotting_data$logenr > -Inf]
    )
    neg_enr_replacement <- 1.25 * min_negat_enr
  } else {
    min_negat_enr <- 0
    neg_enr_replacement <- -0.5
  }
  sanity_check <- length(which(plotting_data$logenr > 0 & 
                                 !(is.na(plotting_data$logenr)) &
                                 plotting_data$logenr < Inf)) > 0
  if (sanity_check) {
    max_posit_enr <- max(
      plotting_data$logenr[plotting_data$logenr > 0 & 
                             !(is.na(plotting_data$logenr)) &
                             plotting_data$logenr < Inf]
    )
    pos_enr_replacement <- 1.25 * max_posit_enr
  } else {
    max_posit_enr <- 0
    pos_enr_replacement <- 0.5
  }
  plotting_data$logenr[
    !(is.na(plotting_data$logenr)) & plotting_data$logenr == -Inf] <- 
    neg_enr_replacement
  plotting_data$logenr[
    !(is.na(plotting_data$logenr)) & plotting_data$logenr == Inf] <- 
    pos_enr_replacement
  
  # Add a binary classification highlighting significant enrichments.
  plotting_data$signif <- plotting_data$Pval_adj_BH < 0.05
  
  # Determine step size on the y and y axes.
  x_extr <- max(-min_negat_enr, max_posit_enr)
  if (x_extr <= 0.2) {
    x_step <- 0.025
  } else if (x_extr <= 0.5) {
    x_step <- 0.1
  } else if (x_extr <= 1.5) {
    x_step <- 0.25
  } else if (x_extr <= 3) {
    x_step <- 0.5
  } else {
    x_step <- 1
  }
  if (max_noninf_p  <= 0.1) {
    y_step <- 0.025
  } else if (max_noninf_p <= 0.5) {
    y_step <- 0.1
  } else if (max_noninf_p <= 1) {
    y_step <- 0.25
  } else if (max_noninf_p <= 2) {
    y_step <- 0.5
  } else if (max_noninf_p <= 10) {
    y_step <- 1
  } else if (max_noninf_p <= 25) {
    y_step <- 5
  } else {
    y_step <- 10
  }
  
  # Prepare data to label the top-10 enriched TEs. If the name or family 
  # granularity is plotted, use a different approach: If there are more than 
  # 10 significantly enriched TEs, label the top-10 of them according to 
  # enrichment score. If there are fewer, only label the significant ones.
  # If the plot shows the "name" category of TEs, make sure that LTR45B
  # and HERVIP10F are labelled.
  if (granularity %in% c("family", "name")) {
    signif_table <- plotting_data[plotting_data$signif, ]
    top_inds <- order(abs(signif_table$OR), decreasing = T)[
      1:min(10, nrow(signif_table))]
    text_data <- signif_table[top_inds, ]
    if (granularity == "name") {
      additional_tes <- setdiff(c("LTR45B", "HERVIP10F-int"), text_data$ID)
      additional_df <- plotting_data[plotting_data$ID %in% additional_tes, ]
      text_data <- rbind(text_data, additional_df)
    }
  } else {
    top_inds <- order(plotting_data$OR, decreasing = T)[
      1:min(10, nrow(plotting_data))]
    text_data <- plotting_data[top_inds, ]
  }
  
  # Generate a plot showing the enrichments.
  x_breaks_inner <- unique(c(rev(-(seq(0, -min_negat_enr, x_step))),
                             seq(0, max_posit_enr, x_step)))
  x_breaks <- c(neg_enr_replacement, x_breaks_inner, pos_enr_replacement)
  plotting_data$signif <- factor(plotting_data$signif, levels = c(T, F))
  enrichment_plot <- ggplot(plotting_data) +
    aes(x = logenr, y = logp, colour = signif) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = c("-Inf", as.character(x_breaks_inner), "Inf"),
      name = "Enrichment [log2(Odds ratio)]\n") +
    scale_y_continuous(breaks = c(seq(0, max_noninf_p, y_step), 
                                  pval_replacement),
                       labels = as.character(
                         c(seq(0, max_noninf_p, y_step), 
                           "Inf")),
                       name = "-log10(q)",
                       expand = expansion(mult = c(0.02, 0.05))) +
    scale_colour_manual(breaks = c(T, F), values = c("red", "grey50"),
                        labels = c("q < 0.05", "q >= 0.05"),
                        name = "Significance",
                        drop = F) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_text_repel(data = text_data, mapping = aes(label = ID),
                    colour = "black", force = 100, min.segment.length = 0,
                    max.time = 10, max.iter = 1000000,
                    max.overlaps = 1000, seed = 2023, size = 3) +
    ggtitle(granularity) +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"))
  if (max(plotting_data$logp) == pval_replacement & nrow(plotting_data) > 1) {
    enrichment_plot <- enrichment_plot +
      geom_hline(yintercept = 1.2 * max_noninf_p, linetype = "dotted")
  }
  
  # Return the generated plot.
  return(enrichment_plot)
  
})

# Save the plots.
outfile_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_enrichm_blood_ccr8_pos_atac_skin_treg_hyperacc_against_rem_reg.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_enrichm_blood_ccr8_pos_atac_skin_treg_hyperacc_against_rem_reg.rds"
)
pdf(outfile_pdf, width = 7, height = 6)
print(plot_list)
dev.off()
saveRDS(plot_list, file = outfile_rds)
