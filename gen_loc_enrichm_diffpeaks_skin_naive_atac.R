# This script analyses enrichments of TE genomic loci (RepeatMasker) in 
#   differential peaks (ATAC level) between skin Treg cells and blood naive Treg 
#   cells.
# Author: Niklas Beumer
# GenomeInfoDb and BSgenome.Hsapiens.UCSC.hg19 must be installed but do 
#   not need to be loaded. bedtools/2.24.0 must be installed.



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
# Modify seqlevels styles.
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
# Generate randomly shuffled versions of the peaks.
#############################################################################

shuffle_genomic_regions_hg19 <- function(regions_gr, chrs_to_use, 
                                         temp_file_location, 
                                         seed_to_use,
                                         shuffled_regions_outfile,
                                         keep_shuffled_bed_on_disk = F) {
  # This function randomly shuffles genomic regions on the hg19 genome using 
  #   bedtools, meaning it returns randomly distributed, non-overlapping regions 
  #   of the same lengths as in the input.
  # regions_gr: GRanges object containing the regions to shuffle.
  # chrs_to_use: Character; The names of the chromosomes (without "chr") to use.
  # temp_file_location: Character; Path to a directory where temporary files 
  #   will be saved. Note, this directory should not exist yet because it will
  #   be deleted afterwards.
  # seed_to_use: Numeric; The random seed to pass on to bedtools.
  # shuffled_regions_outfile: Character; Path to a BED file where the shuffled 
  #   regions will be saved. This may be temporary (see 
  #   keep_shuffled_bed_on_disk).
  # keep_shuffled_bed_on_disk: Logical; Whether to keep the BED file with 
  #   shuffled intervals on disk, if FALSE, the function will delete 
  #   shuffled_regions_outfile.
  # Dependencies: GenomicRanges. GenomeInfoDb and 
  #   BSgenome.Hsapiens.UCSC.hg19 must be installed but do not need to be 
  #   loaded. bedtools/2.29.2 must be installed.
  # Value: A list. The first element is a GRanges object with shuffled genomic 
  #   regions. The second element is a character vector contaning log messages
  #   from bedtools shuffle. Furthermore, the function generates a text file 
  #   containing the random regions (if keep_shuffled_bed_on_disk = T).
  # Author: Niklas Beumer
  
  # Get lengths of the chromosomes in the reference genome.
  # Remove "chr" from chromosome names.
  chr_lengths_all <- GenomeInfoDb::seqlengths(
    BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
  names(chr_lengths_all) <- gsub("chr", "", names(chr_lengths_all))
  chr_lengths <- chr_lengths_all[chrs_to_use]
  
  # Generate the temporary directory.
  dir.create(temp_file_location)
  
  # Save a temporary genome file containing chromosome lengths.
  genome_df <- data.frame(chr = names(chr_lengths),
                          length = chr_lengths)
  genome_tempfile <- paste0(temp_file_location, "/genome_temp.txt")
  write.table(genome_df, file = genome_tempfile, sep = "\t", row.names = F, 
              col.names = F, quote = F)
  
  # Generate and save a temporary BED file containing the regions to shuffle.
  # Use 0-based, half-open, intervals as is the standard for BED files.
  bed_df <- as.data.frame(regions_gr)[, 1:3]
  bed_df$start <- bed_df$start - 1
  bed_tempfile <- paste0(temp_file_location, "/bed_temp.bed")
  write.table(bed_df, file = bed_tempfile, sep = "\t", row.names = F, 
              col.names = F, quote = F)
  
  # Generate the command string for random region generation.
  command_string <- paste0(
    #-- Load bedtools.
    "module load bedtools/2.24.0; ",
    #-- bedtools shuffle command.
    "bedtools shuffle ",
    "-seed ",
    seed_to_use,
    " -noOverlapping ",
    "-i ",
    bed_tempfile,
    " -g ",
    genome_tempfile,
    #-- Re-direct sterr to stout.
    " 2>&1",
    #-- Direct output to the specified file.
    " > ",
    shuffled_regions_outfile,
    #-- Delete temporary input files.
    "; rm ",
    bed_tempfile,
    "; rm ",
    genome_tempfile)
  
  # Run the command for random region generation.
  log_messages <- system(command_string, intern = T)
  
  # Read in the random regions and turn them into a GRanges object with 1-based, 
  # closed, intervals.
  random_regions <- read.table(shuffled_regions_outfile, stringsAsFactors = F)
  random_regions_gr <- makeGRangesFromDataFrame(random_regions, 
                                                seqnames.field = "V1", 
                                                start.field = "V2", 
                                                end.field = "V3",
                                                starts.in.df.are.0based = T)
  
  # If the user did not request to keep the shuffled BED file on disk, delete 
  # it.
  if (!(keep_shuffled_bed_on_disk)) {
    file.remove(shuffled_regions_outfile)
  }
  
  # Remove the temporary directory.
  unlink(temp_file_location, recursive = T)
  
  # Return the GRanges object with the shuffled regions together with the 
  # log messages.
  return(list(random_regions_gr, log_messages))
  
}

# Print a status message.
print("Shuffling genomic intervals.")

# Specify how many shufflings are desired.
num_shufflings <- 10000
# num_shufflings <- 10 # For debugging purposes

# Iterate over the two region classes.
region_classes <- unique(sig_gr$Signature_category)
shuffled_intervals <- lapply(region_classes, FUN = function(x){
  print(x)
  
  # Extract corresponding genomic intervals.
  intervals_gr <- sig_gr[sig_gr$Signature_category == x]
  # intervals_gr <- intervals_gr[1:100] # For debugging purposes.
  
  # Generate random shufflings of the corresponding intervals.
  # Note: In a preliminary analysis, I found out that I can use the
  # default random number generator despite multi-threading, likely because
  # the seed is set each time on each thread.
  seeds <- 1:num_shufflings
  shufflings_all <- mclapply(seeds, FUN = function(y) {
    shuffle_genomic_regions_hg19(
      regions_gr = intervals_gr, 
      chrs_to_use = c(as.character(1:22), "X"), 
      temp_file_location = paste0(location, 
                                  "/te_analysis/tempdir_skinnaive_atac_", 
                                  y), 
      seed_to_use = y, 
      shuffled_regions_outfile = paste0(
        location, "/te_analysis/bed_skinnaive_atac_", y, ".bed")
    )
  }, mc.cores = 10)
  shufflings <- mclapply(shufflings_all, FUN = function(y) {
    y[[1]]
  }, mc.cores = 10)
  logs <- unlist(mclapply(shufflings_all, FUN = function(y) {
    y[[2]]
  }, mc.cores = 10))
  print(logs)
  
  # Return the results from random shuffling.
  return(shufflings)
  
})
names(shuffled_intervals) <- region_classes










#############################################################################
# Investigate enrichments of TEs in the peaks.
#############################################################################

compute_overl_enrichment_wrt_permutations <- function(
    query_regions_gr, subject_regions_gr, subject_regions_granularity_col, 
    shuffled_regions_gr, ncores) {
  # This function computes enrichments and enrichment P values for overlaps
  #   between query regions with subject regions based on a permutation test
  #   using previously shuffled query regions.
  # query_regions_gr: GRanges; The regions to assess enrichment for.
  # subject_regions_gr: GRanges; The regions whose enrichment should be 
  #   assessed.
  # subject_regions_granularity_col: Character; Name of a meta data column in
  #   subject_regions_gr that contains the grouping of subject regions 
  #   according to which the enrichment will be computed.
  # shuffled_regions_gr: List of GRanges objects containing shuffled query
  #   intervals.
  # ncores: Numeric; The number of cores to use for parallel computing.
  # Value: Data frame with results from enrichment analysis
  # Dependencies: GenomicRanges, parallel
  # Author: Niklas Beumer
  
  # Extract the column containing the relevant values for granularity.
  granularity_vals <- subject_regions_gr@elementMetadata[
    , subject_regions_granularity_col]
  
  # Identify overlaps between query intervals and subject intervals.
  query_subj_overl <- findOverlaps(query_regions_gr, subject_regions_gr, 
                                   ignore.strand = T)
  
  query_subj_overl_qhits <- from(query_subj_overl)
  query_subj_overl_shits <- to(query_subj_overl)
  
  # Identify overlaps between shuffled query intervals and subject intervals.
  shuf_subj_overl <- mclapply(1:length(shuffled_regions_gr), FUN = function(x) {
    findOverlaps(shuffled_regions_gr[[x]], subject_regions_gr, 
                 ignore.strand = T)
  }, mc.cores = ncores)
  shuf_subj_overl_qhits <- lapply(shuf_subj_overl, FUN = from)
  shuf_subj_overl_shits <- lapply(shuf_subj_overl, FUN = to)
  
  # Get the region number in the query set and for the shuffled query sets.
  reg_num_query <- length(query_regions_gr)
  reg_num_shuf <- sapply(shuffled_regions_gr, FUN = length)
  
  # Iterate over all groups of subject intervals according to the chosen 
  # granularity.
  unique_groups <- unique(granularity_vals)
  # unique_groups <-
  #   unique_groups[1:min(50, length(unique_groups))] # For debugging purposes
  enrichment_dat <- do.call(rbind, lapply(unique_groups, FUN = function(y) {
    
    # Identify corresponding subject loci.
    corresp_subj_ind <- which(granularity_vals == y)
    
    # Compute the proportion of corresponding query intervals overlapping with 
    # corresponding subject intervals.
    prop_overl_query <- length(unique(query_subj_overl_qhits[
      query_subj_overl_shits %in% corresp_subj_ind])) /
      reg_num_query
    
    # For each shuffled query set, repeat the overlap computation.
    prop_overl_shuffled <- unlist(mclapply(1:length(shuffled_regions_gr), 
                                           FUN = function(x) {
      length(unique(shuf_subj_overl_qhits[[x]][
        shuf_subj_overl_shits[[x]] %in% corresp_subj_ind])) /
                                               reg_num_shuf[x]
    }, mc.cores = ncores))
    
    # Compute an enrichment statistic.
    # This is defined as the proportion of query intervals overlapping with a 
    # subject interval divided by the mean corresponding proportion across the 
    # shuffflings.
    mean_shuf <- mean(prop_overl_shuffled)
    if (mean_shuf == 0) {
      if (prop_overl_query != 0) {
        enrichment_stat <- Inf
      } else {
        enrichment_stat <- prop_overl_query / mean_shuf
      }
    } else {
      enrichment_stat <- prop_overl_query / mean_shuf
    }
    
    # Compute the proportion of shufflings that yielded the same or a larger 
    # proportion of query intervals overlapping with corresponding subject 
    # intervals.
    # This is the enrichment P value.
    pval <- length(which(prop_overl_shuffled >= prop_overl_query)) /
      length(shuffled_regions_gr)
    
    # Return a data frame containing the results.
    temp_df <- data.frame(ID = y, 
                          Overl_prop_query = prop_overl_query,
                          Mean_overl_prop_shuf = mean_shuf,
                          Enrichment = enrichment_stat, 
                          Pval = pval)
    
    # Return this data frame.
    return(temp_df)
    
  }))
  
  # Correct P values for multiple hypothesis testing using the Benjamini-
  # Hochberg method.
  enrichment_dat$Pval_adj_BH <- p.adjust(enrichment_dat$Pval, method = "BH")
  
  # Return value.
  return(enrichment_dat)
  
}


# Print a status message.
print("Computing enrichments")

# Specify the meta data columns by which TEs will be grouped.
# For ma analysis. I will group TEs into four different levels:
# 1. All TEs together (regardless of their class etc.)
# 2. Repeat classes (column "V12")
# 3. Repeat family (Column "V13")
# 4. Repeat name (Column "V11")
grouping_cols <- c("all_tes", "V12", "V13", "V11")
# grouping_cols <- c("all_tes", "V12") # For debugging purposes
granularity_names <- c("all", "class", "family", "name")
# granularity_names <- c("all", "class") # For debugging purposes

# Iterate over the region classes.
enrichment_results_complete <- lapply(region_classes, FUN = function(i) {
  print(i)
  
  # Identify corresponding peaks.
  corresp_peaks_gr <- sig_gr[sig_gr$Signature_category == i]
  # corresp_peaks_gr <- corresp_peaks_gr[1:100] # For debugging purposes.
  
  # Extract corresponding peak shufflings.
  corresp_shufflings <- shuffled_intervals[[i]]
  
  # Iterate over the different levels of granularity with respect to the TEs.
  enrichment_dat_all_granularities <- lapply(grouping_cols, FUN = function(j) {
    print(j)
    compute_overl_enrichment_wrt_permutations(
      query_regions_gr = corresp_peaks_gr, 
      subject_regions_gr = rmsk_filt_gr, 
      subject_regions_granularity_col = j, 
      shuffled_regions_gr = corresp_shufflings, 
      ncores = 10
    )
  })
  names(enrichment_dat_all_granularities) <- granularity_names
  
  # Return the computed enrichment data.
  return(enrichment_dat_all_granularities)
  
})
names(enrichment_results_complete) <- region_classes

# Save the computed enrichment data.
enrichment_results_complete_outfile <- paste0(
  location, 
  "/te_analysis/gen_loc_enrichm_diffpeaks_skin_naive_atac.rds"
)
saveRDS(enrichment_results_complete, file = enrichment_results_complete_outfile)









#############################################################################
# Visualise enrichment results.
#############################################################################

# Iterate over the different levels of enrichment results.
plot_list <- do.call(c, lapply(region_classes, FUN = function(region_class) {
  lapply(granularity_names, FUN = function(granularity) {
    
    # Extract corresponding enrichment data.
    plotting_data <- enrichment_results_complete[[region_class]][[granularity]]
    
    # Compute logarithms of enrichment and adjusted P value.
    plotting_data$logenr <- log2(plotting_data$Enrichment)
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
    } else {
      y_step <- 1
    }
    
    # Prepare data to label the top-10 enriched TEs.
    top_inds <- order(plotting_data$Enrichment, decreasing = T)[
      1:min(10, nrow(plotting_data))]
    text_data <- plotting_data[top_inds, ]
    
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
        name = "Enrichment\n[log2(Prop. overlaps div. by mean prop.\noverlaps across rand. shufflings)]") +
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
      ggtitle(paste0(region_class, "; ", granularity)) +
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
}))

# Save the plots.
outfile_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_enrichm_diffpeaks_skin_naive_atac.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_enrichm_diffpeaks_skin_naive_atac.rds"
)
pdf(outfile_pdf, width = 7, height = 6)
print(plot_list)
dev.off()
saveRDS(plot_list, file = outfile_rds)

