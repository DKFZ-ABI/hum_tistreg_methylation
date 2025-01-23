# This script analyses correlation between methylation changes and chromatin
#   accessibility changes in DMR-peak pairs regarding differentiality between
#   skin Treg cells and blood naive Treg cells (preparatory analysis for
#   DMR-peak-gene links).
# Author: Niklas Beumer



# Load required packages.
library(bsseq)
library(Signac)
library(GenomicRanges)
library(ggplot2)
library(testit)
library(qpdf)
library(parallel)
library(cowplot)
library(rtracklayer)


# Define a location on /yyy.
b330_space <- "/yyy/"
location <- paste0(b330_space, "yyy/hm_treg_bs_rgnsbg")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis",
                     format(Sys.time(), "%m-%d-%y"),
                     sep = "/")
if (!dir.exists(plot_outdir)) {
  dir.create(plot_outdir)
}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the BSseq object containing the smoothed methylation data
# and restrict to the relevant cell types.
meth_data_file <- paste(
  location,
  "/preprocessing_etc/quality_control_after_alignment/2022-03-14_bsseq_object_combined_all_samples_smoothed.rds",
  sep = "/"
)
meth_data <- readRDS(meth_data_file)
relevant_cell_types <- c("Blood naive Treg", "Skin Treg")
meth_data <- meth_data[, pData(meth_data)$Cell_type %in% relevant_cell_types]

# Specify the colour palette for the cell types.
bsseq_metadata <- pData(meth_data)
cell_types_all <- unique(bsseq_metadata$Cell_type)
cell_type_col_palette <- c("blue", "darkorchid1")
names(cell_type_col_palette) <- cell_types_all

# Append this colour_palette to the meta data of the BSseq object.
bsseq_metadata$col <- cell_type_col_palette[bsseq_metadata$Cell_type]
pData(meth_data) <- bsseq_metadata

# Read in the list of differentially methylated regions between skin Tregs and
# blood naive Tregs.
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

# Read in the list of differential peaks between skin Tregs and blood naive
# Tregs.
atac_reg_file <- paste0(
  location,
  "/differential_accessibility/diff_acc_analysis_Skin_TregBlood_naive_Treg_results_with_significance.txt"
)
atac_reg <- read.table(atac_reg_file,
                       header = T,
                       stringsAsFactors = F)
atac_reg <- atac_reg[atac_reg$significant, ]
atac_reg$signature_class <- sapply(
  atac_reg$avg_log2FC,
  FUN = function(x) {
    paste0(ifelse(x < 0, yes = "Skin", no = "Blood_naive"),
           "_Treg__hyperaccessibility")
  }
)

# Specify the prefix of the files containing genomic positions of exons,
# transcription start sites, genes etc.
annotation_pref <-
  "/yyy/hg19_GRCh37_1000genomes/databases/gencode/gencode19/GencodeV19_"

# Read in the files containing gene and exon annotation for normal chromosomes.
genes <- read.table(
  paste0(annotation_pref, "Genes_plain.bed.gz"),
  header = T,
  stringsAsFactors = F,
  comment.char = ""
)
exons <- read.table(
  paste0(annotation_pref, "Exons_plain.bed.gz"),
  header = T,
  stringsAsFactors = F,
  comment.char = ""
)

# Increase start and end positions by 1 nucleotide in order to convert intervals
# to 1-based closed intervals.
genes$chromStart <- genes$chromStart + 1
exons$start <- exons$start + 1

# Generate a GRanges object out of the gene body data.
genes_gr <- GenomicRanges::makeGRangesFromDataFrame(
  genes,
  seqnames.field = "X.chrom",
  start.field = "chromStart",
  end.field = "chromEnd",
  keep.extra.columns = T
)


# Infer the location of the transcription start sites and prepare the 
# corresponding data to appear in plots.
tss_pos <- sapply(
  1:nrow(genes),
  FUN = function(x) {
    ifelse(genes$strand[x] == "+",
           yes = genes$chromStart[x],
           no = genes$chromEnd[x])
  }
)
tss_df <- data.frame(
  chr = genes$X.chrom,
  pos = tss_pos,
  gene = genes$name,
  strand = genes$strand
)

# Read in the information on Malte's scATAC-seq data for the cell types.
atac_seq_data_file <- paste0(
  location,
  "/materials_for_scATAC_seq_tracks/Overview_Maltes_Bam_files.csv"
)
atac_seq_data_overview <- read.csv(atac_seq_data_file, stringsAsFactors = F)

# Read in the sample mapping file.
rna_sample_mapping_path <- paste0(location, 
                                  "/sample_mapping_rnaseq_only_biol_rep.txt")
rna_sample_mapping <- read.table(
  rna_sample_mapping_path,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)





#####################################################
# Identify DMR-peak links.
#####################################################

# Turn both region sets into GRanges objects.
meth_reg_gr <- makeGRangesFromDataFrame(meth_reg, keep.extra.columns = T)
atac_reg_gr <- StringToGRanges(rownames(atac_reg), starts.in.df.are.0based = T)
seqlevelsStyle(atac_reg_gr) <- "NCBI"
atac_reg_gr$peak_ID <- rownames(atac_reg)
atac_reg_gr$signature_class <- atac_reg$signature_class

# Identify overlaps between DMRs and differential peaks.
meth_atac_reg_overl <- findOverlaps(meth_reg_gr, atac_reg_gr, ignore.strand = T)
meth_reg_overlapping_gr <- meth_reg_gr[from(meth_atac_reg_overl)]
atac_reg_overlapping_gr <- atac_reg_gr[to(meth_atac_reg_overl)]




#######################################################################
# Assign the intersecting scATAC-seq peaks to genes.
#######################################################################

# Identify 2kb-regions in front of transcription start sites.
genes$tss_loc <- sapply(
  1:nrow(genes),
  FUN = function(x) {
    ifelse(genes$strand[x] == "+",
           yes = genes$chromStart[x],
           no = genes$chromEnd[x])
  }
)
regions_2kb_start <- sapply(
  1:nrow(genes),
  FUN = function(x) {
    ifelse(genes$strand[x] == "+",
           yes = genes$tss_loc[x] - 2000,
           no = genes$tss_loc[x] + 1)
  }
)
regions_2kb_end <- sapply(
  1:nrow(genes),
  FUN = function(x) {
    ifelse(genes$strand[x] == "+",
           yes = genes$tss_loc[x] - 1,
           no = genes$tss_loc[x] + 2000)
  }
)
regions_2kb <- data.frame(
  chr = genes$X.chrom,
  start = regions_2kb_start,
  end = regions_2kb_end,
  gene = genes$name,
  strand = genes$strand
)
regions_2kb_gr <- makeGRangesFromDataFrame(regions_2kb, keep.extra.columns = T)

# Find overlaps between intersecting peaks and gene bodies or 2kb-regions in
# front of transcription start sites.
atac_genes_overlaps <- findOverlaps(atac_reg_overlapping_gr, genes_gr, ignore.strand = T)
genes_query_hits <- from(atac_genes_overlaps)
genes_subject_hits <- to(atac_genes_overlaps)
atac_2kb_overlaps <- findOverlaps(atac_reg_overlapping_gr, regions_2kb_gr, ignore.strand = T)
twokb_query_hits <- from(atac_2kb_overlaps)
twokb_subject_hits <- to(atac_2kb_overlaps)

# For each overlapping peak, collect information on what gene it is assigned to
# and by what overlap type.
# A peak is assigned to a gene if it overlaps with the gene body or with the
# 2kb-region in front of the transcription start site.
gene_names <- genes_gr$name
gene_names_2kb <- regions_2kb_gr$gene
gene_assignments <- c()
overlap_types <- c()
for (i in 1:length(atac_reg_overlapping_gr)) {
  this_reg_genes_subject_hits <- genes_subject_hits[genes_query_hits == i]
  this_reg_genes_subject_names <- gene_names[this_reg_genes_subject_hits]
  this_reg_2kb_subject_hits <- twokb_subject_hits[twokb_query_hits == i]
  this_reg_2kb_subject_names <- gene_names_2kb[this_reg_2kb_subject_hits]
  gene_names_unique <- unique(c(this_reg_2kb_subject_names, this_reg_genes_subject_names))
  if (length(gene_names_unique) > 0) {
    overlap_type_str <- sapply(
      gene_names_unique,
      FUN = function(x) {
        if (x %in% this_reg_genes_subject_names) {
          if (x %in% this_reg_2kb_subject_names) {
            str <- "gene_body_and_2kb_bef_tss"
          } else {
            str <- "gene_body"
          }
        } else {
          str <- "2kb_bef_tss"
        }
        return(str)
      }
    )
  } else {
    overlap_type_str <- "Not assigned"
    gene_names_unique <- "Not assigned"
  }
  gene_assignments <- c(gene_assignments, paste(gene_names_unique, collapse = ", "))
  overlap_types <- c(overlap_types, paste(overlap_type_str, collapse = ", "))
}
atac_reg_overlapping_gr$gene_assignments <- gene_assignments
atac_reg_overlapping_gr$overlap_types <- overlap_types





########################################################################
# Compute differences of cell-type-wise average raw methylation values
# in the relevant regions.
########################################################################

# Identify average raw methylation values (on the sample level).
avg_meth_sample_level <- getMeth(meth_data,
                                 regions = meth_reg_overlapping_gr,
                                 type = "raw",
                                 what = "perRegion")

# Compute cell-type-level averages of these sample-level methylation values.
avg_meth_celltype_level <- sapply(
  relevant_cell_types,
  FUN = function(x) {
    celltype_subs <- avg_meth_sample_level[, grep(x, colnames(avg_meth_sample_level), fixed = T)]
    return(rowMeans(celltype_subs, na.rm = T))
  }
)

# Compute differences of cell-type-level methylation values.
avg_meth_celltype_level_diffs <-
  avg_meth_celltype_level[, "Skin Treg"] -
  avg_meth_celltype_level[, "Blood naive Treg"]




##########################################################################
# Identify the effect sizes of accessibility changes for the DMR-peak
# pairs.
##########################################################################

# Extract information for the differential peaks that are in DMR-peak
# pairs.
atac_reg_overlapping <- atac_reg[atac_reg_overlapping_gr$peak_ID, ]
atac_reg_overlapping$peak_ID <- rownames(atac_reg_overlapping)
atac_reg_overlapping$gene_assignments <-
  atac_reg_overlapping_gr$gene_assignments
atac_reg_overlapping$overlap_types <- atac_reg_overlapping_gr$overlap_types

# Invert effect sizes to make the vlues reflect a comparison with blood
# naive Tregs as the base line. Currently, skin Tregs are the base line.
atac_reg_overlapping$avg_log2FC_bl_blood_naive <-
  (-1) * atac_reg_overlapping$avg_log2FC





##########################################################################
# Collect and save all results.
##########################################################################

# Prepare and concatenate results.
meth_reg_overlapping <- as.data.frame(meth_reg_overlapping_gr)
avg_meth_df <- as.data.frame(avg_meth_celltype_level)
colnames(avg_meth_df) <- paste0("Avg_raw_meth_", gsub(" ", "_", colnames(avg_meth_df)))
meth_diff_df <- as.data.frame(avg_meth_celltype_level_diffs)
colnames(meth_diff_df) <- "Avg_raw_meth_diff_bl_blood_naive"
meth_reg_overlapping_w_vals <- cbind(meth_reg_overlapping, avg_meth_df, meth_diff_df)
colnames(meth_reg_overlapping_w_vals) <- paste0("Methylation__", colnames(meth_reg_overlapping_w_vals))
colnames(atac_reg_overlapping) <- paste0("Accessibility__", colnames(atac_reg_overlapping))
combined_results <- cbind(meth_reg_overlapping_w_vals, atac_reg_overlapping)

# Save the DMR-peak links with their effect sizes.
combined_results_outfile <- paste0(
  location,
  "/treg_hierarchies/diff_meth_diff_acc_overlapping_dmr_peak_pairs_meth_acc_corr.txt"
)
write.table(
  combined_results,
  file = combined_results_outfile,
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)






######################################################################
# Generate a correlation plot.
######################################################################

# Compute Pearson's correlation coefficient and test this using a two-tailed
# test.
corr_test_res <- cor.test(
  combined_results$Methylation__Avg_raw_meth_diff_bl_blood_naive,
  combined_results$Accessibility__avg_log2FC_bl_blood_naive,
  method = "pearson",
  alternative = "two.sided"
)
exact_coef <- corr_test_res$estimate
exact_pval <- corr_test_res$p.value

# Save the results from the correlation test.
corr_test_res_outfile <- paste0(
  location,
  "/treg_hierarchies/diff_meth_diff_acc_overlapping_dmr_peak_pairs_meth_acc_corr_test_res.rds"
)
saveRDS(corr_test_res, file = corr_test_res_outfile)

# Check that the P value is below 1e-300.
if (!(exact_pval < 1e-300)) {
  stop("P value is not below 1e-300.")
}

# Generate the plot.
corr_plot <- ggplot(combined_results) +
  aes(x = Methylation__Avg_raw_meth_diff_bl_blood_naive, y = Accessibility__avg_log2FC_bl_blood_naive) +
  geom_point(size = 0.75) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_text(
    x = 0.8 *
      max(
        combined_results$Methylation__Avg_raw_meth_diff_bl_blood_naive
      ),
    y = 0.9 * max(combined_results$Accessibility__avg_log2FC_bl_blood_naive),
    label = paste0("r = ", round(exact_coef, 3), "\n", "P < 1e-300")
  ) +
  xlab("Raw methylation difference\n(Blood naive Tregs vs. Skin Tregs)") +
  ylab("Chr. accessibility log2(FC)\n(Blood naive Tregs vs. Skin Tregs)") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))

# Save the plot.
corr_plot_pdf <- paste0(
  plot_outdir,
  "/treg_hierarchies_diff_meth_diff_acc_dmr_peak_pairs_meth_acc_corr.pdf"
)
corr_plot_rds <- paste0(
  plot_rds_outdir,
  "/treg_hierarchies_diff_meth_diff_acc_dmr_peak_pairs_meth_acc_corr.rds"
)
pdf(corr_plot_pdf, width = 4, height = 4)
print(corr_plot)
dev.off()
saveRDS(corr_plot, file = corr_plot_rds)




############################################################################
# Generate track plots fo the DMR-peak-gene links.
############################################################################

get_cpg_pos_in_region <- function(bsseq_obj, region) {
  # This function extracts the CpG sites overlapping with a specified region.
  # bsseq_obj: The BSseq object to extract CpG sites from.
  # region: A GRanges object specifying the region to look at.
  # Dependencies: bsseq, GenomicRanges (need not be loaded but must be installed)
  # Value: A GRanges object with the CpG sites.
  # Author: Niklas Beumer
  
  # Find overlaps between the region and all CpGs in the BSseq object.
  region_overlaps <- GenomicRanges::findOverlaps(query = region, subject = bsseq_obj@rowRanges)
  
  # Index the CpG positions in the BSseq object by the indices of the subject hits.
  cpgs_in_region <- bsseq_obj@rowRanges[to(region_overlaps), ]
  
  # Return the extracted CpG sites.
  return(cpgs_in_region)
  
}

read_bigwig_file_cell_type <- function(paths_table,
                                       cell_type_col,
                                       bw_path_col,
                                       paths_prefix = "",
                                       cell_type) {
  # This function reads in a normalised bigwig file for a cell type.
  # paths_table: Data frame containing paths to bigwig files.
  # cell_type_col: Character; The name of the column in paths_table that
  #   contains cell type information.
  # bw_path_col: Character; The name of the column in paths_table that contains
  #   bigwig file paths.
  # paths_prefix: Character; A prefix to put in front of every path. If
  #   specified, must end with a "/".
  # cell_type: Character; The cell type to read the bigwig file for.
  # Dependencies: rtracklayer, testit
  # Value: A GRanges object containing the data from the bigwig file.
  # Author: Niklas Beumer
  
  # Specify the path to the bigwig file.
  bigwig_path_prelim <- unique(paths_table[paths_table[, cell_type_col] == cell_type, bw_path_col])
  bigwig_path <- paste0(paths_prefix, bigwig_path_prelim)
  if (length(grep(".bw$", bigwig_path, perl = T)) == 0) {
    # If the path only leads to a directory but not
    # to the file itself, look for the file name.
    bigwig_file_name <- list.files(bigwig_path, pattern = ".bw")
    bigwig_path <- paste0(bigwig_path, "/", bigwig_file_name)
  }
  
  # Read in the bigwig file.
  bw_data <- import(bigwig_path)
  
  # Check that bins whose length is not a multiple of 50 occur only once per chromosome.
  # This behaviour would be expected because the bin at the very end of a chromosome will not be
  # an exact multiple of 50.
  assert(all(table(as.character(
    seqnames(bw_data)
  )[(which(width(ranges(bw_data)) %% 50 != 0))]) == 1))
  
  # Synchronise the chromosome names with the naming convention used in my other work.
  seqlevelsStyle(bw_data) <- "NCBI"
  
  # Return the GRanges object containing the data from the bigwig files.
  return(bw_data)
  
}

plot_meth_and_other_signal_many_regions_custom_titles_custom_highlight_many_cell_types <-
  function(bsseq_obj,
           regions_gr,
           extend = rep(5000, length(regions_gr)),
           genes,
           exons,
           tss_pos,
           atac_bigwig_path_table,
           atac_bigwig_prefix,
           rna_bigwig_path_table,
           add_track_bws = NA,
           add_track_names = NA,
           meth_highlight = NA,
           atac_highlight = NA,
           rna_highlight = NA,
           gene_highlight = NA,
           add_track_highlight = NA,
           plot_titles,
           outfile_pdf,
           ncores) {
    # This function plots smoothed methalytion values, scATAC-seq signals and RNA-seq signals
    #   in several regions using custom plot titles (for more than two cell types.
    # bsseq_obj: The BSseq object to retrieve methylation data from. In the pData, this object must contain
    #   the columns "Cell_type" and "col". All cell types present in the BSseq object will be
    #   plotted.
    # regions_gr: A GRanges object containing the regions to plot.
    # extend: Numeric: The number of bases by which the regions wil be extended on either side.
    #   (The plot will range from the first to the last CpG inside this extended region).
    # genes: A GRanges object containing information on human genes. Gene names have the
    #   meta data column name "name".
    # exons: A data.frame containing data on human exons. Must contain the columns
    #   "gene_name", "start" and "end".
    # tss_pos: A data.frame containing data on human transcription start sites. Must contain the columns
    #   "pos", "gene" and "strand".
    # atac_bigwig_path_table: A data.frame containing information on my bigwig files with Malte's scATAC-seq data
    #   (The table follows my own convention, using the columns "Cell_type", "Bam_file_path" and "Norm_bigwig_dir").
    # atac_bigwig_prefix: Character; A prefix to put in front of all bigwig file paths mentioned in bigwig_path_table.
    #   Must end with "/".
    # rna_bigwig_path_table: A data.frame containing paths to bigwig files with RNA signals.
    # add_track_bws: Character; Vector of paths to BigWig files whose contents will be plotted as additional tracks.
    #   If no additional tracks should be plotted, set this to NA.
    # add_track_names: Character; A vector of names that should appear for the additional tracks specified
    #   in "add_track_bws". If no additional tracks should be plotted, set this to NA. If "add_track_bws" is
    #   set, this argument must contain a value.
    # meth_highlight: GRanges object containing the regions to highlight in the methylation track.
    #   May contain regions outside of the actual plotting window. Set this to NA if no highlighting is desired.
    # atac_highlight: GRanges object containing the regions to highlight in the ATAC track.
    #   May contain regions outside of the actual plotting window. Set this to NA if no highlighting is desired.
    # rna_highlight: GRanges object containing the regions to highlight in the RNA track.
    #   May contain regions outside of the actual plotting window. Set this to NA if no highlighting is desired.
    # gene_highlight: GRanges object containing the regions to highlight in the gene track.
    #   May contain regions outside of the actual plotting window. Set this to NA if no highlighting is desired.
    # add_track_highlight: GRanges object containing the regions to highlight in the additional tracks.
    #   May contain regions outside of the actual plotting window. Set this to NA if no highlighting is desired
    #   or no additional tracks should be plotted.
    # plot_titles: Character; The titles to use for the single plots. These should be in an order corresponding
    #   to the order of regions_gr.
    # outfile_pdf: Character; Path to a PDF file where the plot will be saved.
    # ncores: Numeric; The number of cores to use.
    # Dependencies: bsseq, ggplot2, cowplot, GenomicRanges, testit, rtracklayer, qpdf, parallel
    # Value: None, the function just generates and saves the plots.
    # Author: Niklas Beumer
    
    # Generate a GRanges object that contains the regions of interest extended by the specified ranges.
    regions_extend_gr <- regions_gr + extend
    
    # Extract the cell types, sample names and colours corresponding to the
    # samples in the BSseq object.
    cell_types_for_plot <- pData(bsseq_obj)$Cell_type
    cell_types_for_plot_unique <- unique(cell_types_for_plot)
    colours_for_plot <- unique(pData(bsseq_obj)$col)
    sample_names_for_plot <- rownames(pData(bsseq_obj))
    
    # Extract smoothed methylation values for the extended regions.
    smoothed_meth_in_regions <- getMeth(bsseq_obj,
                                        regions = regions_extend_gr,
                                        type = "smooth",
                                        what = "perBase")
    
    # Read in Malte's scATAC-seq data.
    cell_type_atac_data <- lapply(
      cell_types_for_plot_unique,
      FUN = function(x) {
        read_bigwig_file_cell_type(
          paths_table = atac_bigwig_path_table,
          cell_type_col = "Cell_type",
          bw_path_col = "Norm_bigwig_dir",
          paths_prefix = atac_bigwig_prefix,
          cell_type = x
        )
      }
    )
    names(cell_type_atac_data) <- c(cell_types_for_plot_unique)
    
    # Read in the RNA track data for the two cell types.
    # Also specify colours to use in the RNA plot.
    rna_data_list <- lapply(
      cell_types_for_plot_unique,
      FUN = function(x) {
        read_bigwig_file_cell_type(
          paths_table = rna_sample_mapping,
          cell_type_col = "Cell_type",
          bw_path_col = "Bigwig_path",
          cell_type = x
        )
      }
    )
    names(rna_data_list) <- cell_types_for_plot_unique
    rna_colours <- colours_for_plot
    
    # If additional tracks should be plotted, read in the corresponding BigWig files.
    # If a bin is only a single nucleotide long, increase its end position by 1 so
    # that it can be displayed in the plot.
    if (!(all(is.na(add_track_bws)))) {
      additional_tracks <- lapply(
        add_track_bws,
        FUN = function(x) {
          print(paste0("Reading ", x))
          bw_data <- import(x)
          seqlevelsStyle(bw_data) <- "NCBI"
          length_1_inds <- which(width(bw_data) == 1)
          end(bw_data)[length_1_inds] <- end(bw_data)[length_1_inds] + 1
          return(bw_data)
        }
      )
    }
    
    
    # Iterate over all specified regions and save a plot for each of them.
    void <- mclapply(
      1:length(regions_extend_gr),
      FUN = function(l) {
        #-------------- Preliminary steps -------------------
        
        # Extract the CpG positions inside the extended region from the BSseq object.
        region_to_plot <- regions_extend_gr[l]
        cpgs_in_extended_region <- get_cpg_pos_in_region(bsseq_obj, region_to_plot)
        
        # Extract the chromosome on which the region is located.
        chromosome_for_plot <- as.character(seqnames(region_to_plot))
        
        
        #-------------- Methylation data -------------------
        
        # Prepare the methylation data for plotting.
        meth_plotting_data <- data.frame(
          Pos = rep(
            start(cpgs_in_extended_region),
            ncol(smoothed_meth_in_regions[[l]])
          ),
          Meth = c(smoothed_meth_in_regions[[l]]),
          Sample = rep(sample_names_for_plot, each = nrow(smoothed_meth_in_regions[[l]])),
          Cell_type = rep(cell_types_for_plot, each = nrow(smoothed_meth_in_regions[[l]]))
        )
        meth_plotting_data$Cell_type <- factor(meth_plotting_data$Cell_type, levels = names(rna_data_list))
        
        # Identify all highlight regions that overlap with the region to plot.
        if (all(!(is.na(meth_highlight)))) {
          overlaps <- findOverlaps(query = meth_highlight,
                                   subject = region_to_plot,
                                   ignore.strand = T)
          regions_to_highlight <- as.data.frame(meth_highlight[from(overlaps)])
          regions_to_highlight$start <- sapply(
            regions_to_highlight$start,
            FUN = function(start) {
              max(c(start, min(meth_plotting_data$Pos)))
            }
          )
          regions_to_highlight$end <- sapply(
            regions_to_highlight$end,
            FUN = function(end) {
              min(c(end, max(meth_plotting_data$Pos)))
            }
          )
          regions_to_highlight <- regions_to_highlight[regions_to_highlight$end > min(meth_plotting_data$Pos) &
                                                         regions_to_highlight$start < max(meth_plotting_data$Pos), ]
        }
        
        # Generate the plot showing smoothed methylation values.
        meth_plot <- ggplot(meth_plotting_data) +
          scale_x_continuous(
            name = paste0("Position (Chr. ", chromosome_for_plot, ")"),
            limits = c(
              min(meth_plotting_data$Pos),
              max(meth_plotting_data$Pos)
            ),
            expand = c(0, 0)
          ) +
          scale_y_continuous(
            name = "Smoothed\nmethylation",
            limits = c(0, 1),
            expand = expansion(mult = c(0.15, 0))
          ) +
          aes(
            x = Pos,
            y = Meth,
            colour = Cell_type,
            group = Sample
          ) +
          scale_colour_manual(
            breaks = names(rna_data_list),
            values = rna_colours,
            name = "Cell type",
            drop = F,
            guide = guide_legend(
              override.aes = list(size = 6),
              # Make line legend appear as rectangles
              nrow = ifelse(length(rna_data_list) > 2, yes = 2, no = 1)
            )
          )
        if (all(!(is.na(meth_highlight)))) {
          meth_plot <- meth_plot +
            geom_rect(
              data = regions_to_highlight,
              mapping = aes(xmin = start, xmax = end),
              ymin = -0.15,
              ymax = 1,
              fill = "yellow",
              alpha = 0.5,
              inherit.aes = F
            )
        }
        meth_plot <- meth_plot +
          geom_line() +
          geom_rug(
            y = 0,
            sides = "b",
            colour = "black",
            length = unit(0.1, "npc")
          ) +
          ggtitle(plot_titles[l]) +
          theme_classic() +
          theme(
            axis.text.y = element_text(colour = "black"),
            axis.ticks.y = element_line(colour = "black"),
            legend.position = "top",
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5)
          )
        
        
        #-------------- scATAC-seq data -------------------
        
        # Identify all highlight regions that overlap with the region to plot.
        if (all(!(is.na(atac_highlight)))) {
          overlaps <- findOverlaps(query = atac_highlight,
                                   subject = region_to_plot,
                                   ignore.strand = T)
          regions_to_highlight <- as.data.frame(atac_highlight[from(overlaps)])
          regions_to_highlight$start <- sapply(
            regions_to_highlight$start,
            FUN = function(start) {
              max(c(start, min(meth_plotting_data$Pos)))
            }
          )
          regions_to_highlight$end <- sapply(
            regions_to_highlight$end,
            FUN = function(end) {
              min(c(end, max(meth_plotting_data$Pos)))
            }
          )
          regions_to_highlight <- regions_to_highlight[regions_to_highlight$end > min(meth_plotting_data$Pos) &
                                                         regions_to_highlight$start < max(meth_plotting_data$Pos), ]
        }
        
        # Restrict the data to those bins that overlap with CpGs in the extended plotting region.
        atac_subset_region <- GRanges(seqnames = chromosome_for_plot, ranges = IRanges(
          min(meth_plotting_data$Pos),
          max(meth_plotting_data$Pos)
        ))
        atac_data_relevant <- lapply(
          cell_type_atac_data,
          FUN = function(x) {
            x[from(findOverlaps(x, atac_subset_region))]
          }
        )
        
        # Prepare the data to plot scATAC-seq tracks.
        atac_plotting_data <- do.call(rbind, lapply(
          1:length(atac_data_relevant),
          FUN = function(x) {
            temp_df <- as.data.frame(atac_data_relevant[[x]])
            temp_df$Cell_type <- names(atac_data_relevant)[x]
            return(temp_df)
          }
        ))
        atac_plotting_data$Cell_type <- factor(atac_plotting_data$Cell_type, levels = cell_types_for_plot_unique)
        atac_plotting_data$start <- sapply(
          atac_plotting_data$start,
          FUN = function(x) {
            max(x, min(meth_plotting_data$Pos))
          }
        )
        atac_plotting_data$end <- sapply(
          atac_plotting_data$end,
          FUN = function(x) {
            min(x, max(meth_plotting_data$Pos))
          }
        )
        
        # Generate the ATAC-seq signal plot.
        atac_plot <- ggplot(atac_plotting_data) +
          aes(
            xmin = start,
            xmax = end,
            ymax = score,
            fill = Cell_type
          ) +
          scale_x_continuous(
            name = paste0("Position (Chr. ", chromosome_for_plot, ")"),
            limits = c(
              min(meth_plotting_data$Pos),
              max(meth_plotting_data$Pos)
            ),
            expand = c(0, 0)
          ) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.05)), name = "Binned scATAC-seq\nread count [RPKM]") +
          scale_fill_manual(
            breaks = levels(atac_plotting_data$Cell_type),
            values = cell_type_col_palette,
            name = "Cell type"
          )
        if (all(!(is.na(atac_highlight)))) {
          atac_plot <- atac_plot +
            geom_rect(
              data = regions_to_highlight,
              mapping = aes(xmin = start, xmax = end),
              ymin = 0,
              ymax = 1.1 * max(atac_plotting_data$score),
              fill = "yellow",
              alpha = 0.5,
              inherit.aes = F
            )
        }
        atac_plot <- atac_plot +
          geom_rect(ymin = 0) +
          facet_wrap( ~ Cell_type, ncol = 1) +
          theme_classic() +
          theme(
            axis.text.y = element_text(colour = "black"),
            axis.ticks.y = element_line(colour = "black"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            strip.text = element_blank(),
            strip.background = element_blank(),
            legend.position = "none"
          )
        
        
        #-------------- RNA track -------------------
        
        # Identify all highlight regions that overlap with the region to plot.
        if (all(!(is.na(rna_highlight)))) {
          overlaps <- findOverlaps(query = rna_highlight,
                                   subject = region_to_plot,
                                   ignore.strand = T)
          regions_to_highlight <- as.data.frame(rna_highlight[from(overlaps)])
          regions_to_highlight$start <- sapply(
            regions_to_highlight$start,
            FUN = function(start) {
              max(c(start, min(meth_plotting_data$Pos)))
            }
          )
          regions_to_highlight$end <- sapply(
            regions_to_highlight$end,
            FUN = function(end) {
              min(c(end, max(meth_plotting_data$Pos)))
            }
          )
          regions_to_highlight <- regions_to_highlight[regions_to_highlight$end > min(meth_plotting_data$Pos) &
                                                         regions_to_highlight$start < max(meth_plotting_data$Pos), ]
        }
        
        # Prepare the data to plot RNA tracks.
        rna_df_list <- lapply(
          1:length(rna_data_list),
          FUN = function(m) {
            region_data <- rna_data_list[[m]][from(findOverlaps(rna_data_list[[m]], atac_subset_region))]
            df_prelim <- as.data.frame(region_data)
            df_prelim$Cell_type <- names(rna_data_list)[m]
            return(df_prelim)
          }
        )
        rna_plotting_data <- do.call(rbind, rna_df_list)
        rna_plotting_data$Cell_type <- factor(rna_plotting_data$Cell_type, levels = names(rna_data_list))
        rna_plotting_data$start <- sapply(
          rna_plotting_data$start,
          FUN = function(x) {
            max(x, min(meth_plotting_data$Pos))
          }
        )
        rna_plotting_data$end <- sapply(
          rna_plotting_data$end,
          FUN = function(x) {
            min(x, max(meth_plotting_data$Pos))
          }
        )
        
        # Generate the RNA track plot.
        rna_plot <- ggplot(rna_plotting_data) +
          aes(
            xmin = start,
            xmax = end,
            ymax = score,
            fill = Cell_type
          ) +
          scale_x_continuous(
            name = paste0("Position (Chr. ", chromosome_for_plot, ")"),
            limits = c(
              min(meth_plotting_data$Pos),
              max(meth_plotting_data$Pos)
            ),
            expand = c(0, 0)
          ) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.05)), name = "Binned RNA-seq\nread count [RPKM]") +
          scale_fill_manual(
            breaks = levels(rna_plotting_data$Cell_type),
            values = rna_colours,
            name = "Cell type"
          )
        if (all(!(is.na(rna_highlight)))) {
          rna_plot <- rna_plot +
            geom_rect(
              data = regions_to_highlight,
              mapping = aes(xmin = start, xmax = end),
              ymin = 0,
              ymax = 1.1 * max(rna_plotting_data$score),
              fill = "yellow",
              alpha = 0.5,
              inherit.aes = F
            )
        }
        rna_plot <- rna_plot +
          geom_rect(ymin = 0) +
          facet_wrap( ~ Cell_type, ncol = 1) +
          theme_classic() +
          theme(
            axis.text.y = element_text(colour = "black"),
            axis.ticks.y = element_line(colour = "black"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            strip.text = element_blank(),
            strip.background = element_blank(),
            legend.position = "none"
          )
        
        
        #-------------- Custom tracks (if requested) ----------------
        
        if (!(all(is.na(add_track_bws)))) {
          if (all(!(is.na(add_track_highlight)))) {
            overlaps <- findOverlaps(query = add_track_highlight,
                                     subject = region_to_plot,
                                     ignore.strand = T)
            regions_to_highlight <- as.data.frame(add_track_highlight[from(overlaps)])
            regions_to_highlight$start <- sapply(
              regions_to_highlight$start,
              FUN = function(start) {
                max(c(start, min(meth_plotting_data$Pos)))
              }
            )
            regions_to_highlight$end <- sapply(
              regions_to_highlight$end,
              FUN = function(end) {
                min(c(end, max(meth_plotting_data$Pos)))
              }
            )
            regions_to_highlight <- regions_to_highlight[regions_to_highlight$end > min(meth_plotting_data$Pos) &
                                                           regions_to_highlight$start < max(meth_plotting_data$Pos), ]
          }
          
          # Restrict the data to those bins that overlap with CpGs in the extended plotting region.
          add_track_subset_region <- GRanges(seqnames = chromosome_for_plot, ranges = IRanges(
            min(meth_plotting_data$Pos),
            max(meth_plotting_data$Pos)
          ))
          add_track_data <- lapply(
            additional_tracks,
            FUN = function(x) {
              x_restr <- x[from(findOverlaps(x, add_track_subset_region))]
              if (length(x_restr) == 0) {
                x_restr <- add_track_subset_region
                x_restr$score <- 0
              }
              return(x_restr)
            }
          )
          
          # Prepare the data to plot additional tracks.
          # If no bins are present for the plotting region, use a placeholder with a score equal to 0.
          add_track_plotting_data <- as.data.frame(add_track_data[[1]])
          add_track_plotting_data$Name <- add_track_names[1]
          if (length(additional_tracks) > 1) {
            for (new_index in 2:length(additional_tracks)) {
              add_track_plotting_data_temp <- as.data.frame(add_track_data[[new_index]])
              add_track_plotting_data_temp$Name <- add_track_names[new_index]
              add_track_plotting_data <- rbind(add_track_plotting_data,
                                               add_track_plotting_data_temp)
            }
          }
          add_track_plotting_data$Name <- factor(add_track_plotting_data$Name, levels = add_track_names)
          add_track_plotting_data$start <- sapply(
            add_track_plotting_data$start,
            FUN = function(x) {
              max(x, min(meth_plotting_data$Pos))
            }
          )
          add_track_plotting_data$end <- sapply(
            add_track_plotting_data$end,
            FUN = function(x) {
              min(x, max(meth_plotting_data$Pos))
            }
          )
          
          # Generate the  plot.
          add_track_plot <- ggplot(add_track_plotting_data) +
            aes(xmin = start,
                xmax = end,
                ymax = score) +
            scale_x_continuous(
              name = paste0("Position (Chr. ", chromosome_for_plot, ")"),
              limits = c(
                min(meth_plotting_data$Pos),
                max(meth_plotting_data$Pos)
              ),
              expand = c(0, 0)
            ) +
            scale_y_continuous(
              expand = expansion(mult = c(0, 0.05)),
              name = "Signal",
              limits = c(0, max(add_track_plotting_data$score))
            )
          if (all(!(is.na(add_track_highlight)))) {
            add_track_plot <- add_track_plot +
              geom_rect(
                data = regions_to_highlight,
                mapping = aes(xmin = start, xmax = end),
                ymin = 0,
                ymax = 1.1 * max(add_track_plotting_data$score),
                fill = "yellow",
                alpha = 0.5,
                inherit.aes = F
              )
          }
          add_track_plot <- add_track_plot +
            geom_rect(ymin = 0, fill = "black") +
            facet_wrap( ~ Name, ncol = 1, strip.position = "right") +
            theme_classic() +
            theme(
              axis.text.y = element_text(colour = "black"),
              axis.ticks.y = element_line(colour = "black"),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              strip.text = element_text(size = 6)
            )
          
        }
        
        
        #-------------- Gene track -------------------
        
        # Identify all highlight regions that overlap with the region to plot.
        if (all(!(is.na(gene_highlight)))) {
          overlaps <- findOverlaps(query = gene_highlight,
                                   subject = region_to_plot,
                                   ignore.strand = T)
          regions_to_highlight <- as.data.frame(gene_highlight[from(overlaps)])
          regions_to_highlight$start <- sapply(
            regions_to_highlight$start,
            FUN = function(start) {
              max(c(start, min(meth_plotting_data$Pos)))
            }
          )
          regions_to_highlight$end <- sapply(
            regions_to_highlight$end,
            FUN = function(end) {
              min(c(end, max(meth_plotting_data$Pos)))
            }
          )
          regions_to_highlight <- regions_to_highlight[regions_to_highlight$end > min(meth_plotting_data$Pos) &
                                                         regions_to_highlight$start < max(meth_plotting_data$Pos), ]
        }
        
        # Identify genes that overlap with the extended region of interest.
        genes_in_extend_region_gr <- genes_gr[to(findOverlaps(regions_extend_gr[l], genes, type = "any"))]
        genes_in_extend_region <- genes_in_extend_region_gr$name
        
        # Prepare data on the genes and exons in the extended region of interest for plotting.
        if (length(genes_in_extend_region_gr) > 0) {
          y_pos_genes <- sapply(
            as.character(strand(genes_in_extend_region_gr)),
            FUN = function(x) {
              ifelse(x == "+", yes = 0.5, no = -0.5)
            }
          )
          
          region_genes_df <- data.frame(
            x = sapply(
              1:length(genes_in_extend_region_gr),
              FUN = function(y) {
                max(start(genes_in_extend_region_gr[y]),
                    min(meth_plotting_data$Pos))
              }
            ),
            xend = sapply(
              1:length(genes_in_extend_region_gr),
              FUN = function(y) {
                min(end(genes_in_extend_region_gr[y]),
                    max(meth_plotting_data$Pos))
              }
            ),
            y = y_pos_genes,
            yend = y_pos_genes
          )
          
          exons_in_extend_region <- exons[exons$gene_name %in% genes_in_extend_region_gr$name &
                                            exons$end > min(meth_plotting_data$Pos) &
                                            exons$start < max(meth_plotting_data$Pos), ]
          
          if (nrow(exons_in_extend_region) > 0) {
            region_exons_df <- data.frame(
              xmin = sapply(
                exons_in_extend_region$start,
                FUN = function(y) {
                  max(y, min(meth_plotting_data$Pos))
                }
              ),
              xmax = sapply(
                exons_in_extend_region$end,
                FUN = function(y) {
                  min(y, max(meth_plotting_data$Pos))
                }
              ),
              ymin = sapply(
                exons_in_extend_region$strand,
                FUN = function(x) {
                  ifelse(x == "+", yes = 0.3, no = -0.7)
                }
              ),
              ymax = sapply(
                exons_in_extend_region$strand,
                FUN = function(x) {
                  ifelse(x == "+", yes = 0.7, no = -0.3)
                }
              )
            )
          }
          gene_names_df <- data.frame(
            label = genes_in_extend_region_gr$name,
            x = sapply(
              1:nrow(region_genes_df),
              FUN = function(j) {
                mean(c(region_genes_df$x[j], region_genes_df$xend[j]))
              }
            ),
            y = sapply(
              as.character(strand(genes_in_extend_region_gr)),
              FUN = function(x) {
                ifelse(x == "+", yes = 0.15, no = -0.15)
              }
            )
          )
          tss_in_extend_region <- tss_pos[tss_pos$gene %in% genes_in_extend_region_gr$name &
                                            tss_pos$pos >= min(meth_plotting_data$Pos) &
                                            tss_pos$pos <= max(meth_plotting_data$Pos), ]
          if (nrow(tss_in_extend_region) > 0) {
            region_tss_plotting_data_1 <- data.frame(
              x = tss_in_extend_region$pos,
              xend = tss_in_extend_region$pos,
              y = sapply(
                1:nrow(tss_in_extend_region),
                FUN = function(x) {
                  ifelse(tss_in_extend_region$strand[x] == "+",
                         yes = 0.5,
                         no = -0.5)
                }
              ),
              yend = sapply(
                1:nrow(tss_in_extend_region),
                FUN = function(x) {
                  ifelse(tss_in_extend_region$strand[x] == "+",
                         yes = 0.85,
                         no = -0.85)
                }
              )
            )
            region_tss_plotting_data_2 <- data.frame(
              x = tss_in_extend_region$pos,
              xend = sapply(
                1:nrow(tss_in_extend_region),
                FUN = function(x) {
                  ifelse(
                    tss_in_extend_region$strand[x] == "+",
                    yes = tss_in_extend_region$pos[x] + 0.075 * diff(range(meth_plotting_data$Pos)),
                    no = tss_in_extend_region$pos[x] - 0.075 * diff(range(meth_plotting_data$Pos))
                  )
                }
              ),
              y = sapply(
                1:nrow(tss_in_extend_region),
                FUN = function(x) {
                  ifelse(tss_in_extend_region$strand[x] == "+",
                         yes = 0.85,
                         no = -0.85)
                }
              ),
              yend = sapply(
                1:nrow(tss_in_extend_region),
                FUN = function(x) {
                  ifelse(tss_in_extend_region$strand[x] == "+",
                         yes = 0.85,
                         no = -0.85)
                }
              )
            )
          }
        }
        
        # Generate a gene annotaton track.
        gene_plot <- ggplot() +
          scale_y_continuous(limits = c(-1, 1)) +
          scale_x_continuous(
            name = paste0("Position (Chr. ", chromosome_for_plot, ")"),
            limits = c(
              min(meth_plotting_data$Pos),
              max(meth_plotting_data$Pos)
            ),
            expand = c(0, 0)
          )
        if (all(!(is.na(gene_highlight)))) {
          gene_plot <- gene_plot +
            geom_rect(
              data = regions_to_highlight,
              mapping = aes(xmin = start, xmax = end),
              ymin = -1.2,
              ymax = 1.2,
              fill = "yellow",
              alpha = 0.5
            )
        }
        if (length(genes_in_extend_region_gr) > 0) {
          gene_plot <- gene_plot +
            geom_segment(
              data = region_genes_df,
              mapping = aes(
                x = x,
                y = y,
                xend = xend,
                yend = yend
              ),
              colour = "black"
            ) +
            geom_text(
              data = gene_names_df,
              mapping = aes(x = x, y = y, label = label),
              size = 3
            )
          if (nrow(exons_in_extend_region) > 0) {
            gene_plot <- gene_plot +
              geom_rect(
                data = region_exons_df,
                mapping = aes(
                  xmin = xmin,
                  xmax = xmax,
                  ymin = ymin,
                  ymax = ymax
                ),
                fill = "grey50"
              )
          }
          if (nrow(tss_in_extend_region) > 0) {
            gene_plot <- gene_plot +
              geom_segment(
                data = region_tss_plotting_data_1,
                mapping = aes(
                  x = x,
                  y = y,
                  xend = xend,
                  yend = yend
                ),
                size = 1,
                lineend = "square"
              ) +
              geom_segment(
                data = region_tss_plotting_data_2,
                mapping = aes(
                  x = x,
                  y = y,
                  xend = xend,
                  yend = yend
                ),
                arrow = arrow(
                  type = "closed",
                  angle = 20,
                  length = unit(0.1, "inches")
                ),
                size = 1
              )
          }
        }
        gene_plot <- gene_plot +
          theme_classic() +
          theme(
            axis.text.x = element_text(colour = "black"),
            axis.ticks.x = element_line(colour = "black"),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank()
          )
        
        
        #-------------- Combination of plots --------------------
        
        # Combine the plots into one and return the combined plot.
        if (!(all(is.na(add_track_bws)))) {
          combined_plot <- plot_grid(
            meth_plot,
            atac_plot,
            rna_plot,
            add_track_plot,
            gene_plot,
            ncol = 1,
            axis = "lr",
            align = "v",
            rel_heights = c(2.5, 1.5, 1.5, 3, 1.5)
          )
        } else {
          combined_plot <- plot_grid(
            meth_plot,
            atac_plot,
            rna_plot,
            gene_plot,
            ncol = 1,
            axis = "l",
            align = "v",
            rel_heights = c(2.5, 1.75, 1.75, 1.5)
          )
        }
        
        # Save the combined plot at a preliminary location.
        pdf(
          paste0(outfile_pdf, "_", l, ".pdf"),
          width = 5.5,
          height = ifelse(!(all(
            is.na(add_track_bws)
          )), yes = 10, no = 8)
        )
        print(combined_plot)
        dev.off()
        
        # Return an empty string so as to not inflate memory.
        return("")
        
      },
      mc.cores = ncores
    )
    
    # Combine the single PDF files into one. Do so in chunks of 10,000 plots.
    # Remove temporary files.
    if (length(regions_gr) > 1) {
      counter <- 0
      for (l in seq(1, length(regions_gr), 10000)) {
        counter <- counter + 1
        files_to_merge <- paste0(outfile_pdf, "_", l:min((l + 9999), length(regions_gr)), ".pdf")
        chunk_outfile_name <- paste0(outfile_pdf, "_intermediate_", counter, ".pdf")
        pdf_combine(files_to_merge, output = chunk_outfile_name)
        file.remove(files_to_merge)
      }
      comb_files_to_merge <- paste0(outfile_pdf, "_intermediate_", 1:counter, ".pdf")
      pdf_combine(comb_files_to_merge, output = outfile_pdf)
      file.remove(comb_files_to_merge)
    } else {
      file.rename(from = paste0(outfile_pdf, "_1.pdf"), to = outfile_pdf)
    }
    
  }





# Generate a GRanges object containing the overlapping regions between
# the DMR and the peak in each DMR-peak link.
meth_gr <- makeGRangesFromDataFrame(do.call(rbind, lapply(
  1:nrow(combined_results),
  FUN = function(x) {
    corresp_reg_id <- combined_results$Methylation__region_ID[x]
    return(meth_reg[meth_reg$region_ID == corresp_reg_id, 1:3])
  }
)))
acc_gr <- StringToGRanges(combined_results$Accessibility__peak_ID,
                          starts.in.df.are.0based = T)
seqlevelsStyle(acc_gr) <- "NCBI"
overl_gr <- pintersect(meth_gr, acc_gr)
assert(all(width(overl_gr) >= 1))

# Generate track plots for the DMR-peak links. In the methylation track,
# highlight all DMRs in DMR-peak links. In the ATAC track, highlight all
# peaks in DMR-peak links. In the gene track, highlight the overlap
# between the corresponding DMR and the corresponding peak.
plot_meth_and_other_signal_many_regions_custom_titles_custom_highlight_many_cell_types(
  bsseq_obj = meth_data,
  regions_gr = overl_gr,
  extend = 3000,
  genes = genes_gr,
  exons = exons,
  tss_pos = tss_df,
  atac_bigwig_path_table = atac_seq_data_overview,
  atac_bigwig_prefix = b330_space,
  rna_bigwig_path_table = rna_sample_mapping,
  add_track_bws = NA,
  add_track_names = NA,
  meth_highlight = meth_gr,
  atac_highlight = acc_gr,
  rna_highlight = NA,
  gene_highlight = overl_gr,
  add_track_highlight = NA,
  plot_titles = paste(
    combined_results$Methylation__region_ID,
    combined_results$Accessibility__peak_ID,
    sep = ";\n"
  ),
  outfile_pdf = paste0(
    location,
    "/treg_hierarchies/diff_meth_diff_acc_overlapping_dmr_peak_pairs_meth_acc_corr_track_plots.pdf"
  ),
  ncores = 10
)
