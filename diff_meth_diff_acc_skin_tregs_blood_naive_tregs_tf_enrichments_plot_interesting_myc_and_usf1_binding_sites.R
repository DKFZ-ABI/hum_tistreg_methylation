# This R script looks at interesting predicted Myc and USF1 binding sites,
#   generates track plots and checks whether the associated genes are
#   differentially expressed. Interesting predicted binding sites are those
#   binding sites that overlap with a DMR and a peak that was not reported to be
#   differential.
# Author: Niklas Beumer



# Load required packages.
library(GenomicRanges)
library(bsseq)
library(ggplot2)
library(testit)
library(qpdf)
library(parallel)
library(cowplot)
library(rtracklayer)
library(Seurat)
library(Signac)
library(BSgenome)
library(BSgenome.Hsapiens.1000genomes.hs37d5)


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

# Read in the Seurat object with motif annotations.
sc_data_cd4_file <- paste0(
  location,
  "/treg_hierarchies/seurat_obj_cd4_skin_treg_blood_naive_treg_updated_w_motifs_and_footprints_and_chromvar-f_homer_relevant_TFs.rds"
)
sc_data_cd4 <- readRDS(sc_data_cd4_file)

# Read in the regions that are differentially methylated between skin Tregs and 
# blood naive Tregs.
dmr_table_file <- paste0(
  location,
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_gap_0.15_sterrratio_0.5_signature_regions_filtered.txt"
)
dmr_table <- read.table(
  dmr_table_file,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)

# Read in the table with results from differential accessibility analysis
# between skin Tregs and blood naive Tregs.
diff_acc_file <- paste0(
  location,
  "/differential_accessibility/diff_acc_analysis_Skin_TregBlood_naive_Treg_results_with_significance.txt"
)
diff_acc <- read.table(diff_acc_file,
                       header = T,
                       stringsAsFactors = F)

# Read in the results from differential expression analysis.
deseq2_res_file <- paste0(
  location,
  "/RNASeq/analysis_results/2022-01-14_diff_gene_expr_DESEq2_Skin_TregBlood_naive_Treg_results_filtered_with_significance.txt"
)
deseq2_res <- read.table(deseq2_res_file,
                         header = T,
                         stringsAsFactors = F)

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

# Read in the BSseq object containing the smoothed methylation data.
meth_data_file <- paste0(
  location,
  "/preprocessing_etc/quality_control_after_alignment/2022-03-14_bsseq_object_combined_all_samples_smoothed.rds"
)
meth_data <- readRDS(meth_data_file)

# Read in the information on Malte's scATAC-seq data for the cell types.
atac_seq_data_file <- paste0(
  location,
  "/materials_for_scATAC_seq_tracks/Overview_Maltes_Bam_files.csv"
)
atac_seq_data_overview <- read.csv(atac_seq_data_file, stringsAsFactors = F)

# Read in the sample mapping file for the RNA-Seq data.
rna_sample_mapping_file <- paste0(location, 
                                  "/sample_mapping_rnaseq_only_biol_rep.txt")
rna_sample_mapping <- read.table(
  rna_sample_mapping_file,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)

# Specify the paths to relevant ChIP-seq tracks from Encode.
chip_paths <- list(
  myc = paste0(
    location,
    "/external_data/2023-03-09_",
    c(
      "GSM822291_hg19_wgEncodeOpenChromChipHepg2CmycSig",
      "GSM822310_hg19_wgEncodeOpenChromChipK562CmycSig",
      "GSM822290_hg19_wgEncodeOpenChromChipGm12878CmycSig"
    ),
    ".bigWig"
  ),
  usf1 = paste0(
    location,
    "/external_data/2023-03-09_",
    c(
      "GSM803527_hg19_wgEncodeHaibTfbsHepg2Usf1Pcr1xRawRep1",
      "GSM803441_hg19_wgEncodeHaibTfbsK562Usf1V0416101RawRep1",
      "GSM803347_hg19_wgEncodeHaibTfbsGm12878Usf1Pcr2xRawRep1"
    ),
    ".bigWig"
  )
)
cell_lines <- c("Hepg2", "K562", "GM12878")

# Load the hs37d5 genome into the workspace.
genome_obj <- BSgenome.Hsapiens.1000genomes.hs37d5






#########################################################################
# Process methylation data.
#########################################################################

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

# Restrict the data to those samples that contain skin Tregs and blood naive Tregs.
# Also generate a version that additionally contains blood CCR8+ Tregs.
cell_types <- c("Skin Treg", "Blood naive Treg")
cell_types_2 <- c("Skin Treg", "Blood CCR8+ Treg", "Blood naive Treg")
cell_types_by_sample <- pData(meth_data)$Cell_type
samples_to_keep <- which(cell_types_by_sample %in% cell_types)
samples_to_keep_2 <- which(cell_types_by_sample %in% cell_types_2)
meth_data_relevant <- meth_data[, samples_to_keep]
meth_data_relevant_2 <- meth_data[, samples_to_keep_2]





#########################################################################
# Process gene data.
#########################################################################

# Increase start and end positions by 1 nucleotide in order to convert intervals to 1-based closed intervals.
genes$chromStart <- genes$chromStart + 1
exons$start <- exons$start + 1

# Generate a GRanges object out of the gene body data.
genes_gr <- makeGRangesFromDataFrame(
  genes,
  seqnames.field = "X.chrom",
  start.field = "chromStart",
  end.field = "chromEnd",
  keep.extra.columns = T
)


# Infer the location of the transcription start sites and prepare the corresponding data to appear in plots.
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





##################################################################################
# Identify interesting predicted binding sites.
# A binding site is interesting if it overlaps with a DMR of the
# skin-treg__hypomethylation class and a peak that was not reported as
# differential.
##################################################################################

# Turn the DMR locations into a GRanges object.
# Only consider those DMRs that are in the skin_treg__hypomethylation class.
dmrs_gr <- makeGRangesFromDataFrame(
  dmr_table[dmr_table$automatic_annotation == "Skin_Treg__hypomethylation", ], 
  keep.extra.columns = T
)
seqlevelsStyle(dmrs_gr) <- "UCSC"

# Identify peaks that were not reported to be significantly differential.
all_peaks <- rownames(sc_data_cd4)
signif_peaks <- rownames(diff_acc)[diff_acc$significant]
nonsignif_peaks <- setdiff(all_peaks, signif_peaks)

# Extract positions of predicted transcription factor binding sites.
motifs_obj <- Motifs(sc_data_cd4)
motifs_pos <- motifs_obj@positions

# Identify positions of predicted Myc binding sites.
# Note: Since there are two Myc motif sets, I only include regions from the first
# motif set that overlap with a region from the second motif set (to be more confident).
myc_pos_1_gr <- motifs_pos$`c-Myc(bHLH)/LNCAP-cMyc-ChIP-Seq(Unpublished)/Homer`
myc_pos_2_gr <- motifs_pos$`c-Myc(bHLH)/mES-cMyc-ChIP-Seq(GSE11431)/Homer`
myc_myc_overl <- findOverlaps(myc_pos_1_gr, myc_pos_2_gr, ignore.strand = T)
myc_pos_consensus_gr <- myc_pos_1_gr[unique(from(myc_myc_overl))]

# Identify positions of predicted UDF1 binding sites.
usf1_pos_gr <- motifs_pos$`USF1(bHLH)/GM12878-Usf1-ChIP-Seq(GSE32465)/Homer`

# For Mac and USF1, identify predicted binding sites that overlap with a DMR of the
# skin_treg__hypomethylation class and with a peak that was not reported to
# be differential.
pred_sites_relevant <- lapply(
  c("myc", "usf1"),
  FUN = function(x) {
    if (x == "myc") {
      tf_pos_gr <- myc_pos_consensus_gr
    } else {
      tf_pos_gr <- usf1_pos_gr
    }
    tf_dmr_overl <- findOverlaps(tf_pos_gr, dmrs_gr, ignore.strand = T)
    tf_restr_1 <- tf_pos_gr[unique(from(tf_dmr_overl))]
    tf_peak_overl <- findOverlaps(
      tf_restr_1,
      StringToGRanges(nonsignif_peaks, starts.in.df.are.0based = T),
      ignore.strand = T
    )
    tf_restr_2 <- tf_restr_1[unique(from(tf_peak_overl))]
    return(tf_restr_2)
  }
)
names(pred_sites_relevant) <- c("myc", "usf1")

# For each predicted binding site, identify DMRs and peaks that this site overlapped with.
pred_sites_relevant <- lapply(
  pred_sites_relevant,
  FUN = function(x) {
    dmr_overlaps <- findOverlaps(x, dmrs_gr, ignore.strand = T)
    x$overlapping_dmrs <- sapply(
      1:length(x),
      FUN = function(y) {
        paste(dmrs_gr$region_ID[to(dmr_overlaps)[from(dmr_overlaps) == y]], collapse = "; ")
      }
    )
    peak_overlaps <- findOverlaps(x,
                                  StringToGRanges(nonsignif_peaks, starts.in.df.are.0based = T),
                                  ignore.strand = T)
    x$overlapping_peaks <- sapply(
      1:length(x),
      FUN = function(y) {
        paste(nonsignif_peaks[to(peak_overlaps)[from(peak_overlaps) == y]], collapse = "; ")
      }
    )
    return(x)
  }
)


# Save information of interesting predicted binding sites.
pred_sites_relevant_outfile <- paste0(
  location,
  "/treg_hierarchies/diff_meth_diff_acc_skin_vs_naive_pred_myc_usf1_sites_overlapping_w_dmr_and_nondifferential_peak.rds"
)
saveRDS(pred_sites_relevant, file = pred_sites_relevant_outfile)






########################################################################
# Generate track plots for the interesting binding sites.
########################################################################

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
  # cell_type_col: Character; The name of the column in paths_table that contains cell type information.
  # bw_path_col: Character; The name of the column in paths_table that contains bigwig file paths.
  # paths_prefix: Character; A prefix to put in front of every path. If specified, must end with a "/".
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

plot_meth_and_other_signal_many_regions_custom_titles_custom_highlight <-
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
    #   in several regions using custom plot titles.
    # bsseq_obj: The BSseq object to retrieve methylation data from. In the pData, this object must contain
    #   the columns "Cell_type" and "col".
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
    
    # Extract smoothed methylation values for the extended regions.
    smoothed_meth_in_regions <- getMeth(bsseq_obj,
                                        regions = regions_extend_gr,
                                        type = "smooth",
                                        what = "perBase")
    
    # Extract the cell types, sample names and colours corresponding to the samples in the BSseq object.
    cell_types_for_plot <- pData(bsseq_obj)$Cell_type
    cell_type1 <- unique(cell_types_for_plot)[1]
    cell_type2 <- unique(cell_types_for_plot)[2]
    colours_for_plot <- pData(bsseq_obj)$col
    colour_1 <- unique(colours_for_plot)[1]
    colour_2 <- unique(colours_for_plot)[2]
    sample_names_for_plot <- rownames(pData(bsseq_obj))
    
    # Find out whether the scATAC-seq data for the two cell types are identical.
    if (atac_bigwig_path_table$Norm_bigwig_dir[atac_bigwig_path_table$Cell_type == cell_type1] ==
        atac_bigwig_path_table$Norm_bigwig_dir[atac_bigwig_path_table$Cell_type == cell_type2]) {
      atac_data_equal <- T
    } else {
      atac_data_equal <- F
    }
    
    # Read in Malte's scATAC-seq data.
    cell_type_1_atac_data <- read_bigwig_file_cell_type(
      paths_table = atac_bigwig_path_table,
      cell_type_col = "Cell_type",
      bw_path_col = "Norm_bigwig_dir",
      paths_prefix = atac_bigwig_prefix,
      cell_type = cell_type1
    )
    if (!atac_data_equal) {
      cell_type_2_atac_data <- read_bigwig_file_cell_type(
        paths_table = atac_bigwig_path_table,
        cell_type_col = "Cell_type",
        bw_path_col = "Norm_bigwig_dir",
        paths_prefix = atac_bigwig_prefix,
        cell_type = cell_type2
      )
    }
    
    # Read in the RNA track data for the two cell types.
    # Also specify colours to use in the RNA plot. Colours for the RNA plot are specified separately from
    # the colours for the methylation plot and the scATAC-seq track because I have to account for the
    # possibility that Skin CCR8+ Tregs need to be covered.
    cell_type_1_rna_data <- read_bigwig_file_cell_type(
      paths_table = rna_sample_mapping,
      cell_type_col = "Cell_type",
      bw_path_col = "Bigwig_path",
      cell_type = cell_type1
    )
    cell_type_2_rna_data <- read_bigwig_file_cell_type(
      paths_table = rna_sample_mapping,
      cell_type_col = "Cell_type",
      bw_path_col = "Bigwig_path",
      cell_type = cell_type2
    )
    rna_data_list <- list(cell_type_1_rna_data, cell_type_2_rna_data)
    names(rna_data_list) <- c(cell_type1, cell_type2)
    rna_colours <- c(colour_1, colour_2)
    
    # If one of the two cell types comprises Skin Tregs, add RNA data for the
    # corresponding CCR8+ subsets. Also add the corresponding colours.
    if ("Skin Treg" %in% cell_types_for_plot) {
      rna_data_list[["Skin CCR8+ Treg"]] <- read_bigwig_file_cell_type(
        paths_table = rna_sample_mapping,
        cell_type_col = "Cell_type",
        bw_path_col = "Bigwig_path",
        cell_type = "Skin CCR8+ Treg"
      )
      rna_colours <- c(rna_colours, "palevioletred1")
    }
    
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
        cell_type_1_atac_data <- cell_type_1_atac_data[from(findOverlaps(cell_type_1_atac_data, atac_subset_region))]
        if (!atac_data_equal) {
          cell_type_2_atac_data <- cell_type_2_atac_data[from(findOverlaps(cell_type_2_atac_data, atac_subset_region))]
        }
        
        # Prepare the data to plot scATAC-seq tracks.
        atac_plotting_data <- as.data.frame(cell_type_1_atac_data)
        if (!atac_data_equal) {
          atac_plotting_data$Cell_type <- cell_type1
          atac_plotting_data_2 <- as.data.frame(cell_type_2_atac_data)
          atac_plotting_data_2$Cell_type <- cell_type2
          atac_plotting_data <- rbind(atac_plotting_data, atac_plotting_data_2)
        } else {
          atac_plotting_data$Cell_type <- "unknown"
        }
        atac_plotting_data$Cell_type <- factor(atac_plotting_data$Cell_type,
                                               levels = c(cell_type1, cell_type2, "unknown"))
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
            values = c(colour_1, colour_2, "black"),
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
            rel_heights = c(2.5, 1.5, ifelse(
              length(rna_data_list) > 2, yes = 2.25, no = 1.5
            ), 3, 1.5)
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
            rel_heights = c(2.5, 1.5, ifelse(
              length(rna_data_list) > 2, yes = 2.25, no = 1.5
            ), 1.5)
          )
        }
        
        # Save the combined plot at a preliminary location.
        pdf(
          paste0(outfile_pdf, "_", l, ".pdf"),
          width = 5.5,
          height = ifelse(!(all(
            is.na(add_track_bws)
          )), yes = 10, no = 7)
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


# Generate a GRanges object for all DMRs (regardless of their tendency).
all_dmrs_gr <- makeGRangesFromDataFrame(dmr_table)

# Generate a GRanges object for all peaks, regardless of their differentiality.
# Modify the style of seqlevens so that they match with the seqlevels in the bsseq object.
all_peaks_gr <- StringToGRanges(all_peaks, starts.in.df.are.0based = T)
seqlevelsStyle(all_peaks_gr) <- "NCBI"

# Free up some memory.
rm("sc_data_cd4")

# Iterate over the two transcription factors.
for (tf in c("myc", "usf1")) {
  sites_to_use <- pred_sites_relevant[[tf]]
  
  # Modify the style of seqlevens so that they match with the seqlevels in the bsseq object.
  seqlevelsStyle(sites_to_use) <- "NCBI"
  
  # Remove strand information. Otherwise, the function will only plot genes on the same
  # strand as the predicted binding site.
  strand(sites_to_use) <- "*"
  
  # Generate custom plot titles.
  plot_titles <- paste0(
    GRangesToString(sites_to_use),
    ifelse(tf == "myc", yes = " (pred. c-Myc site)", no = " (pred. USF1 site)")
  )
  
  # Generate a GRanges object containing all predicted binding sites, extended by
  # 100 base pairs to either side (for better visibility).
  extended_binding_sites <- sites_to_use + 100
  
  # Specify the ChIP-seq tracks to plot.
  chip_paths_to_use <- chip_paths[[tf]]
  
  # Specify names for the ChIP-seq tracks.
  chip_names <- paste0(ifelse(tf == "myc", yes = "c-Myc", no = "USF1"),
                       "\nChIP-seq\n(",
                       cell_lines,
                       ")")
  
  # Generate a track plot.
  plot_meth_and_other_signal_many_regions_custom_titles_custom_highlight(
    bsseq_obj = meth_data_relevant,
    regions_gr = sites_to_use,
    genes = genes_gr,
    exons = exons,
    tss_pos = tss_df,
    plot_titles = plot_titles,
    atac_bigwig_path_table = atac_seq_data_overview,
    atac_bigwig_prefix = b330_space,
    rna_bigwig_path_table = rna_sample_mapping,
    add_track_bws = chip_paths_to_use,
    add_track_names = chip_names,
    meth_highlight = all_dmrs_gr,
    atac_highlight = all_peaks_gr,
    gene_highlight = extended_binding_sites,
    add_track_highlight = extended_binding_sites,
    outfile_pdf = paste0(
      location,
      "/treg_hierarchies/diff_meth_diff_acc_skin_vs_naive_track_plots_pred_",
      tf,
      "_sites_overlapping_w_dmr_and_nondifferential_peak.pdf"
    ),
    ncores = 5
  )
  
  # Generate track plots that also show blood CCR8+ Tregs.
  plot_meth_and_other_signal_many_regions_custom_titles_custom_highlight_many_cell_types(
    bsseq_obj = meth_data_relevant_2,
    regions_gr = sites_to_use,
    genes = genes_gr,
    exons = exons,
    tss_pos = tss_df,
    plot_titles = plot_titles,
    atac_bigwig_path_table = atac_seq_data_overview,
    atac_bigwig_prefix = b330_space,
    rna_bigwig_path_table = rna_sample_mapping,
    add_track_bws = chip_paths_to_use,
    add_track_names = chip_names,
    meth_highlight = all_dmrs_gr,
    atac_highlight = all_peaks_gr,
    gene_highlight = extended_binding_sites,
    add_track_highlight = extended_binding_sites,
    outfile_pdf = paste0(
      location,
      "/treg_hierarchies/diff_meth_diff_acc_skin_vs_naive_track_plots_pred_",
      tf,
      "_sites_overlapping_w_dmr_and_nondifferential_peak_2_bl_ccr8_treg.pdf"
    ),
    ncores = 5
  )
  
}





#######################################################################################
# Plot differential gene expression results for genes in the proximity
# of predicted c-Myc/USF binding sites that display (i) a hypomethylation in skin
# Tregs, (ii) no hyperaccessibility in skin Tregs but a visible accessibility signal 
# in both cell types, (iii) a peak in at least one of the three ChIP-seq tracks and
# (iv) transcriptomic up-regulation in skin Tregs.
########################################################################################

# Specify the genes to plot.
# GNA11 displays only borderline significance after multiple testing correction but its
# epigenetic signals are so pronounced that I still suggest it for validation.
genes_to_plot <- c(
  "MLPH", # corresponds to pred. Myc binding site 2−238466022−238466029 and pred. USF1 binding site 2−238466021−238466030.
  "RAPGEF5", # corresponds to pred. Myc binding site 7−22186887−22186894.
  "GNA11", # corresponds to pred. Myc binding site 19−3108902−3108909 and pred. USF1 binding site 19−3108901−3108910.
  "SHANK3", # corresponds to pred. Myc binding site 22−51170693−51170700.
  "FLNA", # corresponds to pred. Myc binding site X−153579763−153579770.
  "NCALD" # corresponds to pred. USF1 binding site 8−103126303−103126312.
)

# Retrict DESeq2 results to the relevant genes.
diff_exp_restr <- deseq2_res[genes_to_plot, ]
diff_exp_restr$Gene <- factor(rownames(diff_exp_restr), levels = rev(genes_to_plot))

# Compute -log10(q) values.
diff_exp_restr$logq <- (-1) * log10(diff_exp_restr$padj)

# Make sure that genes with negative Log2FCs are oriented to the left.
diff_exp_restr$logq[diff_exp_restr$log2FoldChange < 0] <- 
  (-1) * diff_exp_restr$logq[diff_exp_restr$log2FoldChange < 0]

# Generate a plot showing DESeq2 results.
deseq2_plot <- ggplot(diff_exp_restr) +
  aes(x = Gene, y = logq, fill = log2FoldChange) +
  scale_y_continuous(labels = abs, n.breaks = 10) +
  scale_fill_gradient2(low = "blue", mid = "grey", high = "red", 
                       name = "log2(fold change)") +
  geom_bar(width = 1, stat = "identity") +
  geom_hline(yintercept = c(-1, 1) * -log10(0.05)) +
  xlab("Gene") +
  ylab("-log10(q)") +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.ticks = element_line(colour = "black"))

# Save this plot.
deseq2_plot_pdf <- paste0(
  plot_outdir, 
  "/treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_deseq2_res_for_genes_close_to_high_interest_pred_myc_usf1_sites.pdf"
)
deseq2_plot_rds <- paste0(
  plot_rds_outdir, 
  "/treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_deseq2_res_for_genes_close_to_high_interest_pred_myc_usf1_sites.rds"
)
pdf(deseq2_plot_pdf, width = 5, height = 2)
print(deseq2_plot)
dev.off()
saveRDS(deseq2_plot, file = deseq2_plot_rds)








#########################################################################
# Save sequences for those predicted binding sites that turn out to be 
# of super-high interest in FASTA format.
#########################################################################

# Specify the indices of super-interesting predicted binding sites.
super_interesting_inds <- list(myc = c(61, 139, 317, 355, 359),
                               usf1 = c(74, 160, 316))

# Save the relevant predicted binding sites (+/- 100 base pairs)
# in FASTA format.
#-- Iterate over the transcription factors.
fasta_sequences <- unlist(lapply(names(pred_sites_relevant), FUN = function(tf) {
  #-- Get the relevant binding sites.
  rel_sites_this_tf <- pred_sites_relevant[[tf]]
  #-- Iterate over the corresponding super-interesting predicted binding sites.
  fasta_sequences_this_tf <- unlist(lapply(super_interesting_inds[[tf]], FUN = function(index) {
    #-- Get the location of the corresponding predicted binding site.
    corresp_site <- rel_sites_this_tf[index]
    seqlevelsStyle(corresp_site) <- "NCBI"
    #-- Extend the binding site by 100 base pairs.
    corresp_site_ext <- corresp_site + 100
    #-- Retrieve the sequence.
    sequence <- getSeq(x = genome_obj,
                       names = seqnames(corresp_site_ext),
                       start = start(corresp_site_ext),
                       end = end(corresp_site_ext),
                       as.character = T)
    #-- Mark the sub-sequence corresponding to the predicted binding site.
    sequence_marked <- paste0(substr(sequence, 1, 100),
                              "$",
                              substr(sequence, 101, 100 + width(corresp_site)),
                              "$",
                              substr(sequence, 101 + width(corresp_site), nchar(sequence)))
    #-- Prepare the Fasta header for this sequence.
    binding_site_str <- GRangesToString(corresp_site, sep = c("_", "_"))
    fasta_header <- paste0(">hs37d5_plus_strand_chr", 
                           binding_site_str, 
                           "_pred_", 
                           ifelse(tf == "myc", yes = "cmyc", no = "usf1"), 
                           "_binding_site_extended_by_100_bp_pred_binding_site_marked_by_dollar_signs")
    #-- Prepare the output in FASTA format.
    fasta_output <- c(fasta_header, sequence_marked, "")
    #-- Return value.
    return(fasta_output)
  }))
  #-- Return value.
  return(fasta_sequences_this_tf)
}))

# Save the extracted sequences.
fasta_sequences_outfile <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_diff_acc_skin_vs_naive_pred_myc_usf1_sites_sequences_for_superinteresting_sites.fasta"
)
write(fasta_sequences, file = fasta_sequences_outfile)
