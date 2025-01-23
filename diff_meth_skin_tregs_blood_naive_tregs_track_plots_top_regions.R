# This script generates track plots for the top regions among the DMRs between
#   skin Treg cells and blood naive Treg cells.
# Author: Niklas Beumer



# Load required packages.
library(bsseq)
library(ggplot2)
library(testit)
library(qpdf)
library(parallel)
library(cowplot)
library(rtracklayer)
library(GenomicRanges)


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

# Specify the prefix of the files containing genomic positions of exons, 
# transcription start sites, genes etc.
annotation_pref <- "/yyy/hg19_GRCh37_1000genomes/databases/gencode/gencode19/GencodeV19_"

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

# Read in the BSseq object containing the smoothed methylation data.
meth_data <- readRDS(
  paste(
    location,
    "/preprocessing_etc/quality_control_after_alignment/2022-03-14_bsseq_object_combined_all_samples_smoothed.rds",
    sep = "/"
  )
)

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
cell_types <- c("Skin Treg", "Blood naive Treg")
cell_types_by_sample <- pData(meth_data)$Cell_type
samples_to_keep <- which(cell_types_by_sample %in% cell_types)
meth_data_relevant <- meth_data[, samples_to_keep]

# Read in differentially methylated regions between skin Tregs and blood naive Tregs.
diff_reg_file <- paste0(
  location,
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_gap_0.15_sterrratio_0.5_signature_regions_filtered.txt"
)
diff_reg <- read.table(
  diff_reg_file,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)

# Read in the information on Malte's scATAC-seq data for the cell types.
atac_seq_data_file <- paste0(location,
                             "/materials_for_scATAC_seq_tracks/Overview_Maltes_Bam_files.csv")
atac_seq_data_overview <- read.csv(atac_seq_data_file, stringsAsFactors = F)

# Read in the sample mapping file for the RNA-Seq data.
rna_sample_mapping_file <- paste0(location, "/sample_mapping_rnaseq_only_biol_rep.txt")
rna_sample_mapping <- read.table(
  rna_sample_mapping_file,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)








#########################################################################
# Identify top regions.
# A top region displays a difference in raw methylation of at least 0.5
# (on autosomes) and at least 0.25 (on chromosome X).
#########################################################################

# Make GRanges objects of the DMRs.
skin_hypo_gr <- makeGRangesFromDataFrame(
  diff_reg[diff_reg$automatic_annotation == "Skin_Treg__hypomethylation", ], 
  keep.extra.columns = T
)
naive_hypo_gr <- makeGRangesFromDataFrame(
  diff_reg[diff_reg$automatic_annotation == "Blood_naive_Treg__hypomethylation", ], 
  keep.extra.columns = T
)

# Compute sample-level average raw methylation values for these regions.
skin_hypo_avg_raw_meth <- getMeth(
  meth_data_relevant,
  regions = skin_hypo_gr,
  type = "raw",
  what = "perRegion"
)
naive_hypo_avg_raw_meth <- getMeth(
  meth_data_relevant,
  regions = naive_hypo_gr,
  type = "raw",
  what = "perRegion"
)

# Aggregate these into cell-level values.
skin_hypo_avg_raw_meth_aggr <- sapply(
  cell_types,
  FUN = function(x) {
    celltype_subs <- skin_hypo_avg_raw_meth[, grep(x, colnames(skin_hypo_avg_raw_meth), fixed = T)]
    return(rowMeans(celltype_subs, na.rm = T))
  }
)
naive_hypo_avg_raw_meth_aggr <- sapply(
  cell_types,
  FUN = function(x) {
    celltype_subs <- naive_hypo_avg_raw_meth[, grep(x, colnames(naive_hypo_avg_raw_meth), fixed = T)]
    return(rowMeans(celltype_subs, na.rm = T))
  }
)

# Compute differences in aggregated raw mthylation values betwen the two cell types.
skin_hypo_diffs <- abs(rowDiffs(skin_hypo_avg_raw_meth_aggr))
naive_hypo_diffs <- abs(rowDiffs(naive_hypo_avg_raw_meth_aggr))

# Restrict the DMRs to those regions that display a difference of at least 0.5
# (on autosomes) or at least 0.25 (on chromosome X).
skin_hypo_seqnames <- as.character(seqnames(skin_hypo_gr))
skin_hypo_reg_to_keep <- which(sapply(
  1:length(skin_hypo_gr),
  FUN = function(x) {
    skin_hypo_diffs[x, 1] >= ifelse(skin_hypo_seqnames[x] != "X", yes = 0.5, no = 0.25)
  }
))
skin_hypo_top_regions_gr <- skin_hypo_gr[skin_hypo_reg_to_keep]
naive_hypo_seqnames <- as.character(seqnames(naive_hypo_gr))
naive_hypo_reg_to_keep <- which(sapply(
  1:length(naive_hypo_gr),
  FUN = function(x) {
    naive_hypo_diffs[x, 1] >= ifelse(naive_hypo_seqnames[x] != "X", yes = 0.5, no = 0.25)
  }
))
naive_hypo_top_regions_gr <- naive_hypo_gr[naive_hypo_reg_to_keep]

# Collect and save all data.
skin_hypo_top_regions_df <- cbind(
  as.data.frame(skin_hypo_top_regions_gr),
  data.frame(Abs_raw_meth_diff = skin_hypo_diffs[skin_hypo_reg_to_keep, 1])
)
skin_hypo_top_regions_df_outfile <- paste0(
  location,
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_skin_treg_hypomethylation_top_regions_cutoff_0.5_0.25.txt"
)
write.table(
  skin_hypo_top_regions_df,
  file = skin_hypo_top_regions_df_outfile,
  sep = "\t",
  row.names = F,
  col.names = T,
  quote = F
)
naive_hypo_top_regions_df <- cbind(
  as.data.frame(naive_hypo_top_regions_gr),
  data.frame(Abs_raw_meth_diff = naive_hypo_diffs[naive_hypo_reg_to_keep, 1])
)
naive_hypo_top_regions_df_outfile <- paste0(
  location,
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_blood_naive_treg_hypomethylation_top_regions_cutoff_0.5_0.25.txt"
)
write.table(
  naive_hypo_top_regions_df,
  file = naive_hypo_top_regions_df_outfile,
  sep = "\t",
  row.names = F,
  col.names = T,
  quote = F
)










#####################################################################################
# Generate track plots for the top regions.
#####################################################################################

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



plot_meth_and_other_signal_many_regions_custom_titles <- function(bsseq_obj,
                                                                  regions_gr,
                                                                  extend = rep(5000, length(regions_gr)),
                                                                  genes,
                                                                  exons,
                                                                  tss_pos,
                                                                  atac_bigwig_path_table,
                                                                  atac_bigwig_prefix,
                                                                  rna_bigwig_path_table,
                                                                  highlight_regions = rep(T, length(regions_gr)),
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
  # highlight_regions: Logical, whether to highlight the regions of interest in the plots. Must be of
  #   the same length as regions_gr.
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
  
  
  # Iterate over all specified regions and save a plot for each of them.
  void <- mclapply(
    1:length(regions_extend_gr),
    FUN = function(l) {
      #-------------- Preliminary steps -------------------
      
      # Extract the CpG positions inside the extended region from the BSseq object.
      cpgs_in_extended_region <- get_cpg_pos_in_region(bsseq_obj, regions_extend_gr[l])
      
      # Extract the chromosome on which the region is located.
      chromosome_for_plot <- as.character(seqnames(regions_extend_gr[l]))
      
      
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
      
      # Prepare a data frame to indicate the region of interest in the plot.
      region_plotting_data <- data.frame(start = start(regions_gr[l]),
                                         end = end(regions_gr[l]))
      
      # Identify genes that overlap with the extended region of interest.
      # Prepare a plot title.
      genes_in_extend_region_gr <- genes_gr[to(findOverlaps(regions_extend_gr[l], genes, type = "any"))]
      genes_in_extend_region <- genes_in_extend_region_gr$name
      
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
      if (highlight_regions[l]) {
        meth_plot <- meth_plot +
          geom_rect(
            data = region_plotting_data,
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
      if (highlight_regions[l]) {
        atac_plot <- atac_plot +
          geom_rect(
            data = region_plotting_data,
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
      if (highlight_regions[l]) {
        rna_plot <- rna_plot +
          geom_rect(
            data = region_plotting_data,
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
      
      
      #-------------- Gene track -------------------
      
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
      if (highlight_regions[l]) {
        gene_plot <- gene_plot +
          geom_rect(
            data = region_plotting_data,
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
      
      # Save the combined plot at a preliminary location.
      pdf(paste0(outfile_pdf, "_", l, ".pdf"),
          width = 5.5,
          height = 7)
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


# Generate methylation track plots.
region_set_list <- list(skin_hypo_top_regions_gr, naive_hypo_top_regions_gr)
file_snips <- c("skin_treg_hypomethylation",
                "blood_naive_treg_hypomethylation")
for (i in 1:2) {
  plot_meth_and_other_signal_many_regions_custom_titles(
    bsseq_obj = meth_data_relevant,
    regions_gr = region_set_list[[i]],
    genes = genes_gr,
    exons = exons,
    tss_pos = tss_df,
    plot_titles = region_set_list[[i]]$region_ID,
    atac_bigwig_path_table = atac_seq_data_overview,
    atac_bigwig_prefix = b330_space,
    rna_bigwig_path_table = rna_sample_mapping,
    outfile_pdf = paste0(
      location,
      "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_",
      file_snips[i],
      "_top_regions_cutoff_0.5_0.25_track_plots.pdf"
    ),
    ncores = 10
  )
}
