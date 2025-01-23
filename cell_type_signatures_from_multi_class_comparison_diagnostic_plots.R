# This R script generates some descriptive plots for the cell type signature 
#   regions derived from a multi-class comparison.
# Author: Niklas Beumer



# Load required packages.
library(bsseq)
library(GenomicRanges)
library(ggplot2)
library(viridis)
library(testit)
library(parallel)
library(ComplexHeatmap)
library(circlize)


# Specify a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Specify the prefix of the files containing genomic positions of exons, 
# genes etc.
annotation_pref <- 
  "/yyy/hg19_GRCh37_1000genomes/databases/gencode/gencode19/GencodeV19_"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis",
                     format(Sys.time(), "%m-%d-%y"),
                     sep = "/")
if (!dir.exists(plot_outdir)) {
  dir.create(plot_outdir)
}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the BSseq object containing the smoothed methylation data.
meth_data <- readRDS(
  paste(
    location,
    "/preprocessing_etc",
    "quality_control_after_alignment",
    "2022-03-14_bsseq_object_combined_all_samples_smoothed.rds",
    sep = "/"
  )
)
meth_data_pdata <- pData(meth_data)

# Identify the cell types that are present in the data.
cell_types <- unique(meth_data_pdata$Cell_type)

# Specify the colour code for the cell types.
cell_type_col_palette <- c("blue", "cyan2", "orange", "darkorchid1", "red")
names(cell_type_col_palette) <- unique(meth_data_pdata$Cell_type)





#################################################################################
# Signature region read-in.
#################################################################################

# Read in the signature regions.
sig_reg_file <- paste0(
  location,
  "/differential_methylation/multi_class_signature_extraction_gap_0.15_largestsecondlargest_1.5_sterrratio_0.5_signature_regions_filtered.txt"
)
sig_reg <- read.table(
  sig_reg_file,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)

# Compute centres of signature regions.
# Region centres are defined as the mean between the start and end position of a region.
sig_reg$centre <- (sig_reg$start + sig_reg$end) / 2

# Quantify the numbers of regions in each signature category.
sig_reg_counts <- table(sig_reg$automatic_annotation)

# Identify signature categories with at least 10 regions.
sig_reg_cats_at_least_10 <- names(sig_reg_counts)[sig_reg_counts >= 10]
sig_reg_at_least_10 <- sig_reg[sig_reg$automatic_annotation %in% sig_reg_cats_at_least_10, ]

# Generate a colour palette for the signature categories with at least 10 regions.
set.seed(2022)
sig_reg_cats_at_least_10_col_palette <- 
  sample(rainbow(length(sig_reg_cats_at_least_10)), 
         length(sig_reg_cats_at_least_10))
names(sig_reg_cats_at_least_10_col_palette) <- sig_reg_cats_at_least_10

# Turn the signature region lists into GRanges objects.
sig_reg_gr <- makeGRangesFromDataFrame(sig_reg, keep.extra.columns = T)
sig_reg_at_least_10_gr <- makeGRangesFromDataFrame(sig_reg_at_least_10, keep.extra.columns = T)




#################################################################################
# Gene and exon data read-in and processing
#################################################################################

extract_gene_parts <- function(exons,
                               type,
                               index,
                               starting_from = NA,
                               class_name) {
  # This function extracts gene parts (specific exons or introns) from a list of genes.
  # exons: Data frame containing information on exons. Must contain the columns "X.chrom", "start", "end",
  #   "gene_name" and "strand".
  # type: Character; either "exon" or "intron".
  # index: Numeric; The index of the regions that should be extracted. If this is (e.g.) set to 2,
  #   the second exon/introns of each gene are extracted, if they exist.
  # starting from: Numeric; If not NA, extracts all regions corresponding to "type" starting from the
  #   given number. If (e.g.) set to 3, all exons/introns starting from the third exon/intron are returned,
  #   if they exist. If tghis parameter is not NA, the argument "index" is ignored.
  # class_name: Character: A class name to give all extracted regions in the output.
  # Dependencies: parallel
  # Value: A data frame with the extracted locations.
  # Author: Niklas Beumer
  
  # Iterate over the genes.
  unique_genes <- unique(exons$gene_name)
  locs_df <- do.call(rbind, mclapply(
    unique_genes,
    FUN = function(x) {
      # Extract the strand and chromosome for this gene.
      curr_strand <- unique(exons$strand[exons$gene_name == x])
      curr_chr <- unique(exons$X.chrom[exons$gene_name == x])
      
      # Extract data on all exon for this gene.
      corresp_exons_data <- exons[exons$gene_name == x, ]
      
      #---
      #--- Steps to extract single regions from every gene.
      if (is.na(starting_from)) {
        # Extract a particular exon.
        if (type == "exon") {
          # Check whether the desired region exists for this gene.
          if (index <= nrow(corresp_exons_data)) {
            # Extract the region, taking into account what strand the gene is located on.
            if (curr_strand == "+") {
              region_start <- corresp_exons_data[index, "start"]
              region_end <- corresp_exons_data[index, "end"]
            } else {
              region_start <- corresp_exons_data[nrow(corresp_exons_data) - index + 1, "start"]
              region_end <- corresp_exons_data[nrow(corresp_exons_data) - index + 1, "end"]
            }
            
            # If the region does not exist for this gene, introduce NAs.
          } else {
            region_start = NA
            region_end = NA
          }
          
          # Extract a particular intron.
        } else {
          # Check whether the desired region exists for this gene.
          if (index <= (nrow(corresp_exons_data) - 1)) {
            # The -1 is needed here because the number of introns is always
            # 1 below the number of exons.
            
            # Extract the region, taking into account what strand the gene is located on.
            if (curr_strand == "+") {
              region_start <- corresp_exons_data[index, "end"] + 1
              region_end <- corresp_exons_data[index + 1, "start"] - 1
            } else {
              region_start <- corresp_exons_data[nrow(corresp_exons_data) - index, "end"] + 1
              region_end <- corresp_exons_data[nrow(corresp_exons_data) - index + 1, "start"] - 1
            }
            
            # If the region does not exist for this gene, introduce NAs.
          } else {
            region_start = NA
            region_end = NA
          }
        }
        
        #---
        #--- If the starting_from parameter is specified, extract all relevant region starting
        #--- from the corresponding index.
      } else {
        # Extract a range of exons.
        if (type == "exon") {
          # Check whether the desired regions exist for this gene.
          if (starting_from <= nrow(corresp_exons_data)) {
            # Extract the regions, taking into account what strand the gene is located on.
            if (curr_strand == "+") {
              region_start <- corresp_exons_data[starting_from:nrow(corresp_exons_data), "start"]
              region_end <- corresp_exons_data[starting_from:nrow(corresp_exons_data), "end"]
            } else {
              region_start <- corresp_exons_data[1:(nrow(corresp_exons_data) - starting_from + 1), "start"]
              region_end <- corresp_exons_data[1:(nrow(corresp_exons_data) - starting_from + 1), "end"]
            }
            
            # If the regions do not exist for this gene, introduce NAs.
          } else {
            region_start = NA
            region_end = NA
          }
          
          # Extract a range of introns.
        } else {
          # Check whether the desired regions exist for this gene.
          if (starting_from <= nrow(corresp_exons_data) - 1) {
            # The -1 is needed here because the number of introns is always
            # 1 below the number of exons.
            
            # Extract the regions, taking into account what strand the gene is located on.
            if (curr_strand == "+") {
              region_start <- corresp_exons_data[starting_from:(nrow(corresp_exons_data) - 1), "end"] + 1
              region_end <- corresp_exons_data[(starting_from + 1):nrow(corresp_exons_data), "start"] - 1
            } else {
              region_start <- corresp_exons_data[1:(nrow(corresp_exons_data) - starting_from), "end"] + 1
              region_end <- corresp_exons_data[2:(nrow(corresp_exons_data) - starting_from + 1), "start"] - 1
            }
            
            # If the regions do not exist for this gene, introduce NAs.
          } else {
            region_start = NA
            region_end = NA
          }
        }
      }
      
      # Generate a data frame containing the start and end position as well as some information on
      # the underlying gene.
      region_df <- data.frame(
        chr = curr_chr,
        start = region_start,
        end = region_end,
        class = class_name,
        gene = x,
        strand = curr_strand
      )
      
      # Return this data frame.
      return(region_df)
      
    },
    mc.cores = 10
  ))
  
  # Remove all NA regions.
  locs_df <- locs_df[!(is.na(locs_df$start)), ]
  
  # Return this filtered data frame.
  return(locs_df)
  
}



# Read in the locations of genes and exons.
# Note: I did not discover any signature regions on chromosome MT. Thus, I do not additionally read in
# data regarding genes on this chromosome.
# Assert that gene names are unique.
genes <- read.table(
  paste0(annotation_pref, "Genes_plain.bed.gz"),
  header = T,
  stringsAsFactors = F,
  comment.char = ""
)
assert(length(unique(genes$name)) == nrow(genes))
exons <- read.table(
  paste0(annotation_pref, "Exons_plain.bed.gz"),
  header = T,
  stringsAsFactors = F,
  comment.char = ""
)

# Increase start positions by 1 nucleotide in order to convert intervals to 1-based closed intervals.
genes$chromStart <- genes$chromStart + 1
exons$start <- exons$start + 1

# Extract transcription start sites.
genes$tss_loc <- sapply(
  1:nrow(genes),
  FUN = function(x) {
    ifelse(genes$strand[x] == "+",
           yes = genes$chromStart[x],
           no = genes$chromEnd[x])
  }
)

# Perform some sanity checks before extracting regions.
void <- mclapply(
  1:nrow(genes),
  FUN = function(x) {
    # Get information on the current gene.
    curr_gene <- genes$name[x]
    curr_tss <- genes$tss_loc[x]
    curr_strand <- genes$strand[x]
    curr_exon_num <- genes$exonCount[x]
    corresp_exons_data <- exons[exons$gene_name == curr_gene, ]
    
    # Check that at least one exon exists for this gene
    assert(nrow(corresp_exons_data) > 0)
    
    # Check whether the strand information in the exons data align with the strand information in the genes data.
    assert(all(corresp_exons_data$strand == curr_strand))
    
    # Check whether the number of exons indicated in the genes table is equal to the number of
    # exons shown in the exons table.
    assert(curr_exon_num == nrow(corresp_exons_data))
    
    # Check whether the exons are ordered according to their genomic coordinates.
    assert(all(sort(corresp_exons_data$start) == corresp_exons_data$start))
    
    # Keep memory usage low.
    return("")
    
  },
  mc.cores = 10
)


# Specify regions starting from 2-kb before the transcription start site and ending right before
# the transcription start site.
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
regions_2kb_df <- data.frame(
  chr = genes$X.chrom,
  start = regions_2kb_start,
  end = regions_2kb_end,
  class = "2kb before TSS",
  gene = genes$name,
  strand = genes$strand
)

# For every gene, extract the locations of exons 1, 2 and 3, of introns 1, 2 and 3 as well
# as of the remaining gene body, separately for exons and introns.
exon_1_df <- extract_gene_parts(
  exons = exons,
  type = "exon",
  index = 1,
  class_name = "First exon"
)
exon_2_df <- extract_gene_parts(
  exons = exons,
  type = "exon",
  index = 2,
  class_name = "Second exon"
)
exon_3_df <- extract_gene_parts(
  exons = exons,
  type = "exon",
  index = 3,
  class_name = "Third exon"
)
exon_rem_df <- extract_gene_parts(
  exons = exons,
  type = "exon",
  starting_from = 4,
  class_name = "Remaining exons"
)
intron_1_df <- extract_gene_parts(
  exons = exons,
  type = "intron",
  index = 1,
  class_name = "First intron"
)
intron_2_df <- extract_gene_parts(
  exons = exons,
  type = "intron",
  index = 2,
  class_name = "Second intron"
)
intron_3_df <- extract_gene_parts(
  exons = exons,
  type = "intron",
  index = 3,
  class_name = "Third intron"
)
intron_rem_df <- extract_gene_parts(
  exons = exons,
  type = "intron",
  starting_from = 4,
  class_name = "Remaining introns"
)

# Check that the first exon always co-incides with the transcription start site.
# Furthermore, check that the end of each gene body (as indicated in the genes table)
# is is at the expected place in the regions data.
void <- mclapply(
  1:nrow(genes),
  FUN = function(x) {
    # Get information on the current gene.
    curr_gene <- genes$name[x]
    curr_tss <- genes$tss_loc[x]
    curr_strand <- genes$strand[x]
    genebody_end <- ifelse(curr_strand == "+",
                           yes = genes$chromEnd[x],
                           no = genes$chromStart[x])
    
    # Extract the boundaries of the corresponding first exon.
    first_exon_start <- exon_1_df$start[exon_1_df$gene == curr_gene]
    first_exon_end <- exon_1_df$end[exon_1_df$gene == curr_gene]
    
    # Check that the transcription start site is the same as the start of the first exon.
    assert(curr_tss == ifelse(curr_strand == "+", yes = first_exon_start, no = first_exon_end))
    
    # Check that the end of the gene body is at the expected place in the extracted regions.
    if (curr_strand == "+") {
      if (!(genebody_end %in% exon_rem_df$end)) {
        if (!(genebody_end %in% exon_3_df$end)) {
          if (!(genebody_end %in% exon_2_df$end)) {
            if (!(genebody_end %in% exon_1_df$end)) {
              stop(
                paste(
                  "The gene body end for",
                  curr_gene,
                  "was not found in the extracted regions."
                )
              )
            }
          }
        }
      }
    } else {
      if (!(genebody_end %in% exon_rem_df$start)) {
        if (!(genebody_end %in% exon_3_df$start)) {
          if (!(genebody_end %in% exon_2_df$start)) {
            if (!(genebody_end %in% exon_1_df$start)) {
              stop(
                paste(
                  "The gene body end for",
                  curr_gene,
                  "was not found in the extracted regions."
                )
              )
            }
          }
        }
      }
    }
    
    # Keep memory usage low.
    return("")
    
  },
  mc.cores = 10
)

# Collect all gene-based regions in a common data frame.
gene_regions <- rbind(
  regions_2kb_df,
  exon_1_df,
  intron_1_df,
  exon_2_df,
  intron_2_df,
  exon_3_df,
  intron_3_df,
  exon_rem_df,
  intron_rem_df
)
gene_regions <- gene_regions[order(gene_regions$gene), ]

# Save the relevant regions as a TXT file.
gene_regions_outfile <- paste0(
  location,
  "/differential_methylation/multi_class_signatures_gene_regions_for_region_centre_loc_pie_chart.txt"
)
write.table(
  gene_regions,
  file = gene_regions_outfile,
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)




#####################################################################################
# CpG island location read-in and processing.
#####################################################################################

get_cpg_island_shores_shelves <- function(cgi_gr) {
  # This function computes the location of CpG shores and shelves based on the locations of CpG islands.
  #   According to Yoann's definition, shores are the 2-kb
  #   flanking regions of the CpG islands and shelves are the regions that reach outwards for 2 additional kb.
  # cgi_gr: GRanges object containing the location of CpG islands.
  # Dependency: GenomicRanges
  # Value: A data frame containing CpG islands, shore and shelf regions.
  # Author: Niklas Beumer
  
  # Extend CpG islands by 2 kb to get regions spanning the shores.
  cgi_with_shores_gr <- cgi_gr + 2000
  
  # Extend CpG islands by 4 kb to get regions spanning the shelves.
  cgi_with_shores_and_shelves_gr <- cgi_gr + 4000
  
  # Identify shores and turn them into a data frame.
  shores_gr <- setdiff(cgi_with_shores_gr, cgi_gr)
  shores_df <- as.data.frame(shores_gr)[, 1:3]
  shores_df$class <- "Shore"
  
  # Identify shelves and turn them into a data frame.
  shelves_gr <- setdiff(cgi_with_shores_and_shelves_gr, cgi_with_shores_gr)
  shelves_df <- as.data.frame(shelves_gr)[, 1:3]
  shelves_df$class <- "Shelf"
  
  # Turn the CpG island locations into a data frame.
  cgi_df <- as.data.frame(cgi_gr)[, 1:3]
  cgi_df$class <- "CpG island"
  colnames(cgi_df) <- colnames(shores_df) <- colnames(shelves_df) <- c("chr", "start", "end", "class")
  
  # Return a data frame containing CpG islands, shores and shelves.
  # Harmonise chromosome names with what is used for the signature regions.
  comb_df <- rbind(cgi_df, shores_df, shelves_df)
  comb_df$chr <- gsub("chr", "", comb_df$chr)
  return(comb_df)
  
}

# Read in CpG islands ("cpgIslandExt" data set).
# Increase start positions by 1 in order to make the intervals 1-based and closed.
cpg_islands_file <- paste0(location,
                           "/external_data/2022-03-07_cpg_islands_ext_ucsc_hg19.txt.gz")
cpg_islands <- read.table(cpg_islands_file,
                          stringsAsFactors = F,
                          sep = "\t")
cpg_islands$V3 <- cpg_islands$V3 + 1
cpg_islands_gr <- makeGRangesFromDataFrame(
  cpg_islands,
  seqnames.field = "V2",
  start.field = "V3",
  end.field = "V4"
)

# Identify shore and shelf locations.
cpg_islands_shores_shelves <- get_cpg_island_shores_shelves(cpg_islands_gr)

# Save these regions in text format.
islands_shores_shelves_outfile <- paste0(
  location,
  "/differential_methylation/cpg_island_shore_shelf_regions_from_cpgIslandExt_list_for_region_centre_loc_pie_chart.txt"
)
write.table(
  cpg_islands_shores_shelves,
  file = islands_shores_shelves_outfile,
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)





###########################################################
# Methylation box plots.
###########################################################

# Compute within-cell-type average raw methylation values for all signature regions
# in categories with at least 10 regions. Save these values.
avg_raw_meth <- getMeth(meth_data,
                        regions = sig_reg_at_least_10_gr,
                        type = "raw",
                        what = "perRegion")
avg_raw_meth_aggr <- do.call(rbind, lapply(
  sig_reg_cats_at_least_10,
  FUN = function(x) {
    curr_signature_subs <- avg_raw_meth[sig_reg_at_least_10_gr$automatic_annotation == x, ]
    temp_df <- do.call(rbind, lapply(
      cell_types,
      FUN = function(y) {
        curr_celltype_subs <- curr_signature_subs[, grep(y, colnames(curr_signature_subs), fixed = T)]
        region_means <- apply(curr_celltype_subs, 1, FUN = mean, na.rm = T)
        temp_df_2 <- data.frame(
          Signature = x,
          Cell_type = y,
          Region = sig_reg_at_least_10$region_ID[sig_reg_at_least_10$automatic_annotation == x],
          Mean_raw_meth = region_means
        )
        return(temp_df_2)
      }
    ))
    return(temp_df)
  }
))
avg_raw_meth_aggr_outfile <- paste0(
  location,
  "/differential_methylation/multi_class_signatures_avg_raw_meth_aggr_by_region_and_cell_type.txt"
)
write.table(
  avg_raw_meth_aggr,
  file = avg_raw_meth_aggr_outfile,
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)


# Generate a box plot showing how strongly the signature regions (top-14) are methylated in
# the different cell types (based on raw methylation).
avg_raw_meth_aggr$Cell_type <- factor(avg_raw_meth_aggr$Cell_type, levels = cell_types)
facet_label_vect <- sapply(
  sig_reg_cats_at_least_10,
  FUN = function(x) {
    paste0(x, "\n(", sig_reg_counts[x], " regions)")
  }
)
label_function <- labeller(Signature = facet_label_vect)
signature_box_plot <- ggplot(avg_raw_meth_aggr) +
  aes(x = Cell_type, y = Mean_raw_meth, colour = Signature) +
  scale_y_continuous(limits = c(0, 1), name = "Raw methylation (within-cell-type mean)") +
  scale_colour_manual(values = sig_reg_cats_at_least_10_col_palette) +
  geom_boxplot(show.legend = F) +
  facet_wrap(
    ~ Signature,
    strip.position = "right",
    ncol = 2,
    labeller = label_function
  ) +
  xlab("Cell type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    strip.text.y = element_text(angle = 0),
    panel.spacing.y = unit(0.1, "inches")
  )
signature_box_plot_pdf <- paste0(
  plot_outdir,
  "/multi_class_signatures_meth_boxplots_by_category_cats_w_at_least_10_regions.pdf"
)
signature_box_plot_rds <- paste0(
  plot_rds_outdir,
  "/multi_class_signatures_meth_boxplots_by_category_cats_w_at_least_10_regions.rds"
)
pdf(signature_box_plot_pdf, width = 18, height = 8)
print(signature_box_plot)
dev.off()
saveRDS(signature_box_plot, file = signature_box_plot_rds)







###############################################################################
# Raw methylation heat map.
###############################################################################

# Add information on how many regions constitute each signature.
sig_reg_at_least_10$annotation_with_count <- facet_label_vect[sig_reg_at_least_10$automatic_annotation]

# Down-sample the Blood CCR8+ Treg Skin Tconv Skin Treg signature to 3,000 regions
# (Otherwise, plotting the heat map will not work).
corresp_lines <- which(
  sig_reg_at_least_10_gr$automatic_annotation == "Blood_CCR8+_Treg__Skin_Tconv__Skin_Treg__hypomethylation"
)
set.seed(2022)
lines_to_keep <- sample(corresp_lines, 3000)
lines_to_rm <- setdiff(corresp_lines, lines_to_keep)
sig_reg_at_least_10_downsampled_gr <- sig_reg_at_least_10_gr[-lines_to_rm]

# Generate a column annotation.
cell_types_by_sample <- meth_data_pdata$Cell_type
col_anno_df <- data.frame(
  Cell_type = factor(cell_types_by_sample, levels = cell_types),
  row.names = rownames(meth_data_pdata)
)
col_anno_obj <- HeatmapAnnotation(
  df = col_anno_df,
  col = list(Cell_type = cell_type_col_palette),
  annotation_label = "Cell type",
  annotation_legend_param = list(
    title = "Cell type",
    title_gp = gpar(fontsize = 17, fontface = "bold"),
    labels_gp = gpar(fontsize = 15),
    grid_height = unit(8, "mm"),
    grid_width = unit(8, "mm")
  )
)

# Specify a custom order of the signatures.
sig_reg_cats_at_least_10_sorted <- c(
  "Blood_CCR8+_Treg__Blood_naive_Treg__Skin_Treg__hypomethylation",
  "Blood_naive_Tconv__Skin_Tconv__hypomethylation",
  "Blood_CCR8+_Treg__Skin_Tconv__Skin_Treg__hypomethylation",
  "Blood_naive_Tconv__Blood_naive_Treg__hypomethylation",
  "Skin_Treg__hypomethylation",
  "Blood_CCR8+_Treg__Skin_Treg__hypomethylation",
  "Blood_naive_Tconv__Blood_naive_Treg__Skin_Tconv__hypomethylation",
  "Skin_Tconv__Skin_Treg__hypomethylation",
  "Blood_CCR8+_Treg__Blood_naive_Treg__Skin_Tconv__Skin_Treg__hypomethylation",
  "Blood_naive_Tconv__hypomethylation",
  "Blood_CCR8+_Treg__Blood_naive_Tconv__Skin_Tconv__Skin_Treg__hypomethylation",
  "Blood_naive_Treg__hypomethylation",
  "Blood_CCR8+_Treg__Blood_naive_Tconv__Blood_naive_Treg__Skin_Treg__hypomethylation",
  "Skin_Tconv__hypomethylation",
  "Blood_CCR8+_Treg__hypomethylation"
)

# Iterate over all signatures to generate methylation heat maps.
single_heat_maps <- lapply(
  sig_reg_cats_at_least_10_sorted,
  FUN = function(x) {
    #-- Get a GRanges object for this signature.
    this_signature_gr <- sig_reg_at_least_10_downsampled_gr[sig_reg_at_least_10_downsampled_gr$automatic_annotation == x]
    #-- Extract average raw methylation values for this signature.
    avg_raw_meth_this_signature <- getMeth(
      meth_data,
      regions = this_signature_gr,
      what = "perRegion",
      type = "raw"
    )
    #-- Generate a row annotation colour-coding the current signature.
    row_anno_df <- data.frame(Signature = rep(x, length(this_signature_gr)))
    row_anno_df$Signature <- factor(row_anno_df$Signature, levels = sig_reg_cats_at_least_10_sorted)
    row_anno_obj <- rowAnnotation(
      df = row_anno_df,
      col = list(Signature = sig_reg_cats_at_least_10_col_palette),
      show_annotation_name = ifelse(
        x == sig_reg_cats_at_least_10_sorted[length(sig_reg_cats_at_least_10)],
        yes = T,
        no = F
      ),
      annotation_legend_param = list(
        title_gp = gpar(fontsize = 17, fontface = "bold"),
        labels_gp = gpar(fontsize = 15),
        grid_height = unit(8, "mm"),
        grid_width = unit(8, "mm")
      )
    )
    #-- Generate a bar-row annotation depicting the number of regions in this signature.
    maximum_reg_num <- max(table(sig_reg_at_least_10_gr$automatic_annotation))
    region_num <- length(which(sig_reg_at_least_10_gr$automatic_annotation == x))
    region_num_anno_obj <- rowAnnotation(
      `log10(\nRegion\nnumber)` = anno_barplot(
        x = rep(log10(region_num), nrow(avg_raw_meth_this_signature)),
        bar_width = 1,
        gp = gpar(fill = "black"),
        ylim = c(0, log10(maximum_reg_num) + 0.5),
        axis = x == sig_reg_cats_at_least_10_sorted[length(sig_reg_cats_at_least_10_sorted)],
        width = unit(0.8, "inches")
      ),
      Placeholder = anno_block(
        labels = paste0(region_num, "\nRegions"),
        show_name = F,
        labels_rot = 0,
        labels_gp = gpar(fontsize = 9)
      ),
      show_annotation_name = x == sig_reg_cats_at_least_10_sorted[length(sig_reg_cats_at_least_10_sorted)],
      annotation_name_gp = gpar(fontsize = 12)
    )
    #-- Generate the heat map for this signature.
    this_signature_heatmap <- Heatmap(
      avg_raw_meth_this_signature,
      col = colorRamp2(seq(1, 0, length.out = 200), viridis(200)),
      cluster_rows = F,
      cluster_columns = F,
      show_row_names = F,
      show_column_names = F,
      heatmap_legend_param = list(
        title = "Methylation",
        labels_gp = gpar(fontsize = 15),
        title_gp = gpar(fontsize = 17, fontface = "bold")
      ),
      left_annotation = row_anno_obj,
      right_annotation = region_num_anno_obj,
      show_heatmap_legend = ifelse(
        x == sig_reg_cats_at_least_10_sorted[1],
        yes = T,
        no = F
      ),
      height = 1
    )
    return(this_signature_heatmap)
  }
)

# Merge the single heat maps into one and save this merged heat map.
indices <- 1:length(sig_reg_cats_at_least_10_sorted)
heatmap_drawing_string <- paste0(
  "draw(col_anno_obj %v% ",
  paste0("single_heat_maps[[", indices, "]]", collapse = " %v% "),
  ", merge_legends = T)"
)
signature_heatmap_outfile_pdf <- paste0(
  plot_outdir,
  "/multi_class_signatures_raw_meth_heatmap_cats_w_at_least_10_regions.pdf"
)
signature_heatmap_outfile_rds <- paste0(
  plot_rds_outdir,
  "/multi_class_signatures_raw_meth_heatmap_cats_w_at_least_10_regions.rds"
)
pdf(signature_heatmap_outfile_pdf,
    width = 16,
    height = 10)
eval(parse(text = heatmap_drawing_string))
dev.off()
saveRDS(single_heat_maps, file = signature_heatmap_outfile_rds)





##################################################################
# Histograms showing the distribution of region lengths.
##################################################################

# Combine regions with a length greater than 3,000 into a common bin.
sig_reg_clipped_3000bp <- sig_reg
sig_reg_clipped_3000bp$width <- sig_reg_clipped_3000bp$end - sig_reg_clipped_3000bp$start + 1
sig_reg_clipped_3000bp$width <- sapply(
  sig_reg_clipped_3000bp$width,
  FUN = function(x) {
    if (x > 3000) {
      return(3750)
    } else {
      return(x)
    }
  }
)
sig_reg_at_least_10_clipped_3000bp <- sig_reg_clipped_3000bp[sig_reg_clipped_3000bp$automatic_annotation %in% sig_reg_cats_at_least_10, ]

# Combined for all signature categories.
reg_length_hist_combined <- ggplot(sig_reg_clipped_3000bp) +
  aes(x = width) +
  scale_y_continuous(name = "Count", expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(
    name = "Signature region length [bp]",
    expand = c(0, 0),
    limits = c(-25, 4025),
    breaks = c(0, 1000, 2000, 3000, 3750),
    labels = c("0", "1000", "2000", "3000", ">3000")
  ) +
  geom_histogram(fill = "black", binwidth = 50) +
  geom_vline(xintercept = 3025, linetype = "dotted") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))
reg_length_hist_combined_pdf <- paste0(plot_outdir,
                                       "/multi_class_signatures_region_length_histogram_combined.pdf")
reg_length_hist_combined_rds <- paste0(plot_rds_outdir,
                                       "/multi_class_signatures_region_length_histogram_combined.rds")
pdf(reg_length_hist_combined_pdf,
    width = 3,
    height = 3)
print(reg_length_hist_combined)
dev.off()
saveRDS(reg_length_hist_combined, file = reg_length_hist_combined_rds)

# Faceted by category (only considering categories with at least 10 regions).
label_function <- labeller(automatic_annotation = facet_label_vect)
reg_length_hist_separate <- ggplot(sig_reg_at_least_10_clipped_3000bp) +
  aes(x = width) +
  scale_y_continuous(name = "Count", expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(
    name = "Signature region length [bp]",
    expand = c(0, 0),
    limits = c(-25, 4025),
    breaks = c(0, 1000, 2000, 3000, 3750),
    labels = c("0", "1000", "2000", "3000", ">3000")
  ) +
  geom_histogram(fill = "black", binwidth = 50) +
  geom_vline(xintercept = 3025, linetype = "dotted") +
  facet_wrap(
    ~ automatic_annotation,
    scales = "free_y",
    ncol = 2,
    strip.position = "right",
    labeller = label_function
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    strip.text.y = element_text(angle = 0)
  )
reg_length_hist_separate_pdf <- paste0(
  plot_outdir,
  "/multi_class_signatures_region_length_histogram_separate_cats_w_at_least_10_regions.pdf"
)
reg_length_hist_separate_rds <- paste0(
  plot_rds_outdir,
  "/multi_class_signatures_region_length_histogram_separate_cats_w_at_least_10_regions.pdf"
)
pdf(reg_length_hist_separate_pdf,
    width = 18,
    height = 9)
print(reg_length_hist_separate)
dev.off()
saveRDS(reg_length_hist_separate, file = reg_length_hist_separate_rds)







##################################################################
# Histograms showing the distribution of CpG numbers in regions.
##################################################################

# Combine regions with more than 50 CpGs into a common bin.
sig_reg_clipped_50_cpgs <- sig_reg
sig_reg_clipped_50_cpgs$cpg_num <- sapply(
  sig_reg_clipped_50_cpgs$cpg_num,
  FUN = function(x) {
    if (x > 50) {
      return(55)
    } else {
      return(x)
    }
  }
)
sig_reg_at_least_10_clipped_50_cpgs <- sig_reg_clipped_50_cpgs[sig_reg_clipped_50_cpgs$automatic_annotation %in% sig_reg_cats_at_least_10, ]

# Combined for all signature categories.
cpg_num_hist_combined <- ggplot(sig_reg_clipped_50_cpgs) +
  aes(x = cpg_num) +
  scale_y_continuous(name = "Count", expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(
    name = "Number of CpG sites",
    expand = c(0, 0),
    limits = c(4.5, 60.5),
    breaks = c(10, 20, 30, 40, 50, 55),
    labels = c("10", "20", "30", "40", "50", " >50")
  ) +
  geom_histogram(fill = "black", binwidth = 1) +
  geom_vline(xintercept = 50.5, linetype = "dotted") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))
cpg_num_hist_combined_pdf <- paste0(plot_outdir,
                                    "/multi_class_signatures_cpg_num_histogram_combined.pdf")
cpg_num_hist_combined_rds <- paste0(plot_rds_outdir,
                                    "/multi_class_signatures_cpg_num_histogram_combined.rds")
pdf(cpg_num_hist_combined_pdf,
    width = 3,
    height = 3)
print(cpg_num_hist_combined)
dev.off()
saveRDS(cpg_num_hist_combined, file = cpg_num_hist_combined_rds)

# Faceted by signature category (only categories with at least 10 regions considered).
cpg_num_hist_separate <- ggplot(sig_reg_at_least_10_clipped_50_cpgs) +
  aes(x = cpg_num) +
  scale_y_continuous(name = "Count", expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(
    name = "Number of CpG sites",
    expand = c(0, 0),
    limits = c(4.5, 60.5),
    breaks = c(10, 20, 30, 40, 50, 55),
    labels = c("10", "20", "30", "40", "50", " >50")
  ) +
  geom_histogram(fill = "black", binwidth = 1) +
  geom_vline(xintercept = 50.5, linetype = "dotted") +
  facet_wrap(
    ~ automatic_annotation,
    scales = "free_y",
    ncol = 2,
    strip.position = "right",
    labeller = label_function
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    strip.text.y = element_text(angle = 0)
  )
cpg_num_hist_separate_pdf <- paste0(
  plot_outdir,
  "/multi_class_signatures_cpg_num_histogram_separate_cats_w_at_least_10_regions.pdf"
)
cpg_num_hist_separate_rds <- paste0(
  plot_rds_outdir,
  "/multi_class_signatures_cpg_num_histogram_separate_cats_w_at_least_10_regions"
)
pdf(cpg_num_hist_separate_pdf,
    width = 18,
    height = 9)
print(cpg_num_hist_separate)
dev.off()
saveRDS(cpg_num_hist_separate, file = cpg_num_hist_separate_rds)

# Print the number of regions with more than 50 or more than 100 CpG sites.
num_gt_50 <- length(which(sig_reg$cpg_num > 50))
num_gt_100 <- length(which(sig_reg$cpg_num > 100))
print(paste(
  num_gt_50,
  "out of",
  nrow(sig_reg),
  "signature regions contain more than 50 CpGs"
))
print(paste(
  num_gt_100,
  "out of",
  nrow(sig_reg),
  "signature regions contain more than 100 CpGs"
))





###################################################################################
# Histogram showing the distance of region centres from transcription start sites.
###################################################################################

# For each region, identify the distance to the closest transcription start site.
sig_reg$closest_tss_dist <- sapply(
  1:nrow(sig_reg),
  FUN = function(x) {
    this_chrom <- sig_reg$seqnames[x]
    potential_genes <- genes[genes$X.chrom == this_chrom, ]
    potential_tss_locs <- potential_genes$tss_loc
    potential_tss_locs_dists <- sig_reg$centre[x] - potential_tss_locs
    min_dist_index <- which.min(abs(potential_tss_locs_dists))
    if (potential_genes$strand[min_dist_index] == "+") {
      return(potential_tss_locs_dists[min_dist_index])
    } else {
      return((-1) * potential_tss_locs_dists[min_dist_index])
    }
  }
)

# Remove regions with distances above 10,000 from the visualisation.
# This was also done in the Delacher 2017 paper.
sig_reg_for_dist_hist <- sig_reg[abs(sig_reg$closest_tss_dist) <= 10000, ]
sig_reg_at_least_10_for_dist_hist <- sig_reg_for_dist_hist[sig_reg_for_dist_hist$automatic_annotation %in% sig_reg_cats_at_least_10, ]

# Histogram: combined for all signature categories.
tss_dist_hist_combined <- ggplot(sig_reg_for_dist_hist) +
  aes(x = closest_tss_dist) +
  scale_y_continuous(name = "Count", expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(name = "Distance between region centre\nand closest TSS [bp]", expand = c(0, 0)) +
  geom_histogram(fill = "black", bins = 100) +
  theme_classic() +
  theme(
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.margin = unit(c(0.5, 0.75, 0.5, 0.5), "cm")
  )
tss_dist_hist_combined_pdf <- paste0(plot_outdir,
                                     "/multi_class_signatures_tss_dist_histogram_combined.pdf")
tss_dist_hist_combined_rds <- paste0(plot_rds_outdir,
                                     "/multi_class_signatures_tss_dist_histogram_combined.rds")
pdf(tss_dist_hist_combined_pdf,
    width = 3,
    height = 3)
print(tss_dist_hist_combined)
dev.off()
saveRDS(tss_dist_hist_combined, file = tss_dist_hist_combined_rds)

# Histogram: faceted by signature category.
tss_dist_hist_separate <- ggplot(sig_reg_at_least_10_for_dist_hist) +
  aes(x = closest_tss_dist) +
  scale_y_continuous(name = "Count", expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(name = "Distance between region centre and closest TSS [bp]", expand = c(0, 0)) +
  geom_histogram(fill = "black", bins = 100) +
  facet_wrap(
    ~ automatic_annotation,
    scales = "free_y",
    ncol = 2,
    strip.position = "right",
    labeller = label_function
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    strip.text.y = element_text(angle = 0)
  )
tss_dist_hist_separate_pdf <- paste0(
  plot_outdir,
  "/multi_class_signatures_tss_dist_histogram_separate_cats_w_at_least_10_regions.pdf"
)
tss_dist_hist_separate_rds <- paste0(
  plot_rds_outdir,
  "/multi_class_signatures_tss_dist_histogram_separate_cats_w_at_least_10_regions.rds"
)
pdf(tss_dist_hist_separate_pdf,
    width = 18,
    height = 9)
print(tss_dist_hist_separate)
dev.off()
saveRDS(tss_dist_hist_separate, file = tss_dist_hist_separate_rds)







#############################################################################################
# Pie chart showing how signature region centres distribute across interesting gene regions.
#############################################################################################

assign_region_centre_to_predef_region <- function(regions,
                                                  predef_regions,
                                                  miscellaneous_category,
                                                  ncores) {
  # This function assigns region centres to predefined genomic regions. Region centres overlapping
  #    with predefined regions of several classes are termed "Ambiguous".
  # regions: Data frame containing the regions to assign. Must contain the columns
  #   "seqnames" and "centre".
  # predef_regions: Data frame containing the predefined genomic regions. Must contain
  #   the columns "chr", "start", "end" and "class".
  # miscellaneous_category: Character; A class to assign regionss whose centres do not
  #   overlap with any predefined region.
  # ncores: Numeric; The number of cores to use in parallel computing.
  # Dependency: parellel
  # Value: A character vector containing the assignments.
  # Author: Niklas Beumer
  
  # Iterate over all regions.
  assignments <- unlist(mclapply(
    1:nrow(regions),
    FUN = function(x) {
      # Get the current chromosome and centre.
      curr_chr <- regions$seqnames[x]
      region_centre <- regions$centre[x]
      
      # Identify genomic regions overlapping with the current centre.
      potential_predef_regions <- predef_regions[predef_regions$chr == curr_chr, ]
      overlap_regions <- potential_predef_regions[region_centre >= potential_predef_regions$start &
                                                    region_centre <= potential_predef_regions$end, ]
      
      # Depending on how many classes the overlapping regions have, determine the correct assignment.
      if (nrow(overlap_regions) > 0) {
        annotations <- unique(overlap_regions$class)
        if (length(annotations) == 1) {
          return(annotations)
        } else {
          return("Ambiguous")
        }
      } else {
        return(miscellaneous_category)
      }
    },
    mc.cores = ncores
  ))
  
  # Return the vector of assignments.
  return(assignments)
  
}

# Assign each signature region centre to a particular gene region identified above.
# If it is not located in any of these regions, use the classification "Outside of any region".
# If it overlaps with several regions of different categories (This can happen in the case
# of overlapping or highly proximal genes), classify a signature region centre as "Ambiguous".
sig_reg$region_location <- assign_region_centre_to_predef_region(
  regions = sig_reg,
  predef_regions = gene_regions,
  miscellaneous_category = "Outside of any region",
  ncores = 10
)
sig_reg$region_location <- factor(
  sig_reg$region_location,
  levels = c(
    "2kb before TSS",
    "First exon",
    "Second exon",
    "Third exon",
    "Remaining exons",
    "First intron",
    "Second intron",
    "Third intron",
    "Remaining introns",
    "Ambiguous",
    "Outside of any region"
  )
)
sig_reg_at_least_10 <- sig_reg[sig_reg$automatic_annotation %in% sig_reg_cats_at_least_10, ]

# Pie chart: Combined for all signature categories.
region_loc_pie_chart_combined <- ggplot(sig_reg) +
  aes(x = "A", fill = region_location) +
  scale_fill_manual(
    values = c(
      "mediumseagreen",
      "dodgerblue1",
      "dodgerblue2",
      "dodgerblue3",
      "dodgerblue4",
      "tomato1",
      "tomato2",
      "tomato3",
      "tomato4",
      "lightblue",
      "yellow"
    ),
    name = "Region"
  ) +
  geom_bar(position = "fill", colour = "black") +
  coord_polar(theta = "y", direction = -1) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )
region_loc_pie_chart_combined_pdf <- paste0(
  plot_outdir,
  "/multi_class_signatures_region_loc_pie_chart_gene_regions_combined.pdf"
)
region_loc_pie_chart_combined_rds <- paste0(
  plot_rds_outdir,
  "/multi_class_signatures_region_loc_pie_chart_gene_regions_combined.rds"
)
pdf(region_loc_pie_chart_combined_pdf,
    width = 5,
    height = 3)
print(region_loc_pie_chart_combined)
dev.off()
saveRDS(region_loc_pie_chart_combined, file = region_loc_pie_chart_combined_rds)

# Pie chart: Faceted by signature category
region_loc_pie_chart_separate <- ggplot(sig_reg_at_least_10) +
  aes(x = "A", fill = region_location) +
  scale_fill_manual(
    values = c(
      "mediumseagreen",
      "dodgerblue1",
      "dodgerblue2",
      "dodgerblue3",
      "dodgerblue4",
      "tomato1",
      "tomato2",
      "tomato3",
      "tomato4",
      "lightblue",
      "yellow"
    ),
    name = "Region"
  ) +
  geom_bar(position = "fill", colour = "black") +
  facet_wrap(
    ~ automatic_annotation,
    ncol = 2,
    strip.position = "right",
    labeller = label_function
  ) +
  coord_polar(theta = "y", direction = -1) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.text.y = element_text(angle = 0)
  )
region_loc_pie_chart_separate_pdf <- paste0(
  plot_outdir,
  "/multi_class_signatures_region_loc_pie_chart_gene_regions_separate_cats_w_at_least_10_regions.pdf"
)
region_loc_pie_chart_separate_rds <- paste0(
  plot_rds_outdir,
  "/multi_class_signatures_region_loc_pie_chart_gene_regions_separate_cats_w_at_least_10_regions.rds"
)
pdf(region_loc_pie_chart_separate_pdf,
    width = 18,
    height = 10)
print(region_loc_pie_chart_separate)
dev.off()
saveRDS(region_loc_pie_chart_separate, file = region_loc_pie_chart_separate_rds)






########################################################################################
# Pie chart showing how signature region centres distribute across CpG island regions.
########################################################################################

# Assign each signture region centre to a particular CpG island region determined above.
# If it is not located in any of these regions, use the classification "Outside of any region".
sig_reg$cpg_island_location_cpgIslandExt <- assign_region_centre_to_predef_region(
  regions = sig_reg,
  predef_regions = cpg_islands_shores_shelves,
  miscellaneous_category = "Outside of any region",
  ncores = 10
)

sig_reg$cpg_island_location_cpgIslandExt <- factor(
  sig_reg$cpg_island_location_cpgIslandExt,
  levels = c("CpG island", "Shore", "Shelf", "Outside of any region")
)
sig_reg_at_least_10 <- sig_reg[sig_reg$automatic_annotation %in% sig_reg_cats_at_least_10, ]

# Pie chart: Combined for all signature categories.
region_loc_pie_chart_cgi_combined <- ggplot(sig_reg) +
  aes(x = "A", fill = cpg_island_location_cpgIslandExt) +
  scale_fill_manual(
    values = c(
      "darkolivegreen4",
      "darkolivegreen3",
      "darkolivegreen1",
      "blue"
    ),
    name = "Region"
  ) +
  geom_bar(position = "fill", colour = "black") +
  coord_polar(theta = "y", direction = -1) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )
region_loc_pie_chart_cgi_combined_pdf <- paste0(
  plot_outdir,
  "/multi_class_signatures_region_loc_pie_chart_cpg_islands_from_cpgIslandExt_list_combined.pdf"
)
region_loc_pie_chart_cgi_combined_rds <- paste0(
  plot_rds_outdir,
  "/multi_class_signatures_region_loc_pie_chart_cpg_islands_from_cpgIslandExt_list_combined.rds"
)
pdf(region_loc_pie_chart_cgi_combined_pdf,
    width = 5,
    height = 3)
print(region_loc_pie_chart_cgi_combined)
dev.off()
saveRDS(region_loc_pie_chart_cgi_combined, file = region_loc_pie_chart_cgi_combined_rds)

# Pie chart: Faceted by signature category
region_loc_pie_chart_cgi_separate <- ggplot(sig_reg_at_least_10) +
  aes(x = "A", fill = cpg_island_location_cpgIslandExt) +
  scale_fill_manual(
    values = c(
      "darkolivegreen4",
      "darkolivegreen3",
      "darkolivegreen1",
      "blue"
    ),
    name = "Region"
  ) +
  geom_bar(position = "fill", colour = "black") +
  facet_wrap(
    ~ automatic_annotation,
    ncol = 2,
    strip.position = "right",
    labeller = label_function
  ) +
  coord_polar(theta = "y", direction = -1) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.text.y = element_text(angle = 0)
  )
region_loc_pie_chart_cgi_separate_pdf <- paste0(
  plot_outdir,
  "/multi_class_signatures_region_loc_pie_chart_cpg_islands_from_cpgIslandExt_list_separate_cats_w_at_least_10_regions.pdf"
)
region_loc_pie_chart_cgi_separate_rds <- paste0(
  plot_rds_outdir,
  "/multi_class_signatures_region_loc_pie_chart_cpg_islands_from_cpgIslandExt_list_separate_cats_w_at_least_10_regions.rds"
)
pdf(region_loc_pie_chart_cgi_separate_pdf,
    width = 18,
    height = 10)
print(region_loc_pie_chart_cgi_separate)
dev.off()
saveRDS(region_loc_pie_chart_cgi_separate, file = region_loc_pie_chart_cgi_separate_rds)







###############################################################################################################
# Save the overview over the signature regions together with the statistics and annotations computed for them.
###############################################################################################################

# Save the signature region list.
final_sig_reg_outfile <- paste0(
  location,
  "/differential_methylation/multi_class_signatures_regions_with_additional_statistics_and_annotations.txt"
)
write.table(
  sig_reg,
  file = final_sig_reg_outfile,
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)
