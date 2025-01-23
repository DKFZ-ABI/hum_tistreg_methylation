# This script analyses blood CCR8+ Treg cell positionings, focusing on
#   DMRs (between skin Treg cells and and blood naive Treg cells) that overlap
#   with differential peaks (so-called DMR-peak pairs).
# Author: Niklas Beumer.



# Load required packages.
library(Seurat)
library(Signac)
library(bsseq)
library(GenomicRanges)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(testit)


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

# Specify the naming conventions for the relevant cell types.
relevant_cell_types <- c("Blood naive Treg", "Blood CCR8+ Treg", "Skin Treg")
relevant_cell_types_atac_naming <- c("blood_naive_treg", "blood_ccr8_treg", 
                                     "skin_treg")

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
# Tregs. Add information regarding in which direction each peak was 
# differential.
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

# Read in the BSseq object containing the smoothed methylation data
# and restrict to the relevant cell types.
meth_data <- readRDS(
  paste(
    location,
    "/preprocessing_etc/quality_control_after_alignment/2022-03-14_bsseq_object_combined_all_samples_smoothed.rds",
    sep = "/"
  )
)
meth_data <- meth_data[, pData(meth_data)$Cell_type %in% relevant_cell_types]

# Read in the Seurat object containing scATAC-seq data of CD4+ cells.
sc_data_cd4_file <- paste0(
  b330_space,
  "yyy/scATAC/final/seurat_human_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar.RDS"
)
sc_data_cd4 <- readRDS(sc_data_cd4_file)

# Restrict the scATAC-seq data to the relevant cell types and perform a new 
# TF-IDF.
sc_data_cd4 <- subset(
  sc_data_cd4, subset = treg_tconv_annot %in% relevant_cell_types_atac_naming
)
sc_data_cd4 <- RunTFIDF(sc_data_cd4)

# Specify the prefix of the files containing genomic positions of exons, 
# transcription start sites, genes etc.
annotation_pref <- 
  "/yyy/hg19_GRCh37_1000genomes/databases/gencode/gencode19/GencodeV19_"

# Read in the files containing gene annotations for normal chromosomes.
genes <- read.table(
  paste0(annotation_pref, "Genes_plain.bed.gz"),
  header = T,
  stringsAsFactors = F,
  comment.char = ""
)

# Increase start and end positions by 1 nucleotide in order to convert intervals 
# to 1-based closed intervals.
genes$chromStart <- genes$chromStart + 1

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








############################################################################
# Identify intersections between DMRs and differential peaks.
############################################################################

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

# Find overlaps between intersecting peaks and gene bodies or 2kb-regions in front of transcription start sites.
atac_genes_overlaps <- findOverlaps(atac_reg_overlapping_gr, genes_gr, ignore.strand = T)
genes_query_hits <- from(atac_genes_overlaps)
genes_subject_hits <- to(atac_genes_overlaps)
atac_2kb_overlaps <- findOverlaps(atac_reg_overlapping_gr, regions_2kb_gr, ignore.strand = T)
twokb_query_hits <- from(atac_2kb_overlaps)
twokb_subject_hits <- to(atac_2kb_overlaps)

# For each overlapping peak, collect information on what gene it is assigned to and by what overlap type.
# A peak is assigned to a gene if it overlaps with the gene body or with the 2kb-region
# in front of the transcription start site.
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
# Compute cell-type-wise average raw methylation values in the relevant
# regions.
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






##################################################################
# Compute cell-type-wise chromatin accessibility values for
# the relevant peaks.
##################################################################

# For each cell type and each region, compute the mean accessibility across the cells.
#-- Iterate over the three relevant cell types.
avg_acc_celltype_level <- sapply(
  relevant_cell_types_atac_naming,
  FUN = function(x) {
    #-- Extract normalised accessibility values for the differential peaks and cells from this
    #-- cell type.
    acc_values <- GetAssayData(
      subset(sc_data_cd4, subset = treg_tconv_annot == x),
      assay = "scATAC_raw",
      slot = "data"
    )[atac_reg_overlapping_gr$peak_ID, ]
    #-- Compute mean values across the cells.
    acc_means <- rowMeans(acc_values)
    #-- Return values.
    return(acc_means)
  }
)

# Harmonise column names with the convention I used for methylation data.
colnames(avg_acc_celltype_level) <- relevant_cell_types



##############################################################################
# Identify regions where blood CCR8+ Tregs are closer to skin Tregs
# and regions where they are closer to blood naive Tregs.
##############################################################################

# Methylation level.
meth_blood_ccr8_position <- sapply(
  1:nrow(avg_meth_celltype_level),
  FUN = function(x) {
    naive_ccr8_dist <- abs(avg_meth_celltype_level[x, "Blood CCR8+ Treg"] - 
                             avg_meth_celltype_level[x, "Blood naive Treg"])
    skin_ccr8_dist <- abs(avg_meth_celltype_level[x, "Blood CCR8+ Treg"] - 
                            avg_meth_celltype_level[x, "Skin Treg"])
    if (naive_ccr8_dist > skin_ccr8_dist) {
      output <- "Closer to skin Tregs"
    } else if (naive_ccr8_dist < skin_ccr8_dist) {
      output <- "Closer to blood naive Tregs"
    } else {
      output <- "Exactly in middle"
    }
    return(output)
  }
)

# ATAC level.
acc_blood_ccr8_position <- sapply(
  1:nrow(avg_acc_celltype_level),
  FUN = function(x) {
    naive_ccr8_dist <- abs(avg_acc_celltype_level[x, "Blood CCR8+ Treg"] - 
                             avg_acc_celltype_level[x, "Blood naive Treg"])
    skin_ccr8_dist <- abs(avg_acc_celltype_level[x, "Blood CCR8+ Treg"] - 
                            avg_acc_celltype_level[x, "Skin Treg"])
    if (naive_ccr8_dist > skin_ccr8_dist) {
      output <- "Closer to skin Tregs"
    } else if (naive_ccr8_dist < skin_ccr8_dist) {
      output <- "Closer to blood naive Tregs"
    } else {
      output <- "Exactly in middle"
    }
    return(output)
  }
)



##############################################################################
# Scale the methylation and accessibility data so that they range from 0 to 1
# for each region.
##############################################################################

scale_matr_lines_betw_0_and_1 <- function(matr) {
  # This function scales the lines of a matrix so that values range between
  #   0 and 1.
  # matr: The matrix to scale.
  # Dependencies: none.
  # Value: The scaled matrix.
  # Author: Niklas Beumer
  
  # Iterate over the rows.
  scaled_matr <- t(apply(
    matr,
    1,
    FUN = function(x) {
      # Identify the minimum and maximum value.
      min_val <- min(x)
      max_val <- max(x)
      
      # Subtract the minimum.
      x <- x - min_val
      
      # Divide values by the difference between the maximum and the minimum.
      divisor <- max_val - min_val
      x <- x / divisor
      
      # Return the scaled values.
      return(x)
      
    }
  ))
  
  # Return the scaled matrix.
  return(scaled_matr)
  
}

# Scale the values between 0 and 1.
avg_meth_celltype_level_scaled <-
  scale_matr_lines_betw_0_and_1(avg_meth_celltype_level)
atac_values_celltype_level_scaled <-
  scale_matr_lines_betw_0_and_1(avg_acc_celltype_level)






#########################################################################
# Collect all the computed values.
#########################################################################

# Generate a data frame containing methylation region information.
meth_region_df <- as.data.frame(meth_reg_overlapping_gr)[, c(7, 8, 11, 12)]
colnames(meth_region_df) <- paste0("Methylation__", colnames(meth_region_df))
atac_region_df <- as.data.frame(atac_reg_overlapping_gr)[, 6:9]
colnames(atac_region_df) <- paste0("Accessibility__", colnames(atac_region_df))

# Prepare methylation values.
colnames(avg_meth_celltype_level) <-
  paste0("Methylation__Raw_meth__",
         colnames(avg_meth_celltype_level))
raw_meth_df <- as.data.frame(avg_meth_celltype_level)
colnames(avg_meth_celltype_level_scaled) <-
  paste0("Methylation__Scaled_meth__",
         colnames(avg_meth_celltype_level_scaled))
scaled_meth_df <- as.data.frame(avg_meth_celltype_level_scaled)

# Prepare ATAC values.
colnames(avg_acc_celltype_level) <-
  paste0("Accessibility__Mean_norm_acc__",
         colnames(avg_acc_celltype_level))
mean_acc_df <- as.data.frame(avg_acc_celltype_level)
colnames(atac_values_celltype_level_scaled) <-
  paste0(
    "Accessibility__Scaled_mean_norm_acc__",
    colnames(atac_values_celltype_level_scaled)
  )
scaled_acc_df <- as.data.frame(atac_values_celltype_level_scaled)

# Prepare positionings of blood CCR8+ Tregs.
meth_blood_ccr8_positioning_df <- data.frame(Methylation__Blood_CCR8_Treg_pos = meth_blood_ccr8_position)
acc_blood_ccr8_positioning_df <- data.frame(Accessibility__Blood_CCR8_Treg_pos = acc_blood_ccr8_position)

# Combine all information into a large data frame.
combined_results <- cbind(
  meth_region_df,
  atac_region_df,
  raw_meth_df,
  scaled_meth_df,
  meth_blood_ccr8_positioning_df,
  mean_acc_df,
  scaled_acc_df,
  acc_blood_ccr8_positioning_df
)
colnames(combined_results) <- gsub(" ", "_", colnames(combined_results))

# Save all computed values.
combined_results_outfile <- paste0(
  location,
  "/treg_hierarchies/diff_meth_diff_acc_overlapping_dmr_peak_pairs_w_ccr8_positioning_w_tfidf.txt"
)
write.table(
  combined_results,
  file = combined_results_outfile,
  sep = "\t",
  row.names = F,
  col.names = T,
  quote = F
)






############################################################################
# Plot a heat map showing scaled methylation and scaled accessibility
# in the DMR-peak pairs.
############################################################################

# Iterate over the two methylation tendencies.
meth_tendencies <- unique(combined_results$Methylation__automatic_annotation)
for (tendency in meth_tendencies) {
  # Get the matrices for heatmap plotting.
  rows_to_keep <- which(combined_results$Methylation__automatic_annotation == tendency)
  meth_matr <- as.matrix(combined_results[rows_to_keep, grep("Scaled_meth", colnames(combined_results))])
  colnames(meth_matr) <- gsub("_", " ", gsub("Methylation__Scaled_meth__", "", colnames(meth_matr)))
  acc_matr <- as.matrix(combined_results[rows_to_keep, grep("Scaled_mean_norm_acc", colnames(combined_results))])
  colnames(acc_matr) <- gsub("_",
                             " ",
                             gsub(
                               "Accessibility__Scaled_mean_norm_acc__",
                               "",
                               colnames(acc_matr)
                             ))
  assert(nrow(meth_matr) == nrow(acc_matr))
  rownames(meth_matr) <- rownames(acc_matr) <- 1:nrow(meth_matr)
  
  # Generate a common row dendrogram for the two data levels.
  combined_matr <- cbind(meth_matr, acc_matr)
  combined_matr_dist <- dist(combined_matr)
  combined_matr_hclust <- hclust(combined_matr_dist)
  
  # Generate the colour functions for the heatmaps.
  meth_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                             colors = viridis(200, direction = -1))
  acc_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                            colors = viridis(200, option = "B"))
  
  # Generate a row annotation showing to which cell type blood CCR8+ Tregs are closest
  # for each region.
  row_anno_meth <- rowAnnotation(
    `Blood CCR8+ Treg position` = combined_results$Methylation__Blood_CCR8_Treg_pos[rows_to_keep],
    col = list(
      `Blood CCR8+ Treg position` = c(
        `Closer to skin Tregs` = "blue",
        `Closer to blood naive Tregs` = "darkorchid1",
        `Exactly in middle` = "grey"
      )
    )
  )
  row_anno_acc <- rowAnnotation(
    `Blood CCR8+ Treg position` = combined_results$Accessibility__Blood_CCR8_Treg_pos[rows_to_keep],
    col = list(
      `Blood CCR8+ Treg position` = c(
        `Closer to skin Tregs` = "blue",
        `Closer to blood naive Tregs` = "darkorchid1",
        `Exactly in middle` = "grey"
      )
    )
  )
  
  # Generate single heat maps.
  meth_heatmap <- Heatmap(
    meth_matr,
    cluster_rows = combined_matr_hclust,
    cluster_columns = F,
    col = meth_col_fun,
    right_annotation = row_anno_meth,
    heatmap_legend_param = list(
      title = "Methylation",
      at = c(0, 1),
      labels = c("min", "max")
    )
  )
  acc_heatmap <- Heatmap(
    acc_matr,
    show_row_names = F,
    cluster_columns = F,
    col = acc_col_fun,
    right_annotation = row_anno_acc,
    heatmap_legend_param = list(
      title = "Chrom. accessibility",
      at = c(0, 1),
      labels = c("min", "max")
    )
  )
  
  # Combine the two heat maps into a list.
  # Note: By default, this orders the rows of the second heat map
  # in the same manner as in the first heat map
  # (https://jokergoo.github.io/ComplexHeatmap-reference/book/a-list-of-heatmaps.html;
  # 09 May 2023)
  combined_heatmaps <- meth_heatmap + acc_heatmap
  
  # Save the combined heat map.
  combined_heatmaps_pdf <- paste0(
    plot_outdir,
    "/treg_hierarchies_diff_meth_diff_acc_dmr_peak_pairs_w_ccr8_positioning_heatmap_",
    tendency,
    "_w_tfidf.pdf"
  )
  combined_heatmaps_rds <- paste0(
    plot_rds_outdir,
    "/treg_hierarchies_diff_meth_diff_acc_dmr_peak_pairs_w_ccr8_positioning_heatmap_",
    tendency,
    "_w_tfidf.rds"
  )
  pdf(combined_heatmaps_pdf, width = 8, height = 10)
  draw(combined_heatmaps)
  dev.off()
  saveRDS(combined_heatmaps, file = combined_heatmaps_rds)
  
}