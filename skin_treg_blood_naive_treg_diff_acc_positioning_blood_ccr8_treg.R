# This R script looks at regions that are differentially accessible between skin
#   Treg cells and blood naive Treg cells and analyses where blood CCR8+ Treg
#   cells are positioned in this comparison.
# Author: Niklas Beumer



# Load required packages.
library(Seurat)
library(Signac)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(ggplot2)
library(GenomicRanges)
library(testit)


# Define a location on /yyy.
b330_space <- "/yyy/"
location <- paste0(b330_space, "nbeumer/hm_treg_bs_rgnsbg")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis",
                     format(Sys.time(), "%m-%d-%y"),
                     sep = "/")
if (!dir.exists(plot_outdir)) {
  dir.create(plot_outdir)
}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Specify the relevant cell types.
relevant_cell_types <- c("Blood naive Treg", "Blood CCR8+ Treg", "Skin Treg")
relevant_cell_types_atac_naming <- c("blood_naive_treg", "blood_ccr8_treg", "skin_treg")

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

# Read in the Seurat object containing scATAC-seq data for CD4+ T cells.
sc_data_cd4_file <- paste0(
  b330_space,
  "msimon/scATAC/final/seurat_human_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar.RDS"
)
sc_data_cd4 <- readRDS(sc_data_cd4_file)

# Restrict the scATAC-seq data to the relevant cell types and perform a new
# TF-IDF.
sc_data_cd4_restr <- subset(sc_data_cd4, subset = treg_tconv_annot %in% relevant_cell_types_atac_naming)
sc_data_cd4_restr <- RunTFIDF(sc_data_cd4_restr)

# To enable donor-level heat map visualisation, also subset scATAC-seq data
# to the relevant cell types and cells from the relevant donors.
# Afterwards, perform a new TF-IDF.
blood_donors <- c("1", "2")
skin_donors <- c("4", "5")
sc_data_cd4_restr_f_donorlevel <- subset(
  sc_data_cd4,
  subset = (treg_tconv_annot == "blood_naive_treg" &
              donor %in% blood_donors) |
    (treg_tconv_annot == "blood_ccr8_treg" &
       donor %in% blood_donors) |
    (treg_tconv_annot == "skin_treg" & donor %in% skin_donors)
)
sc_data_cd4_restr_f_donorlevel <- RunTFIDF(sc_data_cd4_restr_f_donorlevel)

# Specify the prefix of the files containing genomic positions of exons,
# transcription start sites, genes etc.
annotation_pref <-
  "/yyy/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/GencodeV19_"

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
genes_gr <- makeGRangesFromDataFrame(
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

# Identify 2kb-regions in front of transcription start sites.
genes$tss_loc <- tss_pos
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







########################################################################
# Assign the differential peaks to genes.
########################################################################

# Generate a GRanges object out of the differential peak locations.
diff_peaks_gr <- StringToGRanges(rownames(atac_reg))
seqlevelsStyle(diff_peaks_gr) <- "NCBI"

# Increase the start positions of the ATAC peaks by 1 nucleotide in order to make
# the corresponding intervals 1-based and closed.
start(diff_peaks_gr) <- start(diff_peaks_gr) + 1

# A region is assigned to a gene if it overlaps with the
# gene body or with the 2kb-region in front of the transcription start site.
regs_genes_overlaps <- findOverlaps(diff_peaks_gr, genes_gr, ignore.strand = T)
genes_query_hits <- from(regs_genes_overlaps)
genes_subject_hits <- to(regs_genes_overlaps)
regs_2kb_overlaps <- findOverlaps(diff_peaks_gr, regions_2kb_gr, ignore.strand = T)
twokb_query_hits <- from(regs_2kb_overlaps)
twokb_subject_hits <- to(regs_2kb_overlaps)
gene_names <- genes_gr$name
gene_names_2kb <- regions_2kb_gr$gene
gene_assignments <- c()
overlap_types <- c()
for (i in 1:length(diff_peaks_gr)) {
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
atac_reg$gene_assignments <- gene_assignments
atac_reg$gene_assignment_overlap_types <- overlap_types





########################################################################
# Compute cell-type-wise average accessibility values in the relevant
# regions.
########################################################################

# For each cell type and each region, compute the mean accessibility across the cells.
#-- Iterate over the three relevant cell types.
avg_acc_celltype_level <- sapply(
  relevant_cell_types_atac_naming,
  FUN = function(x) {
    #-- Extract normalised accessibility values for the differential peaks and cells from this
    #-- cell type.
    acc_values <- GetAssayData(
      subset(sc_data_cd4_restr, subset = treg_tconv_annot == x),
      assay = "scATAC_raw",
      slot = "data"
    )[rownames(atac_reg), ]
    #-- Compute mean values across the cells.
    acc_means <- rowMeans(acc_values)
    #-- Return values.
    return(acc_means)
  }
)

# Harmonise column names with the convention I used for methylation data.
colnames(avg_acc_celltype_level) <- relevant_cell_types







########################################################################
# Foro donor-level heat map visualisation, compute donor-wise average
# accessibility values in the relevant regions.
########################################################################

# For each cell type, each donor and each region, compute the mean accessibility
# across the cells.
#-- Iterate over the three relevant cell types.
avg_acc_donor_level <- do.call(cbind,
                               lapply(
                                 relevant_cell_types_atac_naming,
                                 FUN = function(x) {
                                   sc_data_cd4_this_cell_type <- subset(sc_data_cd4_restr_f_donorlevel, subset = treg_tconv_annot == x)
                                   this_cell_type_donors <- sort(unique(sc_data_cd4_this_cell_type$donor))
                                   norm_acc_data_this_celltype <- sapply(
                                     this_cell_type_donors,
                                     FUN = function(y) {
                                       sc_data_this_donor <- subset(sc_data_cd4_this_cell_type, subset = donor == y)
                                       norm_acc_data <- GetAssayData(sc_data_this_donor, assay = "scATAC_raw", slot = "data")[rownames(atac_reg), ]
                                       print("Matrix extracted")
                                       data_aggr <- rowMeans(norm_acc_data)
                                       return(data_aggr)
                                     }
                                   )
                                   colnames(norm_acc_data_this_celltype) <-
                                     paste0(x, "_donor_", this_cell_type_donors)
                                   return(norm_acc_data_this_celltype)
                                 }
                               ))
assert(all(
  rownames(avg_acc_donor_level) == rownames(avg_acc_celltype_level)
))






##########################################################################################
# Exclude regions for which the differential accessibility tendencies reported by Seurat
# do not match with what I see in my actual accessibility values.
###########################################################################################

# For each region, check whether the differential accessibility tendency reported by Seurat
# matches with what I see in the actual accessibility values.
exclude_bools <- sapply(
  1:nrow(atac_reg),
  FUN = function(x) {
    diff_tend_sign <- sign(atac_reg$avg_log2FC[x])
    skin_acc <- avg_acc_celltype_level[x, "Skin Treg"]
    blood_naive_acc <- avg_acc_celltype_level[x, "Blood naive Treg"]
    return(
      ifelse(
        diff_tend_sign == -1,
        yes = skin_acc < blood_naive_acc,
        no = skin_acc > blood_naive_acc
      )
    )
  }
)

# Identify regions to exclude.
regions_to_exclude <- rownames(atac_reg)[exclude_bools]

# Save information on which regions were excluded.
output_file <- paste0(
  location,
  "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_region_exclusions.txt"
)
if (file.exists(output_file)) {
  file.remove(output_file)
}
capture.output(cat(
  rep(
    "*************************************************************************************************\n",
    2
  )
), file = output_file, append = T)
capture.output(
  print(
    "Regions that were excluded because the diff. acc. tendencies reported by Seurat"
  ),
  file = output_file,
  append = T
)
capture.output(
  print("did not match with what I saw in the actual accessibility:"),
  file = output_file,
  append = T
)
capture.output(cat(
  rep(
    "*************************************************************************************************\n",
    2
  )
), file = output_file, append = T)
capture.output(cat("\n\n"), file = output_file, append = T)
capture.output(print(
  ifelse(length(regions_to_exclude) > 0, yes = regions_to_exclude, no = "No regions were excluded.")
), file = output_file, append = T)

# Perform the exclusions.
avg_acc_celltype_level <- avg_acc_celltype_level[!exclude_bools, ]
avg_acc_donor_level <- avg_acc_donor_level[!exclude_bools, ]
atac_reg <- atac_reg[!exclude_bools, ]







##############################################################################
# Scale the accessibility data so that they range from 0 to 1
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
avg_acc_celltype_level_scaled <-
  scale_matr_lines_betw_0_and_1(avg_acc_celltype_level)
avg_acc_donor_level_scaled <-
  scale_matr_lines_betw_0_and_1(avg_acc_donor_level)






#########################################################################
# For each region, compute whether blood CCR8+ Tregs are closer to
# skin Tregs or to blood naive Tregs.
#########################################################################

# Quantify distances in scaled accessibility between blood CCR8+ Tregs and
# the two extreme cell types.
ccr8_naive_diffs <- avg_acc_celltype_level_scaled[, "Blood naive Treg"] -
  avg_acc_celltype_level_scaled[, "Blood CCR8+ Treg"]
ccr8_naive_dists <- abs(ccr8_naive_diffs)
ccr8_skin_diffs <- avg_acc_celltype_level_scaled[, "Skin Treg"] -
  avg_acc_celltype_level_scaled[, "Blood CCR8+ Treg"]
ccr8_skin_dists <- abs(ccr8_skin_diffs)

# For each region, determine whether blood CCR8+ Tregs are closer to
# blood naive Tregs or closer to skin Tregs.
closest_extrema <- sapply(
  1:nrow(avg_acc_celltype_level_scaled),
  FUN = function(x) {
    if (ccr8_naive_dists[x] < ccr8_skin_dists[x]) {
      return("Blood naive Treg")
    } else if (ccr8_naive_dists[x] > ccr8_skin_dists[x]) {
      return("Skin Treg")
    } else {
      return("Exactly in middle")
    }
  }
)





#########################################################################
# Collect all the computed values and add peak widths.
#########################################################################

# Prepare peak names.
peak_names_df <- data.frame(Peak = rownames(atac_reg))

# Compute peak widths.
widths_df <- data.frame(Width = width(diff_peaks_gr))

# Prepare accessibility values.
colnames(avg_acc_celltype_level) <-
  paste0("Norm_acc__", colnames(avg_acc_celltype_level))
norm_acc_df <- as.data.frame(avg_acc_celltype_level)
colnames(avg_acc_celltype_level_scaled) <-
  paste0("Scaled_norm_acc__",
         colnames(avg_acc_celltype_level_scaled))
scaled_acc_df <- as.data.frame(avg_acc_celltype_level_scaled)

# Prepare results from distance analysis to extreme cell types.
dist_df <- data.frame(
  Scaled_norm_acc_diff_blood_ccr8_treg_blood_naive_treg = ccr8_naive_diffs,
  Scaled_norm_acc_dist_blood_ccr8_treg_blood_naive_treg = ccr8_naive_dists,
  Scaled_norm_acc_diff_blood_ccr8_treg_skin_treg = ccr8_skin_diffs,
  Scaled_norm_acc_dist_blood_ccr8_treg_skin_treg = ccr8_skin_dists
)
closest_extrema_df <- data.frame(Cell_type_closest_to_blood_ccr8_treg  = closest_extrema)

# Combine all information into a large data frame.
combined_results <- cbind(
  peak_names_df,
  atac_reg,
  widths_df,
  norm_acc_df,
  scaled_acc_df,
  dist_df,
  closest_extrema_df
)

# Save all computed values.
combined_results_outfile <- paste0(
  location,
  "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_acc_discrep_regions_info.txt"
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
# Plot heat maps showing scaled accessibility in the differentially
# accessible regions. Once, use cell-type-level values and once use
# donor-level values.
############################################################################

# Generate a version of the accessibility information, ordered by the difference between
# blood CCR8+ Tregs and skin Tregs.
order_to_use <- order(
  combined_results$Scaled_norm_acc_diff_blood_ccr8_treg_skin_treg *
    sapply(
      combined_results$avg_log2FC,
      FUN = function(x) {
        ifelse(x < 0, yes = -1, no = 1)
      }
    )
)
combined_results_ordered <- combined_results[order_to_use, ]
avg_acc_donor_level_scaled <- avg_acc_donor_level_scaled[order_to_use, ]


###############
# Heat map on the cell type level.

# Iterate over two large categories (decreasing accessibility during differentiation,
# increasing accessibility during differentiation).
for (category in c("<", ">")) {
  # Extract values for this category.
  combined_results_cat <- combined_results_ordered[eval(parse(
    text = paste0("combined_results_ordered$avg_log2FC", category, " 0")
  )), ]
  
  # Get the matrix for heatmap plotting.
  acc_matr <- as.matrix(combined_results_cat[, grep("Scaled_norm_acc__", colnames(combined_results_cat))])
  colnames(acc_matr) <- gsub("Scaled_norm_acc__", "", colnames(acc_matr))
  rownames(acc_matr) <- combined_results_cat$region_ID
  
  # Generate the colour function for the heatmaps.
  acc_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                            colors = viridis(200, option = "B"))
  
  # Generate a row annotation showing to which cell type blood CCR8+ Tregs are closest
  # for each region.
  combined_results_cat$anno_values <- sapply(
    combined_results_cat$Cell_type_closest_to_blood_ccr8_treg,
    FUN = function(x) {
      ifelse(x == "Skin Treg", yes = "Closer to skin Tregs", no = "Closer to blood naive Tregs")
    }
  )
  row_anno <- rowAnnotation(
    `Blood CCR8+ Treg position` =
      combined_results_cat$anno_values,
    col = list(
      `Blood CCR8+ Treg position` = c(
        `Closer to skin Tregs` = "blue",
        `Closer to blood naive Tregs` = "darkorchid1",
        `Exactly in middle` = "grey"
      )
    )
  )
  
  # Generate the heat map.
  acc_heatmap <- Heatmap(
    acc_matr,
    show_row_names = F,
    cluster_rows = F,
    cluster_columns = F,
    col = acc_col_fun,
    right_annotation = row_anno,
    heatmap_legend_param = list(
      title = "Accessibility",
      at = c(0, 1),
      labels = c("min", "max")
    ),
    column_names_max_height = unit(10, "cm"),
  )
  
  
  # Save the heat map.
  file_snip <- ifelse(category == "<", yes = "Skin_Treg__hyperaccessibility", no = "Blood_naive_Treg_hyperaccesibility")
  acc_heatmap_pdf <- paste0(
    plot_outdir,
    "/treg_hierarchies_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_acc_discrep_heatmap_",
    file_snip,
    ".pdf"
  )
  acc_heatmap_rds <- paste0(
    plot_rds_outdir,
    "/treg_hierarchies_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_acc_discrep_heatmap_",
    file_snip,
    ".rds"
  )
  pdf(acc_heatmap_pdf, width = 6, height = 10)
  draw(acc_heatmap)
  dev.off()
  saveRDS(acc_heatmap, file = acc_heatmap_rds)
  
}


###############
# Heat map on the donor level.

# Iterate over two large categories (decreasing accessibility during differentiation,
# increasing accessibility during differentiation).
for (category in c("<", ">")) {
  # Extract values for this category.
  acc_matr <- avg_acc_donor_level_scaled[eval(parse(
    text = paste0("combined_results_ordered$avg_log2FC", category, " 0")
  )), ]
  
  # Generate the colour function for the heatmaps.
  acc_col_fun <- colorRamp2(breaks = seq(0, 1, length.out = 200),
                            colors = viridis(200, option = "B"))
  
  # Generate a row annotation showing to which cell type blood CCR8+ Tregs are closest
  # for each region.
  anno_values <- sapply(
    combined_results_ordered$Cell_type_closest_to_blood_ccr8_treg[eval(parse(
      text = paste0("combined_results_ordered$avg_log2FC", category, " 0")
    ))],
    FUN = function(x) {
      ifelse(x == "Skin Treg", yes = "Closer to skin Tregs", no = "Closer to blood naive Tregs")
    }
  )
  row_anno <- rowAnnotation(
    `Blood CCR8+ Treg position` = anno_values,
    col = list(
      `Blood CCR8+ Treg position` = c(
        `Closer to skin Tregs` = "blue",
        `Closer to blood naive Tregs` = "darkorchid1",
        `Exactly in middle` = "grey"
      )
    )
  )
  
  # Generate the heat map.
  acc_heatmap <- Heatmap(
    acc_matr,
    show_row_names = F,
    cluster_rows = F,
    cluster_columns = F,
    col = acc_col_fun,
    right_annotation = row_anno,
    heatmap_legend_param = list(
      title = "Accessibility",
      at = c(0, 1),
      labels = c("min", "max")
    ),
    column_names_max_height = unit(10, "cm"),
  )
  
  
  # Save the heat map.
  file_snip <- ifelse(category == "<", yes = "Skin_Treg__hyperaccessibility", no = "Blood_naive_Treg_hyperaccesibility")
  acc_heatmap_pdf <- paste0(
    plot_outdir,
    "/treg_hierarchies_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_acc_discrep_heatmap_",
    file_snip,
    "_by_donor.pdf"
  )
  acc_heatmap_rds <- paste0(
    plot_rds_outdir,
    "/treg_hierarchies_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_acc_discrep_heatmap_",
    file_snip,
    "_by_donor.rds"
  )
  pdf(acc_heatmap_pdf, width = 6, height = 10)
  draw(acc_heatmap)
  dev.off()
  saveRDS(acc_heatmap, file = acc_heatmap_rds)
  
}





#####################################################################
# Plot histograms of the differences between blood CCR8+ Tregs
# and skin Tregs.
#####################################################################

# Iterate over two large categories (decreasing accessibility during differentiation,
# increasing accessibility during differentiation).
for (category in c("<", ">")) {
  # Extract differences between skin Tregs and blood CCR8+ Tregs.
  diffs <- combined_results$Scaled_norm_acc_diff_blood_ccr8_treg_skin_treg[eval(parse(text = paste0(
    "combined_results$avg_log2FC", category, " 0"
  )))]
  
  # Generate the histogram.
  histogram <- ggplot() +
    aes(x = diffs) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_histogram(binwidth = 0.02, fill = "black") +
    xlab("Scaled norm. acc. skin Tregs\nminus scaled norm. acc. blood CCR8+ Tregs") +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"))
  
  # Save the histogram.
  file_snip <- ifelse(category == "<", yes = "Skin_Treg__hyperaccessibility", no = "Blood_naive_Treg_hyperaccesibility")
  histogram_pdf <- paste0(
    plot_outdir,
    "/treg_hierarchies_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_acc_discrep_acc_diff_hist_",
    file_snip,
    ".pdf"
  )
  histogram_rds <- paste0(
    plot_rds_outdir,
    "/treg_hierarchies_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_acc_discrep_acc_diff_hist_",
    file_snip,
    ".rds"
  )
  pdf(histogram_pdf, width = 5, height = 4)
  print(histogram)
  dev.off()
  saveRDS(histogram, file = histogram_rds)
}





############################################################################################
# Count in how many regions blood CCR8+ Tregs are closer to skin Tregs / blood naive Tregs.
############################################################################################

# Specify the file name under which output from this section will be saved.
output_file <- paste0(
  location,
  "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_acc_discrep_regions_counts.txt"
)
if (file.exists(output_file)) {
  file.remove(output_file)
}

# Compute peak widths.
peak_widths <- width(diff_peaks_gr)
combined_results$width <- peak_widths

# Print a status message.
capture.output(cat(
  rep(
    "*************************************************************************************************\n",
    2
  )
), file = output_file, append = T)
capture.output(
  print(
    "Numbers of regions in which blood CCR8+ Tregs are closer to skin Tregs / blood naive Tregs:"
  ),
  file = output_file,
  append = T
)
capture.output(cat(
  rep(
    "*************************************************************************************************\n",
    2
  )
), file = output_file, append = T)
capture.output(cat("\n\n"), file = output_file, append = T)

# Iterate over two large categories (decreasing accessibility during differentiation,
# increasing accessibility during differentiation).
for (category in c("<", ">")) {
  # Print what category is currently quantified.
  capture.output(
    cat("############################################\n"),
    file = output_file,
    append = T
  )
  capture.output(print(
    ifelse(category == "<", yes = "Skin_Treg__hyperaccessibility", no = "Blood_naive_Treg_hyperaccesibility")
  ), file = output_file, append = T)
  capture.output(cat("\n"), file = output_file, append = T)
  
  # Restrict the data to those regions that correspond to the current category.
  combined_results_cat <- combined_results[eval(parse(text = paste0(
    "combined_results$avg_log2FC", category, " 0"
  ))), ]
  
  # Compute some statistics for all regions.
  total_reg_count <- nrow(combined_results_cat)
  total_base_count <- sum(combined_results_cat$Width)
  
  # Compute statistics for the identified region categories.
  closer_to_naive_reg_count <- length(
    which(
      combined_results_cat$Cell_type_closest_to_blood_ccr8_treg == "Blood naive Treg"
    )
  )
  closer_to_naive_reg_perc <- round(100 * closer_to_naive_reg_count / total_reg_count, 2)
  closer_to_naive_base_count <- sum(combined_results_cat$Width[combined_results_cat$Cell_type_closest_to_blood_ccr8_treg == "Blood naive Treg"])
  closer_to_naive_base_perc <- round(100 * closer_to_naive_base_count / total_base_count, 2)
  closer_to_skin_reg_count <- length(which(
    combined_results_cat$Cell_type_closest_to_blood_ccr8_treg == "Skin Treg"
  ))
  closer_to_skin_reg_perc <- round(100 * closer_to_skin_reg_count / total_reg_count, 2)
  closer_to_skin_base_count <- sum(combined_results_cat$Width[combined_results_cat$Cell_type_closest_to_blood_ccr8_treg == "Skin Treg"])
  closer_to_skin_base_perc <- round(100 * closer_to_skin_base_count / total_base_count, 2)
  
  # Print the results.
  capture.output(print(
    paste0(
      "Closer to blood naive Tregs: ",
      closer_to_naive_reg_count,
      " regions (",
      closer_to_naive_reg_perc,
      "%), ",
      closer_to_naive_base_count,
      " bases (",
      closer_to_naive_base_perc,
      "%)."
    )
  ), file = output_file, append = T)
  capture.output(print(
    paste0(
      "Closer to skin Tregs: ",
      closer_to_skin_reg_count,
      " regions (",
      closer_to_skin_reg_perc,
      "%), ",
      closer_to_skin_base_count,
      " bases (",
      closer_to_skin_base_perc,
      "%)."
    )
  ), file = output_file, append = T)
  capture.output(cat("\n\n"), file = output_file, append = T)
  
}
