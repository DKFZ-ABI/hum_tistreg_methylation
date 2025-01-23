# This R script identifies cell type signature regions from a multi-class
#   comparisons.
# Author: Niklas Beumer



# Load required packages.
library(bsseq)
library(ggplot2)
library(testit)
library(qpdf)
library(parallel)
library(cowplot)



############################################################
# Perform necessary steps at the beginning of the analysis.
############################################################

# Define a location on /yyy.
b330_space <- "/yyy/"
location <- paste0(b330_space, "yyy/hm_treg_bs_rgnsbg")

# Define the paths to the locations of input and output.
input_location <- paste0(location,
                         "/preprocessing_etc/quality_control_after_alignment")
output_location <- paste0(location, "/differential_methylation")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis",
                     format(Sys.time(), "%m-%d-%y"),
                     sep = "/")
if (!dir.exists(plot_outdir)) {
  dir.create(plot_outdir)
}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Specify the prefix of the files containing genomic positions of exons, transcription start sites, genes etc.
annotation_pref <- "yyy/hg19_GRCh37_1000genomes/databases/gencode/gencode19/GencodeV19_"

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

# Read in the BSseq object containing the smoothed methylation data.
meth_data <- readRDS(
  paste(
    input_location,
    "/2022-03-14_bsseq_object_combined_all_samples_smoothed.rds",
    sep = "/"
  )
)
# meth_data <- meth_data[c(1:1000000, 10000000:11000000), ] # For debugging purposes

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

# Identify the cell types that are present in the data.
cell_types <- unique(pData(meth_data)$Cell_type)





##########################################
# Find cell type signature regions.
##########################################

compute_cpg_stats <- function(bsseq_obj, ncores) {
  # This function computes CpG-level statistics that can be used for an automatic CpG annotation.
  #   The computed statistics are (i) the largest gap in average (across replicates) smoothed
  #   methylation between two adjacent cell types; (ii) the ratio between the largest gap
  #   and the correponding second-to-largest gap; (iii) the ratio between the maximum (across
  #   cell types) within-cell-type standard eror of smoothed methylation and the largest gap
  #   mentioned under i; the cell types present in the "lower-methylation group" (i.e. the cell
  #   types below the largst gap).
  # ncores: Numeric; The number of cores to use in parallel computing.
  # Dependencies: bsseq, GenomicRanges, parallel
  # Value: A data frame containing the CpG locations together with their computed statistics.
  # Author: Niklas Beumer
  
  # Extract sample-level smoothed methylation values in the CpGs.
  smoothed_meth <- getMeth(bsseq_obj, type = "smooth", what = "perBase")
  
  # For each CpG and each cell type, compute average smoothed methylation estimates across the replicates.
  sample_metadata <- pData(bsseq_obj)
  cell_types_to_use <- unique(sample_metadata$Cell_type)
  avg_meth_celltypes <- sapply(
    cell_types_to_use,
    FUN = function(x) {
      celltype_subs <- smoothed_meth[, rownames(sample_metadata)[sample_metadata$Cell_type == x]]
      celltype_means <- unlist(mclapply(
        1:nrow(celltype_subs),
        FUN = function(y) {
          mean(celltype_subs[y, ], na.rm = T)
        },
        mc.cores = ncores
      ))
      return(celltype_means)
    }
  )
  
  # For each CpG, compute the maximum within-cell-type standard errors of smoothed methylation values.
  sterrs_celltypes <- sapply(
    cell_types_to_use,
    FUN = function(x) {
      celltype_subs <- smoothed_meth[, rownames(sample_metadata)[sample_metadata$Cell_type == x]]
      celltype_sterrs <- unlist(mclapply(
        1:nrow(celltype_subs),
        FUN = function(y) {
          values <- celltype_subs[y, ]
          nona_values <- values[!is.na(values)]
          sterr <- sd(nona_values) / sqrt(length(nona_values))
          return(sterr)
        },
        mc.cores = ncores
      ))
      return(celltype_sterrs)
    }
  )
  max_sterrs <- rowMaxs(sterrs_celltypes)
  
  # Compute CpG statistics.
  #-- Iterate over the CpGs.
  cpg_stats <- do.call(rbind, mclapply(
    1:nrow(bsseq_obj),
    FUN = function(x) {
      #-- Extract cell-type-wise smoothed methylation values for this region.
      cpg_values <- avg_meth_celltypes[x, ]
      
      #-- Extract the maximum within-cell type standard error of smoothed methylation.
      cpg_max_sterr <- max_sterrs[x]
      
      #-- Sort the cell types by their methylation.
      cpg_values_sorted <- sort(cpg_values)
      
      #-- Compute methylation difference between adjacent cell types.
      cpg_values_diffs <- diff(cpg_values_sorted)
      
      #-- Identify the largest difference between adjacent cell types.
      cutoff <- which.max(cpg_values_diffs)
      max_diff <- cpg_values_diffs[cutoff]
      
      #-- Identify the ratio between the largest difference and the second-to-largest difference.
      differences_sorted <- sort(cpg_values_diffs, decreasing = T)
      largest_secondlargest_ratio <- differences_sorted[1] / differences_sorted[2]
      
      #-- Identify the ratio between the maximum within-cell-type standard error and the
      #-- largst difference between adjacent cell types.
      max_sterr_max_diff_ratio <- cpg_max_sterr / max_diff
      
      #-- Identify the cell types that would be in the lower-methylation group.
      lower_meth_group <- cpg_values_sorted[1:cutoff]
      
      #-- Collect the annotation and the computed statistics.
      this_cpg_stats_df <- data.frame(
        max_gap = max_diff,
        largest_secondlargest_ratio = largest_secondlargest_ratio,
        max_sterr_max_gap_ratio = max_sterr_max_diff_ratio,
        lower_meth_group = paste(gsub(" ", "_", sort(
          names(lower_meth_group)
        )), collapse = "__")
      )
      return(this_cpg_stats_df)
      
    },
    mc.cores = 10
  ))
  
  # Extract locations of CpG sites.
  cpg_locs <- rowRanges(bsseq_obj)
  cpg_locs_df <- as.data.frame(cpg_locs)[, 1:3]
  
  # Concatenate all relevant data.
  colnames(avg_meth_celltypes) <- gsub(" ", "_", colnames(avg_meth_celltypes))
  colnames(avg_meth_celltypes) <- paste0("Avg_smoothed_meth__", colnames(avg_meth_celltypes))
  return_df <- cbind(cpg_locs_df, avg_meth_celltypes, cpg_stats)
  
  # Return the collected data.
  return(return_df)
  
  
}


aggregate_cpgs_into_regions_by_anno <- function(cpg_info,
                                                annotation_column,
                                                max_cpg_dist,
                                                ncores) {
  # This function aggregates CpGs into regions according to a discrete annotation. Neighbouring
  #    CpGs with the same annotation that are closer than a certain distance will be aggregated into regions.
  # cpg_info: Data frame containing information on the CpGs. Must contain the columns
  #   "seqnames", "start", "end" and the value of annotation_column.
  # annotation_column: Character; the name of a discrete column in cpg_info according
  #   to which aggregation shall be performed. CpGs that should not be part of regions
  #   can be marked by the value "None" in this column.
  # max_cpg_dist: Numeric; The maximum difference between CpG start positions two CpGs can
  #   have to be part of the same aggregated region. If they are further apart, they are
  #   assigned to different regions.
  # ncores: Numeric; The number of cores to use in parallel computing.
  # Dependencies: parallel, testit
  # Value: A data frame containing information on the aggregated regions.
  # Author: Niklas Beumer
  
  # Identify all chromosomes that are present in the data set.
  all_chroms <- unique(cpg_info$seqnames)
  
  # Iterate over all chromosomes.
  detected_regions <- do.call(rbind, mclapply(
    all_chroms,
    FUN = function(x) {
      # Extract data for the current chromosome.
      cpg_info_curr_chrom <- cpg_info[cpg_info$seqnames == x, ]
      
      # Assert that the CpG sites are ordered correctly.
      assert(all(
        cpg_info_curr_chrom$start == sort(cpg_info_curr_chrom$start)
      ))
      
      # Compute distances between neighbouring CpGs.
      distances_to_next <- diff(cpg_info_curr_chrom$start)
      
      # Identify CpGs that are not sufficiently close together to be in a common signature region.
      distances_to_next_too_large <- distances_to_next > max_cpg_dist
      
      # Initialise the data frame that will later contain the extracted regions.
      detected_regions_curr_chrom <- data.frame(
        seqnames = NA,
        start = NA,
        end = NA,
        cpg_num = NA,
        annotation = NA
      )
      
      # Start at the first CpG.
      curr_index <- 1
      
      # Iterate over all CpGs.
      num_cpgs_on_chrom <- nrow(cpg_info_curr_chrom)
      starts <- cpg_info_curr_chrom$start
      ends <- cpg_info_curr_chrom$end
      annotations <- cpg_info_curr_chrom[, annotation_column]
      while (curr_index <= nrow(cpg_info_curr_chrom)) {
        #-- Check whether this CpG belongs to a signature.
        if (annotations[curr_index] != "None") {
          #-- Identify the signature that this CpG belongs to.
          curr_region_anno <- annotations[curr_index]
          
          #-- Memorise the start position of the CpG.
          curr_region_start <- starts[curr_index]
          
          #-- Start a nested loop that is terminated once the end of the current region is reached.
          curr_index_internal <- curr_index
          curr_region_cpg_num <- 1
          while (curr_index_internal <= nrow(cpg_info_curr_chrom)) {
            #-- Check whether this CpG is not the last one on the current chromosome.
            not_last_on_cur_chrom <- curr_index_internal != num_cpgs_on_chrom
            
            #-- If this is the case, perform two additional checks.
            if (not_last_on_cur_chrom) {
              #-- Check whether the next CpG is more than the minimum distance away.
              next_cpg_outside_mindist <- distances_to_next_too_large[curr_index_internal]
              
              #-- Check whether the next CpG is of another signature as the current region.
              next_cpg_other_signature <- annotations[curr_index_internal + 1] != curr_region_anno
              
              #-- If this is the last CpG on the current chromosome, set the next_cpg_outside_mindist variable
              #-- to TRUE so that the next conditional statement will evaulate to FALSE.
            } else {
              next_cpg_outside_mindist <- T
              
            }
            
            
            #-- If one of the requirements are not met, perform the necessary steps to terminate this
            #-- region (Note: If the current CpG is te last one of the current chromosome,
            #-- next_cpg_outside_mindist will automatically be TRUE).
            if (next_cpg_outside_mindist | next_cpg_other_signature) {
              #-- set the end of the current region to this CpG
              curr_region_end <- ends[curr_index_internal]
              
              #-- Collect information on the current region.
              temp_df <- data.frame(
                seqnames = x,
                start = curr_region_start,
                end = curr_region_end,
                cpg_num = curr_region_cpg_num,
                annotation = curr_region_anno
              )
              detected_regions_curr_chrom <- rbind(detected_regions_curr_chrom, temp_df)
              
              #-- Set the running variable of the outer loop to the first CpG after the current region.
              curr_index <- curr_index_internal + 1
              
              #-- Break the loop.
              break()
              
              #-- If all requirements are met, keep the nested loop going.
            } else {
              #-- increment the running variable for the internal loop.
              curr_index_internal <- curr_index_internal + 1
              
              #-- Increase the CpG number by 1.
              curr_region_cpg_num <- curr_region_cpg_num + 1
              
            }
            
          }
          
          
          #-- If the current CpG does not belong to a signature, only perform the necessary steps to keep the loop going.
        } else {
          #-- Increment the running variable.
          curr_index <- curr_index + 1
          
        }
        
      }
      
      # Remove the placeholder line.
      detected_regions_curr_chrom <- detected_regions_curr_chrom[-1, ]
      
      # Return the list of detected regions on the current chromosome.
      return(detected_regions_curr_chrom)
      
    },
    mc.cores = 10
  ))
  
  # Restore the original name of the annotation column.
  colnames(detected_regions)[5] <- annotation_column
  
  # Return the list of aggregated regions.
  return(detected_regions)
  
}


log_reg_leave_1_out_cv_on_raw_meth <- function(bsseq_obj,
                                               regions_info,
                                               regions_gr,
                                               regions_subset = 1:length(regions_gr),
                                               sig_col,
                                               ncores) {
  # This function assesses how well logistic regression models can predict membership
  #   of a methylation sample to a certain group of cell types using leave-1-out cross-
  #   validation.
  # bsseq_obj: The BSseq object to retrieve methylation data from.
  # regions_info: Dataframe where each line corresponds to one region. Must contain at least
  #   one column containing the name of a signature. Entries in this column must follow the convention
  #   "Cell_type_1__Cell_type_2__...__hypomethylation" (Blank spaces in cell type names replaced by
  #   an underscore, cell types separated by two underscores each).
  # regions_gr: GRanges object containing the locations of regions. Must be in the same order
  #   as regions_info.
  # regions_subset: Numeric: Indices of regions that should only be considered.
  # sig_col: Character; Name of the column in regions_info containing the signature names. Entries in
  #   this column must follow the convention "Cell_type_1__Cell_type_2__...__hypomethylation"
  #   (Blank spaces in cell type names replaced by an underscore, cell types separated by two
  #   underscores each).
  # ncores: Numeric; The number of cores to use in parallel computing.
  # Dependencies: bsseq, parallel
  # Value: A data frame with one line for each round of leave-1-out cross-validation and each
  #   region. The data include details about the fitted models, predictions and assessments of
  #   whether the predictions are correct and how many samples are correctly predicted by each
  #   region.
  # Author: Niklas Beumer
  
  # Extract raw methylation values (region averages).
  avg_raw_meth_regions <- getMeth(bsseq_obj,
                                  regions = regions_gr,
                                  what = "perRegion",
                                  type = "raw")
  
  # Iterate over the regions.
  log_reg_res <- do.call(rbind,
                         mclapply(
                           regions_subset,
                           FUN = function(x) {
                             #-- Extract methylation values for this region.
                             avg_raw_meth_region <- avg_raw_meth_regions[x, ]
                             
                             #-- Identify the cell types that are hypomethylated in this region.
                             sig <- regions_info[x, sig_col]
                             sig_split <- strsplit(sig, split = "__")[[1]]
                             hypomethylated_cells <- gsub("_", " ", sig_split[1:(length(sig_split) - 1)])
                             group_1_samples <- unlist(lapply(
                               hypomethylated_cells,
                               FUN = function(y) {
                                 grep(y,
                                      names(avg_raw_meth_region),
                                      fixed = T,
                                      value = T)
                               }
                             ))
                             group_2_samples <- setdiff(names(avg_raw_meth_region), group_1_samples)
                             
                             #-- Iterate over all samples to perform leave-1-out cross-validation.
                             temp_df <- do.call(rbind, lapply(
                               names(avg_raw_meth_region),
                               FUN = function(y) {
                                 #-- Generate the training set.
                                 training_meth_vals <- avg_raw_meth_region[names(avg_raw_meth_region) != y]
                                 training_responses <- sapply(
                                   names(training_meth_vals),
                                   FUN = function(z) {
                                     ifelse(z %in% group_1_samples, yes = 0, no = 1)
                                   }
                                 )
                                 training_df <- data.frame(Meth = training_meth_vals, Resp = training_responses)
                                 
                                 #-- Fit the logistic model.
                                 log_model <- glm(Resp ~ Meth, data = training_df, family = "binomial")
                                 
                                 #-- Prepare the validation data set.
                                 validation_df <- data.frame(Meth = avg_raw_meth_region[y])
                                 
                                 #-- Predict group membership for the validation sample.
                                 validation_prediction <- predict(log_model, newdata = validation_df, type = "response")
                                 validation_outcome <- ifelse(
                                   validation_prediction < 0.5,
                                   yes = "group_1",
                                   no = ifelse(
                                     validation_prediction > 0.5,
                                     yes = "group_2",
                                     no = "Ambiguous"
                                   )
                                 )
                                 
                                 #-- Compare with the ground truth for this sample.
                                 validation_ground_truth <- ifelse(y %in% group_1_samples, yes = "group_1", no = "group_2")
                                 validation_result <- ifelse(
                                   validation_outcome == validation_ground_truth,
                                   yes = "Correct prediction",
                                   no = ifelse(
                                     validation_outcome == "Ambiguous",
                                     yes = "Ambiguous prediction",
                                     no = "Incorrect prediction"
                                   )
                                 )
                                 
                                 #-- Collect all the results from this round of leave-1-out cross-validation.
                                 temp_df_2 <- data.frame(
                                   seqnames = seqnames(regions_gr)[x],
                                   start = start(regions_gr)[x],
                                   end = end(regions_gr)[x],
                                   sig_name = sig,
                                   left_out_sample = gsub("\n", " ", y),
                                   model_intercept = log_model$coefficients["(Intercept)"],
                                   model_slope = log_model$coefficients["Meth"],
                                   predicted_value = validation_prediction,
                                   prediction_outcome = validation_outcome,
                                   ground_truth = validation_ground_truth,
                                   verdict = validation_result
                                 )
                                 colnames(temp_df_2)[4] <- sig_col
                                 
                                 #-- Return the collected results for this round of the cross-validation.
                                 return(temp_df_2)
                                 
                               }
                             ))
                             
                             #-- Quantify how many samples are correctly classified by this region.
                             temp_df$Number_correctly_classified_samples_this_region <- length(which(temp_df$verdict == "Correct prediction"))
                             
                             #-- Return the results for this region.
                             
                             return(temp_df)
                             
                           },
                           mc.cores = ncores
                         ))
  
}


multi_class_signature_extraction <- function(bsseq_obj,
                                             genes_gr,
                                             ncores,
                                             consider_cpgs,
                                             cpg_end_offset = 0,
                                             gap_cutoffs,
                                             largest_secondlargest_ratio_cutoffs,
                                             max_sderr_gap_ratio_cutoffs,
                                             cpg_stat_outfile_snip,
                                             max_cpg_dist,
                                             min_cpg_num,
                                             min_mean_cov,
                                             log_reg_max_num_misclassified,
                                             log_reg_outfile_snip,
                                             regions_outfile_snip,
                                             cats_for_dependencies_plots,
                                             dependencies_pdf_snip,
                                             dependencies_rds_snip) {
  # This function identifies cell type signatures using a multi-class comparison.
  #### General:
  # bsseq:obj: The BSseq object to use.
  # genes_gr: GRanges object containing locations of genes. Must contain the meta data
  #   column "name".
  # ncores: Numeric; The number of cores to use in parallel computing.
  #### CpG annotation:
  # consider_cpgs: Logical; Vector of values indicating whether to consider / not consider
  #   a CpG for the annotation. Must be in an order corresponding to the order of CpGs in
  #   bsseq_obj.
  # cpg_end_offset: Numeric; Number of bases to add to the "end" position of each CpG. This
  #   is useful if two complementary Cs of the same CpG were aggregated but the position of
  #   the CpG only includes the first C.
  # gap_cutoffs: Numeric; Vector of gap cut-offs to use. For each CpG, within-cell-type average
  #   smooth methylation values are computed and cell types are ordered. The maximum methylation
  #   difference that occurs between two adjacent cell types is compared with these thresholds.
  #   CpGs above this threshold are considered for signature regions.
  # largest_secondlargest_ratio_cutoffs: Numeric; Vector of cut-offs imposed on the ratio
  #   between the largest gap and the second-to-largest gap between two adjacent cell types.
  # max_sderr_gap_ratio_cutoffs: Numeric; Vector of cut-offs imposed on the standard-error-based
  #   CpG statistic. For each within-cell-type mean (see above), a standard error is computed.
  #   The largest standard error across the cell types is divided by the largest gap between two
  #   adjacent cell types. This statistic is then used for a filter using the specified cut-offs.
  # cpg_stat_outfile_snip: Character; Snippet that will be used in a text file where computed
  #   CpG statistics and annotations will be saved. The file name will be
  #   >>snippet<<_>>filtering_params<<_cpg_annotation.txt
  #### Aggregation into regions and region filtering:
  # max_cpg_dist. Numeric: The maximum distance between two adjacent CpGs so that they can be
  #   aggregated into the same region.
  # min_cpg_num: Numeric; The minimum number of CpGs a region must contain to be kept.
  # min_mean_cov: Numeric; The minimum mean coverage (across CpGs) a region must display in every
  #   sample to be kept. This will consider all CpGs, not only those specified by
  #   consider_cpgs.
  # log_reg_max_num_misclassified: Numeric:The maximum number of samples that are mis-classified
  #   during leave-1-out cross validation using logistic regression.
  # log_reg_outfile_snip: Character; A snippet used in file names for results from leave-1-out
  #   cross validation using logistic regression. File names will be
  #   >>snippet<<_>>filtering_params<<_log_reg_results_detailed.txt
  # regions_filtered_outfile_snip: Character; A snippet used in file names for (un)filtered
  #   region tables. File names will be >>snippet<<_>>filtering_params<<_signature_regions_(un)filtered.txt
  # cats_for_dependencies_plots: Character; Vector of signature categories for which this function will
  #   plot how the number of regions, the number of covered CpGs and the number of covered bases depends
  #   on the filtering parameters. Category names follow the convention
  #   "Cell_type_1__Cell_type_2__...__hypomethylation" (Blank spaces in cell type names replaced by
  #   an underscore, cell types separated by two underscores each).
  # dependencies_pdf_snip: Character; A snippet that will be used in file names for plots showing
  #   how the number of regions, covered CpG sites and covered based depende on the filtering parameters.
  #   File names will be >>snippet<<_>>metric<<_by_filtering_prms.pdf
  # dependencies_rds_snip: Character; A snippet that will be used in file names for plots showing
  #   how the number of regions, covered CpG sites and covered based depende on the filtering parameters.
  #   File names will be >>snippet<<_>>metric<<_by_filtering_prms.rds
  #### Dependencies etc.
  # Dependencies: bsseq, GenomicRanges (need not be loaded but must be installed), parallel,
  #   testit.
  # Value: A list of data frames containing filtered signature regions. Each element of the list
  #   corresponds to a specific combination of filtering parameters. The names of list elements
  #   indicate what parameters were used. In addition, the function saves text files showing the
  #   CpG annotations, the results of leave-1-out cross-validation using logistic regression and
  #   filtered and unfiltered regions. It also generates plots showing how the number of regions,
  #   the number of CpGs covered and the number of bases covered depends on the filtering
  #   parameters.
  # Author: Niklas Beumer
  
  # Compute CpG-level statistics for the annotation.
  cpg_statistics <- compute_cpg_stats(bsseq_obj = bsseq_obj[consider_cpgs, ], ncores = ncores)
  # cpg_statistics <- cpg_statistics[1:500000, ] # For debugging purposes.
  
  # Increase the CpG end positions by the specified numbers of nucleotides.
  cpg_statistics$end <- cpg_statistics$end + cpg_end_offset
  
  # Iterate over all conditions that should be tested.
  # Collect filtered signature regions in a list.
  filtered_regions <- list()
  for (gap_cutoff in gap_cutoffs) {
    for (largest_secondlargest_ratio_cutoff in largest_secondlargest_ratio_cutoffs) {
      for (max_sderr_gap_ratio_cutoff in max_sderr_gap_ratio_cutoffs) {
        # Generate strings identifying the current parameter set.
        param_str <- paste0(
          "gap_",
          gap_cutoff,
          "_largestsecondlargest_",
          largest_secondlargest_ratio_cutoff,
          "_sterrratio_",
          max_sderr_gap_ratio_cutoff
        )
        param_str_condensed <- paste0(
          "gp",
          gap_cutoff,
          "_lrgstscndlrgst",
          largest_secondlargest_ratio_cutoff,
          "_strrrtio",
          max_sderr_gap_ratio_cutoff
        )
        cat("\n")
        print(param_str)
        
        # Perform the automatic annotation of CpGs.
        cpg_anno <- unlist(mclapply(
          1:nrow(cpg_statistics),
          FUN = function(x) {
            if (cpg_statistics[x, "max_gap"] < gap_cutoff |
                cpg_statistics[x, "largest_secondlargest_ratio"] < largest_secondlargest_ratio_cutoff |
                cpg_statistics[x, "max_sterr_max_gap_ratio"] > max_sderr_gap_ratio_cutoff) {
              annotation_string <- "None"
            } else {
              annotation_string <- paste0(cpg_statistics[x, "lower_meth_group"], "__hypomethylation")
            }
          },
          mc.cores = ncores
        ))
        cpg_statistics$automatic_annotation <- cpg_anno
        
        # Save the annotations for the single CpGs.
        cpg_stat_outfile <- paste0(cpg_stat_outfile_snip,
                                   "_",
                                   param_str,
                                   "_cpg_annotation.txt")
        write.table(
          cpg_statistics,
          file = cpg_stat_outfile,
          sep = "\t",
          row.names = F,
          col.names = T,
          quote = F
        )
        
        # Aggregate the CpGs into signature regions.
        signature_regions <- aggregate_cpgs_into_regions_by_anno(
          cpg_info = cpg_statistics,
          annotation_column = "automatic_annotation",
          max_cpg_dist = max_cpg_dist,
          ncores = ncores
        )
        
        # Generate region IDs.
        signature_regions$region_ID <- paste0(param_str_condensed,
                                              "_region_",
                                              1:nrow(signature_regions))
        
        # Compute region widths.
        signature_regions$width <- signature_regions$end - signature_regions$start + 1
        
        # For each signature region, find the nearest gene and the distance to this gene.
        signature_regions_gr <- GenomicRanges::makeGRangesFromDataFrame(signature_regions)
        nearest_genes <- GenomicRanges::distanceToNearest(x = signature_regions_gr,
                                                          subject = genes_gr,
                                                          ignore.strand = T)
        signature_regions_hits <- from(nearest_genes)
        gene_hits <- to(nearest_genes)
        gene_names <- genes_gr$name
        distances <- nearest_genes@elementMetadata$distance
        signature_regions$nearest_gene <- unlist(mclapply(
          1:nrow(signature_regions),
          FUN = function(x) {
            genes <- gene_names[gene_hits[signature_regions_hits == x]]
            paste(genes, collapse = ", ")
          },
          mc.cores = ncores
        ))
        signature_regions$distance_to_nearest_gene <- unlist(mclapply(
          1:nrow(signature_regions),
          FUN = function(x) {
            value <- unique(distances[signature_regions_hits == x])
            assert(length(value) == 1)
            value
          },
          mc.cores = ncores
        ))
        
        # Annotate regions to genes. A region is assigned to a gene if it overlaps with the
        # gene body or with the 2kb-region in front of the transcription start site.
        regs_genes_overlaps <- findOverlaps(signature_regions_gr, genes_gr, ignore.strand = T)
        genes_query_hits <- from(regs_genes_overlaps)
        genes_subject_hits <- to(regs_genes_overlaps)
        regs_2kb_overlaps <- findOverlaps(signature_regions_gr, regions_2kb_gr, ignore.strand = T)
        twokb_query_hits <- from(regs_2kb_overlaps)
        twokb_subject_hits <- to(regs_2kb_overlaps)
        gene_names <- genes_gr$name
        gene_names_2kb <- regions_2kb_gr$gene
        gene_assignments <- c()
        overlap_types <- c()
        for (i in 1:length(signature_regions_gr)) {
          this_reg_genes_subject_hits <- genes_subject_hits[genes_query_hits == i]
          this_reg_genes_subject_names <- gene_names[this_reg_genes_subject_hits]
          this_reg_2kb_subject_hits <- twokb_subject_hits[twokb_query_hits == i]
          this_reg_2kb_subject_names <- gene_names_2kb[this_reg_2kb_subject_hits]
          gene_names_unique <- unique(c(
            this_reg_2kb_subject_names,
            this_reg_genes_subject_names
          ))
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
          gene_assignments <- c(gene_assignments,
                                paste(gene_names_unique, collapse = ", "))
          overlap_types <- c(overlap_types,
                             paste(overlap_type_str, collapse = ", "))
        }
        signature_regions$gene_assignments <- gene_assignments
        signature_regions$gene_assignment_overlap_types <- overlap_types
        
        # Filter the signature regions so that they contain at least the specified number of CpGs.
        cpg_num_filt_colname <- paste0("at_least_", min_cpg_num, "_cpgs")
        signature_regions[, cpg_num_filt_colname] <- signature_regions$cpg_num >= min_cpg_num
        signature_regions$retained <- signature_regions[, cpg_num_filt_colname]
        print(
          paste0(
            length(which(!(
              signature_regions$retained
            ))),
            " out of ",
            nrow(signature_regions),
            " (",
            round(100 * length(which(
              !(signature_regions$retained)
            )) / nrow(signature_regions), 2),
            "%) regions are removed because they contain less than ",
            min_cpg_num,
            " CpGs."
          )
        )
        
        # Keep regions that display an average (across CpGs) coverage of at least the
        # specified value in all samples. This is inspired by Delacher et al. (2017) Nat.
        # Immunol. 18, 1160-1172. # Use all CpGs in the data, not only those filtered
        # by consider_cpgs.
        current_regions <- which(signature_regions$retained)
        signature_granges <- GenomicRanges::makeGRangesFromDataFrame(signature_regions[current_regions, ])
        signature_coverages <- getCoverage(bsseq_obj,
                                           regions = signature_granges,
                                           type = "Cov",
                                           what = "perRegionAverage")
        insufficient_cov <- which(apply(
          signature_coverages,
          1,
          FUN = function(x) {
            !(all(x >= min_mean_cov))
          }
        ))
        regions_to_remove <- current_regions[insufficient_cov]
        print(
          paste0(
            length(regions_to_remove),
            " out of ",
            length(which(signature_regions$retained)),
            " (",
            round(100 * length(regions_to_remove) / length(
              which(signature_regions$retained)
            ), 2),
            "%) regions are removed due to the coverage>=",
            min_mean_cov,
            "-criterion."
          )
        )
        cov_filt_colname <- paste0("cov_at_least_", min_mean_cov, "_all_samples")
        signature_regions[, cov_filt_colname] <- "not analysed"
        signature_regions[current_regions, cov_filt_colname] <- T
        signature_regions[regions_to_remove, cov_filt_colname] <- F
        signature_regions$retained[regions_to_remove] <- F
        
        # Keep regions where computations on raw methylation would identify the same
        # lower-methylation group.
        current_regions <- which(signature_regions$retained)
        signature_granges <- GenomicRanges::makeGRangesFromDataFrame(signature_regions[current_regions, ])
        avg_raw_meth_samples <- getMeth(
          bsseq_obj,
          regions = signature_granges,
          type = "raw",
          what = c("perRegion")
        )
        sample_metadata <- pData(bsseq_obj)
        cell_types_to_use <- unique(sample_metadata$Cell_type)
        avg_raw_meth_celltypes <- sapply(
          cell_types_to_use,
          FUN = function(x) {
            celltype_subs <- avg_raw_meth_samples[, rownames(sample_metadata)[sample_metadata$Cell_type == x]]
            celltype_means <- unlist(mclapply(
              1:nrow(celltype_subs),
              FUN = function(y) {
                mean(celltype_subs[y, ], na.rm = T)
              },
              mc.cores = ncores
            ))
            return(celltype_means)
          }
        )
        raw_lower_meth_groups <- unlist(mclapply(
          1:nrow(avg_raw_meth_celltypes),
          FUN = function(x) {
            region_values <- avg_raw_meth_celltypes[x, ]
            region_values_sorted <- sort(region_values)
            region_values_diffs <- diff(region_values_sorted)
            cutoff <- which.max(region_values_diffs)
            lower_meth_group <- region_values_sorted[1:cutoff]
            return(paste(gsub(" ", "_", sort(
              names(lower_meth_group)
            )), collapse = "__"))
          },
          mc.cores = ncores
        ))
        regions_to_remove <- current_regions[paste0(raw_lower_meth_groups, "__hypomethylation") !=
                                               signature_regions$automatic_annotation[current_regions]]
        print(
          paste0(
            length(regions_to_remove),
            " out of ",
            length(which(signature_regions$retained)),
            " (",
            round(100 * length(regions_to_remove) / length(
              which(signature_regions$retained)
            ), 2),
            "%) regions are removed due to the criterion according to which raw methylation must identify the same lower-methylation group."
          )
        )
        raw_comp_colname <- paste0("raw_meth_identifies_same_lower_meth_group")
        signature_regions[, raw_comp_colname] <- "not analysed"
        signature_regions[current_regions, raw_comp_colname] <- T
        signature_regions[regions_to_remove, raw_comp_colname] <- F
        
        # Keep regions whose raw methylation correctly classifies at least the specified
        # number of samples in a leave-1-out cross-validation using logistic regression.
        #-- Perform the cross-validation.
        current_regions <- which(signature_regions$retained)
        log_reg_res <- log_reg_leave_1_out_cv_on_raw_meth(
          bsseq_obj = bsseq_obj[consider_cpgs, ],
          regions_info = signature_regions,
          regions_gr = signature_regions_gr,
          regions_subset = current_regions,
          sig_col = "automatic_annotation",
          ncores = ncores
        )
        #-- Save the detailed results from cross-validation.
        log_reg_res_outfile <- paste0(log_reg_outfile_snip,
                                      "_",
                                      param_str,
                                      "_log_reg_results_detailed.txt")
        write.table(
          log_reg_res,
          file = log_reg_res_outfile,
          row.names = F,
          col.names = T,
          sep = "\t",
          quote = F
        )
        #-- Collect information on which region classifies a sufficient number of samples correctly.
        correct_classes_by_region <- log_reg_res$Number_correctly_classified_samples_this_region[seq(1, nrow(log_reg_res), by = ncol(bsseq_obj))]
        #-- Remove regions with insufficient numbers of correctly classified samples.
        regions_to_remove <- current_regions[which(
          correct_classes_by_region < ncol(bsseq_obj) - log_reg_max_num_misclassified
        )]
        print(
          paste0(
            length(regions_to_remove),
            " out of ",
            length(which(signature_regions$retained)),
            " (",
            round(100 * length(regions_to_remove) / length(
              which(signature_regions$retained)
            ), 2),
            "%) regions are removed because logistic regression mis-classified more than ",
            log_reg_max_num_misclassified,
            " sample(s) incorrectly during leave-1-out CV."
          )
        )
        log_reg_filt_colname <- paste0(
          "correctly_classified_samples_log_reg_at_least_",
          ncol(bsseq_obj) - log_reg_max_num_misclassified
        )
        signature_regions[, log_reg_filt_colname] <- "not analysed"
        signature_regions[current_regions, log_reg_filt_colname] <- T
        signature_regions[regions_to_remove, log_reg_filt_colname] <- F
        signature_regions$retained[regions_to_remove] <- F
        
        # Save the information on signature regions with all statistics and the information
        # on whether they are retained during filtering.
        regions_unfiltered_outfile <- paste0(regions_outfile_snip,
                                             "_",
                                             param_str,
                                             "_signature_regions_unfiltered.txt")
        write.table(
          signature_regions,
          file = regions_unfiltered_outfile,
          sep = "\t",
          row.names = F,
          col.names = T,
          quote = F
        )
        
        # Save the information on signature regions passing all filtering steps.
        signature_regions_filtered <- signature_regions[signature_regions$retained, ]
        regions_filtered_outfile <- paste0(regions_outfile_snip,
                                           "_",
                                           param_str,
                                           "_signature_regions_filtered.txt")
        write.table(
          signature_regions_filtered,
          file = regions_filtered_outfile,
          sep = "\t",
          row.names = F,
          col.names = T,
          quote = F
        )
        
        # Add the filtered list of regions to the results list.
        filtered_regions[[param_str]] <- signature_regions_filtered
      }
    }
  }
  
  # Collect the numbers of regions, numbers of covered CpGs and numbers of covered bases
  # in the most interesting signature categories.
  region_nums_cpg_counts_base_counts <- do.call(rbind, lapply(
    gap_cutoffs,
    FUN = function(x) {
      do.call(rbind,
              lapply(
                largest_secondlargest_ratio_cutoffs,
                FUN = function(y) {
                  do.call(rbind,
                          lapply(
                            max_sderr_gap_ratio_cutoffs,
                            FUN = function(z) {
                              param_str <- param_str <- paste0("gap_",
                                                               x,
                                                               "_largestsecondlargest_",
                                                               y,
                                                               "_sterrratio_",
                                                               z)
                              results <- filtered_regions[[param_str]]
                              temp_df <- do.call(rbind,
                                                 lapply(
                                                   cats_for_dependencies_plots,
                                                   FUN = function(i) {
                                                     results_subs <- results[results$automatic_annotation == i, ]
                                                     region_num <- nrow(results_subs)
                                                     cpg_num <- sum(results_subs$cpg_num)
                                                     base_num <- sum(results_subs$end - results_subs$start + 1)
                                                     temp_df_2 <- data.frame(
                                                       gap_cutoff = x,
                                                       largest_secondlargest_ratio_cutoff = y,
                                                       max_sderr_gap_ratio_cutoff = z,
                                                       automatic_annotation = i,
                                                       region_count = region_num,
                                                       cpg_count = cpg_num,
                                                       base_count = base_num
                                                     )
                                                     return(temp_df_2)
                                                   }
                                                 ))
                              return(temp_df)
                            }
                          ))
                }
              ))
    }
  ))
  region_nums_cpg_counts_base_counts$gap_cutoff <- factor(region_nums_cpg_counts_base_counts$gap_cutoff)
  region_nums_cpg_counts_base_counts$largest_secondlargest_ratio_cutoff <- factor(region_nums_cpg_counts_base_counts$largest_secondlargest_ratio_cutoff)
  region_nums_cpg_counts_base_counts$max_sderr_gap_ratio_cutoff <- factor(region_nums_cpg_counts_base_counts$max_sderr_gap_ratio_cutoff)
  
  
  # Plot and save how the number of identified regions, the number of CpGs in these regions and the number of
  # bases in these regions depends on the parameters.
  y_vals <- paste0(c("region", "cpg", "base"), "_count")
  y_titles <- c("Region number",
                "CpG number in signature regions",
                "Base number in signature regions")
  single_file_snips <- paste0(c("region", "cpg", "base"), "_number")
  for (i in 1:length(y_vals)) {
    dependencies_plot <- ggplot(region_nums_cpg_counts_base_counts) +
      aes_string(
        x = "gap_cutoff",
        y = y_vals[i],
        colour = "largest_secondlargest_ratio_cutoff",
        shape = "max_sderr_gap_ratio_cutoff"
      ) +
      geom_point(position = position_jitter(
        width = 0.2,
        height = 0,
        seed = 2022
      )) +
      facet_wrap( ~ automatic_annotation,
                  ncol = 2,
                  scales = "free_y") +
      ylab(y_titles[i]) +
      theme_classic() +
      theme(
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.position = "top"
      )
    dependencies_plot_pdf <- paste0(dependencies_pdf_snip,
                                    "_",
                                    single_file_snips[i],
                                    "_by_filtering_prms.pdf")
    dependencies_plot_rds <- paste0(dependencies_rds_snip,
                                    "_",
                                    single_file_snips[i],
                                    "_by_filtering_prms.pdf")
    pdf(
      dependencies_plot_pdf,
      width = 12,
      height = 2.5 * ceiling(length(cats_for_dependencies_plots) / 2)
    )
    print(dependencies_plot)
    dev.off()
    saveRDS(dependencies_plot, file = dependencies_plot_rds)
  }
  
  # Return the list containing filtered regions for several filtering parameter combinations.
  return(filtered_regions)
  
}



# Extract coverages.
coverages <- getCoverage(meth_data, what = "perBase", type = "Cov")

# Keep CpGs that have a coverage of at least two in at least two samples of each cell type.
# This is the criterion applied to a pairwise comparison in the bsseq tutorial
# (https://www.bioconductor.org/packages/release/bioc/vignettes/bsseq/inst/doc/bsseq_analysis.html; 19 Sep 2022).
# 25,823,027 out of 27,999,538 CpG sites (92.2%) fulfil this criterion.
cell_type_wise_bools <- sapply(
  cell_types,
  FUN = function(x) {
    rowSums(coverages[, meth_data$Cell_type == x] >= 2) >= 2
  }
)
aggregated_bools <- apply(cell_type_wise_bools, 1, FUN = all)

# Perform the identification of signature regions using a multi-class comparison.
interesting_categories <- c(
  "Blood_CCR8+_Treg__Blood_naive_Treg__Skin_Treg__hypomethylation",
  "Blood_CCR8+_Treg__Skin_Treg__hypomethylation",
  "Skin_Treg__hypomethylation",
  "Blood_CCR8+_Treg__hypomethylation"
)
cell_type_signature_regions <- multi_class_signature_extraction(
  bsseq_obj = meth_data,
  genes_gr = genes_gr,
  ncores = 10,
  consider_cpgs = aggregated_bools,
  # Increase CpG end positions by 1 nucleotide to account for the fact that
  # I previously aggregated the two Cs on the complementary strands.
  cpg_end_offset = 1,
  gap_cutoffs = c(0.15, 0.2, 0.25),
  largest_secondlargest_ratio_cutoffs = c(1.25, 1.5, 2),
  max_sderr_gap_ratio_cutoffs = c(0.5, 0.25, 0.1),
  cpg_stat_outfile_snip = paste0(output_location, "/multi_class_signature_extraction"),
  max_cpg_dist = 300,
  min_cpg_num = 3,
  min_mean_cov = 2,
  log_reg_max_num_misclassified = 1,
  log_reg_outfile_snip = paste0(output_location, "/multi_class_signature_extraction"),
  regions_outfile_snip = paste0(output_location, "/multi_class_signature_extraction"),
  cats_for_dependencies_plots = interesting_categories,
  dependencies_pdf_snip = paste0(plot_outdir, "/multi_class_signature_extraction"),
  dependencies_rds_snip = paste0(plot_rds_outdir, "/multi_class_signature_extraction")
)









######################################################################
# Visualise signature regions in categories that are particularly
# interesting.
######################################################################


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

plot_meth_many_regions_cell_type_average_custom_titles <- function(bsseq_obj,
                                                                   regions_gr,
                                                                   extend = rep(5000, length(regions_gr)),
                                                                   genes,
                                                                   exons,
                                                                   tss_pos,
                                                                   highlight_regions = rep(T, length(regions_gr)),
                                                                   plot_titles,
                                                                   outfile_pdf,
                                                                   ncores) {
  # This function plots smoothed methylation values in several regions. The displayed values are within-cell-type averages.
  # bsseq_obj: The BSseq object to retrieve methylation data from. In the pData, this object must contain
  #   the columns "Cell_type" and "col".
  # regions_gr: A GRanges object containing the regions to plot.
  # extend: Numeric: The number of bases by which the regions will be extended on either side.
  #   (The plot will range from the first to the last CpG inside this extended region).
  # genes: A GRanges object containing information on human genes. Gene names have the
  #   meta data column name "name".
  # exons: A data.frame containing data on human exons. Must contain the columns
  #   "gene_name", "start" and "end".
  # tss_pos: A data.frame containing data on human transcription start sites. Must contain the columns
  #   "pos", "gene" and "strand".
  # highlight_regions: Logical, whether to highlight the regions of interest in the plots. Must be of
  #   the same length as regions_gr.
  # plot_titles: Character; The titles to use for the single plots. These should be in an order corresponding
  #   to the order of regions_gr.
  # outfile_pdf: Character; Path to a PDF file where the plot will be saved.
  # ncores: Numeric; The number of cores to use.
  # Dependencies: bsseq, ggplot2, cowplot, GenomicRanges, testit, qpdf, parallel
  # Value: None, the function just generates and saves the plots.
  # Author: Niklas Beumer
  
  
  # Generate a GRanges object that contains the regions of interest extended by the specified ranges.
  regions_extend_gr <- regions_gr + extend
  
  # Extract the cell types, sample names and colours corresponding to the samples in the BSseq object.
  cell_types_for_plot <- pData(bsseq_obj)$Cell_type
  cell_types_for_plot_unique <- unique(cell_types_for_plot)
  colours_for_plot <- unique(pData(bsseq_obj)$col)
  
  # Extract smoothed methylation values for the extended regions.
  smoothed_meth_in_regions <- getMeth(bsseq_obj,
                                      regions = regions_extend_gr,
                                      type = "smooth",
                                      what = "perBase")
  
  # Compute within-cell-type averages of smoothed methylation values.
  avg_smoothed_meth_in_regions <- lapply(
    smoothed_meth_in_regions,
    FUN = function(l) {
      avg_values <- unlist(lapply(
        cell_types_for_plot_unique,
        FUN = function(m) {
          corresp_values <- l[, cell_types_for_plot == m]
          apply(corresp_values, 1, FUN = mean, na.rm = T)
        }
      ))
      matrix(avg_values, ncol = length(cell_types_for_plot_unique))
    }
  )
  
  
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
          ncol(avg_smoothed_meth_in_regions[[l]])
        ),
        Meth = c(avg_smoothed_meth_in_regions[[l]]),
        Cell_type = rep(
          cell_types_for_plot_unique,
          each = nrow(avg_smoothed_meth_in_regions[[l]])
        )
      )
      meth_plotting_data$Cell_type <- factor(meth_plotting_data$Cell_type, levels = cell_types_for_plot_unique)
      
      # Prepare a data frame to indicate the region of interest in the plot.
      region_plotting_data <- data.frame(start = start(regions_gr[l]),
                                         end = end(regions_gr[l]))
      
      # Identify genes that overlap with the extended region of interest.
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
          name = "Smoothed\nmethylation\n(cell type average)",
          limits = c(0, 1),
          expand = expansion(mult = c(0.15, 0))
        ) +
        aes(x = Pos, y = Meth, colour = Cell_type) +
        scale_colour_manual(
          breaks = cell_types_for_plot_unique,
          values = colours_for_plot,
          name = "Cell type",
          drop = F,
          guide = guide_legend(
            override.aes = list(size = 6),
            # Make line legend appear as rectangles
            nrow = ifelse(
              length(cell_types_for_plot_unique) > 3,
              yes = 2,
              no = 1
            )
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
        gene_plot,
        ncol = 1,
        axis = "l",
        align = "v",
        rel_heights = c(2.5, 1.5)
      )
      
      # Save the combined plot at a preliminary location.
      pdf(paste0(outfile_pdf, "_", l, ".pdf"),
          width = 5.5,
          height = 5)
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


# Extract the set of regions that I will use later.
# I decided on using a gap cutoff of 0.15, a largest-secondlargest-ratio cutoff of 1.5
# and a sterr-ratio cutoff of 0.5.
final_parameter_set <- "gap_0.15_largestsecondlargest_1.5_sterrratio_0.5"
signature_regions_final <- cell_type_signature_regions[[final_parameter_set]]

# For each of the interesting categories, generate methylation track plots.
for (cat in interesting_categories) {
  cat_subs <- signature_regions_final[signature_regions_final$automatic_annotation == cat, ]
  cat_subs_gr <- makeGRangesFromDataFrame(cat_subs)
  plot_meth_many_regions_cell_type_average_custom_titles(
    bsseq_obj = meth_data,
    regions_gr = cat_subs_gr,
    genes = genes_gr,
    exons = exons,
    tss_pos = tss_df,
    plot_titles = cat_subs$region_ID,
    outfile_pdf = paste0(
      output_location,
      "/multi_class_signature_extraction_",
      final_parameter_set,
      "_",
      cat,
      ".pdf"
    ),
    ncores = 10
  )
}
