# This script adds clonotype information to Seurat objects with scRNA/TCR-seq
#   data. It also performs the relevant clonotype matches between samples.
# Author: Niklas Beumer
# Run this script with 10 cores!!!



# Load required package(s).
library(Seurat)
library(testit)
library(parallel)


# Define a location on /xxx.
b330_space <- "/xxx/"
location <- paste0(b330_space, "nbeumer/hm_treg_bs_rgnsbg")
location_output <- "/xxx/data/nbeumer"

# Read in the Seurat objects of scRNA/TCR-seq data with cell type annotation.
sc_data_files <- paste0(
  location, 
  c("/treg_hierarchies/scrnatcr_donor6_seur_obj_w_anno_for_tcr_matching.rds",
    "/treg_hierarchies/scrnatcr_donor7_seur_obj_w_anno_for_tcr_matching.rds")
)
sc_data <- lapply(sc_data_files, FUN = readRDS)
names(sc_data) <- c("Donor_6", "Donor_7")

# Generate meta data slots in the Seurat objects highlighting skin Tregs.
# These have previously been identified in the Delacher 2021 paper.
# Donor 6: Cells from cluster 10 that originate from sample R3 are skin Tregs.
# Donor 7: Cells from cluster 8 that originate from sample R7 are skin Tregs.
sc_data[[1]]$is_skin_treg <- sc_data[[1]]$seurat_clusters == 10 & 
  sc_data[[1]]$orig.ident == "MD-scRNA-TCR_R3"
sc_data[[2]]$is_skin_treg <- sc_data[[2]]$seurat_clusters == 8 & 
  sc_data[[2]]$orig.ident == "MD-scRNA-TCR_R7"

# Specify the prefix of files that contain TCR clonotype information.
clonotype_pref <- paste0(b330_space,
                         "Regensburg/10X_scTCR_human_charles/results/merged")







#############################################################################
# Add available clonotype information to the Seurat objects.
#############################################################################

# Iterate over the two Seurat objects.
sc_data <- lapply(sc_data, FUN = function(seur_obj) {
  
  # Initialise the relevant meta data slots in the Seurat object.
  seur_obj$sample_specific_clonotype_id <- NA
  seur_obj$cdr3s_aa <- NA
  seur_obj$cdr3s_nt <- NA
  
  # Identify all the samples present in this Seurat object.
  all_samples <- unique(seur_obj$orig.ident)
  
  # Iterate over the samples to read in TCR information.
  for (sample in all_samples) {
    
    # Print a status message.
    print(paste0("Adding TCR data for sample ", sample, "."))
    
    # Read in the clonotypes.csv file and the filtered_contig_annotations.csv
    # file from the corresponding Cell Ranger output.
    cellranger_output_loc <- paste0(
      clonotype_pref,
      "/",
      sub("scRNA-TCR", "scTCR", sub("scRNA_", "scTCR_", sample)),
      "/outs"
    )
    clonotypes_file <- paste0(cellranger_output_loc, "/clonotypes.csv")
    clonotypes <- read.csv(clonotypes_file)
    contig_anno_file <- paste0(cellranger_output_loc, 
                               "/filtered_contig_annotations.csv")
    contig_anno <- read.csv(contig_anno_file)
    
    # Assert that each clonotype appears once in the clonotype table.
    assert(length(unique(clonotypes$clonotype_id)) == nrow(clonotypes))
    
    # Assert that each clonotype mentioned in the contig annotations table
    # is present in the clonotype table.
    assert(all(contig_anno$raw_clonotype_id %in% clonotypes$clonotype_id))
    
    # Assert that each cell barcode is only associated with a single clonotype
    # ID.
    barcodes_in_tab <- unique(contig_anno$barcode)
    assert(all(sapply(barcodes_in_tab, FUN = function(x) {
      length(
        unique(contig_anno$raw_clonotype_id[contig_anno$barcode == x])
      ) == 1
    })))
    
    # Assert that all cell barcodes in the contig annotation file end with "-1".
    assert(all(sapply(contig_anno$barcode, FUN = function(x) {
      substr(x, nchar(x) - 1, nchar(x)) == "-1"
    })))
    
    # For each barcode that is present in the contig annotation table and in 
    # the Seurat object, add relevant information to the Seurat object.
    # Pay attention to the fact that barcodes in the Seurat object use a 
    # slightly different naming convention compared to those in the contig 
    # annotation table.
    barcodes_seur_this_sample <- colnames(seur_obj)[
      seur_obj$orig.ident == sample
    ]
    barcodes_seur_this_sample_table_conv <- sapply(
      barcodes_seur_this_sample, 
      FUN = function(x) {
        paste0(strsplit(x, split = "_")[[1]][1], "-1")
      })
    barcodes_for_loop <- barcodes_seur_this_sample_table_conv[
      barcodes_seur_this_sample_table_conv %in% contig_anno$barcode
    ]
    barcodes_for_loop_seur <- barcodes_seur_this_sample[
      barcodes_seur_this_sample_table_conv %in% contig_anno$barcode
    ]
    clonotype_ids <- unlist(mclapply(barcodes_for_loop, 
                                     FUN = function(barcode) {
      unique(
        contig_anno$raw_clonotype_id[contig_anno$barcode == barcode]
      ) 
    }, mc.cores = 10))
    seur_obj$sample_specific_clonotype_id[barcodes_for_loop_seur] <- 
      clonotype_ids
    aa_seqs <- unlist(mclapply(barcodes_for_loop, FUN = function(barcode) {
      corresp_clonotype_id <- unique(
        contig_anno$raw_clonotype_id[contig_anno$barcode == barcode]
      )  
      corresp_aa_seq <- clonotypes$cdr3s_aa[
        clonotypes$clonotype_id == corresp_clonotype_id
      ]
      return(corresp_aa_seq)
    }, mc.cores = 10))
    seur_obj$cdr3s_aa[barcodes_for_loop_seur] <- aa_seqs
    nt_seqs <- unlist(mclapply(barcodes_for_loop, FUN = function(barcode) {
      corresp_clonotype_id <- unique(
        contig_anno$raw_clonotype_id[contig_anno$barcode == barcode]
      )  
      corresp_nt_seq <- clonotypes$cdr3s_nt[
        clonotypes$clonotype_id == corresp_clonotype_id
      ]
      return(corresp_nt_seq)
    }, mc.cores = 10))
    seur_obj$cdr3s_nt[barcodes_for_loop_seur] <- nt_seqs
    
  }
  
  # Return the updated Seurat object.
  return(seur_obj)
  
})





#############################################################################
# Add information regarding which cells are clonotype-matched with 
# (i) blood naive Tregs, (ii) blood CCR8+ Tregs and (iii) skin Tregs.
# Two cells are considered clonotype-matched if they have exactly the same
# nucleotide sequence string and this sequence string contains at least 
# one alpha chain and one beta chain.
#############################################################################

# Add information on which cell contains a clonotype that also appears in 
# skin Tregs.
sc_data <- lapply(sc_data, FUN = function(seur_obj) {
  skin_clonotypes <- unique(seur_obj$cdr3s_nt[
    seur_obj$is_skin_treg & !(is.na(seur_obj$cdr3s_nt))
  ])
  skin_clonotypes_w_a_and_b <- grep("TRA", 
                                    grep("TRB", 
                                         skin_clonotypes, 
                                         value = T), 
                                    value = T)
  seur_obj$has_clonotype_appearing_in_skin_treg <- 
    seur_obj$cdr3s_nt %in% skin_clonotypes_w_a_and_b
  return(seur_obj)
})

# Add information on which cell contains a clonotype that also appears in 
# blood CCR8+ Tregs.
sc_data <- lapply(sc_data, FUN = function(seur_obj) {
  ccr8_clonotypes <- unique(seur_obj$cdr3s_nt[
    seur_obj$Anno_Niklas_for_TCR_matching == "Blood_CCR8_Tregs" & 
      !(is.na(seur_obj$cdr3s_nt))
  ])
  ccr8_clonotypes_w_a_and_b <- grep("TRA", 
                                    grep("TRB", 
                                         ccr8_clonotypes, 
                                         value = T), 
                                    value = T)
  seur_obj$has_clonotype_appearing_in_blood_ccr8_treg <- 
    seur_obj$cdr3s_nt %in% ccr8_clonotypes_w_a_and_b
  return(seur_obj)
})

# Add information on which cell contains a clonotype that also appears in 
# blood naive Tregs.
sc_data <- lapply(sc_data, FUN = function(seur_obj) {
  naive_clonotypes <- unique(seur_obj$cdr3s_nt[
    seur_obj$Anno_Niklas_for_TCR_matching == "Blood_RA_Tregs" & 
      !(is.na(seur_obj$cdr3s_nt))
  ])
  naive_clonotypes_w_a_and_b <- grep("TRA", 
                                     grep("TRB", 
                                          naive_clonotypes, 
                                          value = T), 
                                     value = T)
  seur_obj$has_clonotype_appearing_in_blood_ra_treg <- 
    seur_obj$cdr3s_nt %in% naive_clonotypes_w_a_and_b
  return(seur_obj)
})






########################################################################
# Add meta data slots highlighting (i) skin Tregs sharing a clonotype
# with blood CCR8+ Tregs and (ii) blood CCR8+ Tregs sharing a clonotype
# with skin Tregs.
########################################################################

# Add the aforementioned meta data slots.
sc_data <- lapply(sc_data, FUN = function(seur_obj) {
  seur_obj$is_skintreg_and_shares_clonotype_w_bloodccr8treg <-
    seur_obj$is_skin_treg & seur_obj$has_clonotype_appearing_in_blood_ccr8_treg
  seur_obj$is_bloodccr8treg_and_shares_clonotype_w_skintreg <-
    seur_obj$Anno_Niklas_for_TCR_matching == "Blood_CCR8_Tregs" & 
    seur_obj$has_clonotype_appearing_in_skin_treg
  return(seur_obj)
})






########################################################################
# Save the Seurat objects with the new information.
########################################################################

# Save Seurat objects.
void <- lapply(names(sc_data), FUN = function(x) {
  print(x)
  file_name <- paste0(
    location_output, 
    "/scrnatcr_",
    tolower(sub("_", "", x)),
    "_seur_obj_w_clonotype_info_and_matching.rds"
  )
  saveRDS(sc_data[[x]], file = file_name)
})
