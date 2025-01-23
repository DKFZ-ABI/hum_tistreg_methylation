# This script generates Seurat objects with scRNA/TCR-seq data including my
#   cell type annotation. These objects can be used for TCR matching.
# Author: Niklas Beummer



# Load required package(s).
library(Seurat)


# Specify a location on /xxx.
b330_space <- "/xxx"

# Read in the Seurat objects with all cells.
seur_obj_donor6 <- readRDS(paste0(
  b330_space, 
  "/imbusch/d100/10x_scRNA_human/niklas/donor6_CD4_all_samples.rds"
))
seur_obj_donor7 <- readRDS(paste0(
  b330_space, 
  "/imbusch/d100/10x_scRNA_human/scRNA-seq_analysis_donor7_CD4_all_samples.rds"
))

# Read in the Seurat objects with annotations.
seur_objs_anno_donor6_snips <- c(
  "hm_treg_bs_rgnsbgscrnatcr_donor6_sample28_seur_obj_umap_anno.rds",
  "hm_treg_bs_rgnsbgscrnatcr_donor6_sample29_seur_obj_umap_anno.rds"
)
seur_objs_anno_donor6 <- lapply(seur_objs_anno_donor6_snips, FUN = function(x) {
  seur_file <- paste0(b330_space, "/nbeumer/", x)
  return(readRDS(seur_file))
})
seur_objs_anno_donor7_snips <- c(
  "hm_treg_bs_rgnsbgscrnatcr_donor7_sample24_seur_obj_umap_anno.rds",
  "hm_treg_bs_rgnsbgscrnatcr_donor7_sample25_seur_obj_umap_anno.rds"
)
seur_objs_anno_donor7 <- lapply(seur_objs_anno_donor7_snips, FUN = function(x) {
  seur_file <- paste0(b330_space, "/nbeumer/", x)
  return(readRDS(seur_file))
})

# Add my cell type annnotation to the larger Seurat objects for the two donors.
anno_slots_donor6 <- c("Cell_type_annotation_10_res_0.7", 
                       "Cell_type_annotation_20_res_0.9")
seur_obj_donor6$Anno_Niklas_for_TCR_matching <- "Other"
for (i in 1:length(seur_objs_anno_donor6)) {
  seur_obj <- seur_objs_anno_donor6[[i]]
  annos <- seur_obj[[anno_slots_donor6[i]]]
  seur_obj_donor6$Anno_Niklas_for_TCR_matching[rownames(annos)] <- annos[, 1]
} 
anno_slots_donor7 <- c("Cell_type_annotation_10_res_0.7", 
                       "Cell_type_annotation_30_res_0.8")
seur_obj_donor7$Anno_Niklas_for_TCR_matching <- "Other"
for (i in 1:length(seur_objs_anno_donor7)) {
  seur_obj <- seur_objs_anno_donor7[[i]]
  annos <- seur_obj[[anno_slots_donor7[i]]]
  seur_obj_donor7$Anno_Niklas_for_TCR_matching[rownames(annos)] <- annos[, 1]
} 

# From donor 6, there is a sample (sample 30) that was sorted for blood CCR8+
# Tregs. Add this information to the corresponding cell type annotation.
seur_obj_donor6$Anno_Niklas_for_TCR_matching[
  seur_obj_donor6$orig.ident == "MD-scRNA_R30"
] <- "Blood_CCR8_Tregs"

# Save the Seurat objects with the added cell type annotation.
donor6_outfile <- paste0(
  b330_space, 
  "/nbeumer/hm_treg_bs_rgnsbg/treg_hierarchies/scrnatcr_donor6_seur_obj_w_anno_for_tcr_matching.rds"
)
saveRDS(seur_obj_donor6, file = donor6_outfile)
donor7_outfile <- paste0(
  b330_space, 
  "/nbeumer/hm_treg_bs_rgnsbg/treg_hierarchies/scrnatcr_donor7_seur_obj_w_anno_for_tcr_matching.rds"
)
saveRDS(seur_obj_donor7, file = donor7_outfile)
