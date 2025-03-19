# This script examines scRNA/TCR-seq data of clone matched skin Treg cells, 
#   blood CCR8+ Treg cells and blood RA+ Treg cells and validates that the skin 
#   Treg cell gene expression signature gets partly lost in blood CCR8+ Treg 
#   cells.
# Author: Niklas Beumer



# Load required packages.
library(Seurat)
library(ggplot2)


# Define a location on /xxx.
b330_space <- "/xxx/internal/"
location <- paste0(b330_space, "nbeumer/hm_treg_bs_rgnsbg")
location_2 <- "/xxx/data/nbeumer"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/yyyy/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the Seurat objects containing scRNA-seq data with clone information.
# There is a meta data slot called 
# "is_skintreg_and_shares_clonotype_w_bloodccr8treg". If this is TRUE, the cell 
# is a skin Treg that is clonotype-matched with blood CCR8+ Tregs.
# There is a meta data slot called 
# "is_bloodccr8treg_and_shares_clonotype_w_skintreg". If this is TRUE, the cell
# is a blood CCR8+ Treg that is clonotype-matched with skin Tregs.
sc_data_files <- paste0(
  location_2, 
  c("/scrnatcr_donor6_seur_obj_w_clonotype_info_and_matching.rds",
    "/scrnatcr_donor7_seur_obj_w_clonotype_info_and_matching.rds")
)
sc_data <- lapply(sc_data_files, FUN = readRDS)
names(sc_data) <- c("Donor_6", "Donor_7")

# Read in the differential expression table between skin Tregs and blood naive
# Tregs together with blood CCR8+ Treg positionings.
diff_exp_file <- paste0(
  location, 
  "/treg_hierarchies/diff_expr_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_expr_discrep_genes_info.txt"
)
diff_exp <- read.table(diff_exp_file, header = T, stringsAsFactors = F, 
                       sep = "\t")






############################################################################
# Validate a loss of the skin Treg gene expression signature in blood CCR8+
# Tregs.
############################################################################

# Restrict single-cell data to clone-matched skin Tregs and blood CCR8+ Tregs.
sc_data_relev <- lapply(sc_data, FUN = function(seur_obj) {
  seur_subs <- subset(
    seur_obj, 
    subset = is_skintreg_and_shares_clonotype_w_bloodccr8treg | 
      is_bloodccr8treg_and_shares_clonotype_w_skintreg
  )
  skin_annos <- seur_subs$is_skin_treg
  seur_subs$Cell_type <- sapply(1:ncol(seur_subs), FUN = function(x) {
    ifelse(skin_annos[x], yes = "Skin Treg", no = "Blood CCR8+ Treg")
  })
  return(seur_subs)
})

# Add bloodd naive Tregs that bear clonotypes that are shared between skin Tregs
# and blood CCR8+ Tregs.
sc_data_relev <- lapply(1:length(sc_data_relev), FUN = function(x) {
  all_relev_clonotypes <- unique(sc_data_relev[[x]]$cdr3s_nt)
  matching_blood_naive <- which(
    sc_data[[x]]$Anno_Niklas_for_TCR_matching == "Blood_RA_Tregs" & 
      sc_data[[x]]$cdr3s_nt %in% all_relev_clonotypes
  )
  matching_blood_naive_seur <- sc_data[[x]][, matching_blood_naive]
  matching_blood_naive_seur$Cell_type <- "Blood CD45RA+ Treg"
  new_seur_obj <- merge(sc_data_relev[[x]], matching_blood_naive_seur)
  return(new_seur_obj)
})

# Extract skin Treg hyperexpression genes and blood naive Treg hyperexpression 
# genes.
skin_hyperexpr_genes <- diff_exp$Gene[diff_exp$log2FoldChange < 0]
blood_naive_hyperexpr_genes <- diff_exp$Gene[diff_exp$log2FoldChange > 0]

# Compute module scores for the skin Treg hyperexpression and blood naive Treg 
# hyperexpression genes in the clone-matched cells.
sc_data_relev <- lapply(sc_data_relev, FUN = function(seur_obj) {
  seur_obj <- AddModuleScore(seur_obj, 
                             features = list(skin_hyperexpr_genes),
                             name = "Skin_Treg_hyperexpr")
  seur_obj <- AddModuleScore(seur_obj, 
                             features = list(blood_naive_hyperexpr_genes),
                             name = "Blood_naive_Treg_hyperexpr")
  return(seur_obj)
})
names(sc_data_relev) <- names(sc_data)

# Save the Seurat objects with computed module scores.
void <- lapply(names(sc_data_relev), FUN = function(x) {
  file_name <- paste0(
    location_2, 
    "/scrnatcr_",
    tolower(sub("_", "", x)),
    "_seur_obj_w_clonotype_info_and_matching_w_skinnaivediffexp_mod_scores.rds"
  )
  saveRDS(sc_data_relev[[x]], file = file_name)
})

# Collect module scores.
module_scores_collected <- do.call(
  rbind, 
  lapply(names(sc_data_relev), FUN = function(x) {
    seur_obj <- sc_data_relev[[x]]
    mod_score_df <- data.frame(
      Donor = factor(rep(x, 2 * ncol(seur_obj)), levels = names(sc_data_relev)),
      Cell = rep(colnames(seur_obj), 2),
      Cell_type = rep(seur_obj$Cell_type, 2),
      Module_score_type = rep(
        c("Skin Treg hyperexpression", "Blood naive Treg hyperexpression"),
        each = ncol(seur_obj)
      ),
      Module_score = c(seur_obj$Skin_Treg_hyperexpr1, 
                       seur_obj$Blood_naive_Treg_hyperexpr1)
    )
    return(mod_score_df)
  })
)

# Visualise the computed module scores and save the resulting plot.
module_scores_collected$Cell_type <- factor(
  module_scores_collected$Cell_type,
  levels = c("Blood CD45RA+ Treg", "Skin Treg", "Blood CCR8+ Treg")
)
mod_score_plot <- ggplot() +
  scale_y_continuous(name = "Module score") +
  scale_fill_manual(
    breaks = c("Blood CD45RA+ Treg", "Skin Treg", "Blood CCR8+ Treg"),
    values = c("darkorchid1", "blue", "orange"),
                    name = "Cell type") +
  geom_violin(data = module_scores_collected,
              mapping = aes(x = Donor, y = Module_score, fill = Cell_type),
              draw_quantiles = c(0.25, 0.5, 0.75), 
              scale = "width",
              width = 0.7) +
  facet_wrap(~Module_score_type, nrow = 1, scales = "free_y") +
  xlab("Donor") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))
mod_score_pdf <- paste0(
  plot_outdir, 
  "/tcr_matched_skin_treg_blood_ccr8_blood_naive_mod_scores_skin_treg_expr_signature.pdf"
)
mod_score_rds <- paste0(
  plot_rds_outdir, 
  "/tcr_matched_skin_treg_blood_ccr8_blood_naive_mod_scores_skin_treg_expr_signature.rds"
)
pdf(mod_score_pdf, width = 8, height = 4)
print(mod_score_plot)
dev.off()
saveRDS(mod_score_plot, file = mod_score_rds)
