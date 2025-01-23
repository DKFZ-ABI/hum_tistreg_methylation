# This script examines scRNA/TCR-seq data of clone-matched skin Treg cells and 
#   blood CCR8+ Treg cells and validates that the skin Treg cell gene expression 
#   signature gets partly lost in blood CCR8+ Treg cells.
# Author: Niklas Beumer



# Load required packages.
library(Seurat)
library(ggplot2)


# Define a location on /xxx.
b330_space <- "/xxx/"
location <- paste0(b330_space, "nbeumer/hm_treg_bs_rgnsbg")
location_2 <- "/xxx/data/nbeumer"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/yyy/hm_treg_bs_rgnsbg/analysis", 
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
    return(mod_score_df )
  })
)

# Perform statistical tests (Two-tailed Wilcoxon rank sum tests) assessing 
# whether skin Tregs and blood CCR8+ Tregs display differences in the module 
# scores.
score_types <- c("Skin Treg hyperexpression", 
                 "Blood naive Treg hyperexpression")
wilcox_test_res <- lapply(names(sc_data_relev), FUN = function(x) {
  wilcox_test_res_this_donor <- lapply(score_types, FUN = function(y) {
    rel_mod_scores <- module_scores_collected[
      module_scores_collected$Donor == x &
        module_scores_collected$Module_score_type == y, 
    ]
    wilcox.test(
      rel_mod_scores$Module_score[
        rel_mod_scores$Cell_type == "Skin Treg"
      ],
      rel_mod_scores$Module_score[
        rel_mod_scores$Cell_type == "Blood CCR8+ Treg"
      ]
    )
  })
  names(wilcox_test_res_this_donor) <- score_types
  return(wilcox_test_res_this_donor)
})
names(wilcox_test_res) <- names(sc_data_relev)


# Save the results of the statistical tests.
wilcox_outfile <- paste0(
  location, 
  "/treg_hierarchies/clone_matched_skin_treg_blod_ccr8_treg_wilcox_test_res_on_skinnaivediffexp_mod_scores.rds"
)
saveRDS(wilcox_test_res, file = wilcox_outfile)

# Prepare the test results for visualisation.
test_res_for_plot <- data.frame(
  Module_score_type = rep(score_types, each = 2),
  x_pos = rep(c(1, 2), 2),
  y_pos = unlist(lapply(score_types, FUN = function(x) {
    rep(1.1 * 
      max(
        module_scores_collected$Module_score[
          module_scores_collected$Module_score_type == x
        ]
      ), each = 2)
  })),
  Text = do.call(c, lapply(score_types, FUN = function(x) {
    sapply(names(sc_data_relev), FUN = function(y) {
      paste0("P = ", signif(wilcox_test_res[[y]][[x]]$p.value, 3))
    })
  }))
)

# Visualise the computed module scores and save the resulting plot.
mod_score_plot <- ggplot() +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1)),
                     name = "Module score") +
  scale_fill_manual(breaks = c("Blood CCR8+ Treg", "Skin Treg"),
                    values = c("orange", "blue"),
                    name = "Cell type") +
  geom_violin(data = module_scores_collected,
              mapping = aes(x = Donor, y = Module_score, fill = Cell_type),
              draw_quantiles = c(0.25, 0.5, 0.75), 
              scale = "width",
              width = 0.7) +
  geom_text(data = test_res_for_plot,
            mapping = aes(label = Text, x = x_pos, y = y_pos)) +
  facet_wrap(~Module_score_type, nrow = 1, scales = "free_y") +
  xlab("Donor") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))
mod_score_pdf <- paste0(
  plot_outdir, 
  "/tcr_matched_skin_treg_blod_ccr8_treg_mod_scores_skin_treg_expr_signature.pdf"
)
mod_score_rds <- paste0(
  plot_rds_outdir, 
  "/tcr_matched_skin_treg_blod_ccr8_treg_mod_scores_skin_treg_expr_signature.rds"
)
pdf(mod_score_pdf, width = 8, height = 4)
print(mod_score_plot)
dev.off()
saveRDS(mod_score_plot, file = mod_score_rds)
