# This script optimises annotation plots for scRNA/TCR-seq data for our 
#   manuscript.
# Author: Niklas Beumer



# Load required package(s).
library(Seurat)
library(ggplot2)
library(viridis)


# Define a location on /xxx.
b330_space <- "/xxx"
location <- paste0(b330_space, "/nbeumer/hm_treg_bs_rgnsbg")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/yyy/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the Seurat objects with the annotation.
seur_obj_28 <- readRDS(paste0(
  b330_space,
  "/nbeumer/hm_treg_bs_rgnsbgscrnatcr_donor6_sample28_seur_obj_umap_anno.rds"
))
seur_obj_29 <- readRDS(paste0(
  b330_space,
  "/nbeumer/hm_treg_bs_rgnsbgscrnatcr_donor6_sample29_seur_obj_umap_anno.rds"
))
seur_obj_24 <- readRDS(paste0(
  b330_space,
  "/nbeumer/hm_treg_bs_rgnsbgscrnatcr_donor7_sample24_seur_obj_umap_anno.rds"
))
seur_obj_25 <- readRDS(paste0(
  b330_space,
  "/nbeumer/hm_treg_bs_rgnsbgscrnatcr_donor7_sample25_seur_obj_umap_anno.rds"
))



#############
# Sample 28 (donor 6)

# UMAP plot by cluster.
plot <- DimPlot(seur_obj_28, 
                reduction = "UMAP_10", 
                group.by = "Clusters_10_res_0.7") +
  theme(aspect.ratio = 1) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2)))
plot$labels$title <- NULL
outfile_pdf <- paste0(
  plot_outdir, 
  "/scrnatcr_donor6_sample28_umap_space_10_by_clusters_10_res_0.7_optim_for_publication.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir, 
  "/scrnatcr_donor6_sample28_umap_space_10_by_clusters_10_res_0.7_optim_for_publication.rds"
)
pdf(outfile_pdf, width = 10, height = 3)
print(plot)
dev.off()
saveRDS(plot, file = outfile_rds)

# UMAP plot by annotation.
plot_file <- paste0(
  plot_rds_outdir, 
  "/scrnatcr_donor6_sample28_umap_space_10_by_final_anno_10_res_0.7.rds"
)
plot <- readRDS(plot_file)
plot <- plot + 
  theme(aspect.ratio = 1)
outfile_pdf <- paste0(
  plot_outdir, 
  "/scrnatcr_donor6_sample28_umap_space_10_by_final_anno_10_res_0.7_optim_for_publication.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir, 
  "/scrnatcr_donor6_sample28_umap_space_10_by_final_anno_10_res_0.7_optim_for_publication.rds"
)
pdf(outfile_pdf, width = 10, height = 3)
print(plot)
dev.off()
saveRDS(plot, file = outfile_rds)

# UMAP plots by relevant expression values and module scores.
relevant_plots <- c(
  readRDS(paste0(
    plot_rds_outdir, 
    "/scrnatcr_donor6_sample28_umap_space_10_by_expr_general_markers.rds"
  ))[c(4, 10, 13, 17)],
  list(
    readRDS(paste0(
      plot_rds_outdir, 
      "/scrnatcr_donor6_sample28_umap_space_10_by_expr_bloodnaive_boodccr8_diffgenes_mod_scores.rds"
    ))[[1]],
    readRDS(paste0(
      plot_rds_outdir, 
      "/scrnatcr_donor6_sample28_umap_space_10_by_expr_bloodmem_boodccr8_diffgenes_mod_scores.rds"
    ))[[2]]
  )
)
relevant_plots[[5]]$labels$title <- "Blood CCR8 sig 1"
relevant_plots[[6]]$labels$title <- "Blood CCR8 sig 2"
for(plot in relevant_plots) {
  title <- plot$labels$title
  print(title)
  plot_new <- plot
  plot_new$labels$title <- NULL
  plot_new$labels$x <- "UMAP 1"
  plot_new$labels$y <- "UMAP 2"
  plot_new <- plot_new + 
    scale_colour_viridis_c(option = "G") +
    theme(aspect.ratio = 1)
  outfile_pdf <- paste0(plot_outdir, 
                        "/scrnatcr_donor6_sample28_umap_space_10_by_", 
                        gsub(" ", "_", title), 
                        "_optim_for_publication.pdf")
  outfile_rds <- paste0(plot_rds_outdir, 
                        "/scrnatcr_donor6_sample28_umap_space_10_by_", 
                        gsub(" ", "_", title), 
                        "_optim_for_publication.rds")
  pdf(outfile_pdf, width = 10, height = 3)
  print(plot_new)
  dev.off()
  saveRDS(plot_new, file = outfile_rds)
}



#############
# Sample 29 (donor 6)

# UMAP plot by cluster.
plot <- DimPlot(seur_obj_29, 
                reduction = "UMAP_20", 
                group.by = "Clusters_20_res_0.9") +
  theme(aspect.ratio = 1) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2)))
plot$labels$title <- NULL
outfile_pdf <- paste0(
  plot_outdir, 
  "/scrnatcr_donor6_sample29_umap_space_20_by_clusters_20_res_0.9_optim_for_publication.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir, 
  "/scrnatcr_donor6_sample29_umap_space_20_by_clusters_20_res_0.9_optim_for_publication.rds"
)
pdf(outfile_pdf, width = 10, height = 3)
print(plot)
dev.off()
saveRDS(plot, file = outfile_rds)

# UMAP plot by annotation.
plot_file <- paste0(
  plot_rds_outdir, 
  "/scrnatcr_donor6_sample29_umap_space_20_by_final_anno_20_res_0.9.rds"
)
plot <- readRDS(plot_file)
plot <- plot + 
  theme(aspect.ratio = 1)
outfile_pdf <- paste0(
  plot_outdir, 
  "/scrnatcr_donor6_sample29_umap_space_20_by_final_anno_20_res_0.9_optim_for_publication.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir, 
  "/scrnatcr_donor6_sample29_umap_space_20_by_final_anno_20_res_0.9_optim_for_publication.rds"
)
pdf(outfile_pdf, width = 10, height = 3)
print(plot)
dev.off()
saveRDS(plot, file = outfile_rds)

# UMAP plots by relevant expression values and moddule scores.
relevant_plots <- c(
  readRDS(paste0(
    plot_rds_outdir, 
    "/scrnatcr_donor6_sample29_umap_space_20_by_expr_general_markers.rds"
  ))[10],
  list(
    readRDS(paste0(
      plot_rds_outdir, 
      "/scrnatcr_donor6_sample29_umap_space_20_by_expr_bloodnaive_boodccr8_diffgenes_mod_scores.rds"
    ))[[1]],
    readRDS(paste0(
      plot_rds_outdir, 
      "/scrnatcr_donor6_sample29_umap_space_20_by_expr_bloodmem_boodccr8_diffgenes_mod_scores.rds"
    ))[[2]]
  )
)

# Optimise plot titles.
relevant_plots[[2]]$labels$title <- "Blood CCR8 sig 1"
relevant_plots[[3]]$labels$title <- "Blood CCR8 sig 2"
for(plot in relevant_plots) {
  title <- plot$labels$title
  print(title)
  plot_new <- plot
  plot_new$labels$title <- NULL
  plot_new$labels$x <- "UMAP 1"
  plot_new$labels$y <- "UMAP 2"
  plot_new <- plot_new + 
    scale_colour_viridis_c(option = "G") +
    theme(aspect.ratio = 1)
  outfile_pdf <- paste0(plot_outdir, 
                        "/scrnatcr_donor6_sample29_umap_space_20_by_", 
                        gsub(" ", "_", title), 
                        "_optim_for_publication.pdf")
  outfile_rds <- paste0(plot_rds_outdir, 
                        "/scrnatcr_donor6_sample29_umap_space_20_by_", 
                        gsub(" ", "_", title), 
                        "_optim_for_publication.rds")
  pdf(outfile_pdf, width = 10, height = 3)
  print(plot_new)
  dev.off()
  saveRDS(plot_new, file = outfile_rds)
}



#############
# Sample 24 (donor 7)

# UMAP plot by cluster.
plot <- DimPlot(seur_obj_24, 
                reduction = "UMAP_10", 
                group.by = "Clusters_10_res_0.7") +
  theme(aspect.ratio = 1) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2)))
plot$labels$title <- NULL
outfile_pdf <- paste0(
  plot_outdir, 
  "/scrnatcr_donor7_sample24_umap_space_10_by_clusters_10_res_0.7_optim_for_publication.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir, 
  "/scrnatcr_donor7_sample24_umap_space_10_by_clusters_10_res_0.7_optim_for_publication.rds"
)
pdf(outfile_pdf, width = 10, height = 3)
print(plot)
dev.off()
saveRDS(plot, file = outfile_rds)

# UMAP plot by annotation.
plot_file <- paste0(
  plot_rds_outdir, 
  "/scrnatcr_donor7_sample24_umap_space_10_by_final_anno_10_res_0.7.rds"
)
plot <- readRDS(plot_file)
plot <- plot + 
  theme(aspect.ratio = 1)
outfile_pdf <- paste0(
  plot_outdir, 
  "/scrnatcr_donor7_sample24_umap_space_10_by_final_anno_10_res_0.7_optim_for_publication.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir, 
  "/scrnatcr_donor7_sample24_umap_space_10_by_final_anno_10_res_0.7_optim_for_publication.rds"
)
pdf(outfile_pdf, width = 10, height = 3)
print(plot)
dev.off()
saveRDS(plot, file = outfile_rds)

# UMAP plots by relevant expression values and moddule scores.
relevant_plots <- c(
  readRDS(paste0(
    plot_rds_outdir, 
    "/scrnatcr_donor7_sample24_umap_space_10_by_expr_general_markers.rds"
  ))[c(4, 10, 13, 17)],
  list(
    readRDS(paste0(
      plot_rds_outdir, 
      "/scrnatcr_donor7_sample24_umap_space_10_by_expr_bloodnaive_boodccr8_diffgenes_mod_scores.rds"
    ))[[1]],
    readRDS(paste0(
      plot_rds_outdir, 
      "/scrnatcr_donor7_sample24_umap_space_10_by_expr_bloodmem_boodccr8_diffgenes_mod_scores.rds"
    ))[[2]]
  )
)
relevant_plots[[5]]$labels$title <- "Blood CCR8 sig 1"
relevant_plots[[6]]$labels$title <- "Blood CCR8 sig 2"
for(plot in relevant_plots) {
  title <- plot$labels$title
  print(title)
  plot_new <- plot
  plot_new$labels$title <- NULL
  plot_new$labels$x <- "UMAP 1"
  plot_new$labels$y <- "UMAP 2"
  plot_new <- plot_new + 
    scale_colour_viridis_c(option = "G") +
    theme(aspect.ratio = 1)
  outfile_pdf <- paste0(plot_outdir, 
                        "/scrnatcr_donor7_sample24_umap_space_10_by_", 
                        gsub(" ", "_", title), 
                        "_optim_for_publication.pdf")
  outfile_rds <- paste0(plot_rds_outdir, 
                        "/scrnatcr_donor7_sample24_umap_space_10_by_", 
                        gsub(" ", "_", title), 
                        "_optim_for_publication.rds")
  pdf(outfile_pdf, width = 10, height = 3)
  print(plot_new)
  dev.off()
  saveRDS(plot_new, file = outfile_rds)
}



#############
# Sample 25 (donor 7)

# UMAP plot by cluster.
plot <- DimPlot(seur_obj_25, 
                reduction = "UMAP_30", 
                group.by = "Clusters_30_res_0.8") +
  theme(aspect.ratio = 1) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2)))
plot$labels$title <- NULL
outfile_pdf <- paste0(
  plot_outdir, 
  "/scrnatcr_donor7_sample25_umap_space_30_by_clusters_30_res_0.8_optim_for_publication.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir, 
  "/scrnatcr_donor7_sample25_umap_space_30_by_clusters_30_res_0.8_optim_for_publication.rds"
)
pdf(outfile_pdf, width = 10, height = 3)
print(plot)
dev.off()
saveRDS(plot, file = outfile_rds)

# UMAP plot by annotation.
plot_file <- paste0(
  plot_rds_outdir, 
  "/scrnatcr_donor7_sample25_umap_space_20_by_final_anno_30_res_0.8.rds"
)
plot <- readRDS(plot_file)
plot <- plot + 
  theme(aspect.ratio = 1)
outfile_pdf <- paste0(
  plot_outdir, 
  "/scrnatcr_donor7_sample25_umap_space_30_by_final_anno_30_res_0.8_optim_for_publication.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir, 
  "/scrnatcr_donor7_sample25_umap_space_30_by_final_anno_30_res_0.8_optim_for_publication.rds"
)
pdf(outfile_pdf, width = 10, height = 3)
print(plot)
dev.off()
saveRDS(plot, file = outfile_rds)

# UMAP plots by relevant expression values and moddule scores.
relevant_plots <- c(
  readRDS(paste0(
    plot_rds_outdir, 
    "/scrnatcr_donor7_sample25_umap_space_30_by_expr_general_markers.rds"
  ))[10],
  list(
    readRDS(paste0(
      plot_rds_outdir, 
      "/scrnatcr_donor7_sample25_umap_space_30_by_expr_bloodnaive_boodccr8_diffgenes_mod_scores.rds"
    ))[[1]],
    readRDS(paste0(
      plot_rds_outdir, 
      "/scrnatcr_donor7_sample25_umap_space_30_by_expr_bloodmem_boodccr8_diffgenes_mod_scores.rds"
    ))[[2]]
  )
)
relevant_plots[[2]]$labels$title <- "Blood CCR8 sig 1"
relevant_plots[[3]]$labels$title <- "Blood CCR8 sig 2"
for(plot in relevant_plots) {
  title <- plot$labels$title
  print(title)
  plot_new <- plot
  plot_new$labels$title <- NULL
  plot_new$labels$x <- "UMAP 1"
  plot_new$labels$y <- "UMAP 2"
  plot_new <- plot_new + 
    scale_colour_viridis_c(option = "G") +
    theme(aspect.ratio = 1)
  outfile_pdf <- paste0(plot_outdir, 
                        "/scrnatcr_donor7_sample25_umap_space_30_by_", 
                        gsub(" ", "_", title), 
                        "_optim_for_publication.pdf")
  outfile_rds <- paste0(plot_rds_outdir, 
                        "/scrnatcr_donor7_sample25_umap_space_30_by_", 
                        gsub(" ", "_", title), 
                        "_optim_for_publication.rds")
  pdf(outfile_pdf, width = 10, height = 3)
  print(plot_new)
  dev.off()
  saveRDS(plot_new, file = outfile_rds)
}
