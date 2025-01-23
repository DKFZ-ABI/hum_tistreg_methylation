# This script generates an scATAC-seq UMAP space containing only the considered 
#   cell types.
# Author: Niklas Beumer



# Load required packages.
library(Seurat)
library(Signac)
library(ggplot2)
library(harmony)


# Define a location on /yyy.
b330_space <- "/yyy/OE0436/internal/"
location <- paste0(b330_space, "nbeumer/hm_treg_bs_rgnsbg")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the Seurat object containing scATAC-seq data of CD4+ cells.
sc_data_cd4_file <- paste0(
  b330_space, 
  "yyy/scATAC/final/seurat_human_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar.RDS"
)
sc_data_cd4 <- readRDS(sc_data_cd4_file)






############################################################################
# Subset the Seurat object to the relevant cell types and generate a new
# UMAP representation.
############################################################################

# Specify the relevant cell types .
cell_types <- c("skin_treg", "skin_tconv", "blood_ccr8_treg", 
                "blood_naive_treg", "blood_naive_tconv")
cell_types_names_f_plot <- c("Skin Treg", "Skin Tconv", "Blood CCR8+ Treg", 
                             "Blood naive Treg", "Blood naive Tconv")

# Define colours for the cell types.
cell_type_col_palette <- c("blue", "cyan2", "orange", "darkorchid1", "red")

# Restrict scATAC-seq data to the relevant cell types.
sc_data_cd4_restr <- subset(
  sc_data_cd4, 
  subset = treg_tconv_annot %in% cell_types
)
# set.seed(1000) # For debugging purposes
# sc_data_cd4_restr <- sc_data_cd4_restr[,
#   sample(1:ncol(sc_data_cd4_restr), 1000)
# ] # For debugging purposes

# Run a new TF-IDF.
sc_data_cd4_restr <- RunTFIDF(sc_data_cd4_restr)

# Identify top features.
sc_data_cd4_restr <- FindTopFeatures(sc_data_cd4_restr)

# Run singular value decomposition.
sc_data_cd4_restr <- RunSVD(sc_data_cd4_restr, 
                            reduction.key = "RelevantcellsLSI_",
                            reduction.name = "Relevant_cells_LSI")

# Correct donor-specific effects by Harmony.
# project.dim = F is required to prevent an error during the execution 
# (according to preliminary testing).
sc_data_cd4_restr <- RunHarmony(sc_data_cd4_restr, group.by.vars = "donor", 
                                reduction = "Relevant_cells_LSI", 
                                plot_convergence = T, 
                                reduction.save = "Relevant_cells_Harmony",
                                assay.use = "scATAC_raw", project.dim = F) 

# Assess correlation of SVD components with sequencing depth.
depth_cor_plot <- DepthCor(sc_data_cd4_restr, 
                           reduction = "Relevant_cells_Harmony")
depth_cor_pdf <- paste0(
  plot_outdir, 
  "/scatacseq_data_relevant_celltypes_svd_harmony_depthcor.pdf"
)
depth_cor_rds <- paste0(
  plot_rds_outdir, 
  "/scatacseq_data_relevant_celltypes_svd_harmony_depthcor.rds"
)
pdf(depth_cor_pdf, width = 4, height = 4)
print(depth_cor_plot)
dev.off()
saveRDS(depth_cor_plot, file = depth_cor_rds)

# Perform UMAP.
# Ignore the first SVD component as it is strongly correlated with sequencing 
# depth.
sc_data_cd4_restr <- RunUMAP(sc_data_cd4_restr, 
                             dims = 2:50, 
                             reduction = "Relevant_cells_Harmony", 
                             reduction.name = "Relevant_cells_UMAP", 
                             reduction.key = "RelevantcellsUMAP_")

# Plot cell types in UMAP space and save this plot.
celltype_plot <- DimPlot(sc_data_cd4_restr, group.by = "treg_tconv_annot",
                         reduction = "Relevant_cells_UMAP") +
  scale_colour_manual(breaks = cell_types, labels = cell_types_names_f_plot,
                      values = cell_type_col_palette, name = "Cell type") +
  theme(aspect.ratio = 1)
celltype_plot_pdf <- paste0(
  plot_outdir, 
  "/scatacseq_data_relevant_celltypes_umap_by_celltype.pdf"
)
celltype_plot_rds <- paste0(
  plot_rds_outdir, 
  "/scatacseq_data_relevant_celltypes_umap_by_celltype.rds"
)
pdf(celltype_plot_pdf, width = 5, height = 5)
print(celltype_plot)
dev.off()
saveRDS(celltype_plot, file = celltype_plot_rds)

# Plot donors in UMAP space and save this plot.
donor_plot <- DimPlot(sc_data_cd4_restr, group.by = "donor",
                      reduction = "Relevant_cells_UMAP") +
  theme(aspect.ratio = 1)
donor_plot_pdf <- paste0(
  plot_outdir, 
  "/scatacseq_data_relevant_celltypes_umap_by_donor.pdf"
)
donor_plot_rds <- paste0(
  plot_rds_outdir, 
  "/scatacseq_data_relevant_celltypes_umap_by_donor.rds"
)
pdf(donor_plot_pdf, width = 5, height = 5)
print(donor_plot)
dev.off()
saveRDS(donor_plot, file = donor_plot_rds)

# Save the Seurat object containing the UMAP space of relevant cell types.
seur_outfile <- paste0(
  location, 
  "/differential_accessibility/scatacseq_data_relevant_celltypes_seur_obj_w_umap.rds"
)
saveRDS(sc_data_cd4_restr, file = seur_outfile)

