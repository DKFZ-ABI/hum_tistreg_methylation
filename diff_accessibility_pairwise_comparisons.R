# This script analyses differential chromatin accessibility between any pair of
#   cell types that were also included in my methylation analysis.
# Author: Niklas Beumer



# Load required packages.
library(Seurat)
library(Signac)
library(ggplot2)
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

# Read in the Seurat object containing scATAC-seq data of CD4+ cells.
sc_data_cd4_file <- paste0(
  b330_space,
  "yyy/scATAC/final/seurat_human_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar.RDS"
)
sc_data_cd4 <- readRDS(sc_data_cd4_file)






##############################################################
# Differential analysis.
##############################################################

volcano_plot_seur <- function(diff_res, outfile_pdf, outfile_rds) {
  # This function generates a volcano plot based on results from Seurat's
  #   FindMarkers function.
  # diff_res: Data frame from the FindMarkers function. Must contain the
  #   columns "p_val_adj" and "avg_log2FC".
  # outfile_pdf: Character; Path to a PDF file where the plot will be
  #   saved.
  # outfile_rds: Character; Path to an RDS file where the plot will be
  #   saved.
  # Dependencies: ggplot2
  # Value: None, the function just saves the plot.
  # Author: Niklas Beumer
  
  #-- Compute logs of P values.
  diff_res$logp <- -log10(diff_res$p_val_adj)
  
  #-- Take care of the situation in which P values equal to 0 are reported.
  max_logp <- max(diff_res$logp[diff_res$logp < Inf])
  replacement_val <- 1.5 * max_logp
  diff_res$logp[diff_res$logp == Inf] <- replacement_val
  
  #-- Add a column showing whether a transcription factor displayed statistical significance.
  diff_res$signif <- as.character(diff_res$p_val_adj < 0.05)
  
  #-- Generate the plot.
  step_size <- 10 ^ floor(log10(max_logp) - 0.3)
  y_breaks <- seq(0, max_logp, step_size)
  volcano_plot <- ggplot(diff_res) +
    aes(x = avg_log2FC, y = logp, colour = signif) +
    scale_y_continuous(
      breaks = c(y_breaks, replacement_val),
      labels = c(y_breaks, "Inf"),
      expand = expansion(mult = c(0.01, 0.05)),
      name = "-log10(P)"
    ) +
    scale_colour_manual(breaks = c("TRUE", "FALSE"),
                        values = c("red", "grey50")) +
    geom_point(size = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    xlab("avg. Log2(FC)") +
    theme_classic() +
    theme(
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      legend.position = "none"
    )
  
  #-- Save the volcano plot.
  pdf(outfile_pdf, width = 2, height = 2)
  print(volcano_plot)
  dev.off()
  saveRDS(volcano_plot, file = outfile_rds)
  
}

diff_peak_analysis_seur <- function(seur_obj,
                                    baseline_cond,
                                    contrast_cond,
                                    outfile_txt,
                                    outfile_pdf,
                                    outfile_rds) {
  # This function performs differential peak analysis (scATAC-seq) using
  #   functionality from Seurat. The function will use the logistic regression
  #   test for differential analysis and will use the meta data slot
  #   "peak_region_fragments" as the latent variable, as suggested here:
  #   https://stuartlab.org/signac/articles/pbmc_vignette.html; 03-11-2022
  # seur_obj: The Seuratobject to use. The function will use the default assay.
  #   The Idents() must be set to a slot where baseline_cond and contrast_cond
  #   can be found.
  # baseline_cond: Character; Value of Idents(seur_obj) that marks the base line
  #   group. Peaks with higher accessibility in this group will be indicated by
  #   negative Log2FCs.
  # contrast_cond: Character; Value of Idents(seur_obj) that marks the contrast
  #   group. Peaks with higher accessibility in this group will be indicated
  #   by positive Log2FCs.
  # outfile_txt: Character; Path to a TXT file where the results from differential
  #   analysis will be saved.
  # outfile_pdf: Character; Path to a PDF file where a volcano plot will be saved.
  # outfile_rds: Character; Path to an RDS file where a volcano plot will be saved.
  # Value: None, the function just saves the results.
  # Dependencies: Seurat, Signac, ggplot2.
  # Author: Niklas Beumer
  
  # Subset the Seurat object to the categories that shall be compared.
  seur_obj <- seur_obj[, Idents(seur_obj) %in% c(baseline_cond, contrast_cond)]
  
  # Re-run TF-IDF.
  seur_obj <- RunTFIDF(seur_obj)
  
  # Perform differential peak analysis.
  # In preliminary analysis, I found that the category specified under "ident.2"
  # will be used as the base line.
  diff_peak_res <- FindMarkers(
    seur_obj,
    ident.1 = contrast_cond,
    ident.2 = baseline_cond,
    test.use = "LR",
    latent.vars = "peak_region_fragments"
  )
  
  # Add information on whether a peak is significantly differential (adj.
  # P value < 0.05)
  diff_peak_res$significant <- diff_peak_res$p_val_adj < 0.05
  
  # Save the results from differential peak analysis.
  write.table(
    diff_peak_res,
    file = outfile_txt,
    sep = "\t",
    row.names = T,
    col.names = T,
    quote = F
  )
  
  # Generate and save a volcano plot for the differential accessibility
  # results.
  volcano_plot_seur(diff_res = diff_peak_res,
                    outfile_pdf = outfile_pdf,
                    outfile_rds = outfile_rds)
  
}



# Specify the cell types that were also used in my methylation analysis.
cell_types_my_spelling <- c(
  "Skin_Treg",
  "Skin_Tconv",
  "Blood_CCR8+_Treg",
  "Blood_naive_Treg",
  "Blood_naive_Tconv"
)
cell_types_atac_spelling <- c(
  "skin_treg",
  "skin_tconv",
  "blood_ccr8_treg",
  "blood_naive_treg",
  "blood_naive_tconv"
)

# Restrict scATAC-seq data to the cell types that were also included in my
# methylation analysis.
sc_data_cd4_restr <- subset(sc_data_cd4,
                            subset = treg_tconv_annot %in% cell_types_atac_spelling)
# set.seed(2022) # For debugging purposes.
# sc_data_cd4_restr <- sc_data_cd4_restr[sample(rownames(sc_data_cd4_restr), 10000),
#                                        sample(colnames(sc_data_cd4_restr), 1000)] # For debugging purposes.

# Set the Idents() of theSeurat object to the cell type classification.
Idents(sc_data_cd4_restr) <- "treg_tconv_annot"

# Iterate over the pairs of cell types.
for (i in 1:(length(cell_types_my_spelling) - 1)) {
  for (j in (i + 1):length(cell_types_my_spelling)) {
    # Specifythe cell types that are currently used.
    cell_type_1_my_spelling <- cell_types_my_spelling[i]
    cell_type_1_atac_spelling <- cell_types_atac_spelling[i]
    cell_type_2_my_spelling <- cell_types_my_spelling[j]
    cell_type_2_atac_spelling <- cell_types_atac_spelling[j]
    
    # Perform differential accessibility analysis using cell type 2 as the base line.
    diff_peak_res_outfile <- paste0(
      location,
      "/differential_accessibility/diff_acc_analysis_",
      cell_type_1_my_spelling,
      cell_type_2_my_spelling,
      "_results_with_significance.txt"
    )
    volc_file_snip <- paste0("/diff_acc_",
                             cell_type_1_my_spelling,
                             cell_type_2_my_spelling,
                             "_volc")
    diff_peak_analysis_seur(
      seur_obj = sc_data_cd4_restr,
      baseline_cond = cell_type_1_atac_spelling,
      contrast_cond = cell_type_2_atac_spelling,
      outfile_txt = diff_peak_res_outfile,
      outfile_pdf = paste0(plot_outdir, volc_file_snip, ".pdf"),
      outfile_rds = paste0(plot_rds_outdir, volc_file_snip, ".rds")
    )
    
  }
}
