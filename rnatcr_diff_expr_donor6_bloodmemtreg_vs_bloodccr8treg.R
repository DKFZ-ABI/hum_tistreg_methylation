# This script analyses differential gene expression between blood CD45RO+ Treg
#   cells and blood CD45RO+CCR8+ Treg cells in scRNA/scTCR data from donor 6.
# Author: Niklas Beumer



# Load required packages.
library(Seurat)
library(ggplot2)


# Define a location on /yyy.
b330_space <- "/yyy/"
location <- paste0(b330_space, "nbeumer/hm_treg_bs_rgnsbg")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the Seurat object containing scRNA-seq data from donor 6. 
# The sample MD-scRNA_R26 contains blood CD4+ T cells.
# The sample MD-scRNA_R28 contains blood CD4+CD25 T cells.
# The sample MD-scRNA_R29 contains blood CD4+CD25+CD45RO+ T cells.
# The sample MD-scRNA_R30 contains blood CD4+CD25+CD45RO+CCR8+ T cells.
sc_data_file <- paste0(
  b330_space, 
  "imbusch/d100/10x_scRNA_human/niklas/donor6_CD4_all_samples.rds"
)
sc_data <- readRDS(sc_data_file)

# Restrict single-cell data to the relevant samples.
sc_data_restr <- subset(
  sc_data, subset = orig.ident %in% c("MD-scRNA_R29", "MD-scRNA_R30")
)

# Analyse differential gene expression between blood CD4+CD25+CD45RO+ T cells
# and blood CD4+CD25+CD45RO+CCR8+ T cells.
Idents(sc_data_restr) <- sc_data_restr$orig.ident
diff_expr <- FindMarkers(sc_data_restr,
                         ident.1 = "MD-scRNA_R30",
                         ident.2 = "MD-scRNA_R29")

# Add information on whether a gene is significantly differential (adj. 
# P value < 0.05)
diff_expr$significant <- diff_expr$p_val_adj < 0.05

# Save the results from differential expression analysis.
outfile_txt <- paste0(
  location, 
  "/treg_hierarchies/scrnatcr_diff_expr_bloodmemtreg_vs_bloodccr8treg_bl_bloodmemtreg.txt"
)
write.table(diff_expr, file = outfile_txt, sep = "\t",
            row.names = T, col.names = T, quote = F)

# Generate a volcano plot for the differential analysis.
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
  
  #-- Add a column showing whether a transcription factor displayed statistical 
  #-- significance.
  diff_res$signif <- as.character(diff_res$p_val_adj < 0.05)
  
  #-- Generate the plot.
  step_size <- 10 ^ floor(log10(max_logp) - 0.3)
  y_breaks <- seq(0, max_logp, step_size)
  volcano_plot <- ggplot(diff_res) +
    aes(x = avg_log2FC, y = logp, colour = signif) +
    scale_y_continuous(breaks = c(y_breaks, replacement_val),
                       labels = c(y_breaks, "Inf"),
                       expand = expansion(mult = c(0.01, 0.05)),
                       name = "-log10(P)") +
    scale_colour_manual(breaks = c("TRUE", "FALSE"), 
                        values = c("red", "grey50")) +
    geom_point(size = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    xlab("avg. Log2(FC)") +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "none")
  
  #-- Save the volcano plot.
  pdf(outfile_pdf, width = 2, height = 2)
  print(volcano_plot)
  dev.off()
  saveRDS(volcano_plot, file = outfile_rds)
  
}
void <- volcano_plot_seur(
  diff_res = diff_expr, 
  outfile_pdf = paste0(
    plot_outdir, 
    "/scrnatcr_diff_expr_bloodmemtreg_vs_bloodccr8treg_bl_bloodmemtreg_volc.pdf"
  ),
  outfile_rds = paste0(
    plot_rds_outdir, 
    "/scrnatcr_diff_expr_bloodmemtreg_vs_bloodccr8treg_bl_bloodmemtreg_volc.rds"
  )
)
