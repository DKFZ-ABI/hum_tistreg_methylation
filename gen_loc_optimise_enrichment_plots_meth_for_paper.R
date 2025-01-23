# This script optimises the TE enrichment (methylation, skin Treg hypometh.,
#   subfamily level) plot for our publication.
# Author: Niklas Beumer



# Load required package(s)
library(ggplot2)
library(ggrepel)

# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")





# Read in the list of plots.
plot_list <- readRDS(paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_enrichm_dmrs_skin_naive.rds"
))

# Optimise the plot.
curr_plot <- plot_list[[4]]
curr_plot_data <- curr_plot$data
curr_plot_data_signif <- curr_plot_data[curr_plot_data$Pval_adj_BH < 0.05, ]
names_to_show <- c(
  curr_plot_data_signif$ID[
    order(curr_plot_data_signif$Enrichment, decreasing = T)[1:10]
  ], 
  "LTR45B", 
  "HERVIP10F-int",
  "MIR3",
  "L3",
  "L2b",
  "MLT1D",
  "MIR",
  "MIRc",
  "MER117",
  "MIRb",
  "MER20",
  "MER94",
  "L2a",
  "MER5B",
  "MLT1J",
  "MLT1H"
)
curr_plot$layers[[3]]$data <- 
  curr_plot_data_signif[curr_plot_data_signif$ID %in% names_to_show, ]
curr_plot_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_enrichm_dmrs_skin_naive_skin_side_name_optim_f_publication.pdf"
)
curr_plot_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_enrichm_dmrs_skin_naive_skin_side_name_optim_f_publication.rds"
)
pdf(curr_plot_pdf, width = 7, height = 6)
print(curr_plot)
dev.off()
saveRDS(curr_plot, file = curr_plot_rds)
