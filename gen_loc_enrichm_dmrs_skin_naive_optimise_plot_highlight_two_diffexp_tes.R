# This R script omtimises the enrichment plot of single TE elements in the skin 
#   Treg hypomethylation class of DMRs. The plot will now highlight the two TEs
#   that are also differentially expressed between the two cell types.
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

# Read in the original plot.
orig_plots_list <- readRDS(paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_enrichm_dmrs_skin_naive.rds"
))
orig_plot <- orig_plots_list[[4]]

# Specify the elements that should be highlighted in addition.
elements_to_add <- c("HERVIP10F-int", "LTR45B")

# Add these elements to the text plotting data.
text_plotting_data <- orig_plot$layers[[3]]$data
all_plotting_data <- orig_plot$data
add_data <- all_plotting_data[all_plotting_data$ID %in% elements_to_add, ]
text_plotting_data_new <- rbind(text_plotting_data, add_data)

# Update the text plotting data.
new_plot <- orig_plot
new_plot$layers[[3]]$data <- text_plotting_data_new

# Save the updated_plot.
new_plot_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_enrichm_dmrs_skin_naive_skin_hypo_name_w_diffexp_tes.pdf"
)
new_plot_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_enrichm_dmrs_skin_naive_skin_hypo_name_w_diffexp_tes.rds"
)
pdf(new_plot_pdf, width = 7, height = 6)
print(new_plot)
dev.off()
saveRDS(new_plot, file = new_plot_rds)


