# This script optimises an enrichment plot of TEs in skin Treg hypomethylation
#   regions stratified by blood CCR8+ Treg cell positionings (against remaining 
#   regions) for the publication. All significant TEs are labelled now.
# Author: Niklas Beumer



# Load required packages.
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

# Read in the list of corresponding plots and extract the relevant plot.
plot_list_file <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_enrichm_blood_ccr8_pos_meth_skin_treg_hypometh_against_rem_reg.rds"
)
plot_list <- readRDS(plot_list_file)
relevant_plot <- plot_list[[4]]

# Modify the plot so that all significantly enriched TEs are labelled.
point_data <- relevant_plot$data
text_data_old <- relevant_plot$layers[[3]]$data
all_signif_tes <- point_data$ID[point_data$Pval_adj_BH < 0.05]
tes_to_add <- setdiff(all_signif_tes, text_data_old$ID)
text_data_new <- rbind(text_data_old, 
                       point_data[point_data$ID %in% tes_to_add, ])
relevant_plot$layers[[3]]$data <- text_data_new

# Save the optimised plot.
outfile_pdf <- paste0(
  plot_outdir, 
  "/te_analysis_gen_loc_enrichm_blood_ccr8_pos_meth_skin_treg_hypometh_against_rem_reg_name_optim_f_publication.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir, 
  "/te_analysis_gen_loc_enrichm_blood_ccr8_pos_meth_skin_treg_hypometh_against_rem_reg_name_optim_f_publication.rds"
)
pdf(outfile_pdf, width = 7, height = 6)
print(relevant_plot)
dev.off()
saveRDS(relevant_plot, file = outfile_rds)
