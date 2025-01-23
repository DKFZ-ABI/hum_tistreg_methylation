# This script visualises FACS results regarding CCR8 expression on skin Treg
#   cells.
# Author: Niklas Beumer



# Load required package(s).
library(ggplot2)


# Define a location on /yyy.
location <- "/yyy/nbeumer/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Hard-code percentages of CCR8-positive cells among skin Tregs and blood
# Tregs, as determined by FACS.  The values correspond to the percentages
# of CD45RA-CCR8+ cells among sorted Tregs.
# These data are paired (meaning that we have one blood value and one skin 
# value from each donor).
ccr8_percentages <- data.frame(
  Donor = rep(1:5, 2),
  Cell_type = rep(c("Blood Treg", "Skin Treg"), each = 5),
  Percentage = c(29.4, 8.34, 15.6, 5.08, 22.5,
                 80.5, 36.3, 83.7, 51.4, 85.0)
)

# Perform a two-tailed paired t test test assessing whether skin
# contains a higher percentage of CD45RA-CCR8+ Tregs than blood.
t_test_res <- t.test(ccr8_percentages$Percentage[1:5],
                     ccr8_percentages$Percentage[6:10],
                     paired = T,
                     alternative = "two.sided")

# Save the test results.
test_outfile <- paste0(
  location, 
  "/facs_data_ccr8_expression_skin_bloood_tregs_paired_t_teest_res.rds"
)
saveRDS(t_test_res, file = test_outfile)

# Visualise the percentages.
segments_df <- data.frame(x = rep(1, 5),
                          xend = rep(2, 5),
                          y = ccr8_percentages$Percentage[1:5],
                          yend = ccr8_percentages$Percentage[6:10])
ccr8_perc_plot <- ggplot(ccr8_percentages) +
  aes(x = Cell_type, y = Percentage) +
  scale_y_continuous(
    limits = c(0, 100), 
    expand = c(0, 0),
    name = "Percentage of CD45RA-CCR8+ cells\namong Treg cells [%]"
  ) +
  geom_point(shape = 21, fill = "grey", size = 2) +
  geom_segment(data = segments_df,
               mapping = aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(x = 1.5, y = 1.1 * max(ccr8_percentages$Percentage),
            label = paste0("P = ", signif(t_test_res$p.value, 3))) +
  geom_segment(x = 1, 
               xend = 2, 
               y = 1.05 * max(ccr8_percentages$Percentage),
               yend = 1.05 * max(ccr8_percentages$Percentage)) +
  xlab("Cell type") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))

# Save the plot.
outfile_pdf <- paste0(plot_outdir, 
                      "/facs_data_ccr8_expression_skin_bloood_tregs.pdf")
outfile_rds <- paste0(plot_rds_outdir, 
                      "/facs_data_ccr8_expression_skin_bloood_tregs.rds")
pdf(outfile_pdf, width = 3, height = 3)
print(ccr8_perc_plot)
dev.off()
saveRDS(ccr8_perc_plot, file = outfile_rds)
