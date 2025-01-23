# This script generates common Homer plots for differentially methylated or 
#   differentially accessible regions defined by specific blood CCR8+ Treg 
#   cell positionings.
# Author: Niklas Beumer



# Load required package(s)
library(ggplot2)


# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                      format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in Homer results (against the remaining regions) for DMRs defined
# by specific blood CCR8+ Treg positionings.
file_snips <- c("skin_treg_hypomethylation_closer_to_blood_naive",
                "skin_treg_hypomethylation_closer_to_skin")
homer_results_meth <- lapply(file_snips, FUN = function(x) {
  homer_file <- paste0(
    location, 
    "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_", 
    x, 
    "_homer_against_remaining_regions/knownResults.txt"
  )
  read.table(homer_file, sep = "\t", header = T, stringsAsFactors = F, 
             comment.char = "")
})
names(homer_results_meth) <- c(
  "Skin_Treg__hypomethylation\n(blood CCR8+ Tregs are closer to blood naive Tregs)",
  "Skin_Treg__hypomethylation\n(blood CCR8+ Tregs are closer to skin Tregs)"
)

# Read in Homer results (against the remaining regions) for diff. peaks defined
# by specific blood CCR8+ Treg positionings.
file_snips <- c("skin_treg_hyperaccessibility_closer_to_blood_naive",
                "skin_treg_hyperaccessibility_closer_to_skin")
homer_results_acc <- lapply(file_snips, FUN = function(x) {
  homer_file <- paste0(
    location, 
    "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_", 
    x, 
    "_homer_against_remaining_regions/knownResults.txt"
  )
  read.table(homer_file, sep = "\t", header = T, stringsAsFactors = F, 
             comment.char = "")
})
names(homer_results_acc) <- c(
  "Skin_Treg__hyperaccessibility\n(blood CCR8+ Tregs are closer to blood naive Tregs)",
  "Skin_Treg__hyperaccessibility\n(blood CCR8+ Tregs are closer to skin Tregs)"
)

# Specify the motifs that should be plotted.
# These comprise (i) a collection of bZIP motifs and (ii) a collection of 
# bHLH motifs with CACGT sub-motif.
tfs_to_show_bhlh <- c(
  "NPAS(bHLH)/Liver-NPAS-ChIP-Seq(GSE39860)/Homer",
  "BMAL1(bHLH)/Liver-Bmal1-ChIP-Seq(GSE39860)/Homer",
  "bHLHE41(bHLH)/proB-Bhlhe41-ChIP-Seq(GSE93764)/Homer",
  "NPAS2(bHLH)/Liver-NPAS2-ChIP-Seq(GSE39860)/Homer",
  "MNT(bHLH)/HepG2-MNT-ChIP-Seq(Encode)/Homer",
  "n-Myc(bHLH)/mES-nMyc-ChIP-Seq(GSE11431)/Homer",
  "c-Myc(bHLH)/mES-cMyc-ChIP-Seq(GSE11431)/Homer",
  "c-Myc(bHLH)/LNCAP-cMyc-ChIP-Seq(Unpublished)/Homer",
  "CLOCK(bHLH)/Liver-Clock-ChIP-Seq(GSE39860)/Homer",
  "USF1(bHLH)/GM12878-Usf1-ChIP-Seq(GSE32465)/Homer"
)

tfs_to_show_bzip <- c(
  "BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer",
  "Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer",
  "AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer",
  "Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer",
  "JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer",
  "Fosl2(bZIP)/3T3L1-Fosl2-ChIP-Seq(GSE56872)/Homer",
  "Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer",
  "Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer",
  "Bach2(bZIP)/OCILy7-Bach2-ChIP-Seq(GSE44420)/Homer",
  "HLF(bZIP)/HSC-HLF.Flag-ChIP-Seq(GSE69817)/Homer",
  "Atf1(bZIP)/K562-ATF1-ChIP-Seq(GSE31477)/Homer",
  "NFIL3(bZIP)/HepG2-NFIL3-ChIP-Seq(Encode)/Homer",
  "bZIP:IRF(bZIP,IRF)/Th17-BatF-ChIP-Seq(GSE39756)/Homer",
  "MafA(bZIP)/Islet-MafA-ChIP-Seq(GSE30298)/Homer"
)






########################################################################
# Generate a common dot plot for the Homer results.
########################################################################

# Concatenate Homer results into a single data frame.
names_combined <- c(names(homer_results_meth), names(homer_results_acc))
comb_list <- c(homer_results_meth, homer_results_acc)
homer_results_combined <- do.call(
  rbind, 
  lapply(1:length(comb_list), FUN = function(x) {
    curr_df <- comb_list[[x]]
    colnames(curr_df)[6] <- "Number..of.Target.Sequences.with.Motif"
    colnames(curr_df)[8] <- "Number..of.Background.Sequences.with.Motif"
    curr_df$Region_class <- names_combined[x]
    return(curr_df)
}))

# Compute enrichment scores.
target_perc <- as.numeric(gsub(
  "%", 
  "", 
  homer_results_combined$X..of.Target.Sequences.with.Motif
))
background_perc <- as.numeric(gsub(
  "%", 
  "", 
  homer_results_combined$X..of.Background.Sequences.with.Motif
))
homer_results_combined$Enrichment_score <- target_perc / background_perc

# Collect values for these motifs.
relevant_values_bhlh <- homer_results_combined[
  homer_results_combined$Motif.Name %in% tfs_to_show_bhlh, ]
relevant_values_bzip <- homer_results_combined[
  homer_results_combined$Motif.Name %in% tfs_to_show_bzip, ]

# Restrict bHLH data to the methylation level (as these motifs had been 
# shown to be primarily affected by methylation).
relevant_values_bhlh <- relevant_values_bhlh[
  relevant_values_bhlh$Region_class %in% names_combined[1:2], ]

# Clip Benjamini-Hochberg-adjusted P values.
relevant_values_bhlh$qval_clipped <- 
  sapply(relevant_values_bhlh$q.value..Benjamini., FUN = function(x) {
    max(x, 0.001)
})
relevant_values_bzip$qval_clipped <- 
  sapply(relevant_values_bzip$q.value..Benjamini., FUN = function(x) {
    max(x, 0.001)
})

# Compute log-transformed q values.
relevant_values_bhlh$logq <- -log10(relevant_values_bhlh$qval_clipped) + 0.75
relevant_values_bzip$logq <- -log10(relevant_values_bzip$qval_clipped) + 0.75

# Add a column indicating statistical significance.
relevant_values_bhlh$significant <- 
  relevant_values_bhlh$q.value..Benjamini. < 0.05  
relevant_values_bzip$significant <- 
  relevant_values_bzip$q.value..Benjamini. < 0.05  

# Introduce a column showing a binary enrichment classification.
relevant_values_bhlh$greater1 <- sapply(
  relevant_values_bhlh$Enrichment_score, FUN = function(x) {
    ifelse(x > 1, yes = "enriched", no = "depleted")
})
relevant_values_bzip$greater1 <- sapply(
  relevant_values_bzip$Enrichment_score, FUN = function(x) {
    ifelse(x > 1, yes = "enriched", no = "depleted")
})

# Generate and save dot plots.
# Generate a dot plot showing homer results.
relevant_values_bhlh$Motif.Name <- factor(relevant_values_bhlh$Motif.Name, 
                                          levels = rev(tfs_to_show_bhlh))
relevant_values_bzip$Motif.Name <- factor(relevant_values_bzip$Motif.Name, 
                                          levels = rev(tfs_to_show_bzip))
relevant_values_bhlh$Region_class <- factor(relevant_values_bhlh$Region_class, 
                                            levels = names_combined)
relevant_values_bzip$Region_class <- factor(relevant_values_bzip$Region_class, 
                                            levels = names_combined)
dot_plot_bhlh <- ggplot(relevant_values_bhlh) +
  aes(x = Region_class, y = Motif.Name, size = logq, 
      fill = greater1, colour = significant) +
  scale_size_area(breaks = 0.75:3.75, labels = c("0", "1", "2", ">=3"),
                  name = "-log10(q)",
                  na.value = 1) +
  scale_fill_manual(name = "Enrichment",
                    breaks = c("depleted", "enriched"),
                    values = c("blue", "red"),
                    na.value = "black") +
  scale_colour_manual(breaks = c(T, F), values = c("cyan", "grey"),
                      labels = c("q < 0.05", "q >= 0.05"), 
                      name = "Significance",
                      na.value = "black") +
  geom_point(stroke = 1, shape = 21) +
  xlab("Region class") +
  ylab("Motif") +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, 
                                   vjust = 0.5),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))
dot_plot_bhlh_pdf <- paste0(
  plot_outdir, 
  "/treg_hierarchies_diff_meth_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_homer_against_rem_reg_common_dot_plot_bHLH.pdf"
)
dot_plot_bhlh_rds <- paste0(
plot_rds_outdir, 
"/treg_hierarchies_diff_meth_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_homer_against_rem_reg_common_dot_plot_bHLH.rds"
)
pdf(dot_plot_bhlh_pdf, width = 6, height = 7)
print(dot_plot_bhlh)
dev.off()
saveRDS(dot_plot_bhlh, file = dot_plot_bhlh_rds)
dot_plot_bzip <- ggplot(relevant_values_bzip) +
  aes(x = Region_class, y = Motif.Name, size = logq, 
      fill = greater1, colour = significant) +
  scale_size_area(breaks = 0.75:3.75, labels = c("0", "1", "2", ">=3"),
                  name = "-log10(q)",
                  na.value = 1) +
  scale_fill_manual(name = "Enrichment",
                    breaks = c("depleted", "enriched"),
                    values = c("blue", "red"),
                    na.value = "black") +
  scale_colour_manual(breaks = c(T, F), values = c("cyan", "grey"),
                      labels = c("q < 0.05", "q >= 0.05"), 
                      name = "Significance",
                      na.value = "black") +
  geom_point(stroke = 1, shape = 21) +
  xlab("Region class") +
  ylab("Motif") +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, 
                                   vjust = 0.5),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))
dot_plot_bzip_pdf <- paste0(
  plot_outdir, 
  "/treg_hierarchies_diff_meth_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_homer_against_rem_reg_common_dot_plot_bZIP.pdf"
)
dot_plot_bzip_rds <- paste0(
  plot_rds_outdir, 
  "/treg_hierarchies_diff_meth_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_homer_against_rem_reg_common_dot_plot_bZIP.rds"
)
pdf(dot_plot_bzip_pdf, width = 7, height = 7)
print(dot_plot_bzip)
dev.off()
saveRDS(dot_plot_bzip, file = dot_plot_bzip_rds)

