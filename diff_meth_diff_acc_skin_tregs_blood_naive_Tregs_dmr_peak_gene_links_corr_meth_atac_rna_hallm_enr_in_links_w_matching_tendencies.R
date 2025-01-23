# This script analyses enrichment of Hallmark gene sets in DMR-peak-gene links 
#   that display matching differential tendencies on the methylation, 
#   accessibility and expression level.
# Author: Niklas Beumer




# Load required package(s)
library(ggplot2)
library(testit)


# Define a location on /xxx.
location <- "/xxxx/nbeumer/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/yyy/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the list of DMR-peak-gene links with differentiality effect sizes.
links_file <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_diff_acc_diff_gex_dmr_peak_gene_links_meth_acc_gex_corr.txt"
)
links <- read.table(links_file, header = T, stringsAsFactors = F, sep = "\t")

# Read in the file containing gene annotations for normal chromosomes and 
# extract all considered gene names.
annotation_pref <- 
  "/xxx/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/GencodeV19_"
genes <- read.table(paste0(annotation_pref, "Genes_plain.bed.gz"), header = T, 
                    stringsAsFactors = F, comment.char = "")
all_genes <- genes$name

# Read in the Hallmark gene sets from MSigDB.
hallmark_genesets_file <- paste0(
  location, 
  "/external_data/2024-05-16_snapshot_msigdb_hallmark_gene_sets_H_human.rds"
)
hallmark_genesets <- readRDS(hallmark_genesets_file)




###############################################################################
# Extract region sets defined by specific differential tendencies on the 
# three lvels. Specifically, I focus on links that are hypomethylated,
# hypereaccessible and hyperexpressed in blood naive Tregs as well as on links
# that are hypomethylated, hyperaccessible and hyperexpressed in skin Tregs.
###############################################################################

# Extract the relevant links.
meth_target <- sapply(
  links$Methylation__automatic_annotation, 
  FUN = function(x) {
    strsplit(x, split = "__")[[1]][1]
  })
acc_target <- sapply(
  links$Accessibility__signature_class, 
  FUN = function(x) {
    strsplit(x, split = "__")[[1]][1]
  })
gex_target <- sapply(
  links$Expression__signature_class, 
  FUN = function(x) {
    strsplit(x, split = "__")[[1]][1]
  })
blood_naive_side <- links[meth_target == acc_target & 
                            acc_target == gex_target & 
                            gex_target == "Blood_naive_Treg", ]
skin_side <- links[meth_target == acc_target & 
                     acc_target == gex_target & 
                     gex_target == "Skin_Treg", ]
region_sets <- list(blood_naive_side, skin_side)

# Generate file name snippets for each of these region sets.
file_snips <- paste0("hypometh_hyperacc_hyperexpr_", 
                     c("blood_naive",
                       "skin"))




###############################################################################
# Analyse enrichment of Hallmark gene sets (compared to the genome).
###############################################################################

# Identify unique names of Hallmark gene sets.
hallmark_names <- unique(hallmark_genesets$gs_name)

# Iterate over the region sets.
hallm_enr_res <- lapply(1:length(region_sets), FUN = function(i) {
  
  # Extract the set of regions to use.
  region_set_to_use <- region_sets[[i]]
  
  # Identify genes that are part of at least one DMR-peak-gene link.
  regions_genes <- unique(region_set_to_use$Expression__gene_name)
  
  # Identify all other genes in the Gencode annotation.
  all_other_genes <- setdiff(all_genes, regions_genes)
  
  # Analyse enrichment of Hallmark gene sets in the extracted genes.
  # (Fisher's exact test)
  #-- Iterate over the gene sets.
  fisher_test_res <- do.call(rbind, lapply(hallmark_names, FUN = function(x) {
    #-- Extract genes for this gene set.
    corresp_genes <- hallmark_genesets$gene_symbol[
      hallmark_genesets$gs_name == x]
    #-- Identify and quantify overlaps between the genes associated with 
    #-- extracted regions and genes from this gene set.
    n_regions_overlap <- length(intersect(regions_genes, corresp_genes))
    n_regions_not_overlap <- length(setdiff(regions_genes, corresp_genes))
    #-- Identify and quantify overlaps between the remaining regions and the
    #-- gene set.
    n_complement_overlap <- length(intersect(all_other_genes, corresp_genes))
    n_complement_not_overlap <- length(setdiff(all_other_genes, corresp_genes))
    #-- Perform a one-tailed Fisher's exact test to assess whether the current 
    #-- gene set is enriched among the extracted regions.
    fisher_input <- matrix(c(n_regions_overlap, n_regions_not_overlap, 
                             n_complement_overlap, n_complement_not_overlap),
                           ncol = 2, byrow = T)
    fisher_test_res <- fisher.test(fisher_input, alternative = "greater")
    or <- fisher_test_res$estimate
    pval <- fisher_test_res$p.value
    #-- Collect results.
    temp_df <- data.frame(
      Hallmark_set = x,
      n_in_class_overlap_hallm_set = n_regions_overlap,
      n_in_class_not_overlap_hallm_set = n_regions_not_overlap,
      n_not_in_class_overlap_hallm_set = n_complement_overlap,
      n_not_in_class_not_overlap_hallm_set = n_complement_not_overlap,
      OR = or,
      Fisher_pval = pval
    )
    #-- Return value.
    return(temp_df)
  }))
  
  
  # Correct P values for multiple hypothesis testing (Benjamini-Hochberg 
  # correction).
  fisher_test_res$Fisher_pval_adj_BH <- p.adjust(fisher_test_res$Fisher_pval, 
                                                 method = "BH")

  # Save the results from Fisher's exact test.
  fisher_test_res_outfile <- paste0(
    location, 
    "/treg_hierarchies/multiom_skintregsig_",
    file_snips[i], 
    "_hallm_enr.txt")
  write.table(fisher_test_res, file = fisher_test_res_outfile, sep = "\t",
              row.names = F, col.names = T, quote = F)
  
  # Generate and save a volcano plot.
  fisher_test_res$logp <- -log10(fisher_test_res$Fisher_pval_adj_BH)
  max_logp <- max(fisher_test_res$logp[fisher_test_res$logp < Inf])
  logp_replacement_val <- 1.5 * max_logp
  fisher_test_res$logp[fisher_test_res$logp == Inf] <- logp_replacement_val
  fisher_test_res$logor <- log10(fisher_test_res$OR)
  max_logor <- max(fisher_test_res$logor[fisher_test_res$logor < Inf])
  min_logor <- min(fisher_test_res$logor[fisher_test_res$logor > -Inf])
  max_logor_replacement_val <- 1.5 * max_logor
  if (min_logor < 0) {
    min_logor_replacement_val <- 1.5 * min_logor
  } else {
    min_logor_replacement_val <- -1.5
  }
  fisher_test_res$logor[fisher_test_res$logor == Inf] <- 
    max_logor_replacement_val
  fisher_test_res$logor[fisher_test_res$logor == -Inf] <- 
    min_logor_replacement_val
  fisher_test_res$significant <- fisher_test_res$Fisher_pval_adj_BH < 0.05
  step_size_y <- 10 ^ floor(log10(max_logp) - 0.3)
  y_breaks <- seq(0, max_logp, step_size_y)
  step_size_x <- 10 ^ floor(abs(log10(max(c(max_logor, min_logor))) - 0.3))
  x_breaks <- c(rev(-seq(0, 
                         ifelse(min_logor < 0, yes = -min_logor, no = min_logor), 
                         step_size_x)),
                seq(step_size_x, 
                    max_logor, 
                    step_size_x))
  volcano_plot <- ggplot(fisher_test_res) +
    aes(x = logor, y = logp, colour = significant) +
    scale_x_continuous(breaks = c(min_logor_replacement_val, x_breaks, 
                                  max_logor_replacement_val),
                       labels = c("-Inf", x_breaks, "Inf"),
                       name = "log10(OR)") +
    scale_y_continuous(breaks = c(y_breaks, logp_replacement_val),
                       labels = c(y_breaks, "Inf"),
                       expand = expansion(mult = c(0.01, 0.05)),
                       name = "-log10(q)") +
    scale_colour_manual(breaks = c("TRUE", "FALSE"), 
                        values = c("red", "grey50")) +
    geom_point(size = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "none")
  volcano_plot_pdf <- paste0(
    plot_outdir, 
    "/multiom_skintregsig_",
    file_snips[i], 
    "_hallm_enr.pdf"
  )
  volcano_plot_rds <- paste0(
    plot_rds_outdir, 
    "/multiom_skintregsig_",
    file_snips[i], 
    "_hallm_enr.rds"
  )
  pdf(volcano_plot_pdf, width = 3, height = 3)
  print(volcano_plot)
  dev.off()
  saveRDS(volcano_plot, file = volcano_plot_rds)
  
  # Return the results from Fisher's exact tests.
  return(fisher_test_res)
  
})

# Generate names for the collected results.
names(hallm_enr_res) <- file_snips




###########################################################################
# Visualise results for all hallmark gene sets.
###########################################################################

# Concatenate enrichment data.
data_to_plot <- do.call(rbind, lapply(file_snips, FUN = function(x) {
  fisher_res_to_use <- hallm_enr_res[[x]]
  fisher_res_to_use$link_set <- x
  return(fisher_res_to_use)
}))

# Assert that there were no odds ratios equal to Inf or -Inf (because these 
# might have been clipped in the logor column previously).
assert(all(data_to_plot$OR > -Inf))
assert(all(data_to_plot$OR < Inf))

# Clip log-transformed P values.
data_to_plot$logp_clipped <- sapply(data_to_plot$logp, FUN = function(x) {
  min(x, 4)
})

# Visualise enrichments.
data_to_plot$logp_clipped_w_offset <- data_to_plot$logp_clipped + 0.75
hallm_plot <- ggplot(data_to_plot) +
  aes(x = link_set, y = Hallmark_set, fill = logor, 
      size = logp_clipped_w_offset, colour = significant) +
  scale_x_discrete(breaks = file_snips, 
                   labels = c("Blood naive\nTreg", "Skin\nTreg"),
                   name = "Hypometh.,\nhyperacc and\nhyperexpr. in") +
  scale_fill_gradient2(low = "blue", mid = "grey70", high = "red",
                       name = "log10(OR)") +
  scale_size_area(breaks = 0.75:4.75, labels = c("0", "1", "2", "3", ">=4"),
                  name = "-log10(q)") +
  scale_colour_manual(breaks = c(T, F), values = c("cyan", "grey"),
                      labels = c("q < 0.05", "q >= 0.05"), 
                      name = "Significance") +
  geom_point(shape = 21, stroke = 1) +
  ylab("Hallmark gene set") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))

# Save the plot showing enrichments.
outfile_pdf <- paste0(
  plot_outdir, 
  "/multiom_skintregsig_hallm_enr.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir, 
  "/multiom_skintregsig_hallm_enr.rds"
)
pdf(outfile_pdf, width = 6, height = 13)
print(hallm_plot)
dev.off()
saveRDS(hallm_plot, file = outfile_rds)






###########################################################################
# Visualise results for the hallmark gene sets that are enriched on the 
# skin Treg side with q < 0.01. Only visualise enrichments on the 
# skin Treg side.
###########################################################################

# Restrict data to those gene sets that are enriched on the skin Treg side with
# q < 0.01 and only keep enrichment results for the skin Treg side.
data_to_plot_2 <- data_to_plot[
  data_to_plot$Fisher_pval_adj_BH < 0.01 & 
    data_to_plot$link_set == "hypometh_hyperacc_hyperexpr_skin", 
]

# Generate the plot.
# Order the hallmark terms once by P value and once by odds ratio.
order_cols <- c("logp", "logor")
hallm_plot_2 <- lapply(c("logp", "logor"), FUN = function(x) {
  data_ordered <- data_to_plot_2
  data_ordered$Hallmark_set <- factor(
    data_ordered$Hallmark_set, 
    levels = data_ordered$Hallmark_set[order(data_ordered[, x])]
  )
  ggplot(data_ordered) +
    aes(x = link_set, y = Hallmark_set, fill = logor, 
        size = logp_clipped_w_offset) +
    scale_x_discrete(breaks = "hypometh_hyperacc_hyperexpr_skin", 
                     labels = "Skin\nTreg",
                     name = "Hypometh.,\nhyperacc and\nhyperexpr. in") +
    scale_fill_viridis_c(option = "E", name = "log10(OR)") +
    scale_size_area(breaks = 0.75:4.75, labels = c("0", "1", "2", "3", ">=4"),
                    limits = c(0.75, 4.75), name = "-log10(q)") +
    geom_point(shape = 21, stroke = 1) +
    ylab("Hallmark gene set") +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"))
})
names(hallm_plot_2) <- order_cols

# Save the modified plots showing enrichments.
outfile_pdf <- paste0(
  plot_outdir, 
  "/multiom_skintregsig_hallm_enr_q_below_1perc_skin_side.pdf"
)
outfile_rds <- paste0(
  plot_rds_outdir, 
  "/multiom_skintregsig_hallm_enr_q_below_1perc_skin_side.rds"
)
pdf(outfile_pdf, width = 6, height = 8)
print(hallm_plot_2)
dev.off()
saveRDS(hallm_plot_2, file = outfile_rds)




#############################################################################
# Generate a table with the descriptions of the Hallmark gene sets.
#############################################################################

# Collect the descriptions of the gene sets.
gene_set_descriptions <- do.call(rbind, 
                                 lapply(hallmark_names, FUN = function(x) {
  description <- unique(
    hallmark_genesets$gs_description[hallmark_genesets$gs_name == x]
  )
  assert(length(description) == 1)
  temp_df <- data.frame(Hallmark_set = x,
                        Description_according_to_source = description)
  return(temp_df)
}))

# Save the collected descriptions.
descriptions_outfile <- paste0(
  location, 
  "/treg_hierarchies/multiom_skintregsig_hallm_enr_geneset_descr.txt"
)
write.table(gene_set_descriptions, file = descriptions_outfile, sep = "\t", 
            col.names = T, row.names = F, quote = F)
