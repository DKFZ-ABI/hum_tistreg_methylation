# This script analyses differential gene expression between blood CD45RA+
#   Treg cells and blood CD45RA+ Tconv cells using donor ID as an additional
#   covariate.
# Author: Niklas Beumer



# Load required packages.
library(DESeq2)
library(testit)
library(ggplot2)


# Define a location on /xxx.
location <- "/xxx/nbeumer/hm_treg_bs_rgnsbg"

# Read in the sample mapping files.
rna_sample_mapping_path <- paste0(location, 
                                  "/sample_mapping_rnaseq_only_biol_rep.txt")
rna_sample_mapping <- read.table(rna_sample_mapping_path, header = T, 
                                 stringsAsFactors = F, sep = "\t")

# Specify the prefix of the path to the directory where PCA results will be saved.
output_pref <- paste(location, 
                     "RNASeq/analysis_results", 
                     format(Sys.time(), "%Y-%m-%d"), 
                     sep = "/")

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/yyy/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the two cell types to compare.
cell_type_1 <- "Blood naive Treg"
cell_type_2 <- "Blood naive Tconv"







#################################################################
# Prepare raw gene expression counts.
#################################################################

# Extract the paths to raw gene expression counts.
count_paths <- unique(rna_sample_mapping$Raw_expr_path[
  rna_sample_mapping$Cell_type %in% c(cell_type_1, cell_type_2)
])

# Read in all files with raw count data and merge them into a common data frame.
counts_single <- lapply(count_paths, FUN = function(x) {
  read.table(x, header = T, stringsAsFactors = F)
})
for (i in 1:(length(counts_single)-1)) {
  assert(all(counts_single[[i]]$EnsemblID == counts_single[[i+1]]$EnsemblID))
  assert(all(counts_single[[i]]$GeneSymbol == counts_single[[i+1]]$GeneSymbol))
}
gene_name_cols <- counts_single[[1]][, c("EnsemblID", "GeneSymbol")]
counts <- do.call(cbind, lapply(counts_single, FUN = function(x) {
  x[, -which(colnames(x) %in% c("EnsemblID", "GeneSymbol"))]
}))
counts <- cbind(gene_name_cols, counts)

# Remove all genes that have a non-unique gene symbol.
# This removes 1843 genes mapping to 209 gene symbols.
print("Removing genes that have a non-unique gene symbol.")
nrow_before <- nrow(counts)
gene_symb_table <- table(counts$GeneSymbol)
unique_symbs <- names(gene_symb_table)[gene_symb_table == 1]
rem_gene_symbs <- names(gene_symb_table)[gene_symb_table > 1]
counts_unique <- counts[counts$GeneSymbol %in% unique_symbs, ]
nrow_after <- nrow(counts_unique)
print(paste(nrow_after, "out of", nrow_before, "genes are retained in this process."))
print(paste("The removed gene symbols are", paste(rem_gene_symbs, collapse = ", ")))
cat("\n")






########################################################################
# Perform differential gene expression analysis.
########################################################################

# Find the sample names corresponding to the two cell types.
# Also keep track of which sample belongs to which cell type.
lines_to_consider <- rna_sample_mapping$Cell_type %in% 
  c(cell_type_1, cell_type_2)
samples_to_consider <- rna_sample_mapping$Sample[lines_to_consider]
corresp_cell_types <- rna_sample_mapping$Cell_type_wo_blank[lines_to_consider]
corresp_donors <- rna_sample_mapping$Patient_ID_readable[lines_to_consider]
names(corresp_donors) <- samples_to_consider

# Add donor IDs that were missing from the table.
corresp_donors["p024_RNAseq_Naive_Treg_Green1_S1_R1"] <- "Donor A"
corresp_donors["p024_RNAseq_Naive_Treg_Red1_S11_R1"] <- "Donor B"

# Prepare the data to be used by DESeq2.
counts_samples <- as.matrix(counts_unique[, samples_to_consider])
rownames(counts_samples) <- counts_unique$GeneSymbol
# counts_samples <- counts_samples[1:1000, ] # For debugging purposes.
covariate_df <- data.frame(Cell_type = factor(corresp_cell_types), 
                           Donor = factor(gsub(" ", "_", corresp_donors)),
                           row.names = samples_to_consider)
Deseq2_input <- DESeqDataSetFromMatrix(countData = counts_samples,
                                       colData = covariate_df,
                                       design = ~ Donor + Cell_type)

# Specify which one of the two cell types will be used as the "reference".
Deseq2_input$Cell_type <- relevel(Deseq2_input$Cell_type, 
                                  ref = gsub(" ", "_", cell_type_1))

# Run DESeq2 analysis.
Deseq2_input_analysed <- DESeq(Deseq2_input)

# Extract the results of the differential analysis.
# Set the significance threshold used in computing expression level cut-offs to 0.05.
diff_gene_expr <- results(Deseq2_input_analysed, alpha = 0.05)

# Save the object that contains the DESeq2 results as an RDS file.
results_outfile <- paste0(output_pref, 
                          "_diff_gene_expr_DESEq2_", 
                          gsub(" ", "_", cell_type_1), 
                          gsub(" ", "_", cell_type_2), 
                          "_w_donor_as_covariate_results_object.rds")
saveRDS(diff_gene_expr, file = results_outfile)

# Generate an MA plot.
ma_plot_outfile <- paste0(plot_outdir, 
                          "/diff_gene_expr_DESeq2_", 
                          gsub(" ", "_", cell_type_1), 
                          gsub(" ", "_", cell_type_2), 
                          "_w_donor_as_covariate_MA_plot.pdf")
pdf(ma_plot_outfile, width = 3, height = 3)
plotMA(diff_gene_expr, alpha = 0.05)
dev.off()

# Extract the table containing the differential gene expression results and 
# remove genes for which the adjusted P value or the Log2FC is NA. Then, order 
# by Log2FC
diff_gene_expr_df <- as.data.frame(diff_gene_expr)
diff_gene_expr_df <- diff_gene_expr_df[
  !(is.na(diff_gene_expr_df$log2FoldChange) | is.na(diff_gene_expr_df$padj)), 
]
diff_gene_expr_df <- diff_gene_expr_df[
  order(diff_gene_expr_df$log2FoldChange), 
]

# Add information on whether a gene is significantly differentially expressed.
# Differential expression is considered significant if the BH-adjusted P value 
# is below 0.05 and the absolute log2FC is above 0.5.
diff_gene_expr_df$significant <- abs(diff_gene_expr_df$log2FoldChange) > 0.5 & 
  diff_gene_expr_df$padj < 0.05

# Print a message stating how many genes are found to display statistical 
# significance.
num_signif <- length(which(diff_gene_expr_df$significant))
print(paste(num_signif, "genes display statistical significance"))
cat("\n")

# Save the table containing genes with available log2FCs and P values together 
# with the significance classification in text format.
filtered_diff_outfile <- paste0(
  output_pref, 
  "_diff_gene_expr_DESEq2_",
  gsub(" ", "_", cell_type_1), 
  gsub(" ", "_", cell_type_2), 
  "_w_donor_as_covariate_results_filtered_with_significance.txt"
)
write.table(diff_gene_expr_df, file = filtered_diff_outfile, quote = F)

# Add a column showing the negative decadic logarithm of the BH-adjusted P 
# value.
diff_gene_expr_df$logp <- (-1) * log10(diff_gene_expr_df$padj)

# Generate a volcano plot displaying the results of the differential gene 
# expression analysis.
volcano_plot <- ggplot(diff_gene_expr_df) +
  aes(x = log2FoldChange, y = logp, colour = significant) +
  scale_colour_manual(breaks = c(T, F), values = c("red", "grey")) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05)), 
                     name = "-log10(q)") +
  geom_point(size = 0.5) +
  geom_hline(yintercept = (-1) * log10(0.05), linetype = "dotted") +
  geom_vline(xintercept = -0.5, linetype = "dotted") +
  geom_vline(xintercept = 0.5, linetype = "dotted") +
  guides(colour = "none") +
  xlab("log2(fold change)") +
  ggtitle(paste(cell_type_1, "vs.", cell_type_2)) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        plot.title = element_text(size = 10))

# Save the volcano plot.
volc_outfile_pdf <- paste0(plot_outdir, 
                           "/diff_gene_expr_DESeq2_", 
                           gsub(" ", "_", cell_type_1), 
                           gsub(" ", "_", cell_type_2), 
                           "_w_donor_as_covariate_volcano_plot.pdf")
volc_outfile_rds <- paste0(plot_rds_outdir, 
                           "/diff_gene_expr_DESeq2_", 
                           gsub(" ", "_", cell_type_1), 
                           gsub(" ", "_", cell_type_2), 
                           "_w_donor_as_covariate_volcano_plot.rds")
pdf(volc_outfile_pdf, width = 3.5, height = 3.5)
print(volcano_plot)
dev.off()
saveRDS(volcano_plot, file = volc_outfile_rds)