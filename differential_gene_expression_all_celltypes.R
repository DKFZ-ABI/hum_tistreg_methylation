# This script performs differential gene expression analysis for all T cell
#   types that were also subjected to Tagmentation-based whole-genome bisulphite
#   sequencing.
# Author: Niklas Beumer



# Load required packages.
library(DESeq2)
library(testit)
library(ggplot2)


# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Read in the sample mapping files.
rna_sample_mapping_path <- paste0(location, "/sample_mapping_rnaseq_only_biol_rep.txt")
rna_sample_mapping <- read.table(
  rna_sample_mapping_path,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)
wgbs_sample_mapping_path <- paste0(location, "/sample_mapping_twgbs.txt")
wgbs_sample_mapping <- read.table(
  wgbs_sample_mapping_path,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)

# Specify the prefix of the path to the directory where PCA results will be 
# saved.
output_pref <- paste(location,
                     "RNASeq/analysis_results",
                     format(Sys.time(), "%Y-%m-%d"),
                     sep = "/")

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis",
                     format(Sys.time(), "%m-%d-%y"),
                     sep = "/")
if (!dir.exists(plot_outdir)) {
  dir.create(plot_outdir)
}








#################################################################
# Prepare raw gene expression counts.
#################################################################

# Extract the paths to raw gene expression counts.
count_paths <- unique(rna_sample_mapping$Raw_expr_path)

# Read in all files with raw count data and merge them into a common data frame.
counts_single <- lapply(
  count_paths,
  FUN = function(x) {
    read.table(x, header = T, stringsAsFactors = F)
  }
)
for (i in 1:(length(counts_single) - 1)) {
  assert(all(counts_single[[i]]$EnsemblID == counts_single[[i + 1]]$EnsemblID))
  assert(all(counts_single[[i]]$GeneSymbol == counts_single[[i + 1]]$GeneSymbol))
}
gene_name_cols <- counts_single[[1]][, c("EnsemblID", "GeneSymbol")]
counts <- do.call(cbind, lapply(
  counts_single,
  FUN = function(x) {
    x[, -which(colnames(x) %in% c("EnsemblID", "GeneSymbol"))]
  }
))
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
print(paste(
  nrow_after,
  "out of",
  nrow_before,
  "genes are retained in this process."
))
print(paste(
  "The removed gene symbols are",
  paste(rem_gene_symbs, collapse = ", ")
))
cat("\n")

# Correct one column name where "CCR8" was replaced by "CCRR"
current_colnames <- colnames(counts_unique)
current_colnames[current_colnames == "p024_RNAseq_CCRR_Treg_Yellow5_S25_R1"] <- 
  "p024_RNAseq_CCR8_Treg_Yellow5_S25_R1"
colnames(counts_unique) <- current_colnames




########################################################################
# Perform differential gene expression analysis.
########################################################################

# Extract the cell types that were covered by WGBS.
cell_types_to_analyse <- unique(wgbs_sample_mapping$Cell_type)
cell_types_for_outfile <- sapply(
  cell_types_to_analyse,
  FUN = function(x) {
    paste(strsplit(x, split = " ")[[1]], collapse = "_")
  }
)

# Check that these cell types are present in the RNS-seq sample mapping file.
assert(all(cell_types_to_analyse %in% rna_sample_mapping$Cell_type))

# Iterate over all cell type pairs that were also considered for differential 
# methylation analysis.
for (i in 1:(length(cell_types_to_analyse) - 1)) {
  for (j in (i + 1):length(cell_types_to_analyse)) {
    # Identify the cell types to use in this comparison.
    cell_type_1 <- cell_types_to_analyse[i]
    cell_type_1_for_outfile <- cell_types_for_outfile[i]
    cell_type_2 <- cell_types_to_analyse[j]
    cell_type_2_for_outfile <- cell_types_for_outfile[j]
    
    # Print a status message.
    print(
      paste(
        "Analysing differential gene expression between",
        cell_type_1,
        "and",
        cell_type_2
      )
    )
    
    # Find the sample names corresponding to the two cell types.
    # Also keep track of which sample belongs to which cell type.
    lines_to_consider <- rna_sample_mapping$Cell_type %in% c(cell_type_1, cell_type_2)
    samples_to_consider <- rna_sample_mapping$Sample[lines_to_consider]
    corresp_cell_types <- rna_sample_mapping$Cell_type_wo_blank[lines_to_consider]
    
    # Prepare the data to be used by DESeq2.
    counts_samples <- as.matrix(counts_unique[, samples_to_consider])
    rownames(counts_samples) <- counts_unique$GeneSymbol
    # counts_samples <- counts_samples[1:1000, ] # For debugging purposes.
    cell_type_df <- data.frame(Cell_type = factor(corresp_cell_types),
                               row.names = samples_to_consider)
    Deseq2_input <- DESeqDataSetFromMatrix(
      countData = counts_samples,
      colData = cell_type_df,
      design = ~ Cell_type
    )
    
    # Specify which one of the two cell types will be used as the "reference".
    # With respect to the cell type name, account for a disparity with respect to the CCR8+ cell type.
    Deseq2_input$Cell_type <- relevel(
      Deseq2_input$Cell_type,
      ref = sub("CCR8+", "CCR8", cell_type_1_for_outfile, fixed = T)
    )
    
    
    # Run DESeq2 analysis.
    Deseq2_input_analysed <- DESeq(Deseq2_input)
    
    # Extract the results of the differential analysis.
    # Set the significance threshold used in computing expression level cut-offs to 0.05.
    diff_gene_expr <- results(Deseq2_input_analysed, alpha = 0.05)
    
    # Save the object that contains the DESeq2 results as an RDS file.
    results_outfile <- paste0(
      output_pref,
      "_diff_gene_expr_DESEq2_",
      cell_type_1_for_outfile,
      cell_type_2_for_outfile,
      "_results_object.rds"
    )
    saveRDS(diff_gene_expr, file = results_outfile)
    
    # Generate an MA plot.
    ma_plot_outfile <- paste0(
      plot_outdir,
      "/diff_gene_expr_DESeq2_",
      cell_type_1_for_outfile,
      cell_type_2_for_outfile,
      "_MA_plot.pdf"
    )
    pdf(ma_plot_outfile, width = 3, height = 3)
    plotMA(diff_gene_expr, alpha = 0.05)
    dev.off()
    
    # Extract the table cpmtaining the differential gene expression results and remove genes
    # fow shih the adjusted P value or the Log2FC is NA. Then, order by Log2FC
    diff_gene_expr_df <- as.data.frame(diff_gene_expr)
    diff_gene_expr_df <- diff_gene_expr_df[!(is.na(diff_gene_expr_df$log2FoldChange) |
                                               is.na(diff_gene_expr_df$padj)), ]
    diff_gene_expr_df <- diff_gene_expr_df[order(diff_gene_expr_df$log2FoldChange), ]
    
    # Add information on whether a gene is significantly differentially expressed.
    # Differential expression is considered significant if the BH-adjusted P value is
    # below 0.05 and the absolute log2FC is above 0.5.
    diff_gene_expr_df$significant <- abs(diff_gene_expr_df$log2FoldChange) > 0.5 &
      diff_gene_expr_df$padj < 0.05
    
    # Print a message stating how many genes are found to display statistical significance.
    num_signif <- length(which(diff_gene_expr_df$significant))
    print(paste(num_signif, "genes display statistical significance"))
    cat("\n")
    
    # Save the table containing genes with available log2FCs and P values together with the
    # significance classification in text format.
    filtered_diff_outfile <- paste0(
      output_pref,
      "_diff_gene_expr_DESEq2_",
      cell_type_1_for_outfile,
      cell_type_2_for_outfile,
      "_results_filtered_with_significance.txt"
    )
    write.table(diff_gene_expr_df, file = filtered_diff_outfile, quote = F)
    
    # Add a column showing the negative decadic logarithm of the BH-adjusted P value.
    diff_gene_expr_df$logp <- (-1) * log10(diff_gene_expr_df$padj)
    
    # Generate a volcano plot displaying the results of the differential gene expression analysis.
    volcano_plot <- ggplot(diff_gene_expr_df) +
      aes(x = log2FoldChange, y = logp, colour = significant) +
      scale_colour_manual(breaks = c(T, F), values = c("red", "grey")) +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.05)), name = "-log10(q)") +
      geom_point(size = 0.5) +
      geom_hline(yintercept = (-1) * log10(0.05),
                 linetype = "dotted") +
      geom_vline(xintercept = -0.5, linetype = "dotted") +
      geom_vline(xintercept = 0.5, linetype = "dotted") +
      guides(colour = "none") +
      xlab("log2(fold change)") +
      ggtitle(paste(cell_type_1, "vs.", cell_type_2)) +
      theme_classic() +
      theme(
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        plot.title = element_text(size = 10)
      )
    
    # Save the volcano plot.
    volc_outfile_pdf <- paste0(
      plot_outdir,
      "/diff_gene_expr_DESeq2_",
      cell_type_1_for_outfile,
      cell_type_2_for_outfile,
      "_volcano_plot.pdf"
    )
    volc_outfile_rds <- paste0(
      plot_rds_outdir,
      "/diff_gene_expr_DESeq2_",
      cell_type_1_for_outfile,
      cell_type_2_for_outfile,
      "_volcano_plot.rds"
    )
    pdf(volc_outfile_pdf, width = 3.5, height = 3.5)
    print(volcano_plot)
    dev.off()
    saveRDS(volcano_plot, file = volc_outfile_rds)
    
  }
}
