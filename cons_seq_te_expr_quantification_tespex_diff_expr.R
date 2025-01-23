# This script analyses differential TE expression between cell types based on TE 
#   expression counts from TEspeX.
# Author: Niklas Beumer



# Load required packages.
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(testit)


# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the RNA sample mapping file.
sample_mapping_file <- paste0(location, 
                              "/sample_mapping_rnaseq_only_biol_rep.txt")
sample_mapping <- read.table(sample_mapping_file, header = T, 
                             stringsAsFactors = F, sep = "\t")

# Read in the mapping statistics from TEspeX.
mapping_stats_file <- paste0(
  location, 
  "/te_analysis/cons_seq_tespex_results_bulk_rna_all/mapping_stats.txt"
)
mapping_stats <- read.table(mapping_stats_file, header = T, 
                            stringsAsFactors = F)

# Read in the expression counts from TEspeX.
tespex_counts_file <- paste0(
  location, 
  "/te_analysis/cons_seq_tespex_results_bulk_rna_all/outfile.txt"
)
tespex_counts <- read.table(tespex_counts_file, header = T, 
                            stringsAsFactors = F, sep = "\t")

# Define the relevant cell types.
relevant_cell_types <- c("Skin Treg", "Skin Tconv", "Blood CCR8+ Treg", 
                         "Blood naive Treg", "Blood naive Tconv")
cell_types_for_outfile <- gsub(" ", "_", relevant_cell_types)






#############################################################################
# Perform differential TE expression analysis.
#############################################################################

# Check that all TE names are unique.
assert(length(unique(tespex_counts$TE)) == nrow(tespex_counts))

# Correct one Fastq name where "CCR8" was replaced by "CCRR"
colnames(tespex_counts) <- gsub("CCRR", "CCR8", colnames(tespex_counts))
mapping_stats$SRR <- gsub("CCRR", "CCR8", mapping_stats$SRR)

# Iterate over all cell type pairs.
for (i in 1:(length(relevant_cell_types)-1)) {
  for (j in (i+1):length(relevant_cell_types)) {
    
    # Identify the cell types to use in this comparison.
    cell_type_1 <- relevant_cell_types[i]
    cell_type_1_for_outfile <- cell_types_for_outfile[i]
    cell_type_2 <- relevant_cell_types[j]
    cell_type_2_for_outfile <- cell_types_for_outfile[j]
    
    # Print a status message.
    print(paste("Analysing differential TE expression between", 
                cell_type_1, "and", cell_type_2))
    
    # Find the sample names corresponding to the two cell types.
    # Also keep track of which sample belongs to which cell type.
    lines_to_consider <- sample_mapping$Cell_type %in% c(cell_type_1, 
                                                         cell_type_2)
    samples_to_consider <- sample_mapping$Sample[lines_to_consider]
    corresp_cell_types <- sample_mapping$Cell_type_wo_blank[lines_to_consider]
    
    # Prepare the data to be used by DESeq2.
    cols_to_consider <- sapply(samples_to_consider, FUN = function(samplename) {
      grep(samplename, colnames(tespex_counts))
    })
    counts_samples <- tespex_counts[, cols_to_consider]
    rownames(counts_samples) <- tespex_counts$TE
    cell_type_df <- data.frame(Cell_type = factor(corresp_cell_types), 
                               row.names = samples_to_consider)
    Deseq2_input <- DESeqDataSetFromMatrix(countData = counts_samples,
                                           colData = cell_type_df,
                                           design = ~ Cell_type)
    
    # Specify which one of the two cell types will be used as the "reference".
    # With respect to the cell type name, account for a disparity with respect 
    # to the CCR8+ cell type.
    Deseq2_input$Cell_type <- relevel(Deseq2_input$Cell_type, 
                                      ref = sub("CCR8+", 
                                                "CCR8", 
                                                cell_type_1_for_outfile, 
                                                fixed = T))
    
    
    # Run DESeq2 analysis.
    # Use the total number of mapped reads (as returned by TEspeX) as size 
    # factors, as recommended here: https://github.com/fansalon/TEspeX; 
    # 03-Jan-2024. Use factor 10^(-7) to bring values close to 1. Otherwise, 
    # the dispersion estimation will return an error.
    sizeFactors(Deseq2_input) <- 
      sapply(colnames(Deseq2_input), FUN = function(samplename) {
        corresp_mappedreads <- mapping_stats$mapped[
          grep(samplename, mapping_stats$SRR)]
      }) * 10^(-7)
    Deseq2_input_analysed <- estimateDispersions(Deseq2_input)
    Deseq2_input_analysed <- nbinomWaldTest(Deseq2_input_analysed)
    
    # Extract the results of the differential analysis.
    # Set the significance threshold used in computing expression level cut-offs 
    # to 0.05.
    diff_expr <- results(Deseq2_input_analysed, alpha = 0.05)
    
    # Save the object that contains the DESeq2 results as an RDS file.
    results_outfile <- paste0(
      location, 
      "/te_analysis/cons_seq_tespex_diff_expr_DESEq2_", 
      cell_type_1_for_outfile, 
      cell_type_2_for_outfile, 
      "_results_object.rds"
    )
    saveRDS(diff_expr, file = results_outfile)
    
    # Generate an MA plot.
    ma_plot_outfile <- paste0(
      plot_outdir, 
      "/te_analysis_cons_seq_DESeq2_", 
      cell_type_1_for_outfile, 
      cell_type_2_for_outfile, 
      "_MA_plot.pdf"
    )
    pdf(ma_plot_outfile, width = 3, height = 3)
    plotMA(diff_expr, alpha = 0.05)
    dev.off()
    
    # Extract the table cpmtaining the differential expression results and 
    # remove TEs for which the adjusted P value or the Log2FC is NA. Then, order 
    # by Log2FC.
    diff_expr_df <- as.data.frame(diff_expr)
    diff_expr_df <- diff_expr_df[!(is.na(diff_expr_df$log2FoldChange) | 
                                     is.na(diff_expr_df$padj)), ]
    diff_expr_df <- diff_expr_df[order(diff_expr_df$log2FoldChange), ]
    
    # Add information on whether a TE is significantly differentially expressed.
    # Differential expression is considered significant if the BH-adjusted P 
    # value is below 0.05 and the absolute log2FC is above 0.5.
    diff_expr_df$significant <- abs(diff_expr_df$log2FoldChange) > 0.5 & 
      diff_expr_df$padj < 0.05
    
    # Print a message stating how many TEs are found to display statistical 
    # significance.
    num_signif <- length(which(diff_expr_df$significant))
    print(paste(num_signif, "TEs display statistical significance"))
    cat("\n")
    
    # Save the table containing TEs with available log2FCs and P values 
    # together with the significance classification in text format.
    filtered_diff_outfile <- paste0(
      location, 
      "/te_analysis/cons_seq_tespex_diff_expr_DESEq2_", 
      cell_type_1_for_outfile, 
      cell_type_2_for_outfile, 
      "_results_filtered_with_significance.txt"
    )
    write.table(diff_expr_df, file = filtered_diff_outfile, quote = F)
    
    # Add a column showing the negative decadic logarithm of the BH-adjusted 
    # P value.
    diff_expr_df$logp <- (-1) * log10(diff_expr_df$padj)
    
    # Prepare data for labelling all significant TEs in the volcano plot.
    diff_expr_df_signif <- diff_expr_df[diff_expr_df$significant, ]
    diff_expr_df_signif$TE <- rownames(diff_expr_df_signif)
    
    # Generate a volcano plot displaying the results of the differential TE 
    # expression analysis.
    volcano_plot <- ggplot(diff_expr_df) +
      aes(x = log2FoldChange, y = logp, colour = significant) +
      scale_colour_manual(breaks = c(T, F), values = c("red", "grey")) +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.05)), 
                         name = "-log10(q)") +
      geom_point(size = 0.5) +
      geom_hline(yintercept = (-1) * log10(0.05), linetype = "dotted") +
      geom_vline(xintercept = -0.5, linetype = "dotted") +
      geom_vline(xintercept = 0.5, linetype = "dotted") +
      geom_text_repel(data = diff_expr_df_signif, mapping = aes(label = TE),
                      colour = "black", size = 2) +
      guides(colour = "none") +
      xlab("log2(fold change)") +
      ggtitle(paste(cell_type_1, "vs.", cell_type_2)) +
      theme_classic() +
      theme(axis.text = element_text(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            plot.title = element_text(size = 10))
    
    # Save the volcano plot.
    volc_outfile_pdf <- paste0(
      plot_outdir, 
      "/te_analysis_cons_seq_DESeq2_", 
      cell_type_1_for_outfile, 
      cell_type_2_for_outfile, 
      "_volcano_plot.pdf"
    )
    volc_outfile_rds <- paste0(
      plot_rds_outdir, 
      "/te_analysis_cons_seq_DESeq2_", 
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
