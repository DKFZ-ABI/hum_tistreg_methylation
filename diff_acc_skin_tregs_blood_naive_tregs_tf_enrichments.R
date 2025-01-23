# This script analyses what transcription factor binding sites are enriched
#   in regions that are differentially accessible between skin Treg cells and
#   blood naive Treg cells (homer, comparisons against whole genome).
# Author: Niklas Beumer



# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis",
                     format(Sys.time(), "%m-%d-%y"),
                     sep = "/")
if (!dir.exists(plot_outdir)) {
  dir.create(plot_outdir)
}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the list of differentially accessible peaks between skin Tregs and 
# blood naive Tregs.
acc_reg_file <- paste0(
  location,
  "/differential_accessibility/diff_acc_analysis_Skin_TregBlood_naive_Treg_results_with_significance.txt"
)
acc_reg <- read.table(acc_reg_file, header = T, stringsAsFactors = F)

# Restrict the differentially accessible peaks to those that display statistical 
# significance after multiple testing correction.
acc_reg <- acc_reg[acc_reg$p_val_adj < 0.05, ]






##########################################################################
# Analyse enrichment of transcription factor binding sites in the
# peaks, stratefied by accessibility tendency.
##########################################################################

# Define file snippets.
file_snippets <- c("Blood_naive_Treg__hyperaccessibility",
                   "Skin_Treg__hyperaccessibility")
names(file_snippets) <- c(">", "<")

# Iterate over the two methylation tendencies.
for (tendency in names(file_snippets)) {
  # Extract the peak set to use.
  peak_set_to_use <- rownames(acc_reg)[
    eval(parse(text = paste("acc_reg$avg_log2FC", tendency, "0")))
  ]
  
  # Generate a directory to hold homer outputs and temporary inputs
  # (if it doesn't already exist).
  homer_out_dir <- paste0(
    location,
    "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_",
    file_snippets[tendency],
    "_homer_against_genome"
  )
  if (!(dir.exists(homer_out_dir))) {
    dir.create(homer_out_dir)
  }
  
  # Generate a temporary bed file for the peak set.
  peaks_split <- lapply(
    peak_set_to_use,
    FUN = function(x) {
      strsplit(x, split = "-")[[1]]
    }
  )
  bed_prep_df <- data.frame(
    chr = sapply(
      peaks_split,
      FUN = function(x) {
        x[1]
      }
    ),
    start = sapply(
      peaks_split,
      FUN = function(x) {
        x[2]
      }
    ),
    end = sapply(
      peaks_split,
      FUN = function(x) {
        x[3]
      }
    )
  )
  bed_prep_df$Region_ID <- peak_set_to_use
  bed_prep_df$Empty_col_5 <- "mock"
  bed_prep_df$Strand = "."
  # bed_prep_df <- bed_prep_df[1:20, ] # For debugging purposes.
  bed_outfile <- paste0(homer_out_dir, "/temp_region_file.bed")
  write.table(
    bed_prep_df,
    file = bed_outfile,
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
  )
  
  # Generate a command string for homer analysis.
  command_string <- paste0(
    #-- Load homer module
    "module load homer/4.11; ",
    
    #-- Call homer's motif enrichment pipeline.
    "findMotifsGenome.pl ",
    bed_outfile,
    # Path to bed file containing region locations.
    " hg19 ",
    # Genome version
    homer_out_dir,
    # Path to directory where homer outputs will appear.
    " -size given ",
    # Use the full regions, not only parts of them.
    # " -size 50 ", # For debugging purposes.
    "-p 10 ",
    # 10 CPUs
    "-preparsedDir ",
    homer_out_dir,
    # Path to directory where homer will save parsed files.
    
    #-- Capture stderr next to stout.
    " 2>&1; ",
    
    #-- Remove temporary bed file.
    "rm ",
    bed_outfile
  )
  
  # Run homer analysis.
  homer_logs <- system(command_string, intern = T)
  
  # Save the messages from the homer process.
  homer_logs_outfile <- paste0(homer_out_dir, "/homer_logs.txt")
  write(homer_logs, file = homer_logs_outfile)
  
}
