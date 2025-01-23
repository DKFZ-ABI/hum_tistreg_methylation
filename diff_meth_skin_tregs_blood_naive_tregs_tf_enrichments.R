# This script analyses what transcription factor binding sites are enriched
#   in regions that are differentially methylated between skin Treg cells and
#   blood naive Treg cells (homer, comparisons against whole genome).
# Author: Niklas Beumer



# Load required package(s)
library(ggplot2)
library(parallel)
library(testit)


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

# Read in the list of differentially methylated regions between skin Tregs and 
# blood naive Tregs.
meth_reg_file <- paste0(
  location,
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_gap_0.15_sterrratio_0.5_signature_regions_filtered.txt"
)
meth_reg <- read.table(
  meth_reg_file,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)
# meth_reg <- meth_reg[1:5000, ] # For debugging purposes.








##########################################################################
# Analyse enrichment of transcription factor binding sites in the
# regions, stratefied by methylation tendency.
##########################################################################

# Iterate over the two methylation tendencies.
for (tendency in c("Blood_naive_Treg__hypomethylation",
                   "Skin_Treg__hypomethylation")) {
  # Extract the region set to use.
  region_set_to_use <- meth_reg[meth_reg$automatic_annotation == tendency, ]
  
  # Generate a directory to hold homer outputs and temporary inputs
  # (if it doesn't already exist).
  homer_out_dir <- paste0(
    location,
    "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_",
    tendency,
    "_homer_against_genome"
  )
  if (!(dir.exists(homer_out_dir))) {
    dir.create(homer_out_dir)
  }
  
  # Generate a temporary bed file for the peak set.
  # Decrease start positions by 1 nucleotide to make intervals 0-based and half-open.
  bed_prep_df <- region_set_to_use[, 1:3]
  bed_prep_df$seqnames <- paste0("chr", bed_prep_df$seqnames)
  bed_prep_df$start <- bed_prep_df$start - 1
  bed_prep_df$Region_ID <- region_set_to_use$region_ID
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
