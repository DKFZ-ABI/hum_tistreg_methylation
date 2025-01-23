# This script analyses amplicon data for the comparison between blood naive Treg 
#   cells and blood nave Tconv cells and the pilot run therein.
# Author: Niklas Beumer
# Requires the modules bismark/0.20.0, bowtie2/2.3.5.1, samtools/1.11.0 and 
#   bedtools/2.24.0. Requires cutadapt to be installed (This script was
#   developed with cutadapt v4.8 in Python v3.8.3).



# Load required package(s).
library(testit)
library(GenomicRanges)
library(ggplot2)


# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Specify the name of the amplicon experiment that is analysed here.
experiment_name <- "05_2024_pilot_run"

# Specify the location of Fastq files.
fastq_loc <- paste0(location, "/amplicon_analysis/raw_data/", experiment_name)

# Specify the location of data output.
output_loc <- paste0(location, "/amplicon_analysis/analysis_results/", 
                     experiment_name)

# Specify the samples to consider.
samples <- c("Tamara10", "Tamara70")

# Based on the sample names, specify the names of the Fastq files.
fastq_files <- lapply(samples, FUN = function(x) {
  temp_vect <- paste0(fastq_loc, "/", x, "_", c("R1", "R2"), "_001.fastq.gz")
  names(temp_vect) <- c("R1", "R2")
  return(temp_vect)
})
names(fastq_files) <- samples

# Specify the barcodes used for different cell types and the corresponding cell 
# type names. Note that the barcodes will only show up in Read 1.
barcodes <- c("TATATC", "GAGAGA")
cell_types <- c("Tconv", "Treg")
names(barcodes) <- cell_types

# Specify the location of the reference genome prepared for the "Bismark" 
# aligner.
ref_gen_loc <- paste0(
  location, 
  "/external_data/2024-06-06_GRCh37_reference_genome_from_ensembl_75/raw"
)

# Specify the exact locations of the amplicons.
amplicon_names <- c("R8225_TIGIT", "R8226_TIGIT", "R1198_IL2RA", 
                    "R10871#2_LRRN3", "R10425_SYTL3", "R115_TNFRSF8", 
                    "R12143-1_FOXP3", "FOXP3TSDR_FOXP3", "R12143-2_FOXP3", 
                    "R6632_CTLA4", "R132_TNFRSF1B", "R122_TNFRSF1B")
amplicon_locs <- GRanges(seqnames = c("3", "3", "10", "7", "6", "1", "X", "X", 
                                      "X", "2", "1", "1"), 
                         ranges = IRanges(start = c(114012711, 114014310, 
                                                    6102137, 110731169, 
                                                    159083980, 12188066, 
                                                    49119989, 49117054, 
                                                    49120347, 204736407, 
                                                    12263695, 12231169), 
                                          end = c(114013071, 114014730, 
                                                  6102543, 110731507, 
                                                  159084279, 12188381, 
                                                  49120374, 49117406, 
                                                  49120704, 204736768, 
                                                  12264119, 12231538)))
amplicon_locs$Amplicon_name <- amplicon_names

# Specify adaptor sequences.
r1_adapt_5pr <- "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
r1_5pr_treg <- paste0(r1_adapt_5pr, barcodes["Treg"])
r1_5pr_tconv <- paste0(r1_adapt_5pr, barcodes["Tconv"])
r2_adapt_5pr <- "GACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
r1_adapt_3pr <- toupper(paste0(rev(strsplit(
  gsub("A", "t", gsub("C", "g", gsub("G", "c", gsub("T", "a", r2_adapt_5pr)))),
  split = ""
)[[1]]), collapse = ""))
r2_3pr_treg <- toupper(paste0(rev(strsplit(
  gsub("A", "t", gsub("C", "g", gsub("G", "c", gsub("T", "a", r1_5pr_treg)))),
  split = ""
)[[1]]), collapse = ""))
r2_3pr_tconv <- toupper(paste0(rev(strsplit(
  gsub("A", "t", gsub("C", "g", gsub("G", "c", gsub("T", "a", r1_5pr_tconv)))),
  split = ""
)[[1]]), collapse = ""))









#############################################################################
# Demultiplex reads.
#############################################################################

# Print a status message.
cat("####################################################################\n")
cat("Demultiplexing\n")
cat("####################################################################\n\n")

# Read in the Fastq files. Extract reads, read names and Phred scores. While 
# doing so, check that there are equal numbers of reads for R1 and R2 and that
# they are in the same order.
fastq_data <- lapply(samples, FUN = function(x) {
  curr_filenames <- fastq_files[[x]]
  this_sample_reads <- lapply(names(curr_filenames), FUN = function(y) {
    fastq_raw <- readLines(curr_filenames[y])
    # fastq_raw <- fastq_raw[1:4000] # For debugging purposes.
    read_names <- fastq_raw[seq(1, length(fastq_raw), 4)]
    sequences <- fastq_raw[seq(2, length(fastq_raw), 4)]
    phreds <- fastq_raw[seq(4, length(fastq_raw), 4)]
    assert(length(read_names) == length(sequences))
    assert(length(sequences) == length(phreds))
    temp_list <- list(Read_names = read_names,
                      Sequences = sequences,
                      Phred = phreds)
    return(temp_list)
  })
  names(this_sample_reads) <- names(curr_filenames)
  assert(length(this_sample_reads[[1]]$Read_names) == 
           length(this_sample_reads[[2]]$Read_names))
  r1_split <- lapply(this_sample_reads[[1]]$Read_names, FUN = function(z) {
    strsplit(z, split = ":")[[1]]
  })
  r1_split_2 <- lapply(r1_split, FUN = function(z) {
    z[7] <- strsplit(z[7], split = " ")[[1]][1]
    return(z)
  }) # Removes the information whether this is Read 1 or Read 2.
  r2_split <- lapply(this_sample_reads[[2]]$Read_names, FUN = function(z) {
    strsplit(z, split = ":")[[1]]
  })
  r2_split_2 <- lapply(r2_split, FUN = function(z) {
    z[7] <- strsplit(z[7], split = " ")[[1]][1]
    return(z)
  }) # Removes the information whether this is Read 1 or Read 2.
  assert(all(sapply(1:length(r1_split_2), FUN = function(z) {
    all(r1_split_2[[z]] == r2_split_2[[z]])
  })))
  return(this_sample_reads)
})
names(fastq_data) <- samples

# Print some sample statistics.
cat("### Total reads before multiplexing:\n")
void <- sapply(samples, FUN = function(x) {
  read_num <- length(fastq_data[[x]]$R1$Read_names)
  cat(paste0(x, ": ", read_num, " reads\n"))
})
cat("\n\n\n")

# Perform demultiplexing by extracting reads whose first bases exactly match 
# the specified barcodes. Store the matching reads without their barcodes.
# Note that these barcodes will only be present in Read 1. However, I previously
# showed that Read 1 and Read 2 are in the same order in the Fastq files so I
# can just apply the same slicing to Read 1 and Read 2.
# Print resulting statistics to the output.
cat("### Reads extracted during demultiplexing:\n")
void <- lapply(samples, FUN = function(x) {
  all_match_inds <- unlist(lapply(barcodes, FUN = function(y) {
    barcode_len <- nchar(y)
    matching_read_ind <- which(
      sapply(fastq_data[[x]]$R1$Sequences, FUN = function(z) {
        substr(z, 1, barcode_len) == y
      })
    )
    num_match <- length(matching_read_ind)
    cat(paste0(x, ": ", names(barcodes)[barcodes == y], ": ", num_match, 
                 "/", length(fastq_data[[x]]$R1$Sequences), " reads\n"))
    r1_matching_read_names <- fastq_data[[x]]$R1$Read_names[matching_read_ind]
    r1_matching_seqs <- 
      sapply(fastq_data[[x]]$R1$Sequences[matching_read_ind], 
             FUN = function(z) {
               substr(z, barcode_len + 1, nchar(z))
             })
    r1_matching_phreds <- 
      sapply(fastq_data[[x]]$R1$Phred[matching_read_ind],
             FUN = function(z) {
               substr(z, barcode_len + 1, nchar(z))
             })
    r1_lines_to_write <- 
      unlist(lapply(1:length(r1_matching_read_names), FUN = function(i) {
        unlist(lapply(
          list(r1_matching_read_names, 
               r1_matching_seqs, 
               rep("+", length(r1_matching_read_names)),
               r1_matching_phreds), 
          FUN = function(j) {
            j[i]
          }))
      }))
    r1_outfile <- paste0(
      output_loc, "/", x, "_", "R1_001.fastq.gz_demultiplexed_",  
      names(barcodes)[barcodes == y], ".fastq"
    )
    writeLines(r1_lines_to_write, con = r1_outfile, sep = "\n")
    r2_matching_read_names <- fastq_data[[x]]$R2$Read_names[matching_read_ind]
    r2_matching_seqs <- fastq_data[[x]]$R2$Sequences[matching_read_ind]
    r2_matching_phreds <- fastq_data[[x]]$R2$Phred[matching_read_ind]
    r2_lines_to_write <- 
      unlist(lapply(1:length(r2_matching_read_names), FUN = function(i) {
        unlist(lapply(
          list(r2_matching_read_names, 
               r2_matching_seqs, 
               rep("+", length(r2_matching_read_names)),
               r2_matching_phreds), 
          FUN = function(j) {
            j[i]
          }))
      }))
    r2_outfile <- paste0(
      output_loc, "/", x, "_", "R2_001.fastq.gz_demultiplexed_",  
      names(barcodes)[barcodes == y], ".fastq"
    )
    writeLines(r2_lines_to_write, con = r2_outfile, sep = "\n")
    return(matching_read_ind)
  }))
  unmatch_read_inds <- (1:length(fastq_data[[x]]$R1$Sequences))[-all_match_inds]
  num_unmatch <- length(unmatch_read_inds)
  cat(paste0(x, ": Unmatched: ", num_unmatch, 
               "/", length(fastq_data[[x]]$R1$Sequences), " reads\n"))
  r1_unmatch_read_names <- fastq_data[[x]]$R1$Read_names[unmatch_read_inds]
  r1_unmatch_seqs <- fastq_data[[x]]$R1$Sequences[unmatch_read_inds]
  r1_unmatch_phreds <- fastq_data[[x]]$R1$Phred[unmatch_read_inds]
  r1_lines_to_write <- 
    unlist(lapply(1:length(r1_unmatch_read_names), FUN = function(i) {
      unlist(lapply(
        list(r1_unmatch_read_names, 
             r1_unmatch_seqs, 
             rep("+", length(r1_unmatch_read_names)),
             r1_unmatch_phreds), 
        FUN = function(j) {
          j[i]
        }))
    }))
  r1_outfile <- paste0(
    output_loc, "/", x, "_", "R1_001.fastq.gz_demultiplexed_unmatched.fastq"
  )
  writeLines(r1_lines_to_write, con = r1_outfile, sep = "\n")
  r2_unmatch_read_names <- fastq_data[[x]]$R2$Read_names[unmatch_read_inds]
  r2_unmatch_seqs <- fastq_data[[x]]$R2$Sequences[unmatch_read_inds]
  r2_unmatch_phreds <- fastq_data[[x]]$R2$Phred[unmatch_read_inds]
  r2_lines_to_write <- 
    unlist(lapply(1:length(r2_unmatch_read_names), FUN = function(i) {
      unlist(lapply(
        list(r2_unmatch_read_names, 
             r2_unmatch_seqs, 
             rep("+", length(r2_unmatch_read_names)),
             r2_unmatch_phreds), 
        FUN = function(j) {
          j[i]
        }))
    }))
  r2_outfile <- paste0(
    output_loc, "/", x, "_", "R2_001.fastq.gz_demultiplexed_unmatched.fastq"
  )
  writeLines(r2_lines_to_write, con = r2_outfile, sep = "\n")
})






#############################################################################
# Trim adaptors.
#############################################################################

# Print a status message.
cat("####################################################################\n")
cat("Adaptor trimming\n")
cat("####################################################################\n\n")

# Trim adaptors using cutadapt. Depending on whether the corresponding file 
# contains Treg- or Tconv-derived reads, specify the correct 3' adaptor within 
# read 2. Don't trim the 5' end of read 1 because I already trimmed this part 
# during the demultiplexing step.
void <- lapply(samples, FUN = function(x) {
  void2 <- lapply(cell_types, FUN = function(y) {
    cat(paste0("Processing ", x, "; ", y, "\n"))
    message(paste0("Adaptor trimming: Processing ", x, "; ", y))
    trim_command <- paste0(
      "cutadapt ",
      " -a ", r1_adapt_3pr,
      " -G ", r2_adapt_5pr,
      " -A ", ifelse(y == "Tconv", yes = r2_3pr_tconv, no = r2_3pr_treg),
      " -o ",
      output_loc, 
      "/", 
      x, 
      "_R1_001.fastq.gz_demultiplexed_",  
      y, 
      "_trimmed.fastq",
      " -p ",
      output_loc, 
      "/", 
      x, 
      "_R2_001.fastq.gz_demultiplexed_",  
      y, 
      "_trimmed.fastq ",
      output_loc, 
      "/", 
      x, 
      "_R1_001.fastq.gz_demultiplexed_",  
      y, 
      ".fastq ",
      output_loc, 
      "/", 
      x, 
      "_R2_001.fastq.gz_demultiplexed_",  
      y, 
      ".fastq"
    )
    trim_logs <- system(trim_command, intern = T)
    cat(paste(trim_logs, collapse = "\n"))
    cat("\n\n\n")
  })
})






#############################################################################
# Perform alignments using the "Bismark" aligner.
#############################################################################

# Print a status message.
cat(
  "\n\n\n####################################################################\n"
)
cat("Alignment using the 'Bismark' aligner\n")
cat("####################################################################\n\n")

# Align reads and using the "Bismark" aligner. Do so separately for every cell 
# type in every sample. Make "Bismark" save unmapped and  ambiguously mapped 
# reads in separate files.
# I found out that using samtools/1.11.0 instead of samtools/1.3.1 prevents an
# error, which for instance meant that the alignment Bam file was not created.
# This solution was inspired by 
# https://github.com/FelixKrueger/Bismark/issues/353; 07-Jun 2024.
void <- lapply(samples, FUN = function(x) {
  void2 <- lapply(cell_types, FUN = function(y) {
    cat(paste0("Processing ", x, "; ", y, "\n"))
    message(paste0("Bismark aligner: Processing ", x, "; ", y))
    tempdir <- paste0(output_loc,
                      "/bismark_tempdir_",
                      x,
                      y)
    command_string_aln <- paste0(
      "module load bismark/0.20.0; ",
      "module load bowtie2/2.3.5.1; ",
      "module load samtools/1.11.0; ",
      "bismark ",
      #"--non_directional ",
      "--unmapped ",
      "--ambiguous ", 
      "--output_dir ", 
      output_loc, 
      "/", 
      x, 
      "_",  
      y, 
      "_bismark_alignment ",
      "--temp_dir ",
      tempdir,
      " --genome ",
      ref_gen_loc,
      " -1 ", 
      output_loc, 
      "/", 
      x, 
      "_R1_001.fastq.gz_demultiplexed_",  
      y, 
      "_trimmed.fastq ",
      "-2 ", 
      output_loc, 
      "/", 
      x, 
      "_R2_001.fastq.gz_demultiplexed_",  
      y, 
      "_trimmed.fastq; ",
      "rm -r ",
      tempdir
    )
    bismark_logs_aln <- system(command_string_aln, intern = T)
    cat(paste(bismark_logs_aln, collapse = "\n"))
    cat("\n\n\n")
  })
})

# Independently align unaligned R1 reads using the "Bismark" aligner. Do so 
# separately for every cell type in every sample. Make "Bismark" save unmapped 
# and ambiguously mapped reads in separate files.
void <- lapply(samples, FUN = function(x) {
  void2 <- lapply(cell_types, FUN = function(y) {
    cat(paste0("Realigning leftover R1 reads from ", x, "; ", y, "\n"))
    message(paste0("Bismark aligner: Realigning leftover R1 reads from ", x, 
                   "; ", y))
    tempdir <- paste0(output_loc,
                      "/bismark_tempdir_r1realign_",
                      x,
                      y)
    command_string_aln <- paste0(
      "module load bismark/0.20.0; ",
      "module load bowtie2/2.3.5.1; ",
      "module load samtools/1.11.0; ",
      "bismark ",
      "--unmapped ",
      "--ambiguous ", 
      "--output_dir ", 
      output_loc, 
      "/", 
      x, 
      "_",  
      y, 
      "_bismark_alignment_r1_leftover ",
      "--temp_dir ",
      tempdir,
      " --genome ",
      ref_gen_loc,
      " ",
      output_loc, 
      "/", 
      x, 
      "_",
      y, 
      "_bismark_alignment/",
      x,
      "_R1_001.fastq.gz_demultiplexed_",
      y,
      "_trimmed.fastq_unmapped_reads_1.fq.gz; ",
      "rm -r ",
      tempdir
    )
    bismark_logs_aln <- system(command_string_aln, intern = T)
    cat(paste(bismark_logs_aln, collapse = "\n"))
    cat("\n\n\n")
  })
})

# Independently align unaligned R2 reads using the "Bismark" aligner. Do so 
# separately for every cell type in every sample. Make "Bismark" save unmapped 
# and ambiguously mapped reads in separate files.
void <- lapply(samples, FUN = function(x) {
  void2 <- lapply(cell_types, FUN = function(y) {
    cat(paste0("Realigning leftover R2 reads from ", x, "; ", y, "\n"))
    message(paste0("Bismark aligner: Realigning leftover R2 reads from ", x, 
                   "; ", y))
    tempdir <- paste0(output_loc,
                      "/bismark_tempdir_r2realign_",
                      x,
                      y)
    command_string_aln <- paste0(
      "module load bismark/0.20.0; ",
      "module load bowtie2/2.3.5.1; ",
      "module load samtools/1.11.0; ",
      "bismark ",
      "--unmapped ",
      "--ambiguous ",
      "--output_dir ", 
      output_loc, 
      "/", 
      x, 
      "_",  
      y, 
      "_bismark_alignment_r2_leftover ",
      "--temp_dir ",
      tempdir,
      " --genome ",
      ref_gen_loc,
      " ",
      output_loc, 
      "/", 
      x, 
      "_",
      y, 
      "_bismark_alignment/",
      x,
      "_R2_001.fastq.gz_demultiplexed_",
      y,
      "_trimmed.fastq_unmapped_reads_2.fq.gz; ",
      "rm -r ",
      tempdir
    )
    bismark_logs_aln <- system(command_string_aln, intern = T)
    cat(paste(bismark_logs_aln, collapse = "\n"))
    cat("\n\n\n")
  })
})

# Merge Bam files from the different alignment steps. Do so separately for each
# sample and each demultiplexed cell type. Store the resulting Bam file in a
# dedicated location.
void <- lapply(samples, FUN = function(x) {
  void2 <- lapply(cell_types, FUN = function(y) {
    cat(paste0("Merging Bam files for ", x, "; ", y, "\n"))
    message(paste0("After Bismark aligner: Merging Bam files for ", x, "; ", y))
    merge_output_loc <- paste0(
      output_loc,
      "/", 
      x, 
      "_",  
      y, 
      "_bismark_alignment_original_plus_leftover_realignment"
    )
    if (!dir.exists(merge_output_loc)) {
      dir.create(merge_output_loc)
    }
    merge_command <- paste0(
      "module load samtools/1.11.0; ",
      "samtools merge ",
      "-f ",
      merge_output_loc,
      "/merged_bam_file.bam ",
      output_loc, 
      "/", 
      x, 
      "_",  
      y, 
      "_bismark_alignment/",
      x,
      "_R1_001.fastq.gz_demultiplexed_",
      y,
      "_trimmed_bismark_bt2_pe.bam ",
      output_loc, 
      "/", 
      x, 
      "_",  
      y, 
      "_bismark_alignment_r1_leftover/",
      x,
      "_R1_001.fastq.gz_demultiplexed_",
      y,
      "_trimmed.fastq_unmapped_reads_1_bismark_bt2.bam ",
      output_loc, 
      "/", 
      x, 
      "_",  
      y, 
      "_bismark_alignment_r2_leftover/",
      x,
      "_R2_001.fastq.gz_demultiplexed_",
      y,
      "_trimmed.fastq_unmapped_reads_2_bismark_bt2.bam"
    )
    merge_locs <- system(merge_command, intern = T)
    cat(paste(merge_locs, collapse = "\n"))
    cat("\n\n\n")
  })
})





#############################################################################
# Perform methylation calling.
#############################################################################

# Print a status message.
cat(
  "\n\n\n####################################################################\n"
)
cat("Methylation calling\n")
cat("####################################################################\n\n")


# Call methylation values. Do so separately for every cell type in every 
# sample. Combine results from all possible alignment strategies (top/bottom 
# strand, CT or GA converted) into a single output table. According to
# the documentation, the coordinates in the output table are 1-based.
void <- lapply(samples, FUN = function(x) {
  void2 <- lapply(cell_types, FUN = function(y) {
    cat(paste0("Processing ", x, "; ", y, "\n"))
    message(paste0("Methylation calling: Processing ", x, "; ", y))
    command_string_methcall <- paste0(
      "module load bismark/0.20.0; ",
      "module load bowtie2/2.3.5.1; ",
      "module load samtools/1.11.0; ",
      "bismark_methylation_extractor ",
      "--comprehensive ",
      "--output ",
      output_loc,
      "/",
      x,
      "_",
      y,
      "_bismark_methcall ",
      "--cytosine_report ",
      "--genome_folder ",
      ref_gen_loc,
      " ",
      output_loc,
      "/",
      x,
      "_",
      y,
      "_bismark_alignment_original_plus_leftover_realignment/merged_bam_file.bam"
    )
    bismark_logs_methcall <- system(command_string_methcall, intern = T)
    cat(paste(bismark_logs_methcall, collapse = "\n"))
    cat("\n\n\n")
  })
})







############################################################################
# Process methylation data for amplicon regions.
############################################################################

# Print a status message.
cat(
  "\n\n\n####################################################################\n"
)
cat("Processing methylation calls\n")
cat("####################################################################\n\n")

# Use bedtools to restrict methylation calls to CpGs overlapping with the 
# amplicon locations (requires some reformatting of the methylation call files
# so that there is a start and end position). Take into account that all 
# intervals are already 1-based and closed. Thus, I can keep them as they are.
# Afterwards, sum up data for the two complementary strands of the same CpG (if 
# both are part of the amplicon), compute coverage and methylation beta value 
# and stratify data by amplicon. Note that according to the "Bismark" 
# documentation the first count column contains methylated counts while the 
# second count column contains unmethylated counts.
all_rel_calls <- lapply(samples, FUN = function(x) {
  this_sample_res <- do.call(rbind, lapply(cell_types, FUN = function(y) {
    cat(paste0("Processing ", x, "; ", y, "\n"))
    message(paste0("Methylation call processing: Processing ", x, "; ", y))
    meth_call_file <- paste0(output_loc, "/", x, "_", y, "_bismark_methcall/",
                             "merged_bam_file.CpG_report.txt")
    temp_bed <- paste0(output_loc, "/", x, "_", y, "_bismark_methcall/",
                       "merged_bam_file.CpG_report.txt_temp.bed")
    column_dup_command <- paste0(
      "awk -F '\t' 'BEGIN {OFS = FS} {print $1, $2, $2, $3, $4, $5, $6, $7}' ",
      meth_call_file,
      " > ",
      temp_bed)
    dup_logs <- system(column_dup_command, intern = T)
    cat(paste(dup_logs, collapse = "\n"))
    cat("\n")
    amplicon_df <- as.data.frame(amplicon_locs)[, 1:3]
    amplicon_bed_file <- paste0(output_loc, "/", x, "_", y, 
                                "_bismark_methcall/amplicon_temp.bed")
    write.table(amplicon_df, file = amplicon_bed_file, sep = "\t", 
                col.names = F, row.names = F, quote = F)
    intersect_outfile <- paste0(meth_call_file,
                                "_only_rel_cpgs.txt")
    intersect_command <- paste0(
      "module load bedtools/2.24.0; ",
      "bedtools intersect ",
      "-a ", temp_bed,
      " -b ", amplicon_bed_file,
      " -u",
      " > ", intersect_outfile, "; ",
      "rm ", temp_bed, "; ",
      "rm ", amplicon_bed_file
    )
    intersect_logs <- system(intersect_command, intern = T)
    cat(paste(intersect_logs, collapse = "\n"))
    cat("\n\n\n")
    intersect_res_df <- read.table(intersect_outfile, header = F, 
                                   stringsAsFactors = F)
    intersect_res_df_merged <- data.frame(chr = c(), start = c(), end = c(), 
                                          strand = c(), count_methylated = c(), 
                                          count_nonmethylated = c())
    z <- 1
    while(z <= nrow(intersect_res_df)) {
      if (intersect_res_df[z, "V1"] == intersect_res_df[z + 1, "V1"] &
          intersect_res_df[z, "V2"] == intersect_res_df[z + 1, "V2"] - 1 &
          intersect_res_df[z, "V3"] == intersect_res_df[z + 1, "V3"] - 1 &
          intersect_res_df[z, "V4"] != intersect_res_df[z + 1, "V4"] &
          intersect_res_df[z, "V4"] == "+") {
        new_v1 <- intersect_res_df[z, "V1"]
        new_v2 <- intersect_res_df[z, "V2"]
        new_v3 <- intersect_res_df[z + 1, "V3"]
        new_v4 <- "*"
        new_v5 <- sum(intersect_res_df[c(z, z + 1), "V5"])
        new_v6 <- sum(intersect_res_df[c(z, z + 1), "V6"])
        z <- z + 2 # Skip a line since this line was already incorporated.
      } else {
        new_v1 <- intersect_res_df[z, "V1"]
        new_v2 <- intersect_res_df[z, "V2"]
        new_v3 <- intersect_res_df[z, "V3"]
        new_v4 <- "*"
        new_v5 <- intersect_res_df[z, "V5"]
        new_v6 <- intersect_res_df[z, "V6"]
        z <- z + 1 # Go to the next line.
      }
      new_df <- data.frame(chr = new_v1, start = new_v2, end = new_v3, 
                           strand = new_v4, count_methylated = new_v5, 
                           count_nonmethylated = new_v6)
      intersect_res_df_merged <- rbind(intersect_res_df_merged, new_df)
    }
    intersect_res_df_merged$Coverage <- 
      intersect_res_df_merged$count_methylated + 
      intersect_res_df_merged$count_nonmethylated
    intersect_res_df_merged$Methylation_beta <- 
      intersect_res_df_merged$count_methylated / 
      intersect_res_df_merged$Coverage
    intersect_res_df_merged$Sample <- x
    intersect_res_df_merged$Cell_type <- y
    intersect_res_df_merged_gr <- makeGRangesFromDataFrame(
      intersect_res_df_merged,
      keep.extra.columns = T
    )
    intersect_res_amplicon_overl <- findOverlaps(intersect_res_df_merged_gr,
                                                 amplicon_locs)
    amplicon_strat_res <- do.call(
      rbind, 
      lapply(1:length(amplicon_locs), FUN = function(z) {
        temp_df <- intersect_res_df_merged[
          unique(from(intersect_res_amplicon_overl)[
            to(intersect_res_amplicon_overl) == z
          ])
        , ]
        temp_df$Amplicon <- amplicon_locs$Amplicon_name[z]
        return(temp_df)
      }))
    return(amplicon_strat_res)
  }))
  return(this_sample_res)
})
names(all_rel_calls) <- samples

# Save the data.
processed_data_outfile <- paste0(output_loc, 
                                 "/relevant_methylation_results.rds")
write.csv(all_rel_calls, file = processed_data_outfile)







############################################################################
# Visualise the data.
############################################################################

# Print a status message.
cat(
  "\n\n\n####################################################################\n"
)
cat("Visualising results\n")
cat("####################################################################\n\n")

# Plot methylation data in a pearls-on-a-string fashion. The colour represents
# the methylation beta value. The size of the "pearls" represents the coverage 
# of the corresponding CpG. Save these plots.
void <- lapply(samples, FUN = function(x) {
  rel_data <- all_rel_calls[[x]]
  rel_data$logcov <- log10(rel_data$Coverage + 1)
  rel_data$Cell_type <- factor(rel_data$Cell_type, levels = rev(cell_types))
  plot_list <- lapply(amplicon_names, FUN = function(y) {
    amplicon_data <- rel_data[rel_data$Amplicon == y, ]
    amplicon_plot <- ggplot(amplicon_data) +
      aes(x = start, y = Cell_type, fill = Methylation_beta, size = logcov) +
      scale_x_continuous(
        name = paste0("Position (Chr. ", unique(amplicon_data$chr), ")")) +
      scale_fill_viridis_c(limits = c(0, 1), name = "Methylation", 
                           direction = -1) +
      scale_size_area(limits = c(0, 5), name = "log10(\nCoverage + 1)") +
      geom_hline(yintercept = cell_types, colour = "black") +
      geom_point(shape = 21, colour = "black", stroke = 1.5) +
      coord_cartesian(
        xlim = c(start(amplicon_locs)[amplicon_locs$Amplicon_name == y],
                 end(amplicon_locs)[amplicon_locs$Amplicon_name == y]),
        ylim = c(0, 3),
        expand = F,
        default = T,
        clip = "off") +
      ylab("Cell type") +
      ggtitle(paste0("Amplicon: ", gsub("_", " (", y), ")")) +
      theme_classic() +
      theme(axis.text = element_text(colour = "black"),
            axis.ticks.x = element_line(colour = "black"),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank())
    return(amplicon_plot)
  })
  outfile_snip <- paste0("/amplicons_05_2024_pilot_run_results_", x)
  amplicon_plot_pdf <- paste0(plot_outdir, outfile_snip, ".pdf")
  amplicon_plot_rds <- paste0(plot_rds_outdir, outfile_snip, ".rds")
  pdf(amplicon_plot_pdf, width = 8, height = 4)
  print(plot_list)
  dev.off()
  saveRDS(plot_list, file = amplicon_plot_rds)
})
