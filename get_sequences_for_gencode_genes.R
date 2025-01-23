# This script retrieves the full sequences (hs37d5) of gene bodies annotated in 
#   Gencode V19 and saves them in Fasta format.
# Author: Niklas Beumer



# Load required packages.
library(BSgenome)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(testit)


# Define a location on /yyy.
b330_space <- "/yyy/"
location <- paste0(b330_space, "yyy/hm_treg_bs_rgnsbg")

# Specify the prefix of the files containing genomic positions of exons, 
# transcription start sites, genes etc.
annotation_pref <- 
  "/yyy/hg19_GRCh37_1000genomes/databases/gencode/gencode19/GencodeV19_"

# Read in the files containing gene and exon annotation for normal chromosomes.
genes <- read.table(paste0(annotation_pref, "Genes_plain.bed.gz"), header = T, 
                    stringsAsFactors = F, comment.char = "")

# Increase start and end positions by 1 nucleotide in order to convert intervals 
# to 1-based closed intervals.
genes$chromStart <- genes$chromStart + 1

# Check that all gene names are unique. This important for running the TESpex
# tool later.
assert(length(unique(genes$name)) == nrow(genes))

# Get the sequences of full gene bodies.
genome_obj <- BSgenome.Hsapiens.1000genomes.hs37d5
seqs <- getSeq(genome_obj,
               names = genes$X.chrom,
               start = genes$chromStart,
               end = genes$chromEnd,
               as.character = T)

# Save the sequences in Fasta format.
all_gene_names <- genes$name
fasta_strings <- sapply(1:(2*length(seqs)), FUN = function(x) {
  if (x %% 2 != 0) {
    str_to_add <- paste0(">", all_gene_names[x / 2 + 0.5])
  } else {
    str_to_add <- seqs[x / 2]
  }
  return(str_to_add)
})
fasta_outfile <- paste0(
  location, 
  "/te_analysis/sequences_full_gene_bodies_gencode_v19_hs37d5_plus_strand.fa"
)
writeLines(fasta_strings, con = fasta_outfile)


