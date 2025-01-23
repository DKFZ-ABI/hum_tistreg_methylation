# This script computes TPM values from the RPKM values for the RNA-Seq data.
# Author: Niklas Beumer



# Load required packages.
library(testit)

# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Read in the sample mapping file.
sample_mapping_path <- paste0(location, "/sample_mapping_rnaseq.txt")
sample_mapping <- read.table(sample_mapping_path, header = T, 
                             stringsAsFactors = F, sep = "\t")

# Read in the files containing the FPKM values.
rpkm_paths <- unique(sample_mapping$RPKM_path)
rpkm_tables <- lapply(rpkm_paths, FUN = function(x) {
  read.table(x, header = T, stringsAsFactors = F)
})

# Specify the path to the directory where TPM values will be saved.
# Generate this directory if it does not already exist.
tpm_path <- paste0(location, "/RNASeq/tpm")
if (!dir.exists(tpm_path)) {dir.create(tpm_path)}





########################################################
# Collect RPKM values for the relevant samples.
########################################################

# Get the column names representing the Ensembl IDs in each RPKM table.
ensembl_colnames <- sapply(rpkm_tables, FUN = function(x) {
  cols <- colnames(x)
  name <- cols[grep("ensembl", cols, ignore.case = T)]
  return(name)
})

# Get the column names representing the gene symbols in each RPKM table.
gene_symb_colnames <- sapply(rpkm_tables, FUN = function(x) {
  cols <- colnames(x)
  name <- cols[grep("symbol", cols, ignore.case = T)]
  return(name)
})

# Check whether the RPKM tables contain the same Ensembl-IDs in the same order
# and collect the Ensembl-IDs in the data
for (i in 1:(length(rpkm_tables) - 1)) {
  genes_1 <- rpkm_tables[[i]][, ensembl_colnames[i]]
  genes_2 <- rpkm_tables[[i + 1]][, ensembl_colnames[i + 1]]
  assert(all(genes_1 == genes_2))
}
ensembls <- rpkm_tables[[1]][, ensembl_colnames[1]]

# Collect the gene symbols in the data. If there are differences, this is likely 
# because gene symbols were converted to dates by Excel in a part of the data. 
# In this case, use the gene symbol that does not contain a "-", indicating it 
# is not a date.
gene_symbs <- sapply(1:length(rpkm_tables), FUN = function(x) {
  rpkm_tables[[x]][, gene_symb_colnames[x]]
})
gene_symbs_refined <- apply(gene_symbs, 1, FUN = function(x) {
  options <- unique(x)
  if (length(options) == 1) {
    return(options)
  } else {
    non_dates <- options[grep("-", options, invert = T)]
    assert(length(non_dates) == 1)
    return(non_dates)
  }
})

# Concatenate all data into a single data frame.
rpkm_wo_ids <- lapply(1:length(rpkm_tables), FUN = function(x) {
  cols <- colnames(rpkm_tables[[x]])
  cols_to_remove <- which(cols %in% c(ensembl_colnames[x], gene_symb_colnames[x]))
  return(rpkm_tables[[x]][, -cols_to_remove])
})
rpkm_combined <- data.frame(Ensembl_ID = ensembls, Gene_symbol = gene_symbs_refined)
rpkm_combined <- do.call(cbind, c(list(rpkm_combined), rpkm_wo_ids))

# Restrict the RPKM data to the samples that are present in the sample mapping file.
samples_to_consider <- sample_mapping$Sample
match_strings <- sapply(samples_to_consider, FUN = function(x) {
  sub("_R1$", "", sub("p024_RNAseq_", "", x), perl = T)
})
all_cols <- colnames(rpkm_combined)
colnames_to_keep <- unlist(sapply(match_strings, FUN = function(x) {
  all_cols[grep(x, all_cols)]
}))
rpkm_combined_relevant <- rpkm_combined[, c("Ensembl_ID", "Gene_symbol", colnames_to_keep)]
colnames(rpkm_combined_relevant) <- c("Ensembl_ID", "Gene_symbol", samples_to_consider)





################################################################################
# Compute TPM values. The TPM value for a gene can be computed as 10^6 times the 
# RPKM value for that gene divided by the sum of all RPKM values across all 
# genes (see Lena's slides from the sequencing analysis lecture for MoBi 
# Bachelor 5th semester).
################################################################################

# Iterate over the samples and compute TPM values.
tpm_combined_relevant <- rpkm_combined_relevant
for (sample in samples_to_consider) {
  rpkm_sum <- sum(rpkm_combined_relevant[, sample])
  tpm_combined_relevant[, sample] <- sapply(
    rpkm_combined_relevant[, sample], FUN = function(x) {
      10^6 * x / rpkm_sum 
    }
  )
}

# Save the TPM values.
tpm_outfile <- paste0(tpm_path, "/tpm_all.txt")
write.table(tpm_combined_relevant, file = tpm_outfile, col.names = T, 
            row.names = F, sep = "\t", quote = F)
