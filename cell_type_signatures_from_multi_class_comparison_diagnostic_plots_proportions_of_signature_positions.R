# This script computes the proportions of positions (e.g. first exon, first
#   intron, ...) of the multi-class cell type signatures.
# Author: Niklas Beumer.



# Specify a location on yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Read in the signature regions with their respective genome region annotation.
sig_reg_file <- paste0(
  location,
  "/differential_methylation/multi_class_signatures_regions_with_additional_statistics_and_annotations.txt"
)
sig_reg <- read.table(
  sig_reg_file,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)

# Identify all signature categories and all genomic locations.
sig_cats_all <- sort(unique(sig_reg$automatic_annotation))
gene_body_locs_all <- sort(unique(sig_reg$region_location))
cgi_locs_all <- sort(unique(sig_reg$cpg_island_location_cpgIslandExt))

# For each signature category, compute the percentage of signature regions
# centered in each interval w.r.t. gene bodies. Also perform this computation
# for all signature categories together.
# Save this table.
gene_body_quant_1 <- data.frame(
  Signature_category = "All",
  Interval = gene_body_locs_all,
  Percentage = 100 * sapply(
    gene_body_locs_all,
    FUN = function(x) {
      length(which(sig_reg$region_location == x)) / nrow(sig_reg)
    }
  )
)
gene_body_quant_2 <- do.call(rbind, lapply(
  sig_cats_all,
  FUN = function(x) {
    df_to_use <- sig_reg[sig_reg$automatic_annotation == x, ]
    temp_df <- data.frame(
      Signature_category = x,
      Interval = gene_body_locs_all,
      Percentage = 100 * sapply(
        gene_body_locs_all,
        FUN = function(y) {
          length(which(df_to_use$region_location == y)) / nrow(df_to_use)
        }
      )
    )
    return(temp_df)
  }
))
gene_body_quant_complete <- rbind(gene_body_quant_1, gene_body_quant_2)
gene_body_quant_outfile <- paste0(
  location,
  "/differential_methylation/sig_reg_pie_charts_percentages_gene_body.txt"
)
write.table(
  gene_body_quant_complete,
  file = gene_body_quant_outfile,
  sep = "\t",
  row.names = F,
  col.names = T,
  quote = F
)

# For each signature category, compute the percentage of signature regions
# centered in each interval w.r.t. CpG islands. Also perform this computation
# for all signature categories together.
# Save this table.
cgi_quant_1 <- data.frame(
  Signature_category = "All",
  Interval = cgi_locs_all,
  Percentage = 100 * sapply(
    cgi_locs_all,
    FUN = function(x) {
      length(which(sig_reg$cpg_island_location_cpgIslandExt == x)) / nrow(sig_reg)
    }
  )
)
cgi_quant_2 <- do.call(rbind, lapply(
  sig_cats_all,
  FUN = function(x) {
    df_to_use <- sig_reg[sig_reg$automatic_annotation == x, ]
    temp_df <- data.frame(
      Signature_category = x,
      Interval = cgi_locs_all,
      Percentage = 100 * sapply(
        cgi_locs_all,
        FUN = function(y) {
          length(which(df_to_use$cpg_island_location_cpgIslandExt == y)) /
            nrow(df_to_use)
        }
      )
    )
    return(temp_df)
  }
))
cgi_quant_complete <- rbind(cgi_quant_1, cgi_quant_2)
cgi_quant_outfile <- paste0(
  location,
  "/differential_methylation/sig_reg_pie_charts_percentages_cpg_isl.txt"
)
write.table(
  cgi_quant_complete,
  file = cgi_quant_outfile,
  sep = "\t",
  row.names = F,
  col.names = T,
  quote = F
)
