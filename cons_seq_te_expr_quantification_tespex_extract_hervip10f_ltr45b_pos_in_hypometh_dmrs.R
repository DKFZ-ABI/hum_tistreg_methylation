# This script extracts those HERVIP10F-int and LTR45B insertion sites that
#   overlap with Skin Treg hypomethylation DMRs (only insertion sites
#   yielding non-missing methylation levels in my data).
# Author: Niklas Beumer



# Load required packages.
library(GenomicRanges)
library(bsseq)


# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the BSseq object containing the methylation data. Restrict the data
# to the relevant cell types.
meth_data <- readRDS(paste0(
  location, 
  "/preprocessing_etc/quality_control_after_alignment/2022-03-14_bsseq_object_combined_all_samples_smoothed.rds"
))
meth_data <- meth_data[
  , pData(meth_data)$Cell_type %in% c("Skin Treg", "Blood naive Treg")
]


# Read in the regions that are differentially methylated between skin Tregs 
# and blood naive Tregs. Restrict to DMRs of the skin Treg hypomethylation
# class.
dmr_table_file <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_gap_0.15_sterrratio_0.5_signature_regions_filtered.txt"
)
dmr_table <- read.table(dmr_table_file, header = T, stringsAsFactors = F, 
                        sep = "\t")
dmr_table_rel <- dmr_table[
  dmr_table$automatic_annotation == "Skin_Treg__hypomethylation"
, ]

# Turn the DMR locations into a GRanges object.
dmrs_rel_gr <- makeGRangesFromDataFrame(dmr_table_rel, keep.extra.columns = T)

# Read in the locations of repeat elements from RepeatMasker
rmsk_file <- paste0(location, 
                    "/external_data/2023-09-14_rmsk_hg19_from_ucsc.txt.gz")
rmsk <- read.table(rmsk_file, header = F, stringsAsFactors = F)

# Filter the repeat element data so that only those elements remain where I am
# confident that they are transposons.
conf_transp_cats <- c("LINE", "DNA", "SINE", "LTR", "RC", "Retroposon")
rmsk_filt <- rmsk[rmsk$V12 %in% conf_transp_cats, ]

# Generate a GRanges object for the annotated TE insertion sites.
# Take into account that this table originated from UCSC. Intervals are thus 
# very likely 0-based and closed.
# Modify seqlevels styles so that they match with what I used for the signature
# regions.
rmsk_filt_gr <- makeGRangesFromDataFrame(rmsk_filt,
                                         ignore.strand = T,
                                         seqnames.field = "V6",
                                         start.field = "V7",
                                         end.field = "V8",
                                         starts.in.df.are.0based = T,
                                         keep.extra.columns = T)
seqlevelsStyle(rmsk_filt_gr) <- "NCBI"

# Restrict TE data to those TEs that are relevant for my analysis.
relevant_tes <- c("HERVIP10F-int", "LTR45B")
relevant_te_locs <- lapply(relevant_tes, FUN = function(x) {
  rmsk_filt_gr[rmsk_filt_gr$V11 == x]
})
names(relevant_te_locs) <- relevant_tes








################################################################################
# Find those insertion sites of the relevant TEs that overlap with a 
# skin Treg hypomethylation DMR.
################################################################################

# For each of the relevant TE subfamilis, find those that overlap with a 
# skin Treg hypomethylation DMR. For these TEs, collect interesting information.
overl_TE_info <- lapply(relevant_tes, FUN = function(x) {
  te_dmr_overl <- findOverlaps(relevant_te_locs[[x]], dmrs_rel_gr,
                               ignore.strand = T)
  te_dmr_overl_from <- from(te_dmr_overl)
  te_dmr_overl_to <- to(te_dmr_overl)
  te_seqnames <- as.character(seqnames(relevant_te_locs[[x]]))
  te_starts <- start(relevant_te_locs[[x]])
  te_ends <- end(relevant_te_locs[[x]])
  dmr_ID <- dmrs_rel_gr$region_ID
  dmr_gene_assignments <- dmrs_rel_gr$gene_assignments  
  info_df <- do.call(
    rbind, 
    lapply(unique(te_dmr_overl_from), FUN = function(y) {
      te_info <- data.frame(TE_chr = te_seqnames[y],
                            TE_start = te_starts[y],
                            TE_end = te_ends[y])
      corresp_to_vals <- te_dmr_overl_to[te_dmr_overl_from == y]
      corresp_dmr_ids <- paste(sapply(corresp_to_vals, FUN = function(z) {
        dmr_ID[z]
      }), collapse = ", ")
      corresp_gene_assignments <- unique(unlist(
        lapply(corresp_to_vals, FUN = function(z) {
          strsplit(dmr_gene_assignments[z], split = ", ")[[1]]
        })
      ))
      if (!(all(corresp_gene_assignments == "Not assigned"))) {
        corresp_gene_assignments <- corresp_gene_assignments[
          !(corresp_gene_assignments == "Not assigned")
        ]
      }
      
      te_info$DMR_IDs <- corresp_dmr_ids
      te_info$DMR_Gene_assignments <- paste(corresp_gene_assignments,
                                            collapse = ", ")
      return(te_info)
  }))
})
names(overl_TE_info) <- relevant_tes

# Compute methylation differences between blood naive Tregs and skin Tregs
# with respect to the relevant insertion sites. Order the insertion sites by
# the corresponding differences. Only keep those insertion sites that have
# non-missing methylation values.
overl_TE_info_w_methdiff <- lapply(overl_TE_info, FUN = function(x) {
  these_sites_gr <- makeGRangesFromDataFrame(x, 
                                             seqnames.field = "TE_chr",
                                             start.field = "TE_start",
                                             end.field = "TE_end",
                                             keep.extra.columns = T)
  these_sites_meth_vals <- getMeth(meth_data, regions = these_sites_gr, 
                                   type = "raw", what = "perRegion")
  these_sites_celltype_vals <- sapply(
    c("Skin Treg", "Blood naive Treg"), 
    FUN = function(y) {
      celltype_subs <- these_sites_meth_vals[
        , grep(y, colnames(these_sites_meth_vals), fixed = T)
      ]
    return(rowMeans(celltype_subs, na.rm = T))
  })
  colnames(these_sites_celltype_vals) <- paste0(
    "TE_avg_raw_meth_", 
    c("Skin_Treg", "Blood_naive_Treg")
  )
  diffs <- -rowDiffs(these_sites_celltype_vals)
  sites_to_keep <- which(apply(these_sites_celltype_vals, 1, FUN = function(y) {
    all(!is.na(y))
  }))
  updated_df <- cbind(x[sites_to_keep, ],
                      these_sites_celltype_vals[sites_to_keep, ])
  updated_df$TE_meth_diff_naive_vs_skin <- diffs[sites_to_keep]
  updated_df <- updated_df[order(updated_df$TE_meth_diff_naive_vs_skin), ]
  return(updated_df)
})

# Save the tables with relevant TE insertion sites.
outfile_snip_1 <- paste0(location, "/te_analysis/cons_seq_insertion_sites_")
outfile_snip_2 <- "_overl_w_skin_treg_hypometh_dmrs.txt"
void <- lapply(relevant_tes, FUN = function(x) {
  outfile <- paste0(outfile_snip_1, x, outfile_snip_2)
  write.table(overl_TE_info_w_methdiff[[x]], file = outfile, col.names = T, 
              row.names = F, sep = "\t", quote = F)
})
