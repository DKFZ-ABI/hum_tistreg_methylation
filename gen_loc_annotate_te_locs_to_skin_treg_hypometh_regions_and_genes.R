# This script annotates insertion sites of TEs that are enriched in skin Treg
#   hypomethylation regions to skin Treg hypomethylation DMRs and genes.
# Author: Niklas Beumer



# Load required package(s).
library(GenomicRanges)
library(rtracklayer)


# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the results of the enrichment analysis of TEs in skin Treg 
# hypomethylation regions.
enrichment_file <- paste0(location, 
                          "/te_analysis/gen_loc_enrichm_dmrs_skin_naive.rds")
enrichment_res <- readRDS(enrichment_file)
enrichment_res_skin_name <- enrichment_res$Skin_Treg__hypomethylation$name

# Read in the DMRs between skin Tregs and blood naive Tregs. Restrict to skin 
# Treg hypomethylation regions and turn them into a GRanges object.
dmr_file <- paste0(
  location, 
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_gap_0.15_sterrratio_0.5_signature_regions_filtered.txt"
)
dmrs <- read.table(dmr_file, header = T, stringsAsFactors = F, sep = "\t")
dmrs_filt <- dmrs[dmrs$automatic_annotation == "Skin_Treg__hypomethylation", ]
dmrs_gr <- makeGRangesFromDataFrame(dmrs_filt, keep.extra.columns = T)

# Read in the locations of repeat elements from RepeatMasker
rmsk_file <- paste0(location, 
                    "/external_data/2023-09-14_rmsk_hg19_from_ucsc.txt.gz")
rmsk <- read.table(rmsk_file, header = F, stringsAsFactors = F)

# Filter the repeat element data so that only those elements remain where I am
# confident that they are transposons.
conf_transp_cats <- c("LINE", "DNA", "SINE", "LTR", "RC", "Retroposon")
rmsk_filt <- rmsk[rmsk$V12 %in% conf_transp_cats, ]

# Restrict the TE locations to those TEs that were significantly enriched in
# skin Treg hypomethylation regions.
enr_tes <- enrichment_res_skin_name$ID[
  enrichment_res_skin_name$Pval_adj_BH < 0.05]
rmsk_filt_2 <- rmsk_filt[rmsk_filt$V11 %in% enr_tes, ]

# Generate a GRanges object for the annotated TE insertion sites.
# Take into account that this table originated from UCSC. Intervals are thus 
# very likely 0-based and closed.
# Modify seqlevels styles so that they match with what I used for the DMRs.
rmsk_filt_gr <- makeGRangesFromDataFrame(rmsk_filt_2,
                                         ignore.strand = T,
                                         seqnames.field = "V6",
                                         start.field = "V7",
                                         end.field = "V8",
                                         starts.in.df.are.0based = T,
                                         keep.extra.columns = T)
seqlevelsStyle(rmsk_filt_gr) <- "NCBI"








############################################################################
# Annotate TE insertion sites to skin Treg hypomethylation regions and to
# genes.
############################################################################

# Find overlaps between relevant TE insertion sites and skin Treg 
# hypomethylation regions.
te_dmr_overlaps <- findOverlaps(rmsk_filt_gr, dmrs_gr, ignore.strand = T)

# Collect data on TE-DMR overlaps.
query_hits <- from(te_dmr_overlaps)
subject_hits <- to(te_dmr_overlaps)
te_seqnames <- seqnames(rmsk_filt_gr)
te_starts <- start(rmsk_filt_gr)
te_ends <- end(rmsk_filt_gr)
te_names <- rmsk_filt_gr$V11
te_classes <- rmsk_filt_gr$V12
dmr_seqnames <- seqnames(dmrs_gr)
dmr_starts <- start(dmrs_gr)
dmr_ends <- end(dmrs_gr)
dmr_ids <- dmrs_gr$region_ID
dmr_genes <- dmrs_gr$gene_assignments
dmr_overl_types <- dmrs_gr$gene_assignment_overlap_types
overlaps_data <- do.call(rbind, mclapply(1:length(te_dmr_overlaps), 
                                         FUN = function(x) {
  corresp_q_hit <- query_hits[x]
  corresp_s_hit <- subject_hits[x]
  temp_df <- data.frame(
    TE__seqnames = te_seqnames[corresp_q_hit],
    TE__start = te_starts[corresp_q_hit],
    TE_end = te_ends[corresp_q_hit],
    TE__name = te_names[corresp_q_hit],
    TE__class = te_classes[corresp_q_hit],
    DMR__seqnames = dmr_seqnames[corresp_s_hit],
    DMR__start = dmr_starts[corresp_s_hit],
    DMR__end = dmr_ends[corresp_s_hit],
    DMR_ID = dmr_ids[corresp_s_hit],
    DMR__gene_assignments = dmr_genes[corresp_s_hit],
    DMR__gene_assignment_overlap_types = dmr_overl_types[corresp_s_hit]
  )
  return(temp_df)
}, mc.cores = 10))

# Order the overlap data alphabetically be annotated genes.
overlaps_data <- overlaps_data[order(overlaps_data$DMR__gene_assignments), ]

# Save the compiled overlaps data.
overlaps_outfile <- paste0(
  location, 
  "/te_analysis/gen_loc_tes_enriched_in_skin_hypometh_insertion_sites_overl_w_skin_hypometh_dmrs.txt"
)
write.table(overlaps_data, file = overlaps_outfile, sep = "\t", row.names = F, 
            col.names = T, quote = F)

# Generate a BigWig file with the locations of TE regions in the table.
overlaps_table_gr <- makeGRangesFromDataFrame(overlaps_data, 
                                              seqnames.field = "TE__seqnames",
                                              start.field = "TE__start",
                                              end.field = "TE_end",
                                              keep.extra.columns = F)
overlaps_table_gr <- disjoin(overlaps_table_gr)
overlaps_table_gr$score <- 1
overlaps_table_gr <- overlaps_table_gr[
  seqnames(overlaps_table_gr) %in% c(as.character(1:22), "X", "MT")
]
seqlevels(overlaps_table_gr) <- unique(as.character(seqnames(overlaps_table_gr)))
chr_lengths <- GenomeInfoDb::seqlengths(
  BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5)[
    c(as.character(1:22), "X", "MT")
  ]
seqlengths(overlaps_table_gr) <- chr_lengths[seqlevels(overlaps_table_gr)]
bigwig_outfile <- paste0(
  location, 
  "/te_analysis/gen_loc_tes_enriched_in_skin_hypometh_insertion_sites_overl_w_skin_hypometh_dmrs_bigwig.bw"
)
export(overlaps_table_gr, con = bigwig_outfile)
