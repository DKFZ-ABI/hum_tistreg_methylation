# This script delves more into discrepancies between methylation and chromatin
#   accessibility with respect to the enrichment of transcription factors
#   (Comparison between skin Treg cells and blood naive Treg cells).
#   Specifically, this script generates the plts showing accessibility around
#   predicted motif sites.
# Author: Niklas Beumer



# Load required packages.
library(Seurat)
library(Signac)
library(motifmatchr)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(testit)
library(ggplot2)


# Define a location on /yyy.
b330_space <- "/yyy/"
location <- paste0(b330_space, "yyy/hm_treg_bs_rgnsbg")

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/xxx/hm_treg_bs_rgnsbg/analysis",
                     format(Sys.time(), "%m-%d-%y"),
                     sep = "/")
if (!dir.exists(plot_outdir)) {
  dir.create(plot_outdir)
}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the Seurat object containing scATAC-seq data of CD4+ cells.
sc_data_cd4_file <- paste0(
  b330_space,
  "yyy/scATAC/final/seurat_human_CD4_CD25_scATAC_binary_5000_frags_dimred_chromvar.RDS"
)
sc_data_cd4 <- readRDS(sc_data_cd4_file)

# Specify the path to the fragments file.
fragments_path <- paste0(
  b330_space,
  "yyy/10X_scATAC_human/cellranger_aggregated/cellranger_aggr_5000_frags_cd4_cd25/outs/frags_split/fragments_50000_max_fragments.tsv.gz"
)

# Read in the results for known motifs from Homer.
homer_snips <- c(
  "diff_meth_skin_treg_blood_naive_treg_Blood_naive_Treg__hypomethylation",
  "diff_meth_skin_treg_blood_naive_treg_Skin_Treg__hypomethylation",
  "diff_acc_skin_treg_blood_naive_treg_Blood_naive_Treg__hyperaccessibility",
  "diff_acc_skin_treg_blood_naive_treg_Skin_Treg__hyperaccessibility"
)
homer_results <- lapply(
  homer_snips,
  FUN = function(x) {
    homer_file <- paste0(location,
                         "/treg_hierarchies/",
                         x,
                         "_homer_against_genome/knownResults.txt")
    read.table(
      homer_file,
      sep = "\t",
      header = T,
      stringsAsFactors = F,
      comment.char = ""
    )
  }
)
names(homer_results) <- homer_snips






###############################################################################
# Generate chromatin accessibility footprint plots for all interesting
# transcription factor motifs.
###############################################################################

# Update the Seurat object.
sc_data_cd4 <- UpdateSeuratObject(sc_data_cd4)
sc_data_cd4[["scATAC_raw"]] <- as.ChromatinAssay(sc_data_cd4[["scATAC_raw"]])
genome(sc_data_cd4[["scATAC_raw"]]) <- "hg19"
fragments_obj <- CreateFragmentObject(fragments_path, cells = colnames(sc_data_cd4))
Fragments(sc_data_cd4) <- fragments_obj

# Restrict the ATAC-seq data to skin Tregs and blood naive Tregs.
sc_data_cd4 <- subset(sc_data_cd4, subset = treg_tconv_annot %in% c("skin_treg", "blood_naive_treg"))

# Perform a new TF-IDFfor the reduced object.
sc_data_cd4 <- RunTFIDF(sc_data_cd4)

# Specify the motifs that are relevant.
tfs_to_show <- c(
  # Selection of bHLH factors
  "c-Myc(bHLH)/LNCAP-cMyc-ChIP-Seq(Unpublished)/Homer",
  "c-Myc(bHLH)/mES-cMyc-ChIP-Seq(GSE11431)/Homer",
  "USF1(bHLH)/GM12878-Usf1-ChIP-Seq(GSE32465)/Homer",
  "n-Myc(bHLH)/mES-nMyc-ChIP-Seq(GSE11431)/Homer",
  "HIF2a(bHLH)/785_O-HIF2a-ChIP-Seq(GSE34871)/Homer",
  # Selection of bZIP factors.
  "Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer",
  "Fosl2(bZIP)/3T3L1-Fosl2-ChIP-Seq(GSE56872)/Homer",
  "Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer",
  "JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer",
  "Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer",
  "Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer",
  "BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer",
  "Bach2(bZIP)/OCILy7-Bach2-ChIP-Seq(GSE44420)/Homer",
  "AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer",
  "Bach1(bZIP)/K562-Bach1-ChIP-Seq(GSE31477)/Homer"
)

# Read in the motif matrices (position frequency matrices) for the specified motifs.
#-- Iterate over the motifs.
pos_freq_matrs <- do.call(PFMatrixList, lapply(
  tfs_to_show,
  FUN = function(x) {
    #-- Find one homer results directory that contains this motif.
    found <- F
    i <- 0
    while (!found) {
      i <- i + 1
      in_file <- x %in% homer_results[[i]]$Motif.Name
      #-- Identify the index of the motif in that file.
      index <- which(homer_results[[i]]$Motif.Name == x)
      #-- Specify the file containing the corresponding position frequency matrix.
      pfm_file <- paste0(
        location,
        "/treg_hierarchies/",
        homer_snips[i],
        "_homer_against_genome/knownResults/known",
        index,
        ".motif"
      )
      #-- Check whether this file exists.
      if (file.exists(pfm_file)) {
        found <- T
      }
    }
    #-- Read in the motif.
    pfm <- t(read.table(
      pfm_file,
      header = F,
      stringsAsFactors = F,
      comment.char = ">"
    ))
    #-- Add bases as rownames. I previously compared some motif matrices with their consensus
    #-- in homer_results[[i]] and found that bases are indeed in this order.
    rownames(pfm) <- c("A", "C", "G", "T")
    #-- Assert that all entries are rounded to a maximum of 3 decimal digits.
    assert(all(pfm == round(pfm, 3)))
    #-- Generate a PFMatrix object.
    #-- Include a factor of 1,000 because otherwise entries are rounded to 0.
    #-- 1,000 is sufficient here because I previously tested that numbera are present
    #-- with (at most) three decimal digits.
    pfmatr_obj <- PFMatrix(
      ID = x,
      name = x,
      strand = "*",
      profileMatrix = 1000 * pfm
    )
    #-- Return value.
    return(pfmatr_obj)
  }
))
names(pos_freq_matrs) <- tfs_to_show

# Add motifs to the Seurat object.
sc_data_cd4 <- AddMotifs(sc_data_cd4, genome = BSgenome.Hsapiens.UCSC.hg19, pfm = pos_freq_matrs)

# Generate transcription factor footprint information.
sc_data_cd4 <- Footprint(sc_data_cd4, motif.name = tfs_to_show, genome = BSgenome.Hsapiens.UCSC.hg19)

# Plot transcription factor footprints.
footprint_plots <- lapply(
  tfs_to_show,
  FUN = function(x) {
    PlotFootprint(
      sc_data_cd4,
      features = x,
      group.by = "treg_tconv_annot",
      label = F
    ) &
      scale_colour_manual(
        breaks = c("skin_treg", "blood_naive_treg"),
        labels = c("Skin Treg", "Blood naive Treg"),
        values = c("blue", "darkorchid1"),
        name = "Cell type"
      ) &
      ggtitle(gsub("/", "\n", x))
  }
)

# Save these plots.
footprint_plot_pdf <- paste0(
  plot_outdir,
  "/treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_footprint_plots_relevant_tfs.pdf"
)
footprint_plot_rds <- paste0(
  plot_rds_outdir,
  "/treg_hierarchies_diff_meth_diff_acc_skin_vs_naive_tf_homer_against_genome_footprint_plots_relevant_tfs.rds"
)
pdf(footprint_plot_pdf, width = 5, height = 5)
print(footprint_plots)
dev.off()
saveRDS(footprint_plots, file = footprint_plot_rds)