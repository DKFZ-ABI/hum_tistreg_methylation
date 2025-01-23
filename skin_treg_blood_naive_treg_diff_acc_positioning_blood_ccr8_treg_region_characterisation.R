# This script characterises regions where blood CCR8+ Treg cells display a
#   particular positioning (accessibility-wise) with respect to skin Treg cells
#   and blood naive Treg cells.
# Author: Niklas Beumer



# Load required package(s)
library(ggplot2)
library(GenomicRanges)
library(parallel)
library(testit)
library(ggplot2)


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

# Read in the information on differentially accessible regions between
# skin Tregs and blood naive Tregs (with positioning of blood CCR8+ Tregs).
regions_file <- paste0(
  location,
  "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_acc_discrep_regions_info.txt"
)
regions <- read.table(
  regions_file,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)

# Read in the file containing gene annotations for normal chromosomes and 
# extract all considered gene names.
annotation_pref <- "/yyy/hg19_GRCh37_1000genomes/databases/gencode/gencode19/GencodeV19_"
genes <- read.table(
  paste0(annotation_pref, "Genes_plain.bed.gz"),
  header = T,
  stringsAsFactors = F,
  comment.char = ""
)
all_genes <- genes$name

# Read in the ontology gene sets from MSigDB.
ontology_genesets_file <- paste0(
  location,
  "/external_data/2023-01-17_snapshot_msigdb_ontology_gene_sets_c5_human.rds"
)
ontology_genesets <- readRDS(ontology_genesets_file)





###############################################################################
# Extract region sets with interesting positionings of blood CCR8+ Tregs.
###############################################################################

# Blood naive Treg hyperaccessibility regions in which blood CCR8+ Tregs are closer
# to blood naive Tregs.
naive_hyper_ccr8_close_blood_naive <-
  regions[regions$avg_log2FC > 0 &
            regions$Cell_type_closest_to_blood_ccr8_treg == "Blood naive Treg", ]

# Blood naive Treg hyperaccessibility regions in which blood CCR8+ Tregs are closer
# to skin Tregs.
naive_hyper_ccr8_close_skin <-
  regions[regions$avg_log2FC > 0 &
            regions$Cell_type_closest_to_blood_ccr8_treg == "Skin Treg", ]

# Skin Treg hyperaccessibility regions in which blood CCR8+ Tregs are closer
# to blood naive Tregs.
skin_hyper_ccr8_close_blood_naive <-
  regions[regions$avg_log2FC < 0 &
            regions$Cell_type_closest_to_blood_ccr8_treg == "Blood naive Treg", ]

# Skin Treg hyperaccesibility regions in which blood CCR8+ Tregs are closer
# to skin Tregs.
skin_hyper_ccr8_close_skin <-
  regions[regions$avg_log2FC < 0 &
            regions$Cell_type_closest_to_blood_ccr8_treg == "Skin Treg", ]

# Skin Treg hyperaccessibility regions in which blood CCR8+ Tregs are substantially
# hyperaccesible with respect to skin Tregs.
# These are identified as regions where the difference (in scaled accessibility) between
# skin Tregs and blood CCR8+ Tregs is more than one standard deviation below 1.
scaled_acc_diffs_cat <- regions$Scaled_norm_acc_diff_blood_ccr8_treg_skin_treg[regions$avg_log2FC < 0]
scaled_acc_diffs_cat_sd <- sd(scaled_acc_diffs_cat)
skin_hyper_ccr8_hyper_wrt_skin <-
  regions[regions$avg_log2FC < 0 &
            regions$Scaled_norm_acc_diff_blood_ccr8_treg_skin_treg < -scaled_acc_diffs_cat_sd, ]

# Collect these region sets in a list.
region_sets <- list(
  naive_hyper_ccr8_close_blood_naive,
  naive_hyper_ccr8_close_skin,
  skin_hyper_ccr8_close_blood_naive,
  skin_hyper_ccr8_close_skin,
  skin_hyper_ccr8_hyper_wrt_skin
)

# Generate file name snippets for each of these region sets.
file_snips <- c(
  "blood_naive_treg_hyperaccessibility_closer_to_blood_naive",
  "blood_naive_treg_hyperaccessibility_closer_to_skin",
  "skin_treg_hyperaccessibility_closer_to_blood_naive",
  "skin_treg_hyperaccessibility_closer_to_skin",
  "skin_treg_hyperaccessibility_even_larger_in_blood_ccr8_sd_crit"
)

# Generate a list of region sets that constitute complements to the regions in
# the actual region sets. The complement is constituted of the remaining regions
# of the same general accessibility tendency.
complement_sets <- list(
  naive_hyper_ccr8_close_skin,
  naive_hyper_ccr8_close_blood_naive,
  skin_hyper_ccr8_close_skin,
  skin_hyper_ccr8_close_blood_naive,
  regions[regions$avg_log2FC < 0 &
            !(regions$Peak %in% skin_hyper_ccr8_hyper_wrt_skin$Peak), ]
)








#############################################################################
# Analyse enrichment of transcription factors in region sets (compared to
# the whole genome) using Homer
#############################################################################

# Iterate over the region sets.
for (i in 1:length(region_sets)) {
  # Extract the region set to use.
  region_set_to_use <- region_sets[[i]]
  
  # Generate a directory to hold homer outputs and temporary inputs
  # (if it doesn't already exist).
  homer_out_dir <- paste0(
    location,
    "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_",
    file_snips[i],
    "_homer_against_genome"
  )
  if (!(dir.exists(homer_out_dir))) {
    dir.create(homer_out_dir)
  }
  
  # Generate a temporary bed file for the region set.
  bed_prep_df <- do.call(rbind, lapply(
    region_set_to_use$Peak,
    FUN = function(x) {
      x_split <- strsplit(x, split = "-")[[1]]
      return(data.frame(
        Chr = x_split[1],
        Start = x_split[2],
        End = x_split[3]
      ))
    }
  ))
  bed_prep_df$Region_ID <- region_set_to_use$Peak
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




#####################################################################################
# Analyse enrichment of transcription factors in region sets (compared to
# the the remaining regions of the same general accessibility tendency) using Homer
#####################################################################################

# Iterate over the region sets.
for (i in 1:length(region_sets)) {
  # Extract the region set to use.
  region_set_to_use <- region_sets[[i]]
  
  # Extract the complement set to use.
  complement_to_use <- complement_sets[[i]]
  
  # Generate a directory to hold homer outputs and temporary inputs
  # (if it doesn't already exist).
  homer_out_dir <- paste0(
    location,
    "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_",
    file_snips[i],
    "_homer_against_remaining_regions"
  )
  if (!(dir.exists(homer_out_dir))) {
    dir.create(homer_out_dir)
  }
  
  # Generate a temporary bed file for the region set.
  bed_prep_df <- do.call(rbind, lapply(
    region_set_to_use$Peak,
    FUN = function(x) {
      x_split <- strsplit(x, split = "-")[[1]]
      return(data.frame(
        Chr = x_split[1],
        Start = x_split[2],
        End = x_split[3]
      ))
    }
  ))
  bed_prep_df$Region_ID <- region_set_to_use$Peak
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
  
  # Generate a temporary bed file for the complement set.
  bed_prep_df <- do.call(rbind, lapply(
    complement_to_use$Peak,
    FUN = function(x) {
      x_split <- strsplit(x, split = "-")[[1]]
      return(data.frame(
        Chr = x_split[1],
        Start = x_split[2],
        End = x_split[3]
      ))
    }
  ))
  bed_prep_df$Region_ID <- complement_to_use$Peak
  bed_prep_df$Empty_col_5 <- "mock"
  bed_prep_df$Empty_col_5 <- "mock"
  bed_prep_df$Strand = "."
  # bed_prep_df <- bed_prep_df[1:20, ] # For debugging purposes.
  complement_outfile <- paste0(homer_out_dir, "/temp_complement_file.bed")
  write.table(
    bed_prep_df,
    file = complement_outfile,
    sep = "\t",
    row.names = F,
    col.names = F,
    quote = F
  )
  
  # Generate a command string for homer analysis.
  # The complement region set is used as the custom background set.
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
    " -bg ",
    complement_outfile,
    # Path to the custom background file
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
    
    #-- Remove temporary bed files.
    "rm ",
    bed_outfile,
    "; ",
    "rm ",
    complement_outfile
  )
  
  # Run homer analysis.
  homer_logs <- system(command_string, intern = T)
  
  # Save the messages from the homer process.
  homer_logs_outfile <- paste0(homer_out_dir, "/homer_logs.txt")
  write(homer_logs, file = homer_logs_outfile)
  
}




##########################################################################################
# Visualise enrichment of selected bZIP motifs among classes of regions that are
# hyperaccessible in skin Tregs. The motifs plotted here are the same bZIP motifs that
# were plotted in the script treg_hierarchies/skin_treg_blood_naive_treg_diff_meth_
# positioning_blood_ccr8_treg_region_characterisation.R plus some additional ones
# that prominently appeared in the ATAC data.
###########################################################################################

# Read in the results for known motifs from Homer.
homer_results <- lapply(
  file_snips[3:4],
  FUN = function(x) {
    homer_file <- paste0(
      location,
      "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_",
      x,
      "_homer_against_remaining_regions/knownResults.txt"
    )
    read.table(
      homer_file,
      sep = "\t",
      header = T,
      stringsAsFactors = F,
      comment.char = ""
    )
  }
)
names(homer_results) <- c(
  "Skin_Treg__hyperaccessibility\n(blood CCR8+ Tregs are closer to blood naive Tregs)",
  "Skin_Treg__hyperaccessibility\n(blood CCR8+ Tregs are closer to skin Tregs)"
)


# Compute enrichment scores and combine the results for the region classes.
#-- Iterate over homer's results.
homer_results_combined <- do.call(rbind, lapply(
  names(homer_results),
  FUN = function(x) {
    curr_df <- homer_results[[x]]
    #-- Compute enrichment scores.
    target_perc <- as.numeric(gsub("%", "", curr_df$X..of.Target.Sequences.with.Motif))
    background_perc <- as.numeric(gsub("%", "", curr_df$X..of.Background.Sequences.with.Motif))
    curr_df$Enrichment_score <- target_perc / background_perc
    #-- Scale the enrichment scores so that they range from 0 to 1.
    min_score <- min(curr_df$Enrichment_score)
    max_score <- max(curr_df$Enrichment_score)
    curr_df$Enrichment_score_scaled <- (curr_df$Enrichment_score - min_score) / (max_score - min_score)
    #-- Make column names consistent between the different results tables.
    colnames(curr_df)[6] <- "Number..of.Target.Sequences.with.Motif"
    colnames(curr_df)[8] <- "Number..of.Background.Sequences.with.Motif"
    #-- Add information on the region class.
    curr_df$Region_class <- x
    #-- Return value
    return(curr_df)
  }
))

# Save the combined results.
homer_results_combined_outfile <- paste0(
  location,
  "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_skin_treg_hyperaccessibility_homer_against_remaining_regions_combined_results_w_enr_scores.txt"
)
write.table(
  homer_results_combined,
  file = homer_results_combined_outfile,
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

# Specify the transcription factor motifs that should be shown in the plot.
# The motifs to show are the bZIP motifs that were already plotted in the script
# treg_hierarchies/skin_treg_blood_naive_treg_diff_meth_positioning_blood_ccr8_treg_region_characterisation.R
# plus some additional ones that prominently appeared in the ATAC data.
tfs_to_show <- c(
  "HLF(bZIP)/HSC-HLF.Flag-ChIP-Seq(GSE69817)/Homer",
  "Atf1(bZIP)/K562-ATF1-ChIP-Seq(GSE31477)/Homer",
  "BATF(bZIP)/Th17-BATF-ChIP-Seq(GSE39756)/Homer",
  "Atf3(bZIP)/GBM-ATF3-ChIP-Seq(GSE33912)/Homer",
  "NFIL3(bZIP)/HepG2-NFIL3-ChIP-Seq(Encode)/Homer",
  "AP-1(bZIP)/ThioMac-PU.1-ChIP-Seq(GSE21512)/Homer",
  "bZIP:IRF(bZIP,IRF)/Th17-BatF-ChIP-Seq(GSE39756)/Homer",
  "MafA(bZIP)/Islet-MafA-ChIP-Seq(GSE30298)/Homer",
  "Fra1(bZIP)/BT549-Fra1-ChIP-Seq(GSE46166)/Homer",
  "JunB(bZIP)/DendriticCells-Junb-ChIP-Seq(GSE36099)/Homer",
  "Fosl2(bZIP)/3T3L1-Fosl2-ChIP-Seq(GSE56872)/Homer",
  "Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer",
  "Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer",
  "Bach2(bZIP)/OCILy7-Bach2-ChIP-Seq(GSE44420)/Homer"
)

# Collect values for these transcription factors.
relevant_values <- homer_results_combined[homer_results_combined$Motif.Name %in% tfs_to_show, ]

# Clip Benjamini-Hochberg-adjusted P values.
relevant_values$qval_clipped <- sapply(
  relevant_values$q.value..Benjamini.,
  FUN = function(x) {
    max(x, 0.001)
  }
)

# Compute log-transformed q values.
relevant_values$logq <- -log10(relevant_values$qval_clipped) + 0.75

# Add a column indicating statistical significance.
relevant_values$significant <- relevant_values$q.value..Benjamini. < 0.05

# Generate a dot plot showing homer results.
relevant_values$Motif.Name <- factor(relevant_values$Motif.Name, levels = rev(tfs_to_show))
relevant_values$Region_class <- factor(relevant_values$Region_class, levels = names(homer_results))
dot_plot <- ggplot(relevant_values) +
  aes(
    x = Region_class,
    y = Motif.Name,
    size = logq,
    fill = Enrichment_score_scaled,
    colour = significant
  ) +
  scale_size_area(
    breaks = 0.75:3.75,
    labels = c("0", "1", "2", ">=3"),
    name = "-log10(q)",
    limits = c(0.75, 3.75)
  ) +
  scale_fill_viridis_c(
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c("min. f. class", "max. f. class"),
    name = "Enrichment",
    option = "C"
  ) +
  scale_colour_manual(
    breaks = c(T, F),
    values = c("cyan", "grey"),
    labels = c("q < 0.05", "q >= 0.05"),
    name = "Significance"
  ) +
  geom_point(shape = 21, stroke = 1) +
  xlab("Region class") +
  ylab("Motif") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      colour = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

# Save the plot.
dot_plot_pdf <- paste0(
  plot_outdir,
  "/treg_hierarchies_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_skin_treg_hyperaccessibility_homer_against_remaining_regions_dot_plot_selected_motifs.pdf"
)
dot_plot_rds <- paste0(
  plot_rds_outdir,
  "/treg_hierarchies_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_skin_treg_hyperaccessibility_homer_against_remaining_regions_dot_plot_selected_motifs.rds"
)
pdf(dot_plot_pdf, width = 7, height = 7)
print(dot_plot)
dev.off()
saveRDS(dot_plot, file = dot_plot_rds)








###########################################################################################
# Analyse enrichment of gene ontology terms (compared to the genome).
###########################################################################################

# Identify unique gene ontology terms.
ontology_terms <- unique(ontology_genesets$gs_name)


# Iterate over the region sets.
for (i in 1:length(region_sets)) {
  # Extract the set of regions to use.
  region_set_to_use <- region_sets[[i]]
  
  # Identify genes associated with the extracted regions.
  regions_genes <- unique(unlist(
    lapply(
      region_set_to_use$gene_assignments,
      FUN = function(x) {
        strsplit(x, split = ", ")[[1]]
      }
    )
  ))
  regions_genes <- regions_genes[regions_genes != "Not assigned"]
  
  # Identify all other genes in the Gencode annotation.
  all_other_genes <- setdiff(all_genes, regions_genes)
  
  # Analyse enrichment of gene ontology terms in the extracted regions.
  # (Fisher's exact test)
  #-- Iterate over the ontology terms.
  fisher_test_res <- do.call(rbind,
                             mclapply(
                               ontology_terms,
                               FUN = function(x) {
                                 #-- Extract genes for this ontology term.
                                 corresp_genes <- ontology_genesets$gene_symbol[ontology_genesets$gs_name == x]
                                 #-- Identify and quantify overlaps between the genes associated with extracted
                                 #-- regions and genes from this ontology term.
                                 n_regions_overlap <- length(intersect(regions_genes, corresp_genes))
                                 n_regions_not_overlap <- length(setdiff(regions_genes, corresp_genes))
                                 #-- Identify and quantify overlaps between the remaining regions and the
                                 #-- predicted transcription factor binding sites.
                                 n_complement_overlap <- length(intersect(all_other_genes, corresp_genes))
                                 n_complement_not_overlap <- length(setdiff(all_other_genes, corresp_genes))
                                 #-- Perform a one-tailed Fisher's exact test to assess whether the current ontology term is
                                 #-- enriched among the extracted regions.
                                 fisher_input <- matrix(
                                   c(
                                     n_regions_overlap,
                                     n_regions_not_overlap,
                                     n_complement_overlap,
                                     n_complement_not_overlap
                                   ),
                                   ncol = 2,
                                   byrow = T
                                 )
                                 fisher_test_res <- fisher.test(fisher_input, alternative = "greater")
                                 or <- fisher_test_res$estimate
                                 pval <- fisher_test_res$p.value
                                 #-- Collect results.
                                 temp_df <- data.frame(
                                   Ontology_term = x,
                                   n_in_class_overlap_ont_term = n_regions_overlap,
                                   n_in_class_not_overlap_ont_term = n_regions_not_overlap,
                                   n_not_in_class_overlap_ont_term = n_complement_overlap,
                                   n_not_in_class_not_overlap_ont_term = n_complement_not_overlap,
                                   OR = or,
                                   Fisher_pval = pval
                                 )
                                 #-- Return value.
                                 return(temp_df)
                               },
                               mc.cores = 10
                             ))
  
  
  # Correct P values for multiple hypothesis testing (Benjamini-Hochberg correction).
  fisher_test_res$Fisher_pval_adj_BH <- p.adjust(fisher_test_res$Fisher_pval, method = "BH")
  
  # Save the results from Fisher's exact test.
  fisher_test_res_outfile <- paste0(
    location,
    "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_",
    file_snips[i],
    "_go_enr_against_genome.txt"
  )
  
  write.table(
    fisher_test_res,
    file = fisher_test_res_outfile,
    sep = "\t",
    row.names = F,
    col.names = T,
    quote = F
  )
  
  # Generate and save a volcano plot.
  fisher_test_res$logp <- -log10(fisher_test_res$Fisher_pval_adj_BH)
  max_logp <- max(fisher_test_res$logp[fisher_test_res$logp < Inf])
  logp_replacement_val <- 1.5 * max_logp
  fisher_test_res$logp[fisher_test_res$logp == Inf] <- logp_replacement_val
  fisher_test_res$logor <- log10(fisher_test_res$OR)
  max_logor <- max(fisher_test_res$logor[fisher_test_res$logor < Inf])
  min_logor <- min(fisher_test_res$logor[fisher_test_res$logor > -Inf])
  assert(min_logor < 0)
  max_logor_replacement_val <- 1.5 * max_logor
  min_logor_replacement_val <- 1.5 * min_logor
  fisher_test_res$logor[fisher_test_res$logor == Inf] <- max_logor_replacement_val
  fisher_test_res$logor[fisher_test_res$logor == -Inf] <- min_logor_replacement_val
  fisher_test_res$significant <- fisher_test_res$Fisher_pval_adj_BH < 0.05
  step_size_y <- 10 ^ floor(log10(max_logp) - 0.3)
  y_breaks <- seq(0, max_logp, step_size_y)
  step_size_x <- 10 ^ floor(abs(log10(max(
    c(max_logor, min_logor)
  )) - 0.3))
  x_breaks <- c(rev(-seq(0, -min_logor, step_size_x)), seq(step_size_x, max_logor, step_size_x))
  volcano_plot <- ggplot(fisher_test_res) +
    aes(x = logor, y = logp, colour = significant) +
    scale_x_continuous(
      breaks = c(
        min_logor_replacement_val,
        x_breaks,
        max_logor_replacement_val
      ),
      labels = c("-Inf", x_breaks, "Inf"),
      name = "log10(OR)"
    ) +
    scale_y_continuous(
      breaks = c(y_breaks, logp_replacement_val),
      labels = c(y_breaks, "Inf"),
      expand = expansion(mult = c(0.01, 0.05)),
      name = "-log10(q)"
    ) +
    scale_colour_manual(breaks = c("TRUE", "FALSE"),
                        values = c("red", "grey50")) +
    geom_point(size = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic() +
    theme(
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      legend.position = "none"
    )
  volcano_plot_pdf <- paste0(
    plot_outdir,
    "/treg_hierarchies_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_",
    file_snips[i],
    "_go_enr_against_genome.pdf"
  )
  volcano_plot_rds <- paste0(
    plot_rds_outdir,
    "/treg_hierarchies_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_",
    file_snips[i],
    "_go_enr_against_genome.rds"
  )
  pdf(volcano_plot_pdf, width = 3, height = 3)
  print(volcano_plot)
  dev.off()
  saveRDS(volcano_plot, file = volcano_plot_rds)
  
}







###########################################################################################
# Analyse enrichment of gene ontology terms (compared to the remaining regions of the
# same general accessbility tendency).
###########################################################################################


# Iterate over the region sets.
for (i in 1:length(region_sets)) {
  # Extract the set of regions to use.
  region_set_to_use <- region_sets[[i]]
  
  # Extract the complement regions to use.
  complement_to_use <- complement_sets[[i]]
  
  # Identify genes associated with the extracted regions and the complement regions.
  regions_genes <- unique(unlist(
    lapply(
      region_set_to_use$gene_assignments,
      FUN = function(x) {
        strsplit(x, split = ", ")[[1]]
      }
    )
  ))
  regions_genes <- regions_genes[regions_genes != "Not assigned"]
  complement_genes <- unique(unlist(
    lapply(
      complement_to_use$gene_assignments,
      FUN = function(x) {
        strsplit(x, split = ", ")[[1]]
      }
    )
  ))
  complement_genes <- complement_genes[complement_genes != "Not assigned"]
  
  # Analyse enrichment of gene ontology terms in the extracted regions.
  # (Fisher's exact test)
  #-- Iterate over the ontology terms.
  fisher_test_res <- do.call(rbind,
                             mclapply(
                               ontology_terms,
                               FUN = function(x) {
                                 #-- Extract genes for this ontology term.
                                 corresp_genes <- ontology_genesets$gene_symbol[ontology_genesets$gs_name == x]
                                 #-- Extract genes
                                 #-- Identify and quantify overlaps between the genes associated with extracted
                                 #-- regions and genes from this ontology term.
                                 n_regions_overlap <- length(intersect(regions_genes, corresp_genes))
                                 n_regions_not_overlap <- length(setdiff(regions_genes, corresp_genes))
                                 #-- Identify and quantify overlaps between the remaining regions and the
                                 #-- predicted transcription factor binding sites.
                                 n_complement_overlap <- length(intersect(complement_genes, corresp_genes))
                                 n_complement_not_overlap <- length(setdiff(complement_genes, corresp_genes))
                                 #-- Perform a one-tailed Fisher's exact test to assess whether the current ontology term is
                                 #-- enriched among the extracted regions.
                                 fisher_input <- matrix(
                                   c(
                                     n_regions_overlap,
                                     n_regions_not_overlap,
                                     n_complement_overlap,
                                     n_complement_not_overlap
                                   ),
                                   ncol = 2,
                                   byrow = T
                                 )
                                 fisher_test_res <- fisher.test(fisher_input, alternative = "greater")
                                 or <- fisher_test_res$estimate
                                 pval <- fisher_test_res$p.value
                                 #-- Collect results.
                                 temp_df <- data.frame(
                                   Ontology_term = x,
                                   n_in_class_overlap_ont_term = n_regions_overlap,
                                   n_in_class_not_overlap_ont_term = n_regions_not_overlap,
                                   n_not_in_class_overlap_ont_term = n_complement_overlap,
                                   n_not_in_class_not_overlap_ont_term = n_complement_not_overlap,
                                   OR = or,
                                   Fisher_pval = pval
                                 )
                                 #-- Return value.
                                 return(temp_df)
                               },
                               mc.cores = 10
                             ))
  
  
  # Correct P values for multiple hypothesis testing (Benjamini-Hochberg correction).
  fisher_test_res$Fisher_pval_adj_BH <- p.adjust(fisher_test_res$Fisher_pval, method = "BH")
  
  # Save the results from Fisher's exact test.
  fisher_test_res_outfile <- paste0(
    location,
    "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_",
    file_snips[i],
    "_go_enr_against_remaining_regions.txt"
  )
  
  write.table(
    fisher_test_res,
    file = fisher_test_res_outfile,
    sep = "\t",
    row.names = F,
    col.names = T,
    quote = F
  )
  
  # Generate and save a volcano plot.
  fisher_test_res$logp <- -log10(fisher_test_res$Fisher_pval_adj_BH)
  max_logp <- max(fisher_test_res$logp[fisher_test_res$logp < Inf])
  logp_replacement_val <- 1.5 * max_logp
  fisher_test_res$logp[fisher_test_res$logp == Inf] <- logp_replacement_val
  fisher_test_res$logor <- log10(fisher_test_res$OR)
  max_logor <- max(fisher_test_res$logor[fisher_test_res$logor < Inf])
  min_logor <- min(fisher_test_res$logor[fisher_test_res$logor > -Inf])
  assert(min_logor < 0)
  max_logor_replacement_val <- 1.5 * max_logor
  min_logor_replacement_val <- 1.5 * min_logor
  fisher_test_res$logor[fisher_test_res$logor == Inf] <- max_logor_replacement_val
  fisher_test_res$logor[fisher_test_res$logor == -Inf] <- min_logor_replacement_val
  fisher_test_res$significant <- fisher_test_res$Fisher_pval_adj_BH < 0.05
  step_size_y <- 10 ^ floor(log10(max_logp) - 0.3)
  y_breaks <- seq(0, max_logp, step_size_y)
  step_size_x <- 5 * 10 ^ floor(log10(max(c(
    max_logor, min_logor
  ))) - 0.3)
  x_breaks <- c(rev(-seq(0, -min_logor, step_size_x)), seq(step_size_x, max_logor, step_size_x))
  volcano_plot <- ggplot(fisher_test_res) +
    aes(x = logor, y = logp, colour = significant) +
    scale_x_continuous(
      breaks = c(
        min_logor_replacement_val,
        x_breaks,
        max_logor_replacement_val
      ),
      labels = c("-Inf", x_breaks, "Inf"),
      name = "log10(OR)"
    ) +
    scale_y_continuous(
      breaks = c(y_breaks, logp_replacement_val),
      labels = c(y_breaks, "Inf"),
      expand = expansion(mult = c(0.01, 0.05)),
      name = "-log10(q)"
    ) +
    scale_colour_manual(breaks = c("TRUE", "FALSE"),
                        values = c("red", "grey50")) +
    geom_point(size = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic() +
    theme(
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      legend.position = "none"
    )
  volcano_plot_pdf <- paste0(
    plot_outdir,
    "/treg_hierarchies_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_",
    file_snips[i],
    "_go_enr_against_remaining_regions.pdf"
  )
  volcano_plot_rds <- paste0(
    plot_rds_outdir,
    "/treg_hierarchies_diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_",
    file_snips[i],
    "_go_enr_against_remaining_regions.rds"
  )
  pdf(volcano_plot_pdf, width = 3, height = 3)
  print(volcano_plot)
  dev.off()
  saveRDS(volcano_plot, file = volcano_plot_rds)
  
}
