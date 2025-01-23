# This script examines positionings of blood CCR8+ Treg cells and performs
#   Fisher's exact tests between different omics levels.
# Author: Niklas Beumer



# Define a location on /xxx.
location <- "/xxx/nbeumer/hm_treg_bs_rgnsbg"


# Read in the information regarding where blood CCR8+ Tregs are positioned
# in the different omics levels.
ccr8_pos_meth_file <- paste0(
  location,
  "/treg_hierarchies/diff_meth_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_meth_discrep_regions_info.txt"
)
ccr8_pos_meth <- read.table(
  ccr8_pos_meth_file,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)
ccr8_pos_acc_file <- paste0(
  location,
  "/treg_hierarchies/diff_acc_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_acc_discrep_regions_info.txt"
)
ccr8_pos_acc <- read.table(
  ccr8_pos_acc_file,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)
ccr8_pos_gex_file <- paste0(
  location,
  "/treg_hierarchies/diff_expr_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_filtered_f_expr_discrep_genes_info.txt"
)
ccr8_pos_gex <- read.table(
  ccr8_pos_gex_file,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)








################################################################################
# Perform Fisher's exact test to assess whether the poportion of features for
# which blood CCR8+ Tregs are closer to skin Tregs than to blood naive Tregs
# differs between the different omics levels.
################################################################################

# Write the results to a summary file.
summary_file <- paste0(
  location,
  "/treg_hierarchies/diff_meth_acc_expr_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_fisher_test_summary.txt"
)
sink(summary_file)

# Generate a descriptive introduction to the summary file.
cat("Analysis of the proportions of features (DMRs/peaks/genes) where blood CCR8+ Tregs\n")
cat("are closer to skin Tregs or closer to blood naive Tregs. This analysis tests whether the\n")
cat("aforementioned proportions differ between the different omics levels.\n\n\n\n\n\n")

# Prepare a vector that will later contain all P values.
pvals_vect <- c()

# Iterate over the two directions (stronger in blood naive Tregs and stronger in skin Tregs).
directions <- c("Blood_naive_Treg", "Skin_Treg")
for (direction in directions) {
  # Print a section header to the summary file.
  cat(
    "##############################################################################################\n"
  )
  cat(
    "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
  )
  cat(
    paste0(
      "Hypomethylation/hyperaccessibility/hyperexpression in ",
      direction,
      "\n"
    )
  )
  cat(
    "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
  )
  cat(
    "##############################################################################################\n"
  )
  cat("\n\n\n")
  
  # Iterate over all pairs of omics levels.
  all_levels <- c("meth", "acc", "gex")
  for (i in 1:(length(all_levels) - 1)) {
    for (j in (i + 1):length(all_levels)) {
      levels <- c(all_levels[i], all_levels[j])
      
      # Print a section header to the summary file.
      cat(
        "******************************************************************************\n"
      )
      cat(paste0("Comparison: ", levels[1], " vs. ", levels[2], "\n"))
      cat(
        "******************************************************************************\n"
      )
      cat("\n")
      
      # Find the relevant tables.
      rel_tab <- lapply(
        levels,
        FUN = function(x) {
          eval(parse(text = paste0("ccr8_pos_", x)))
        }
      )
      
      # Restrict to those features that correspond to the current direction.
      fc_direction <- ifelse(direction == "Blood_naive_Treg",
                             yes = ">",
                             no = "<")
      rel_tab_restr <- lapply(
        1:length(rel_tab),
        FUN = function(x) {
          tab_temp <- rel_tab[[x]]
          if (levels[x] == "meth") {
            tab_to_use <- tab_temp[tab_temp$automatic_annotation == paste0(direction, "__hypomethylation"), ]
          } else {
            tab_to_use <- tab_temp[eval(parse(text = paste0(
              "tab_temp$",
              ifelse(levels[x] == "acc", yes = "avg_log2FC", no = "log2FoldChange"),
              " ",
              fc_direction,
              " 0"
            ))), ]
          }
          return(tab_to_use)
        }
      )
      
      # Generate the contingency table for Fisher's exact test.
      contingency_tab <- sapply(
        rel_tab_restr,
        FUN = function(x) {
          return(table(x$Cell_type_closest_to_blood_ccr8_treg))
        }
      )
      colnames(contingency_tab) <- paste0("Closer to ", levels)
      
      # Make sure that skin Tregs appear first.
      contingency_tab <- contingency_tab[c("Skin Treg", "Blood naive Treg"), ]
      
      # Perform a two-tailed Fisher's exact test.
      fisher_test_res <- fisher.test(contingency_tab, alternative = "two.sided")
      
      # Save the results from Fisher's exact test together with the contingency table used.
      list_to_save <- list(contingency_tab = contingency_tab,
                           fisher_test_res = fisher_test_res)
      list_to_save_outfile <- paste0(
        location,
        "/treg_hierarchies/diff_meth_acc_expr_skin_treg_blood_naive_treg_positioning_blood_ccr8_treg_fisher_test_stronger_in_",
        direction,
        "_",
        levels[1],
        "_vs_",
        levels[2],
        ".rds"
      )
      saveRDS(list_to_save, file = list_to_save_outfile)
      
      # Print the contingency table, the test result and the exact P value to the summary file.
      cat("Contingency table:\n")
      print(contingency_tab)
      cat("\n\n")
      cat("Results from Fisher's exact test:\n")
      print(fisher_test_res)
      cat("\n")
      cat("Exact P value:\n")
      cat(fisher_test_res$p.value)
      
      # Print a spacer to the summary file.
      cat("\n\n\n\n\n")
      
      # Collect the P value from this comparison in the P values vector.
      pvals_vect <- c(pvals_vect, fisher_test_res$p.value)
      names(pvals_vect)[length(pvals_vect)] <- paste0(direction, "_hyp__", levels[1], "_vs_", levels[2])
      
    }
  }
  
  # Print a spacer to the summary file.
  cat("\n\n\n\n\n\n")
  
}

# Correct P values for multiple hypothesis testing using the Benjamini-Hochberg
# method and print the results to the output file.
cat(
  "******************************************************************************\n"
)
cat("Benjamini-Hochberg-adjusted P values:\n\n")
cat(
  "******************************************************************************\n"
)
cat("\n")
pvals_adj <- p.adjust(pvals_vect, method = "BH")
print(pvals_adj)

# Stop keeping track of all output.
sink()
