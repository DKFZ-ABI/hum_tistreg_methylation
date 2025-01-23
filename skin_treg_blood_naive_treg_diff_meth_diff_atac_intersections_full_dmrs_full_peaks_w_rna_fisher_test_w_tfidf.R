# This script computes proportions and performs Fisher's exact tests for
#   blood CCR8+ Treg cell positionings regarding DMR-peak-gene links
#   (Comparison: Skin Treg cells vs. blood naive Treg cells).
# Author: Niklas Beumer


# Load required package(s).
library(testit)


# Define a location on /yyy.
location <- "/yyy/hm_treg_bs_rgnsbg"

# Read in the information regarding blood CCR8+ Treg positionings w.r.t.
# DMR-peak-gene links.
ccr8_pos_file <- paste0(
  location,
  "/treg_hierarchies/diff_meth_diff_acc_diff_expr_dmr_peak_gene_links_w_ccr8_positioning_w_tfidf.txt"
)
ccr8_pos <- read.table(
  ccr8_pos_file,
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)



#####################################################################################
# Perform Fisher's exact test to assess whether the poportion of features for which
# blood CCR8+ Tregs are closer to skin Tregs than to blood naive Tregs differs
# between the different omics levels.
#####################################################################################

# Write the results to a summary file.
summary_file <- paste0(
  location,
  "/treg_hierarchies/diff_meth_diff_acc_diff_expr_dmr_peak_gene_links_w_ccr8_positioning_fisher_test_summary_w_tfidf.txt"
)
sink(summary_file)

# Generate a descriptive introduction to the summary file.
cat("Analysis of the proportions of features (DMRs/peaks/genes) where blood CCR8+ Tregs\n")
cat("are closer to skin Tregs or closer to blood naive Tregs. This analysis tests whether the\n")
cat("aforementioned proportions differ between the different omics levels.\n")
cat("The analysis is restricted to DMR-peak-gene links.\n\n\n\n\n\n")

# Iterate over the two methylation tendencies (hypomethylation in blood naive
# Tregs and hypomethylation in skin Tregs).
directions <- c("Blood_naive_Treg", "Skin_Treg")
for (direction in directions) {
  # Print a section header to the summary file.
  cat(
    "##############################################################################################\n"
  )
  cat(
    "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
  )
  cat(paste0("Hypomethylation in ", direction, "\n"))
  cat(
    "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
  )
  cat(
    "##############################################################################################\n"
  )
  cat("\n\n\n")
  
  # Restrict to those features that correspond to the current methylation
  # tendency.
  rel_tab_restr <- ccr8_pos[ccr8_pos$Methylation__automatic_annotation ==
                              paste0(direction, "__hypomethylation"), ]
  
  # Iterate over all pairs of omics levels.
  all_levels <- c("meth", "acc", "gex")
  for (i in 1:(length(all_levels) - 1)) {
    for (j in (i + 1):length(all_levels)) {
      levels <- c(all_levels[i], all_levels[j])
      levels_long <- sapply(
        levels,
        FUN = function(x) {
          switch(x, meth = "Methylation", acc = "Accessibility", "Expression")
        }
      )
      
      # Print a section header to the summary file.
      cat(
        "******************************************************************************\n"
      )
      cat(paste0("Comparison: ", levels[1], " vs. ", levels[2], "\n"))
      cat(
        "******************************************************************************\n"
      )
      cat("\n")
      
      # Generate the contingency table for Fisher's exact test.
      contingency_tab <- matrix(c(
        length(which(
          rel_tab_restr[, paste0(levels_long[1], "__Blood_CCR8_Treg_pos")] ==
            "Closer to skin Tregs"
        )),
        length(which(
          rel_tab_restr[, paste0(levels_long[1], "__Blood_CCR8_Treg_pos")] ==
            "Closer to blood naive Tregs"
        )),
        length(which(
          rel_tab_restr[, paste0(levels_long[2], "__Blood_CCR8_Treg_pos")] ==
            "Closer to skin Tregs"
        )),
        length(which(
          rel_tab_restr[, paste0(levels_long[2], "__Blood_CCR8_Treg_pos")] ==
            "Closer to blood naive Tregs"
        ))
      ),
      ncol = 2,
      byrow = T)
      rownames(contingency_tab) <- levels_long
      colnames(contingency_tab) <- c("Closer to skin Tregs", "Closer to blood naive Tregs")
      
      # Compute proportions.
      assert(all(rowSums(contingency_tab) == nrow(rel_tab_restr)))
      proportions <- contingency_tab / nrow(rel_tab_restr)
      
      # Perform a two-tailed Fisher's exact test.
      fisher_test_res <- fisher.test(contingency_tab, alternative = "two.sided")
      
      # Save the results from Fisher's exact test together with the contingency table used.
      list_to_save <- list(
        contingency_tab = contingency_tab,
        proportions = proportions,
        fisher_test_res = fisher_test_res
      )
      list_to_save_outfile <- paste0(
        location,
        "/treg_hierarchies/diff_meth_diff_acc_diff_expr_dmr_peak_gene_links_w_ccr8_positioning_fisher_test_hypometh_in_",
        direction,
        "_",
        levels[1],
        "_vs_",
        levels[2],
        "_w_tfidf.rds"
      )
      saveRDS(list_to_save, file = list_to_save_outfile)
      
      # Print the contingency table, the test result and the exact P value to the summary file.
      cat("Contingency table:\n")
      print(contingency_tab)
      cat("\n\n")
      cat("Proportions:\n")
      print(proportions)
      cat("\n\n")
      cat("Results from Fisher's exact test:\n")
      print(fisher_test_res)
      cat("\n")
      cat("Exact P value:\n")
      cat(fisher_test_res$p.value)
      
      # Print a spacer to the summary file.
      cat("\n\n\n\n\n")
      
    }
  }
  
  # Print a spacer to the summary file.
  cat("\n\n\n\n\n\n\n")
  
}

# Stop keeping track of all output.
sink()
