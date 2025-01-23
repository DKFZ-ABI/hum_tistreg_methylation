# This script performs statistical tests for all box plots that we will show in
#   the publication and where statistics is relevant (exluding those showing
#   fat Treg cells because statistics for these plots are included in different
#   scripts).
# Author: Niklas Beumer



# Load required package(s).
library(testit)


# Define a location on /xxx.
location <- "/xxx/nbeumer/hm_treg_bs_rgnsbg"

# Define the output file that will later contain a summary of the test results.
test_res_summary_file <- paste0(location, 
                                "/te_analysis/all_te_boxplots_stat_tests.txt")

# Start capturing output.
sink(test_res_summary_file)





#############################################################################
# Perform the tests.
#############################################################################

############################
# Skin vs. blood naive; LTR TEs; Methylation
# Two-tailed Student's t tests on the mean methylation across all considered TE 
# insertion sites.

# Read in underlying data.
underl_data_file <- paste0(
  location, 
  "/te_analysis/gen_loc_epigen_validation_ltr_discrep_betw_meth_and_atac_boxplots_data_meth.txt"
)
underl_data_lines <- readLines(underl_data_file)
underl_data_lines_refined <- sapply(seq(1, length(underl_data_lines), 2), 
                                    FUN = function(x) {
  paste(underl_data_lines[x], underl_data_lines[x + 1])
})
underl_data <- read.table(text = underl_data_lines_refined, sep = "\t")
assert(all(underl_data$V14 == "LTR"))
underl_data$region_str <- paste(underl_data$V1, underl_data$V2, underl_data$V3, 
                                sep = "_")

# Extract skin and blood naive values.
skin_vals <- underl_data[grepl("Skin Treg", underl_data$V20), ]
naive_vals <- underl_data[grepl("Blood naive Treg", underl_data$V20), ]

# Aggregate values into one value for each sample by computing the mean across
# the TE insertion sites.
skin_samples <- grep("Skin Treg", unique(underl_data$V20), value = T)
skin_means <- sapply(skin_samples, FUN = function(x) {
  mean(underl_data$V21[underl_data$V20 == x])
})
naive_samples <- grep("Blood naive Treg", unique(underl_data$V20), value = T)
naive_means <- sapply(naive_samples, FUN = function(x) {
  mean(underl_data$V21[underl_data$V20 == x])
})

# Perform a two-tailed Student's t test.
test_res <- t.test(skin_means, naive_means, var.equal = T, 
                   alternative = "two.sided")

# Save the test result.
test_outfile <- paste0(
  location, 
  "/te_analysis/te_boxplots_stat_tests_skin_vs_naive_ltr_meth.rds"
)
saveRDS(test_res, file = test_outfile)

# Print the test results to output.
cat("\n\n\n\n### Skin vs. naive: LTR: Meth\n\n")
print(test_res)
cat(paste0("P = ", test_res$p.value, "\n"))



############################
# Skin vs. blood naive; SINE and DNA TEs; Methylation
# Two-tailed Student's t tests on the mean methylation across all considered TE 
# insertion sites.

# Read in underlying data.
underl_data_file <- paste0(
  location, 
  "/te_analysis/gen_loc_epigen_validation_sine_dna_concordance_betw_meth_and_atac_boxplots_data_meth.txt"
)
underl_data_lines <- readLines(underl_data_file)
underl_data_lines_refined <- sapply(
  seq(1, length(underl_data_lines), 2), 
  FUN = function(x) {
    paste(underl_data_lines[x], underl_data_lines[x + 1])
  })
underl_data <- read.table(text = underl_data_lines_refined, sep = "\t")
assert(all(underl_data$V14 %in% c("SINE", "DNA")))
underl_data$region_str <- paste(underl_data$V1, underl_data$V2, underl_data$V3, 
                                sep = "_")

# Extract values for the different TE classes and the different cell types.
te_classes <- unique(underl_data$V14)
cell_types <- c("Skin Treg", "Blood naive Treg")
separate_vals <- lapply(te_classes, FUN = function(x) {
  restr_table <- underl_data[underl_data$V14 == x, ]
  sublist <- lapply(cell_types, FUN = function(y) {
    restr_table[grepl(y, restr_table$V20, fixed = T), ]
  })
  names(sublist) <- cell_types
  return(sublist)
})
names(separate_vals) <- te_classes

# Aggregate values into sample-level values by computing the mean value across
# the considered TE insertion sites.
separate_vals_aggr <- lapply(separate_vals, FUN = function(x) {
  lapply(x, FUN = function(y) {
    samples <- unique(y$V20)
    means <- sapply(samples, FUN = function(x) {
      mean(y$V21[y$V20 == x])
    })
    return(means)
  })
})

# Iterate over the different TE classes and perform two-tailed Student's t
# tests on the means between the considered insertion sites. Print 
# the results to the output.
cat("\n\n\n\n### Skin vs. naive: SINE and DNA: Meth")
test_res <- lapply(te_classes, FUN = function(x) {
  cat(paste("\n\n##", x))
  single_test_res <- t.test(separate_vals_aggr[[x]][[1]],
                            separate_vals_aggr[[x]][[2]],
                            var.equal = T,
                            alternative = "two.sided")
  print(single_test_res)
  cat(paste0("P = ", single_test_res$p.value, "\n"))
  return(single_test_res)
})
names(test_res) <- te_classes

# Save the test result.
test_outfile <- paste0(
  location, 
  "/te_analysis/te_boxplots_stat_tests_skin_vs_naive_sine_and_dna_meth.rds"
)
saveRDS(test_res, file = test_outfile)




############################
# Skin, blood CCR8 and blood naive; Skin Treg hypometh; Methylation
# Two-tailed Student's t test (pair-wise for each cell type 
# comparison) on means of the considered TE insertion sites.

# Read in underlying data.
underl_data_file <- paste0(
  location, 
  "/te_analysis/gen_loc_epigen_validation_te_discrep_blood_ccr8_pos_boxplots_data_meth.txt"
)
underl_data_lines <- readLines(underl_data_file)
underl_data_lines_refined <- sapply(seq(1, length(underl_data_lines), 2), 
                                    FUN = function(x) {
                                      paste(underl_data_lines[x], underl_data_lines[x + 1])
                                    })
underl_data <- read.table(text = underl_data_lines_refined, sep = "\t")
underl_data$region_str <- paste(underl_data$V1, underl_data$V2, underl_data$V3, 
                                sep = "_")

# Extract values for the different TE classes and the different cell types.
te_classes <- unique(underl_data$V20)
cell_types <- c("Skin Treg", "Blood CCR8+ Treg", "Blood naive Treg")
separate_vals <- lapply(te_classes, FUN = function(x) {
  restr_table <- underl_data[underl_data$V20 == x, ]
  sublist <- lapply(cell_types, FUN = function(y) {
    restr_table[grepl(y, restr_table$V21, fixed = T), ]
  })
  names(sublist) <- cell_types
  return(sublist)
})
names(separate_vals) <- te_classes

# Aggregate values into sample-level values by computing the mean across the
# considered TE insertion sites.
separate_vals_aggr <- lapply(separate_vals, FUN = function(x) {
  lapply(x, FUN = function(y) {
    samples <- unique(y$V21)
    means <- sapply(samples, FUN = function(x) {
      mean(y$V22[y$V21 == x])
    })
    return(means)
  })
})

# Iterate over the different TE classes and cell type comparisons and perform
# two-tailed Student's tests on the means between the considered TE insertion 
# sites. Print the results to the output.
comparisons <- do.call(c, lapply(1:(length(cell_types) - 1), FUN = function(x) {
  lapply((x + 1):length(cell_types), FUN = function(y) {
    c(cell_types[x], cell_types[y])
  })
}))
cat("\n\n\n\n### Skin , CCR8 and naive: Meth")
test_res <- lapply(te_classes, FUN = function(x) {
  cat(paste("\n\n##", x))
  comparisons_res <- lapply(comparisons, FUN = function(y) {
    cat(paste("\n\n#", y[1], "vs.", y[2]))
    single_test_res <- t.test(separate_vals_aggr[[x]][[y[1]]],
                              separate_vals_aggr[[x]][[y[2]]],
                              var.equal = T,
                              alternative = "two.sided")
    print(single_test_res)
    cat(paste0("P = ", single_test_res$p.value, "\n"))
    return(single_test_res)
  })
  names(comparisons_res) <- lapply(comparisons, FUN = function(y) {
    paste(y, collapse = "_vs_")
  })
  pvals <- sapply(comparisons_res, FUN = function(y) {
    y$p.value
  })
  pvals_adj <- p.adjust(pvals, method = "BH")
  cat(paste("\n\n# Benjamini-Hochberg-adjusted P values:\n"))
  print(pvals_adj)
  comparisons_res$Pvals_adj_BH <- pvals_adj
  return(comparisons_res)
})
names(test_res) <- te_classes

# Save the test result.
test_outfile <- paste0(
  location, 
  "/te_analysis/te_boxplots_stat_tests_skin_ccr8_naive_skin_treg_hypometh_dmrs_meth.rds"
)
saveRDS(test_res, file = test_outfile)




############################
# Skin, blood CCR8 and blood naive; LTR45B and HERVIP10F-int; Methylation
# Two-tailed Student's t test (pair-wise for each cell type 
#   comparison) on means across the considered TE insertion sites.

# Read in underlying data.
underl_data_file <- paste0(
  location, 
  "/te_analysis/gen_loc_epigen_validation_tes_w_relevance_on_rna_level_boxplots_data_by_special_elements_meth_w_ccr8.txt"
)
underl_data_lines <- readLines(underl_data_file)
underl_data_lines_refined <- sapply(seq(1, length(underl_data_lines), 2), 
                                    FUN = function(x) {
                                      paste(underl_data_lines[x], 
                                            underl_data_lines[x + 1])
                                    })
underl_data <- read.table(text = underl_data_lines_refined, sep = "\t")
underl_data$region_str <- paste(underl_data$V1, underl_data$V2, underl_data$V3, 
                                sep = "_")

# Extract values for the different TE classes and the different cell types.
te_classes <- unique(underl_data$V13)
cell_types <- c("Skin Treg", "Blood CCR8+ Treg", "Blood naive Treg")
separate_vals <- lapply(te_classes, FUN = function(x) {
  restr_table <- underl_data[underl_data$V13 == x, ]
  sublist <- lapply(cell_types, FUN = function(y) {
    restr_table[grepl(y, restr_table$V20, fixed = T), ]
  })
  names(sublist) <- cell_types
  return(sublist)
})
names(separate_vals) <- te_classes

# Aggregate values into sample-level values by computing the mean across the
# considered TE insertion sites.
separate_vals_aggr <- lapply(separate_vals, FUN = function(x) {
  lapply(x, FUN = function(y) {
    samples <- unique(y$V20)
    means <- sapply(samples, FUN = function(x) {
      mean(y$V21[y$V20 == x])
    })
    return(means)
  })
})

# Iterate over the different TE classes and cell type comparisons and perform
# two-tailed Student's t tests on the means between the considered insertion
# sites. Print the results to the output. Adjust P values by Benjamini-
# Hochberg correction.
comparisons <- do.call(c, lapply(1:(length(cell_types) - 1), FUN = function(x) {
  lapply((x + 1):length(cell_types), FUN = function(y) {
    c(cell_types[x], cell_types[y])
  })
}))
cat("\n\n\n\n### Skin , CCR8 and naive: HERVIP10F and LTR45B: Meth")
test_res <- lapply(te_classes, FUN = function(x) {
  cat(paste("\n\n##", x))
  comparisons_res <- lapply(comparisons, FUN = function(y) {
    cat(paste("\n\n#", y[1], "vs.", y[2]))
    single_test_res <- t.test(separate_vals_aggr[[x]][[y[1]]],
                              separate_vals_aggr[[x]][[y[2]]],
                              var.equal = T,
                              alternative = "two.sided")
    print(single_test_res)
    cat(paste0("P = ", single_test_res$p.value, "\n"))
    return(single_test_res)
  })
  names(comparisons_res) <- lapply(comparisons, FUN = function(y) {
    paste(y, collapse = "_vs_")
  })
  pvals <- sapply(comparisons_res, FUN = function(y) {
    y$p.value
  })
  pvals_adj <- p.adjust(pvals, method = "BH")
  cat(paste("\n\n# Benjamini-Hochberg-adjusted P values:\n"))
  print(pvals_adj)
  comparisons_res$Pvals_adj_BH <- pvals_adj
  return(comparisons_res)
})
names(test_res) <- te_classes

# Save the test result.
test_outfile <- paste0(
  location, 
  "/te_analysis/te_boxplots_stat_tests_skin_ccr8_naive_hervip10f_ltr45b_meth.rds"
)
saveRDS(test_res, file = test_outfile)






# Stop capturing output.
sink()