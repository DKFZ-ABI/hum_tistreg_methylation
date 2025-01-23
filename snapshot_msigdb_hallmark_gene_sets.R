# This R script saves the MSigDB hallmark gene sets to the cluster. 
# Author: Niklas Beumer



# Load required packages.
library(msigdbr)

# Specify a location on /yyy.
location <- 
  "/yyy/hm_treg_bs_rgnsbg/external_data"

# Get the current date.
curr_date <- format(Sys.time(), "%Y-%m-%d")

# Retrieve the gene sets (for humans).
hallm_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")

# Save the gene set collection snapshots.
saveRDS(hallm_gene_sets, 
        file = paste0(location, 
                      "/", 
                      curr_date, 
                      "_snapshot_msigdb_hallmark_gene_sets_H_human.rds"))



