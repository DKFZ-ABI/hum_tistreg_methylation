# This script visualises clonotype sharing between skin Treg cells, blood CCR8+
#   Treg cells and blood CD45RA+ Treg cells.
# Author: Niklas Beumer



# Load required package(s).
library(Seurat)
library(ggplot2)
library(viridis)
library(testit)
library(grid)
library(gridExtra)


# Define a location on /xxx.
b330_space <- "/xxx/"
location <- paste0(b330_space, "nbeumer/hm_treg_bs_rgnsbg")
location_2 <- "/xxx/data/nbeumer"

# Create an output directory for plots, if it doesn't already exist.
plot_outdir <- paste("/yyy/hm_treg_bs_rgnsbg/analysis", 
                     format(Sys.time(), "%m-%d-%y"), sep = "/")
if (!dir.exists(plot_outdir)) {dir.create(plot_outdir)}

# Specify the output directory for the RDS files containing the plots.
plot_rds_outdir <- paste0(location, "/plot_rds_objects")

# Read in the Seurat objects containing scRNA-seq data with clone information.
sc_data_files <- paste0(
  location_2, 
  c("/scrnatcr_donor6_seur_obj_w_clonotype_info_and_matching.rds",
    "/scrnatcr_donor7_seur_obj_w_clonotype_info_and_matching.rds")
)
sc_data <- lapply(sc_data_files, FUN = readRDS)
names(sc_data) <- c("Donor_6", "Donor_7")






#############################################################################
# Visualise clonotype sharing between (i) skin Tregs and blood CCR8+ Tregs 
# and (ii) skin Tregs and blood naive Tregs.
# Two cells are considered clonotype-matched if their clonotype has the
# same nucleotide sequence string and this sequence string contains at least
# one alpha chain sequence and one beta chain sequence.
#############################################################################

# Iterate across the samples to collect clonotype sharing information.
clonotype_sharing <- do.call(rbind, lapply(names(sc_data), FUN = function(x) {
  
  # Identify all clonotypes that appear in skin Tregs and contain at least one
  # alpha chain and at least one beta chain. 
  # Assign more comfortable labels to these clonotypes.
  skin_clonotypes <- unique(sc_data[[x]]$cdr3s_nt[sc_data[[x]]$is_skin_treg])
  skin_clonotypes_filt <- grep("TRA", 
                               grep("TRB", 
                                    skin_clonotypes, 
                                    value = T), 
                               value = T)
  names(skin_clonotypes_filt) <- paste0(x, "_skin_clonotype_", 
                                        1:length(skin_clonotypes_filt))
  
  # Compute the total number of skin Tregs, blood CCR8+ Tregs and blood naive
  # Tregs from this donor.
  # Ignore cells that do not have any clonotype information in order to 
  # avoid biases later.
  total_skintreg_num <- length(
    which(sc_data[[x]]$is_skin_treg & !(is.na(sc_data[[x]]$cdr3s_nt)))
  )
  total_boodccr8treg_num <- length(
    which(sc_data[[x]]$Anno_Niklas_for_TCR_matching == "Blood_CCR8_Tregs" & 
            !(is.na(sc_data[[x]]$cdr3s_nt)))
  )
  total_boodnaivetreg_num <- length(
    which(sc_data[[x]]$Anno_Niklas_for_TCR_matching == "Blood_RA_Tregs" & 
            !(is.na(sc_data[[x]]$cdr3s_nt)))
  )
  
  # For each of the skin Treg clonotypes, compute how many skin Tregs, blood 
  # CCR8+ Tregs and blood naive Tregs belong to this clonotype (and the 
  # corresponding percentages).
  clonotypes_by_cell <- sc_data[[x]]$cdr3s_nt
  skin_anno <- sc_data[[x]]$is_skin_treg
  naive_ccr8_anno <- sc_data[[x]]$Anno_Niklas_for_TCR_matching
  this_donor_results <- do.call(
    rbind, 
    lapply(names(skin_clonotypes_filt), FUN = function(y) {
      skin_num <- length(
        which(clonotypes_by_cell == skin_clonotypes_filt[y] & 
                skin_anno & 
                !(is.na(clonotypes_by_cell)))
      )
      skin_perc <- 100 * skin_num / total_skintreg_num
      bloodccr8_num <- length(
        which(clonotypes_by_cell == skin_clonotypes_filt[y] & 
                naive_ccr8_anno == "Blood_CCR8_Tregs" & 
                !(is.na(clonotypes_by_cell)))
      )
      bloodccr8_perc <- 100 * bloodccr8_num / total_boodccr8treg_num
      bloodnaive_num <- length(
        which(clonotypes_by_cell == skin_clonotypes_filt[y] & 
                naive_ccr8_anno == "Blood_RA_Tregs" & 
                !(is.na(clonotypes_by_cell)))
      )
      bloodnaive_perc <- 100 * bloodnaive_num / total_boodnaivetreg_num
      temp_df <- data.frame(
        Donor = x,
        Clonotype = y,
        Clonotype_seq = skin_clonotypes_filt[y],
        Cell_type = c("Skin Treg", "Blood CCR8+ Treg", "Blood RA+ Treg"),
        Num_cells_from_celltype_w_clonotype = c(skin_num,
                                                bloodccr8_num,
                                                bloodnaive_num),
        Perc_cells_from_celltype_w_clonotype = c(skin_perc,
                                                 bloodccr8_perc,
                                                 bloodnaive_perc)
      )
      return(temp_df)
    }))
  
  # Return the collected information.
  return(this_donor_results)
  
}))

# Save the collected clonotype sharing information.
clonotype_sharing_outfile <- paste0(
  location, 
  "/treg_hierarchies/scrnatcr_clonotypesharing_skintreg_bloodccr8_bloodnaive_centred_on_skintreg.txt"
)
write.table(clonotype_sharing, file = clonotype_sharing_outfile, sep = "\t",
            row.names = F, col.names = T, quote = F)

# Visualise the extent of clonotype sharing.
# To this end, generate a plot for each donor.
# Also include the percentages of shown skin Treg clonotypes that can be
# foundd in the other two cell types.
clonotype_sharing_plots <- do.call(c, lapply(names(sc_data), FUN = function(x) {
  
  # Extract clonotype sharing data for this donor.
  plotting_data <- clonotype_sharing[clonotype_sharing$Donor == x, ]
  plotting_data$Cell_type <- factor(plotting_data$Cell_type, 
                                    levels = c("Blood RA+ Treg", 
                                               "Skin Treg", 
                                               "Blood CCR8+ Treg"))
  
  # Order the clonotypes by the percentage of skin Tregs they belong to.
  assert(length(unique(plotting_data$Clonotype)) == nrow(plotting_data) / 3)
  plotting_data$Clonotype <- factor(
    plotting_data$Clonotype, 
    levels = unique(plotting_data$Clonotype)[
      order(sapply(unique(plotting_data$Clonotype), FUN = function(y) {
        plotting_data$Perc_cells_from_celltype_w_clonotype[
          plotting_data$Clonotype == y & plotting_data$Cell_type == "Skin Treg"
        ]
      }), decreasing = T)
    ]
  )
  
  # Find out which clonotypes are shared between two cell types.
  skin_ccr8_shared_clonotypes <- plotting_data$Clonotype[
    plotting_data$Cell_type == "Blood CCR8+ Treg" & 
      plotting_data$Perc_cells_from_celltype_w_clonotype > 0
  ]
  skin_naive_shared_clonotypes <- plotting_data$Clonotype[
    plotting_data$Cell_type == "Blood RA+ Treg" & 
      plotting_data$Perc_cells_from_celltype_w_clonotype > 0
  ]
  
  # Find out where borders in the stacked bar charts are.
  cumsums_skin <- cumsum(
    sapply(rev(levels(plotting_data$Clonotype)), FUN = function(y) {
      plotting_data$Perc_cells_from_celltype_w_clonotype[
        plotting_data$Cell_type == "Skin Treg" & plotting_data$Clonotype == y
      ]
    })
  )
  cumsums_ccr8 <- cumsum(
    sapply(rev(levels(plotting_data$Clonotype)), FUN = function(y) {
      plotting_data$Perc_cells_from_celltype_w_clonotype[
        plotting_data$Cell_type == "Blood CCR8+ Treg" & 
          plotting_data$Clonotype == y
      ]
    })
  )
  cumsums_naive <- cumsum(
    sapply(rev(levels(plotting_data$Clonotype)), FUN = function(y) {
      plotting_data$Perc_cells_from_celltype_w_clonotype[
        plotting_data$Cell_type == "Blood RA+ Treg" & 
          plotting_data$Clonotype == y
      ]
    })
  )
  
  # Generate data to link fields in the different bars. A link is only plotted
  # if a clonotype appeared in both linked cell types.
  skin_naive_links_x_coord <- rep(c(1.251, 1.749, 1.749, 1.251), 
                                  length(skin_naive_shared_clonotypes))
  skin_naive_links_y_coord <- do.call(
    c, 
    lapply(skin_naive_shared_clonotypes, FUN = function(y) {
      this_clonotype_ind <- which(rev(levels(plotting_data$Clonotype)) == y)
      if (this_clonotype_ind > 1) {
        previous_clonotype_name <- rev(levels(plotting_data$Clonotype))[
          this_clonotype_ind - 1
        ]
      }
      coordinates_vect <- c(
        ifelse(this_clonotype_ind > 1, 
               yes = cumsums_naive[previous_clonotype_name], 
               no = 0),
        ifelse(this_clonotype_ind > 1,
               yes = cumsums_skin[previous_clonotype_name],
               no = 0),
        cumsums_skin[as.character(y)],
        cumsums_naive[as.character(y)]
      )
      return(coordinates_vect)
    })
  )
  skin_ccr8_links_x_coord <- rep(c(2.251, 2.749, 2.749, 2.251), 
                                  length(skin_ccr8_shared_clonotypes))
  skin_ccr8_links_y_coord <- do.call(
    c, 
    lapply(skin_ccr8_shared_clonotypes, FUN = function(y) {
      this_clonotype_ind <- which(rev(levels(plotting_data$Clonotype)) == y)
      if (this_clonotype_ind > 1) {
        previous_clonotype_name <- rev(levels(plotting_data$Clonotype))[
          this_clonotype_ind - 1
        ]
      }
      coordinates_vect <- c(
        ifelse(this_clonotype_ind > 1,
               yes = cumsums_skin[previous_clonotype_name],
               no = 0),
        ifelse(this_clonotype_ind > 1,
               yes = cumsums_ccr8[previous_clonotype_name],
               no = 0),
        cumsums_ccr8[as.character(y)],
        cumsums_skin[as.character(y)]
      )
      return(coordinates_vect)
    })
  )
  polygon_data_skinnaive <- data.frame(
    x = skin_naive_links_x_coord,
    y = skin_naive_links_y_coord,
    Clonotype = rep(skin_naive_shared_clonotypes, each = 4)
  )
  polygon_data_skinccr8 <- data.frame(
    x = skin_ccr8_links_x_coord,
    y = skin_ccr8_links_y_coord,
    Clonotype = rep(skin_ccr8_shared_clonotypes, each = 4)
  )
  
  # Generate the plot showing clonotype abundances and clonotype matches.
  clonotype_sharing_plot <- ggplot(plotting_data) +
    aes(x = Cell_type, y = Perc_cells_from_celltype_w_clonotype, 
        fill = Clonotype) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      name = "Percentage of cells belonging to the clonotype\namong cells with cclonotype information [%]"
    ) +
    scale_fill_manual(values = turbo(length(levels(plotting_data$Clonotype)), 
                                     direction = -1),
                      guide = "none") +
    geom_bar(stat = "identity", colour = "black", size = 0.005, width = 0.5) +
    geom_polygon(data = polygon_data_skinnaive, 
                 mapping = aes(x = x, y = y, fill = Clonotype),
                 alpha = 0.5,
                 inherit.aes = F) +
    geom_polygon(data = polygon_data_skinccr8, 
                 mapping = aes(x = x, y = y, fill = Clonotype),
                 alpha = 0.5,
                 inherit.aes = F) +
    ggtitle(x) +
    xlab("Cell type") +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
           axis.ticks = element_line(colour = "black"))
  
  # Compute the percentages of shown skin Treg clonotypes that are shared with
  # the other two cell types. Prepare these values to be plotted.
  assert(length(skin_naive_shared_clonotypes) == 
           length(unique(skin_naive_shared_clonotypes)))
  assert(length(skin_ccr8_shared_clonotypes) == 
           length(unique(skin_ccr8_shared_clonotypes)))
  perc_of_skin_clonotypes_found_in_naive <- 
    100 * length(skin_naive_shared_clonotypes) / 
    length(unique(plotting_data$Clonotype))
  perc_of_skin_clonotypes_found_in_ccr8 <- 
    100 * length(skin_ccr8_shared_clonotypes) / 
    length(unique(plotting_data$Clonotype))
  clonotype_recall_data <- data.frame(
    Cell_type = factor(c("Blood RA+ Treg", "Blood CCR8+ Treg"),
                       levels = c("Blood RA+ Treg", "Blood CCR8+ Treg")),
    recall = c(perc_of_skin_clonotypes_found_in_naive, 
               perc_of_skin_clonotypes_found_in_ccr8),
    text_y_pos = rep(60, 2),
    text = paste0(round(c(perc_of_skin_clonotypes_found_in_naive, 
                          perc_of_skin_clonotypes_found_in_ccr8),
                        2),
                  "%")
  )
  
  # Generate a plot showing what percentage of skin Treg clonotypes was
  # also found in the other two cell types.
  clonotype_recall_plot <- ggplot(clonotype_recall_data) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       name = "Recall of skin Treg clonotypes") +
    geom_bar(stat = "identity",
             mapping = aes(x = Cell_type, y = recall),
             fill = "black",
             width = 0.5) +
    geom_text(mapping = aes(x = Cell_type, y = text_y_pos, label = text)) +
    ggtitle(x) +
    xlab("Cell type") +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"))
  
  # Return the generated plot.
  return(list(clonotype_sharing_plot, clonotype_recall_plot))
  
}))

# Save the generated plots.
clonotype_sharing_pdf <- paste0(
  plot_outdir, 
  "/scrnatcr_clonotypesharing_skintreg_bloodccr8_bloodnaive_centred_on_skintreg.pdf"
)
clonotype_sharing_rds <- paste0(
  plot_rds_outdir, 
  "/scrnatcr_clonotypesharing_skintreg_bloodccr8_bloodnaive_centred_on_skintreg.rds"
)
pdf(clonotype_sharing_pdf, width = 8, height = 10)
grid.draw(arrangeGrob(grobs = clonotype_sharing_plots[c(1, 3, 2, 4)], ncol = 2))
dev.off()
saveRDS(clonotype_sharing_plots, file = clonotype_sharing_rds)

