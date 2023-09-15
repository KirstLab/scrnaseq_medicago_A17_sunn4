set.seed(1407)

outfile <- "logs/reclustering_roothair.out" # File name of output log
#Check its existence
if ( file.exists(outfile) ) {
    #Delete file if it exists
    file.remove(outfile)
}

my_log <- file(outfile) 
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

suppressMessages( require(monocle3) )
suppressMessages( require(vroom) )
suppressMessages( require(cowplot) )
suppressMessages( require(ggplot2) )
suppressMessages( require(dplyr) )
suppressMessages( require(stringr) )

############################################
## Selecting cells of Epidermis/root hair ##
############################################
cds <- readRDS("rds_files/batched_integrated_clustered_complete_dataset.rds")

# Where to save the results
folder_name <- "RECLUSTERING/Epidermis_roothair/"

( p0 <- plot_cells(cds,
                 color_cells_by="cluster",
                 cell_size = 0.5,
                 group_label_size = 5) )

cluster_number = c(1, 14, 22, 12, 3)

cells_on_clusters <- cds@clusters$UMAP$clusters
cells_on_clusters <- cells_on_clusters[cells_on_clusters %in% cluster_number]
cells_on_clusters <- names(cells_on_clusters)

cds_subset <- cds[, colnames(cds) %in% cells_on_clusters ]

( p1 <- plot_cells(cds_subset,
                   group_label_size = 7,
                   cell_size = 0.75) )

system( paste0( "mkdir -p ", folder_name, "clustering_images") )
ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/epidermis_selected_clusters.svg"),
                p1,
                height=12,
                width=16,
                units="cm",
                bg = "#FFFFFF", 
                dpi = 300)
# Reclustering
cds_subset <- clear_cds_slots(cds_subset)

cds_subset <- monocle3::preprocess_cds(cds_subset, num_dim = 100, verbose = T)
cds_subset <- reduce_dimension(cds_subset, verbose = T)

cds_subset <- cluster_cells(cds_subset,
                            verbose = T,
                            resolution = 1e-03,
                            random_seed = 1407)

( p2 <- plot_cells(cds_subset, 
                   group_label_size = 7,
                   cell_size = 0.75) )

( p3 <- plot_cells(cds_subset,
                   group_label_size = 5,
                   trajectory_graph_segment_size = 1.5,
                   graph_label_size = 5,
                   cell_size = 0.5)  +
        facet_wrap(~Group + timepoint, nrow = 2) +
        theme(strip.text = element_blank(),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_text(size = 40,
                                        hjust = 0.5),
              legend.title=element_blank()) )

( p3.1 <- plot_cells(cds_subset,
                   group_label_size = 5,
                   trajectory_graph_segment_size = 1.5,
                   graph_label_size = 5,
                   cell_size = 0.75)  +
        facet_wrap(~timepoint, nrow = 1) +
        theme(strip.text = element_blank(),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_text(size = 40,
                                        hjust = 0.5),
              legend.title=element_blank()) )

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/epidermis_selected_clusters_RECLUSTERED_by_timepoint.svg"),
                p3.1,
                height=12,
                width=48,
                units="cm",
                bg = "#FFFFFF", 
                dpi = 300)

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/epidermis_selected_clusters_RECLUSTERED_by_time_and_genotype.png"),
                p3,
                height=25,
                width=40,
                units="cm",
                bg = "#FFFFFF",
                dpi = 300)

system( paste0("mkdir -p ", folder_name, "rds_file_subset") )
saveRDS(cds_subset,
        paste0(folder_name,
               "rds_file_subset",
               "/medicago_integrated_subset_epidermis.rds") )

## Identification of markers for each cluster
system( paste0("mkdir -p ",
               folder_name, 
               "top_1000_markers_per_cluster") )

marker_test_res <- top_markers(cds_subset, 
                               group_cells_by="cluster",
                               genes_to_test_per_group = 1000,
                               marker_sig_test = T,
                               cores = 8,
                               verbose = T)

write.csv( marker_test_res, 
           paste0( folder_name, 
                   "top_1000_markers_per_cluster",
                   "/top_1000_markers_per_cluster.csv"),
           row.names = F) 

filtering_criteria = "specificity"

folder_name <- paste0( folder_name, "top_1000_markers_per_cluster/", "by_", filtering_criteria )
system( paste0("mkdir -p ", folder_name) )
for (c in unique(marker_test_res$cell_group) ) {
    
    top_specific_markers <- marker_test_res %>%
        dplyr::filter(marker_test_q_value < 0.05) %>%
        dplyr::filter(fraction_expressing >= 0.20) %>%
        dplyr::filter(cell_group == c) %>%
        dplyr::arrange( desc(specificity) ) 
    
    write.csv(top_specific_markers, 
              paste0(folder_name,
                     "/list_of_top_1000_for_cluster_",
                     c,
                     "_",filtering_criteria,
                     ".csv"),
              row.names = F)
}

closeAllConnections()
