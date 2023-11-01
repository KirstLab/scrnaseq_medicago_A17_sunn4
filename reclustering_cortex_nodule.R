set.seed(1407)

outfile <- "logs/reclustering_cortex_nodule.out"
if ( file.exists(outfile) ) {
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
## Selecting cells from cortex and nodule ##
############################################
cds <- readRDS("rds_files/batched_integrated_clustered_complete_dataset.rds")
folder_name <- "RECLUSTERING/CORTEX_NOD/"

( p0 <- plot_cells(cds,
                   color_cells_by="cluster",
                   cell_size = 1,
                   group_label_size = 5)  +
        facet_wrap(~timepoint) )

colData(cds)$cluster <- monocle3::clusters(cds)

cluster_number = c(1, 3, 4, 6, 9, 15)

cells_on_clusters <- cds@clusters$UMAP$clusters
cells_on_clusters <- cells_on_clusters[cells_on_clusters %in% cluster_number]
cells_on_clusters <- names(cells_on_clusters)

cds_subset <- cds[, colnames(cds) %in% cells_on_clusters ]

( p1 <- plot_cells(cds_subset,
                   cell_size = 1,
                   group_label_size = 8) )

( p1.1 <- plot_cells(cds_subset,
                     cell_size = 1,
                     group_label_size = 8)  +
        facet_wrap(~timepoint, nrow = 1) )

system( paste0( "mkdir -p ", folder_name, "clustering_images") )
ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/cortex_nod_selected_clusters_both_genotypes.png"),
                p1,
                height=12,
                width=12,
                units="cm",
                bg = "#FFFFFF",
                dpi = 300)

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/cortex_nod_selected_clusters_both_genotypes.svg"),
                p1,
                height=12,
                width=12,
                units="cm",
                bg = "#FFFFFF",
                dpi = 300)

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/cortex_nod_selected_clusters_by_time_and_genotype.png"),
                p1.1,
                height=12,
                width=48,
                units="cm",
                bg = "#FFFFFF",
                dpi = 300)

system( paste0("mkdir -p ", folder_name, "rds_file_subset/") )
saveRDS(cds_subset,
        file = paste0(folder_name, "rds_file_subset/", "medicago_integrated_selected_clusters_cortex_nodule.rds") )

# Reclustering
cds_subset <- clear_cds_slots(cds_subset)
cds_subset <- monocle3::preprocess_cds(cds_subset,
                                       num_dim = 100)

cds_subset <- reduce_dimension(cds_subset, verbose = T)
cds_subset = align_cds(cds_subset,
                       num_dim = 100,
                       alignment_group = "batch2", verbose = T)

( p2 <- plot_cells(cds_subset,
                   color_cells_by="batch2",
                   label_cell_groups=FALSE) )

res <- c(1, 0.1, 0.001, 0.0001, 0.00001)
cds_subset <- cluster_cells(cds_subset,
                            verbose = T, 
                            resolution = res, 
                            random_seed = 1407)

( p2 <- plot_cells(cds_subset, 
                   group_label_size = 7,
                   cell_size = 1) )

( p2.1 <- plot_cells(cds_subset,
                     group_label_size = 5,
                     trajectory_graph_segment_size = 1.5,
                     graph_label_size = 5,
                     cell_size = 1)  +
        facet_wrap(~ timepoint, nrow = 1) +
        theme(strip.text = element_blank(),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_text(size = 40,
                                        hjust = 0.5),
              legend.title=element_blank()) )

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/cortex_nod_complete_clustering_before_filtering.png"),
                p2,
                height=12,
                width=12,
                units="cm",
                bg = "#FFFFFF",
                dpi = 300)

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/cortex_nod_complete_clustering_before_filtering.png"),
                p2.1,
                height=12,
                width=12,
                units="cm",
                bg = "#FFFFFF",
                dpi = 300)

saveRDS(cds_subset, "RECLUSTERING/CORTEX_NOD/rds_file_subset/cortex_nod_selected_clusters_RECLUSTERED_before_filterind.rds")

cluster_number = c(7,18,16)
cells_on_clusters <- cds_subset@clusters$UMAP$clusters
cells_on_clusters <- cells_on_clusters[cells_on_clusters %in% cluster_number]
cells_on_clusters <- names(cells_on_clusters)

cds_subset <- cds_subset[, !colnames(cds_subset) %in% cells_on_clusters ]

( p3 <- plot_cells(cds_subset,
                   group_label_size = 5,
                   trajectory_graph_segment_size = 1.5,
                   graph_label_size = 5,
                   cell_size = 0.5)  +
        facet_wrap(~ timepoint, nrow = 1) +
        theme(strip.text = element_blank(),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_text(size = 40,
                                        hjust = 0.5),
              legend.title=element_blank()) )

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/cortex_nod_complete_clustering_after_removal_of_LR.png"),
                p3,
                height=12,
                width=12,
                units="cm",
                bg = "#FFFFFF",
                dpi = 300)

# Reclustering
cds_subset <- clear_cds_slots(cds_subset)
cds_subset <- monocle3::preprocess_cds(cds_subset,
                                       num_dim = 100)
cds_subset <- reduce_dimension(cds_subset, verbose = T)
cds_subset = align_cds(cds_subset,
                       num_dim = 100,
                       alignment_group = "batch2", verbose = T)

( p4 <- plot_cells(cds_subset,
                   color_cells_by="batch2",
                   label_cell_groups=FALSE) )

res <- c(1, 0.1, 0.001, 0.0001, 0.00001)
cds_subset <- cluster_cells(cds_subset,
                            verbose = T, 
                            resolution = res, 
                            random_seed = 1407)

( p5 <- plot_cells(cds_subset, 
                   group_label_size = 8,
                   cell_size = 0.75)  +
        theme(strip.text = element_blank(),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_text(size = 40,
                                        hjust = 0.5),
              legend.title=element_blank()) )

ggplot2::ggsave(
    paste0(folder_name,
           "clustering_images",
           "/reclustered_cortex_nod_edited.svg"),
    p5,
    height=12,
    width=16,
    units="cm",
    bg = "#FFFFFF", 
    dpi = 300)

( p5.1 <- plot_cells(cds_subset,
                     group_label_size = 5,
                     trajectory_graph_segment_size = 1.5,
                     graph_label_size = 5,
                     cell_size = 1)  +
        facet_wrap(~ timepoint, nrow = 1) +
        theme(strip.text = element_blank(),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_text(size = 40,
                                        hjust = 0.5),
              legend.title=element_blank()) )

( p6 <- plot_cells(cds_subset, 
                   group_label_size = 7,
                   cell_size = 0.5) )

( p7 <- plot_cells(cds_subset, 
                   group_label_size = 7,
                   cell_size = 0.5) )
saveRDS(cds_subset, "RECLUSTERING/CORTEX_NOD/rds_file_subset/cortex_nod_selected_clusters_RECLUSTERED.rds")

ggplot2::ggsave(
    paste0(folder_name,
           "clustering_images",
           "/cortex_nod_selected_clusters_RECLUSTERED.png"),
    p6,
    height=12,
    width=16,
    units="cm",
    bg = "#FFFFFF", 
    dpi = 300)

ggplot2::ggsave(
    paste0(folder_name,
           "clustering_images",
           "/cortex_nod_selected_clusters_RECLUSTERED_with_trajectory.png"),
    p6,
    height=12,
    width=16,
    units="cm",
    bg = "#FFFFFF", 
    dpi = 300)

ggplot2::ggsave(
    paste0(folder_name,
           "clustering_images",
           "/cortex_nod_selected_clusters_RECLUSTERED_by_time_and_genotype.png"),
    p3,
    height=12,
    width=40,
    units="cm",
    bg = "#FFFFFF",
    dpi = 300)

( p3.1 <- plot_cells(cds_subset,
                     group_label_size = 5,
                     trajectory_graph_segment_size = 1.5,
                     graph_label_size = 5,
                     cell_size = 0.75)  +
        facet_wrap(~ timepoint, nrow = 1) +
        theme(strip.text = element_blank(),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_text(size = 40,
                                        hjust = 0.5),
              legend.title=element_blank()) )

ggplot2::ggsave(
    paste0(folder_name,
           "clustering_images",
           "/cortex_nod_selected_clusters_RECLUSTERED_by_timepoint.svg"),
    p3.1,
    height=12,
    width=48,
    units="cm",
    bg = "#FFFFFF",
    dpi = 300)

system( paste0("mkdir -p ", folder_name, "rds_file_subset") )
saveRDS(cds_subset,
        paste0(folder_name,
               "rds_file_subset",
               "/medicago_integrated_subset_cortex_nodule.rds") )

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

write.table( marker_test_res, 
             paste0( folder_name, 
                     "top_1000_markers_per_cluster",
                     "/top_1000_markers_per_cluster.tsv"),
             row.names = F, col.names = T, sep = "\t") 

filtering_criteria = "specificity"

folder_name <- paste0(folder_name, "top_1000_markers_per_cluster/", "by_", filtering_criteria )
system( paste0("mkdir -p ", folder_name) )
for (c in unique(marker_test_res$cell_group) ) {
    
    top_specific_markers <- marker_test_res %>%
        dplyr::filter(marker_test_q_value < 0.05) %>%
        dplyr::filter(fraction_expressing >= 0.10) %>%
        dplyr::filter(cell_group == c) %>%
        dplyr::arrange( desc(specificity) )
    
    top_specific_markers_top1000 <- top_specific_markers[1:1000, ]
    
    write.table(top_specific_markers, 
                paste0(folder_name,
                       "/list_of_top_1000_for_cluster_",
                       c,
                       "_",filtering_criteria,
                       ".tsv"),
                row.names = F, col.names = T, sep = "\t") 
}

folder_name <- paste0(folder_name, "/top_100_markers_per_cluster/", "by_", filtering_criteria )
system( paste0("mkdir -p ", folder_name) )
for (c in unique(marker_test_res$cell_group) ) {
    
    top_specific_markers <- marker_test_res %>%
        dplyr::filter(marker_test_q_value < 0.05) %>%
        dplyr::filter(fraction_expressing >= 0.10) %>%
        dplyr::filter(cell_group == c) %>%
        dplyr::arrange( desc(specificity) ) 
    
    top_specific_markers_top100 <- top_specific_markers[1:100, ]
    
    write.table(top_specific_markers_top100, 
                paste0(folder_name,
                       "/list_of_top_100_for_cluster_",
                       c,
                       "_",filtering_criteria,
                       ".tsv"),
                row.names = F, col.names = T, sep = "\t") 
}

closeAllConnections()