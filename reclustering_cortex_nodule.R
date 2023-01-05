set.seed(1407)
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

# Where to save the results
folder_name <- "RECLUSTERING/CORTEX_NOD_sunn4/"

( p0 <- plot_cells(cds,
                   color_cells_by="cluster",
                   cell_size = 0.5,
                   group_label_size = 5) )

colData(cds)$cluster <- monocle3::clusters(cds)

cluster_number = c( 
    2, 16, # Cortex
    6 #nodule
)

cells_on_clusters <- cds@clusters$UMAP$clusters
cells_on_clusters <- cells_on_clusters[cells_on_clusters %in% cluster_number]
cells_on_clusters <- names(cells_on_clusters)

cds_subset <- cds[, colnames(cds) %in% cells_on_clusters ]

( p1_both <- plot_cells(cds_subset,
                        cell_size = 0.75,
                        group_label_size = 8) )

## Select cells from sunn-4
cells_on_sunn <- cds_subset@colData
cells_on_sunn <- cells_on_sunn[cells_on_sunn$Group == "Sunn-4", ]
cells_on_sunn_names <- unique(rownames(cells_on_sunn))

cds_subset <- cds_subset[, colnames(cds_subset) %in% cells_on_sunn_names ]

( p1_sunn4 <- plot_cells(cds_subset,
                   cell_size = 0.75,
                   group_label_size = 8) )

( p1.1 <- plot_cells(cds_subset,
                     cell_size = 0.75,
                     group_label_size = 8)  +
    facet_wrap(~Group + timepoint, nrow = 1) )

system( paste0( "mkdir -p ", folder_name, "clustering_images") )
ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/cortex_nod_selected_clusters_both_genotypes.svg"),
                p1_both,
                height=12,
                width=12,
                units="cm",
                bg = "#FFFFFF", 
                dpi = 300)

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/cortex_nod_selected_clusters_sunn4.svg"),
                p1_sunn4,
                height=12,
                width=16,
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

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/cortex_nod_selected_clusters_both_genotypes_by_timepoint.svg"),
                p1.1,
                height=12,
                width=48,
                units="cm",
                bg = "#FFFFFF", 
                dpi = 300)

system("mkdir -p RECLUSTERING/CORTEX_NOD_sunn4/rds_file_subset/")
saveRDS(cds_subset, file = "RECLUSTERING/CORTEX_NOD_sunn4/rds_file_subset/medicago_integrated_selected_clusters_cortex_nodule.rds")

## Removes cells from lateral root meristem
#cds_subset <- monocle3::choose_cells(cds_subset)
#cells_to_exclude <- as.data.frame( rownames(colData(cds_subset)) )
#write.table(cells_to_exclude,
#            "cells_to_exclude.tsv",
#            col.names = F, row.names = F, quote = F, sep = "\t")

cells_to_exclude <- vroom::vroom(
    "cells_to_exclude.tsv",
    col_names = F,
    delim = "\t")

cds_subset <- cds_subset[, !colnames(cds_subset) %in% cells_to_exclude$X1 ]

( p1 <- plot_cells(cds_subset,
                   cell_size = 0.75,
                   group_label_size = 8) )

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/cortex_nod_selected_clusters_sunn4_after_removing_0h.svg"),
                p1,
                height=12,
                width=16,
                units="cm",
                bg = "#FFFFFF", 
                dpi = 300)

( p1.1 <- plot_cells(cds_subset,
                     cell_size = 0.75,
                     group_label_size = 8)  +
    facet_wrap(~Group + timepoint, nrow = 1) )

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/cortex_nod_selected_clusters_by_time_and_genotype_after_removing_0h.svg"),
                p1.1,
                height=12,
                width=48,
                units="cm",
                bg = "#FFFFFF", 
                dpi = 300)


saveRDS(cds_subset, file = "RECLUSTERING/CORTEX_NOD_sunn4/rds_file_subset/medicago_integrated_selected_clusters_cortex_nodule_after_removing_LR.rds")

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/cortex_nod_selected_clusters_AFTER_REMOVING_LR.png"),
                p1.1,
                height=12,
                width=16,
                units="cm",
                bg = "#FFFFFF", 
                dpi = 300)

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/cortex_nod_selected_clusters_by_time_and_genotype_AFTER_REMOVING_LR.png"),
                p1.1,
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
                   cell_size = 0.5) )

( p3 <- plot_cells(cds_subset,
                   group_label_size = 5,
                   trajectory_graph_segment_size = 1.5,
                   graph_label_size = 5,
                   cell_size = 0.5)  +
        facet_wrap(~Group + timepoint, nrow = 1) +
        theme(strip.text = element_blank(),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_text(size = 40,
                                        hjust = 0.5),
              legend.title=element_blank()) )

ggplot2::ggsave(
    paste0(folder_name,
           "clustering_images",
           "/cortex_nod_selected_clusters_RECLUSTERED.png"),
    p2,
    height=12,
    width=16,
    units="cm",
    bg = "#FFFFFF", 
    dpi = 300)

ggplot2::ggsave(
    paste0(folder_name,
           "clustering_images",
           "/cortex_nod_selected_clusters_RECLUSTERED.svg"),
    p2,
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

write.csv( marker_test_res, 
           paste0( folder_name, 
                   "top_1000_markers_per_cluster",
                   "/top_1000_markers_per_cluster.csv"),
           row.names = F) 

filtering_criteria = "specificity"

folder_name <- paste0(folder_name, "/top_1000_markers_per_cluster/", "by_", filtering_criteria )
system( paste0("mkdir -p ", folder_name) )
for (c in unique(marker_test_res$cell_group) ) {
    
    top_specific_markers <- marker_test_res %>%
        dplyr::filter(marker_test_q_value < 0.05) %>%
        dplyr::filter(fraction_expressing >= 0.10) %>%
        dplyr::filter(cell_group == c) %>%
        dplyr::arrange( desc(specificity) ) 
    
    top_specific_markers_top1000 <- top_specific_markers[1:1000, ]
    
    write.csv(top_specific_markers_top1000, 
              paste0(folder_name,
                     "/list_of_top_1000_for_cluster_",
                     c,
                     "_",filtering_criteria,
                     ".csv"),
              row.names = F)
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
    
    write.csv(top_specific_markers_top100, 
              paste0(folder_name,
                     "/list_of_top_100_for_cluster_",
                     c,
                     "_",filtering_criteria,
                     ".csv"),
              row.names = F)
}