set.seed(1407)

suppressMessages( require(monocle3) )
suppressMessages( require(vroom) )
suppressMessages( require(cowplot) )
suppressMessages( require(ggplot2) )
suppressMessages( require(dplyr) )
suppressMessages( require(stringr) )

cds <- readRDS("rds_files/batched_integrated_clustered_complete_dataset.rds")

( p0 <- plot_cells(cds,
                   color_cells_by="cluster",
                   cell_size = 0.5,
                   group_label_size = 5) )

folder_name <- "RECLUSTERING/PERICYCLE_v2/"

cluster_number = c(  
    17, 19, 5, 18, 13, 7, 16, 14
)

cells_on_clusters <- cds@clusters$UMAP$clusters
cells_on_clusters <- cells_on_clusters[cells_on_clusters %in% cluster_number]
cells_on_clusters <- names(cells_on_clusters)

cds_subset <- cds[, colnames(cds) %in% cells_on_clusters ]

( p1 <- plot_cells(cds_subset,
                   cell_size = 0.75,
                   group_label_size = 8) ) + 
    facet_wrap(~ timepoint, nrow = 1)

system( paste0("mkdir -p ", folder_name, "rds_file_subset") )
saveRDS(cds_subset,
        paste0(folder_name,
               "rds_file_subset",
               "/medicago_integrated_subset_pericycle.rds") )