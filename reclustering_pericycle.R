set.seed(1407)

outfile <- "logs/reclustering_pericycle.out"
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

cds <- readRDS("rds_files/batched_integrated_clustered_complete_dataset.rds")

( p0 <- plot_cells(cds,
                   color_cells_by="cluster",
                   cell_size = 0.5,
                   group_label_size = 5) )

folder_name <- "RECLUSTERING/PERICYCLE/"

cluster_number = c(7)

cells_on_clusters <- cds@clusters$UMAP$clusters
cells_on_clusters <- cells_on_clusters[cells_on_clusters %in% cluster_number]
cells_on_clusters <- names(cells_on_clusters)

cds_subset <- cds[, colnames(cds) %in% cells_on_clusters ]

( p1 <- plot_cells(cds_subset,
                   cell_size = 0.75,
                   group_label_size = 8) ) + 
    facet_wrap(~ timepoint, nrow = 1)

system( paste0( "mkdir -p ", folder_name, "clustering_images") )
ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/pericycle_selected_clusters.svg"),
                p1,
                height=12,
                width=16,
                units="cm",
                bg = "#FFFFFF", 
                dpi = 300)

cds_subset <- clear_cds_slots(cds_subset)
cds_subset <- monocle3::preprocess_cds(cds_subset, num_dim = 100, verbose = T)

cds_subset <- reduce_dimension(cds_subset, verbose = T)

res <- c(1, 0.1, 0.01, 0.001, 0.0001, 0.00001)
cds_subset <- cluster_cells(cds_subset,
                            verbose = T,
                            resolution = 0.001,
                            random_seed = 1407)

( p2 <- plot_cells(cds_subset, 
                   group_label_size = 7,
                   cell_size = 0.75) )

( p3 <- plot_cells(cds_subset,
                   group_label_size = 7,
                   graph_label_size = 7,
                   cell_size = 0.5)  +
        facet_wrap(~Group + timepoint, nrow = 2) +
        theme(strip.text = element_blank(),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_text(size = 40,
                                        hjust = 0.5),
              legend.title=element_blank()) )

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/pericycle_selected_clusters_RECLUSTERED.svg"),
                p2,
                height=12,
                width=16,
                units="cm",
                bg = "#FFFFFF", 
                dpi = 300)

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/pericycle_selected_clusters_RECLUSTERED_by_time_and_genotype.png"),
                p3,
                height=25,
                width=40,
                units="cm",
                bg = "#FFFFFF",
                dpi = 300)

( p1.1 <- plot_cells(cds_subset,
                     cell_size = 0.75,
                     group_label_size = 8) +
        facet_wrap(~ timepoint, nrow = 1) +
        theme(strip.text = element_blank(),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_text(size = 40,
                                        hjust = 0.5),
              legend.title=element_blank()) )

ggplot2::ggsave(filename = paste0(folder_name,
                                  "clustering_images",
                                  "/pericycle_selected_clusters_by_timepoint.svg"),
                p1.1,
                height=12,
                width=48,
                units="cm",
                bg = "#FFFFFF", 
                dpi = 300)

system( paste0("mkdir -p ", folder_name, "rds_file_subset") )
saveRDS(cds_subset,
        paste0(folder_name,
               "rds_file_subset",
               "/medicago_integrated_subset_pericycle.rds") )

closeAllConnections()