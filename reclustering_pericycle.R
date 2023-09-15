set.seed(1407)
suppressMessages( require(monocle3) )
suppressMessages( require(vroom) )
suppressMessages( require(cowplot) )
suppressMessages( require(ggplot2) )
suppressMessages( require(dplyr) )
suppressMessages( require(stringr) )

#####################################
## Selecting cells from pericycle. ##
#####################################
cds <- readRDS("rds_files/batched_integrated_clustered_complete_dataset.rds")

( p0 <- plot_cells(cds,
                   color_cells_by="cluster",
                   cell_size = 0.5,
                   group_label_size = 5) )

# Where to save the results
folder_name <- "RECLUSTERING/PERICYCLE/"

cluster_number = c(  
    8  # Pericycle 
)

cells_on_clusters <- cds@clusters$UMAP$clusters
cells_on_clusters <- cells_on_clusters[cells_on_clusters %in% cluster_number]
cells_on_clusters <- names(cells_on_clusters)

cds_subset <- cds[, colnames(cds) %in% cells_on_clusters ]

( p1 <- plot_cells(cds_subset,
                   cell_size = 0.75,
                   group_label_size = 8) ) + 
    facet_wrap(~ Group + timepoint, nrow = 2)

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

## Differential expression within C8 over time

## Fits the regression
gene_fits <- fit_models(cds_subset,
                        model_formula_str = "~timepoint + 0",
                        cores = 7)

fit_coefs <- coefficient_table(gene_fits)

RH_terms_sig <- fit_coefs %>%
    filter(q_value < 0.05) %>%
    filter(status == "OK") %>%
    group_by(term) %>%
    select(gene_short_name, term, q_value, estimate) %>%
    arrange( desc( term ), desc( estimate ) )

write.table(RH_terms_sig, 
            "DEGs_linear_regression_pericycle_C8_over_time_point.tsv",
            col.names = T,
            row.names = T,
            sep = "\t",
            quote = T)

# Reclustering
cds_subset <- clear_cds_slots(cds_subset)
cds_subset <- monocle3::preprocess_cds(cds_subset, num_dim = 100, verbose = T)

cds_subset <- reduce_dimension(cds_subset, verbose = T)

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

folder_name <- paste0(folder_name, "top_1000_markers_per_cluster/", "by_", filtering_criteria )
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
