set.seed(1407)

outfile <- "logs/trajectory_inference_cortex_nodule.out"
if ( file.exists(outfile) ) {
    file.remove(outfile)
}

my_log <- file(outfile)
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

require(ggplot2)
require(monocle3)
require(tidyverse)
require(paletteer)

## Help functions
theme_umap <- function(base.size = 14) {
    ggplot2::theme_classic(base_size = base.size) + 
        ggplot2::theme(axis.ticks = ggplot2::element_blank(), 
                       axis.text = ggplot2::element_blank(), 
                       plot.subtitle = ggplot2::element_text(face = "italic", size = 11), 
                       plot.caption = ggplot2::element_text(face = "italic", size = 11))
}
guide_umap <- function(key.size = 4) {
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = key.size, alpha = 1)))
}

## Help functions
my_ggsave <- function( name = name,
                       plot = plot,
                       height = 12,
                       width = 16) {
    
    ggplot2::ggsave(
        filename = name,
        plot = plot,
        height = height,
        width = width,
        units = "in",
        bg = "#FFFFFF", 
        dpi = 300)
    
}

my_write_csv <- function( my_obj = my_obj,
                          name = name,
                          col.names = T,
                          row.names = F,
                          quote = F,
                          sep = "\t") {
    
    write.table(my_obj,
                name,
                col.names = T,
                row.names = F,
                quote = F,
                sep = "\t")
}

palette_cluster <- paletteer::paletteer_d("ggsci::default_jama")
palette_celltype <- paletteer::paletteer_d("ggsci::category20_d3")
palette_heatmap <- paletteer::paletteer_d("wesanderson::Zissou1")

cds <- readRDS("RECLUSTERING/CORTEX_NOD/rds_file_subset/medicago_integrated_subset_cortex_nodule.rds")

monocle3::plot_cells(cds,
                     label_cell_groups=T,
                     graph_label_size=1.5,
                     cell_size = 1,
                     group_label_size = 7) 

system("mkdir -p RECLUSTERING/CORTEX_NOD/images/slingshot")

UMAP <- SingleCellExperiment::reducedDims(cds)[["UMAP"]]

colData(cds)$cluster <- monocle3::clusters(cds)

monocle3::plot_cells(cds,
                     label_cell_groups=T,
                     graph_label_size=1.5,
                     cell_size = 1,
                     group_label_size = 7) + 
    facet_wrap(~timepoint)

start_clust <- as.data.frame( colData(cds) ) %>%
    group_by(cluster, timepoint) %>%
    summarise(n_cells = n()) %>%
    filter(timepoint == "0h")

start_clust <- start_clust[start_clust$n_cells == max(start_clust$n_cells ),]

### trajectory inference with slingshot, starting at cluster with more cells from T0h.
sling_res <- slingshot::slingshot( UMAP, 
                                   clusterLabels = colData(cds)$cluster, 
                                   start.clus = start_clust$cluster,  # select starting point
                                   approx_points = 1000,
                                   extend = "n",
                                   omega = TRUE
                                   )

slingshot::slingLineages(sling_res)

### Selecting clusters that are in the lineages of interest
cluster_number = c(1,6,10,9,7)

cells_on_clusters <- cds@clusters$UMAP$clusters
cells_on_clusters <- cells_on_clusters[cells_on_clusters %in% cluster_number]
cells_on_clusters <- names(cells_on_clusters)

cds_subset <- cds[, colnames(cds) %in% cells_on_clusters ]

plot_cells(cds_subset,
           cell_size = 0.75,
           group_label_size = 8)  + 
    facet_wrap(~timepoint)

cells_to_exclude <- vroom::vroom( "RECLUSTERING/CORTEX_NOD/cells_to_exclude.tsv",
                                  col_names = F,
                                  delim = "\t")

cds_subset <- cds_subset[, !colnames(cds_subset) %in% cells_to_exclude$X1 ]

( p1 <- plot_cells(cds_subset,
                   cell_size = 1,
                   group_label_size = 8) )

ggsave(filename = paste("RECLUSTERING/CORTEX_NOD/images/slingshot/selected_clusters_important_lineages.svg"),
       p1,
       width = 8,
       height = 5,
       dpi = 300,
       bg = "white")

( p2 <- plot_cells(cds_subset,
                   cell_size = 0.75,
                   group_label_size = 8)  + 
        facet_wrap(~timepoint) )

ggsave(filename = paste("RECLUSTERING/CORTEX_NOD/images/slingshot/selected_clusters_important_lineages_by_timepoint.svg"),
       p2,
       width = 10,
       height = 6,
       dpi = 300,
       bg = "white")

saveRDS(cds_subset, "RECLUSTERING/CORTEX_NOD/rds_file_subset/clusters_cortex_nodule_selected_for_trajectory.rds")

UMAP_sub <- SingleCellExperiment::reducedDims(cds_subset)[["UMAP"]]
colData(cds_subset)$cluster <- monocle3::clusters(cds_subset)

### trajectory inference with slingshot
sling_res <- slingshot::slingshot( UMAP_sub, 
                                   clusterLabels = colData(cds_subset)$cluster, 
                                   start.clus = start_clust$cluster,
                                   approx_points = 1000,
                                   extend = "n")

slingshot::slingLineages(sling_res)

sling_pt <- slingshot::slingPseudotime(sling_res) %>% 
    as.data.frame() 

col_n <- paste("PT", seq(1:ncol(sling_pt)), sep = "")
sling_pt <- sling_pt %>% 
    magrittr::set_colnames(col_n) 

for (c in 1:ncol(sling_pt)) {
    
    pt_lin <- paste0("PT", c)
    
    ( ps_plot <- UMAP_sub %>% 
            as.data.frame() %>% 
            magrittr::set_colnames(c("UMAP_1", "UMAP_2")) %>% 
            mutate(PT = sling_pt[[pt_lin]]) %>% 
            ggplot(aes(x = UMAP_1, y = UMAP_2, color = PT)) + 
            geom_point(size = 1, alpha = 0.75) + 
            labs(x = "UMAP 1", 
                 y = "UMAP 2", 
                 color = paste("Pseudotime lineage",c )) + 
            scale_color_gradientn(colors = palette_heatmap, 
                                  labels = scales::label_number(accuracy = 0.1),
                                  na.value="white") + 
            theme_umap() )
    
    ggsave(filename = paste("RECLUSTERING/CORTEX_NOD/images/slingshot/Pseudotime_lineage", c, "UMAP.png"),
           ps_plot, width = 12, height = 8, dpi = 300, bg = "white")
    
    ggsave(filename = paste("RECLUSTERING/CORTEX_NOD/images/slingshot/Pseudotime_lineage", c, "UMAP.svg"),
           ps_plot, width = 12, height = 8, dpi = 300, bg = "white")
    
}

saveRDS(sling_res,
        file = paste0("RECLUSTERING/CORTEX_NOD/rds_file_subset/cortex_nodules_slingshot_trajectory.rds")
)

closeAllConnections()