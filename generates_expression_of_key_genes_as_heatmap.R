suppressMessages( require(ComplexHeatmap) )
suppressMessages( require(ggplot2) )
suppressMessages( require(dplyr) )
suppressMessages( require(RColorBrewer) )
suppressMessages( require(circlize) )
suppressMessages( require(monocle3) )
suppressMessages( require(tidyverse) )
suppressMessages( require(vroom) )
suppressMessages( require(viridis) )
suppressMessages( require(docopt) )

set.seed(1407)

my_pheatmap <- function(marker_ids = markers2) {
    
    cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                    cell_group=clusters(cds)[colnames(cds)])
    
    agg_mat <- aggregate_gene_expression(cds, NULL, cell_group_df)
    agg_mat2 <- agg_mat[rownames(agg_mat) %in% marker_ids$gene_id, ]
    
    marker_ids$Gene_ID <- as.character(marker_ids$gene_id)
    marker_ids$synonym <- as.character(marker_ids$acronym)
    
    for(i in 1:nrow(agg_mat2) ) {
        
        rownames(agg_mat2)[i] <- ifelse( rownames(agg_mat2)[i] %in% marker_ids$gene_id,
                                         unique(marker_ids[marker_ids$gene_id %in% rownames(agg_mat2)[i], "acronym"]),
                                         marker_ids$gene_id)
    }
    
    pheatmap::pheatmap(agg_mat2,
                       show_rownames = T,
                       cluster_rows=T,
                       cluster_cols=T,
                       scale="none",
                       color = viridis(20),
                       clustering_method="ward.D2",
                       fontsize=12)
}

my_ggplot <- function(plot_name = plot_name,
                      file_name = file_name,
                      format_to_save = format_to_save,
                      height_n = height_n) { 
    
    ggplot2::ggsave( paste0(file_name, ".", format_to_save),
                     p1,
                     height = height_n,
                     width = 20,
                     units = "cm",
                     dpi = 300,
                     bg = "#FFFFFF")
    
}

cds <- readRDS("rds_files/batched_integrated_clustered_complete_dataset.rds")

( p0 <- plot_cells(cds,
                   color_cells_by="cluster",
                   cell_size = 0.6,
                   group_label_size = 9) )

###################
## Schiessl_2019 ##
###################
system( 'mkdir -p images/heatmaps_Schiessl_2019/' )
markers <- suppressMessages(
    vroom("selected_markers/markers_Schiessl_2019_table_s2_including_annot.tst",
          delim = "\t",
          col_names = T,
          na = "NA") )

markers2 <- markers[markers$Group == "nodule specific", ]

p1 <- my_pheatmap(marker_ids = markers2)

my_ggplot(p1, "images/heatmaps_Schiessl_2019/nodule_specific",
          "svg", round(3 + round(3 + nrow(markers2)*0.4)) )
my_ggplot(p1, "images/heatmaps_Schiessl_2019/nodule_specific",
          "png", round(3 + nrow(markers2)*0.4))

markers2 <- markers[markers$Group == "epidermal infection", ]

p1 <- my_pheatmap(marker_ids = markers2)

my_ggplot(p1, "images/heatmaps_Schiessl_2019/epidermal_infection",
          "svg", round(3 + nrow(markers2)*0.4))
my_ggplot(p1, "images/heatmaps_Schiessl_2019/epidermal_infection",
          "png", round(3 + nrow(markers2)*0.4))

markers2 <- markers[markers$Group == "Lateral root specific", ]

p1 <- my_pheatmap(marker_ids = markers2)

my_ggplot(p1, "images/heatmaps_Schiessl_2019/Lateral_root_specific",
          "svg", round(3 + nrow(markers2)*0.4))
my_ggplot(p1, "images/heatmaps_Schiessl_2019/Lateral_root_specific",
          "png", round(3 + nrow(markers2)*0.4))

markers2 <- markers[markers$Group == "auxin signaling", ]

p1 <- my_pheatmap(marker_ids = markers2)

my_ggplot(p1, "images/heatmaps_Schiessl_2019/auxin_signaling",
          "svg", (round(3 + nrow(markers2)*0.4)) )
my_ggplot(p1, "images/heatmaps_Schiessl_2019/auxin_signaling",
          "png", round(3 + nrow(markers2)*0.4))

markers2 <- markers[markers$Group == "cytokinin signaling", ]

p1 <- my_pheatmap(marker_ids = markers2)

my_ggplot(p1, "images/heatmaps_Schiessl_2019/cytokinin_signaling", 
          "svg", round(3 + nrow(markers2)*0.4))
my_ggplot(p1, "images/heatmaps_Schiessl_2019/cytokinin_signaling",
          "png", round(3 + nrow(markers2)*0.4))

######################
## Roy et al., 2020 ##
######################
system( 'mkdir -p images/heatmaps_roy_et_al_2020/' )
markers <- suppressMessages(
    vroom("selected_markers/figure_2_roy_et_al_2020_including_annot.tsv",
          delim = "\t",
          col_names = T,
          na = "NA") )

markers2 <- markers[markers$Category == "Early signaling", ]

p1 <- my_pheatmap(marker_ids = markers2)

my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Early_signaling",
          "svg", round(3 + nrow(markers2)*0.4))
my_ggplot(p1, "images/heatmaps_roy_et_al_2020/nodule_specific",
          "png", round(3 + nrow(markers2)*0.4))

markers2 <- markers[markers$Category == "Rhizobial Infection", ]

p1 <- my_pheatmap(marker_ids = markers2)

my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Rhizobial_Infection",
          "svg", round(3 + nrow(markers2)*0.4))
my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Rhizobial_Infection",
          "png", round(3 + nrow(markers2)*0.4))

markers2 <- markers[markers$Category == "Nodule Organogenesis", ]

p1 <- my_pheatmap(marker_ids = markers2)

my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Nodule_Organogenesis",
          "svg", round(3 + nrow(markers2)*0.4))
my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Nodule_Organogenesis",
          "png", round(3 + nrow(markers2)*0.4))

markers2 <- markers[markers$Category == "Symbiosome formation", ]

p1 <- my_pheatmap(marker_ids = markers2)

my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Symbiosome_formation",
          "svg", round(3 + nrow(markers2)*0.4))
my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Symbiosome_formation",
          "png", round(3 + nrow(markers2)*0.4))

markers2 <- markers[markers$Category == "Nodule Metabolism and Transport", ]

p1 <- my_pheatmap(marker_ids = markers2)

my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Nodule_Metabolism_and_Transport",
          "svg", round(3 + nrow(markers2)*0.4))
my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Nodule_Metabolism_and_Transport",
          "png", round(3 + nrow(markers2)*0.4))

markers2 <- markers[markers$Category == "Senescence", ]

p1 <- my_pheatmap(marker_ids = markers2)

my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Senescence", "svg", round(3 + nrow(markers2)*0.4))
my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Senescence", "png", round(3 + nrow(markers2)*0.4))

markers2 <- markers[markers$Category == "Defense", ]

p1 <- my_pheatmap(marker_ids = markers2)

my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Defense", "svg", round(3 + nrow(markers2)*0.4))
my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Defense", "png", round(3 + nrow(markers2)*0.4))

markers2 <- markers[markers$Category == "Host Range Restriction", ]

p1 <- my_pheatmap(marker_ids = markers2)

my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Host_Range_Restriction", "svg", round(3 + nrow(markers2)*0.4))
my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Host_Range_Restriction", "png", round(3 + nrow(markers2)*0.4))

markers2 <- markers[markers$Category == "Bacterial Maturation", ]

p1 <- my_pheatmap(marker_ids = markers2)

my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Bacterial_Maturation", "svg", round(3 + nrow(markers2)*0.4))
my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Bacterial_Maturation", "png", round(3 + nrow(markers2)*0.4))

markers2 <- markers[markers$Category == "Autoregulation of Nodule Number", ]

p1 <- my_pheatmap(marker_ids = markers2)

my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Autoregulation_of_Nodule_Number", "svg", round(3 + nrow(markers2)*0.4))
my_ggplot(p1, "images/heatmaps_roy_et_al_2020/Autoregulation_of_Nodule_Number", "png", round(3 + nrow(markers2)*0.4))