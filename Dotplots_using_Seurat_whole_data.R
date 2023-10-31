require(Seurat)
require(vroom)
require(ggplot2)

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

generate_dotplots <- function(file_name = file_name,
                              c_nickname = c_nickname,
                              folder_n = folder_n) {
    
    cds <- readRDS( paste0("rds_files/", file_name, ".rds") )
    
    cds@meta.data$clusters <- paste0(c_nickname, "-",
                                     cds@meta.data$seurat_clusters)
    
    system( paste0("mkdir -p images/dotplots/", folder_n) )
    for( i in unique(markers$Group) ) {
        
        markers2 <- markers[markers$Group == i, ]
        markers2 <- markers2[markers2$gene_id %in% rownames(cds), ]
        
        ( dotp <- DotPlot(cds,
                          dot.scale = 7,
                          features = markers2$gene_id,
                          cluster.idents = T,
                          cols = c("Darkblue", "red"),
                          group.by = "clusters",
                          assay = "RNA") +
                scale_x_discrete(position = "top",
                                 labels = markers2$acronym) +
                xlab("")+
                ylab("") +
                ggtitle("") +
                theme( axis.text.y = element_text(size = rel(1)),
                       axis.text.x = element_text(size = rel(1),vjust = 0.5,
                                                  angle = 90))+
                coord_flip() )
        
        ggsave(filename = paste0("images/dotplots/",folder_n,"/",
                                 file_name, "--",
                                 i, ".svg"),
               plot = dotp,
               height = round(3 + round(3 + nrow(markers2)*0.4)),
               width = 25, units = "cm", bg = "white")
        
        ggsave(filename = paste0("images/dotplots/",folder_n,"/",
                                 file_name, "--",
                                 i, ".png"),
               plot = dotp,
               height = round(3 + round(3 + nrow(markers2)*0.4)),
               width = 25, units = "cm", bg = "white")
    }
}

generate_heatmaps <- function(file_name = file_name,
                              c_nickname = c_nickname,
                              folder_n = folder_n) {
    
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
    
    cds <- readRDS( paste0("rds_files/", file_name, ".rds") )
    
    colData(cds)$clusters <- paste0( c_nickname, "-",
                                     monocle3::clusters(cds) )
    
    system( paste0("mkdir -p images/heatmaps/", folder_n) )
    for( i in unique(markers$Group) ) {
        
        markers2 <- markers[markers$Group == i, ]
        
        p1 <- my_pheatmap(marker_ids = markers2)
        
        my_ggplot(p1, 
                  paste0("images/heatmaps/", folder_n, "/",
                         file_name, "--",
                         i),
                  "svg",
                  round(3 + round(3 + nrow(markers2)*0.4)) )
        
        my_ggplot(p1, paste0("images/heatmaps/", folder_n, "/",
                             file_name, "--",
                             i),
                  "png",
                  round(3 + nrow(markers2)*0.4))
    }
}

### Schiessl_2019
markers <- suppressMessages(
    vroom("selected_markers/markers_Schiessl_2019_table_s2_including_annot.tst",
          delim = "\t",
          col_names = T,
          na = "NA") )

generate_dotplots(file_name = "seurat_formated_whole_dataset",
                  c_nickname = "C",
                  folder_n = "whole_dataset_Schiessl")

generate_heatmaps(file_name = "batched_integrated_clustered_complete_dataset",
                  c_nickname = "C",
                  folder_n = "whole_dataset_Schiessl")

### Roy_2020
markers <- suppressMessages(
    vroom("selected_markers/figure_2_roy_et_al_2020.csv",
          delim = ",",
          col_names = T,
          na = "NA") ) %>%
    dplyr::rename( Group = Category) %>%
    dplyr::rename(acronym = "Gene Symbol")

generate_dotplots(file_name = "seurat_formated_whole_dataset",
                  c_nickname = "C",
                  folder_n = "whole_dataset_Roy")

generate_heatmaps(file_name = "batched_integrated_clustered_complete_dataset",
                  c_nickname = "C",
                  folder_n = "whole_dataset_Roy")
