set.seed(1407)
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
suppressMessages( require(vroom) )
suppressMessages( require(tidyverse) )
suppressMessages( require(monocle3) )

my_pheatmap <- function(marker_ids = markers,
                        c_order = NULL) {
    
    cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                    cell_group=clusters(cds)[colnames(cds)])
    
    agg_mat <- aggregate_gene_expression(cds, NULL, cell_group_df)
    agg_mat2 <- agg_mat[rownames(agg_mat) %in% marker_ids$gene_id, ]
    
    # To keep the order of the input genes
    agg_mat2 <- agg_mat2[ order( match( rownames(agg_mat2),  marker_ids$gene_id)), ]    
    marker_ids$Gene_ID <- as.character(marker_ids$gene_id)
    marker_ids$synonym <- as.character(marker_ids$acronym)
    
    for( i in 1:nrow(agg_mat2) ) {
        
        rownames(agg_mat2)[i] <- ifelse( rownames(agg_mat2)[i] %in% marker_ids$gene_id,
                                         unique(marker_ids[marker_ids$gene_id %in% rownames(agg_mat2)[i], "acronym"]),
                                         marker_ids$gene_id)
    }
    
    if( is.null(c_order) == FALSE) {
                agg_mat2 <- agg_mat2[, c_order]
    }
    
    pheatmap::pheatmap(agg_mat2,
                       show_rownames = T,
                       cluster_rows=F,
                       cluster_cols=F,
                       scale="none",
                       color = viridis(10, option = "D"),
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

markers <- suppressMessages( 
    vroom("selected_markers/including_annot_all_markers_of_interest_selected.tsv",
          delim = "\t",
          col_names = T,
          na = "NA") )

( p1 <- my_pheatmap( marker_ids = markers,
                     c_order = c(
                         1, 3, 14, 22, # Epidermis / Root hair
                         12, # Root hair
                         15, 29, # Lateral root
                         26, # Xylem
                         10, 23, # Stele
                         18, # Cell division vasculature
                         8, # Pericycle
                         7, 24, # Endodermis
                         2, 4, 9, 11, 13, 16, # Cortex
                         6, # Nodule
                         5, 17, 19, 20, 21, 25, 27, 28 # NA
                     ) ) )

ggplot2::ggsave("images/heatmap_of_selected_markers_for_cluster_annotation.png",
                p1,
                height = 30,
                width = 40,
                units = "cm",
                dpi = 300,
                bg = "#FFFFFF")

ggplot2::ggsave("images/heatmap_of_selected_markers_for_cluster_annotation.svg",
                p1,
                height = 30,
                width = 40,
                units = "cm",
                dpi = 300,
                bg = "#FFFFFF")