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
    
    agg_mat2 <- agg_mat2[ order( match( rownames(agg_mat2),  marker_ids$gene_id)), ]    
    marker_ids$Gene_ID <- as.character(marker_ids$gene_id)
    marker_ids$synonym <- as.character(marker_ids$acronym)
    
    for ( i in 1:nrow(agg_mat2) ) {
        
        rownames(agg_mat2)[i] <- ifelse( rownames(agg_mat2)[i] %in% marker_ids$gene_id,
                                         unique(marker_ids[marker_ids$gene_id %in% rownames(agg_mat2)[i], "acronym"]),
                                         marker_ids$gene_id)
    }
    
    if ( is.null(c_order) == FALSE) {
        agg_mat2 <- agg_mat2[, c_order]
    }
    
    NCRs <- agg_mat2[!rowSums(agg_mat2) == 0, ]
    
    pheatmap::pheatmap(NCRs,
                       show_rownames = F,
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

## Heatmap

markers <- vroom("selected_markers/all_NCRs_medicago.tsv")

p1 <- my_pheatmap( marker_ids = markers,
                   c_order = NULL)

ggplot2::ggsave( "images/all_NCRs_medicago.png",
                 p1,
                 height = 30,
                 width = 20,
                 units = "cm",
                 dpi = 300,
                 bg = "#FFFFFF")