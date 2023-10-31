set.seed(1407)

outfile <- "logs/generates_the_heatmap_containing_gene_markers_for_figure1.out"

if ( file.exists(outfile) ) {
    file.remove(outfile)
}

my_log <- file(outfile)
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

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
suppressMessages( require(Seurat) )

my_pheatmap <- function(marker_ids = markers,
                        c_order = NULL) {
    
    cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                    cell_group=clusters(cds)[colnames(cds)])
    
    agg_mat <- aggregate_gene_expression(cds, NULL, cell_group_df)
    agg_mat2 <- agg_mat[rownames(agg_mat) %in% marker_ids$gene_id, ]
    
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

seurat_obj <- readRDS("rds_files/seurat_formated_whole_dataset.rds")

markers_all <- suppressMessages( 
    vroom("selected_markers/including_annot_all_markers_of_interest_selected.tsv",
          delim = "\t",
          col_names = T,
          na = "NA") )

order_levels = c(
    2, 4, # Epidermis
    11, # Root hair
    21,22, # lateral root
    17, 19, #   stele xylem
    5, 18, # Stele
    13, # Stele Phloem
    7, # Pericycle
    16, #vascular bundle
    14, # Stele cell division
    8, #Endodermis
    10, #Endodermis (suberized)
    1,3, #Cortex
    6,9,15, # nodule
    12,20,23,24
)

seurat_obj@meta.data$seurat_clusters2 <- base::factor(x = seurat_obj@meta.data$seurat_clusters, levels = order_levels )

( dotp <- Seurat::DotPlot(seurat_obj,
                          dot.scale = 6,
                          features = rev(markers_all$gene_id),
                          cluster.idents = F,
                          cols = c("Darkblue", "red"),
                          group.by = "seurat_clusters2",
                          assay = "RNA") +
        scale_x_discrete(position = "top") +
        xlab("")+
        ylab("") +
        ggtitle("") +
        theme( axis.text.y = element_text(size = rel(1)) ) +
        coord_flip() )

( p1 <- my_pheatmap( marker_ids = markers_all,
                     c_order = order_levels) )

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
closeAllConnections()