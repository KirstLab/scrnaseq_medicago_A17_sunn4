set.seed(1407)

outfile <- "logs/umap_colored_by_cell_type.out" # File name of output log
#Check its existence
if ( file.exists(outfile) ) {
    #Delete file if it exists
    file.remove(outfile)
}

my_log <- file(outfile)
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

require(monocle3)
require(ggplot2)

cds <- readRDS("rds_files/batched_integrated_clustered_complete_dataset.rds")

colData(cds)$clusters <- monocle3::clusters(cds)

meta_info <- as.data.frame( colData(cds) )

epi_root <- c(1, 3, 14, 22) # Epidermis/Root hair
root_hair <- c(12) # Root hair
lateral_root <- c(15, 29) # Lateral Root
stele <- c(26, 10, 23, 18) # Stele
pericycle <- c(8) # Pericycle 
endodermis <- c(7, 24) # Endodermis
cortex <- c(2,  4, 9, 11, 13, 16) # Cortex
nodule <- c(6) # Nodule related
na <- c(5, 17, 19, 20, 21, 25, 27, 28 ) # NA

colData(cds)$annot <- ifelse (
    colData(cds)$clusters %in% epi_root, "Epidermis/Root hair", 
    ifelse (
        colData(cds)$clusters %in% cortex, "Cortex",
        ifelse (
            colData(cds)$clusters %in% nodule, "Root Nodule Symbiosis", 
            ifelse (
                colData(cds)$clusters %in% endodermis, "Endodermis", 
                ifelse (
                    colData(cds)$clusters %in% stele, "Stele", 
                    ifelse (
                        colData(cds)$clusters %in% pericycle, "Pericycle",
                        ifelse (
                            colData(cds)$clusters %in% root_hair, "Root hair",
                            ifelse (
                                colData(cds)$clusters %in% lateral_root, "Lateral root",
                                ifelse (
                                    colData(cds)$clusters %in% na, "Unidentified", ""
) ) ) ) ) ) ) ) ) 

# to generate transparent colors
t_col <- function(color, percent = 50, name = NULL) {

    rgb.val <- col2rgb(color)
    
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)
    
    invisible(t.col)
}

annot_names <- c("Epidermis/Root hair",
                 "Cortex",
                 "Root Nodule Symbiosis",
                 "Endodermis",
                 "Stele",
                 "Pericycle",
                 "Root hair",
                 "Lateral root",
                 "Unidentified"
)

Tol_muted <- c(
    t_col( '#88CCEE' ),
    t_col( '#882255' ),
    t_col( '#44AA99' ),
    t_col( '#332288' ),
    t_col( '#117733' ),
    t_col( '#AA4499' ),
    t_col( '#999933' ),
    t_col( '#CC6677' ),
    t_col( "gray") 
)

( p1 <- plot_cells(cds,
                   group_cells_by = "cluster",
                   color_cells_by = "annot",
                   label_cell_groups = F,
                   cell_size = 0.65)  +
        theme(legend.title = element_blank(),
              legend.text = element_text(size = 14) ) +
        theme(strip.text = element_text(size = 14),
              plot.title = element_text(size=20,
                                        face = "bold") ) + 
        scale_colour_manual(breaks = annot_names,
                            values = Tol_muted) )

ggplot2::ggsave(filename = "images/clustered_dataset_colored_by_annotation.svg",
                p1,
                height=10,
                width=18,
                units="cm",
                dpi = 300,
                bg = "#FFFFFF")

( p1_without_leg <- plot_cells(cds,
                               group_cells_by = "cluster",
                               color_cells_by = "annot",
                               label_cell_groups = F,
                               cell_size = 0.65)  +
        theme(legend.title = element_blank(),
              legend.text = element_text(size = 14) ) +
        theme(strip.text = element_text(size = 14),
              plot.title = element_text(size=20,
                                        face = "bold"),
              legend.position="none") + 
        scale_colour_manual(breaks = annot_names,
                            values = Tol_muted) )

ggplot2::ggsave(filename = "images/clustered_dataset_colored_by_annotation_without_legend_fixed.svg",
                p1_without_leg,
                height=10,
                width=15,
                units="cm",
                dpi = 300,
                bg = "#FFFFFF")

( p2 <- plot_cells(cds,
                   group_cells_by = "cluster",
                   cell_size = 0.4,
                   group_label_size = 6)  +
        theme(legend.title = element_blank(),
              legend.text = element_text(size = 14) ) +
        theme(strip.text = element_text(size = 14),
              plot.title = element_text(size=20,
                                        face = "bold"),
              legend.position="none") )

ggplot2::ggsave(filename = "images/clustered_dataset_clusters_for_figure1.svg",
                p2,
                height=10,
                width=15,
                units="cm",
                bg = "#FFFFFF")

closeAllConnections()
