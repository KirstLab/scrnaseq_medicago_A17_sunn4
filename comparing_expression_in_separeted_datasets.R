set.seed(1407)

suppressMessages( require(docopt) )
suppressMessages( require(vroom) )
suppressMessages( require(tidyverse) )
suppressMessages( require(monocle3) )
suppressMessages( library(ggplot2) )
suppressMessages( library(cowplot) )

a17_cds <- readRDS("rds_files/batched_integrated_clustered_complete_dataset_A17.rds")
sunn4_cds <- readRDS("rds_files/batched_integrated_clustered_complete_dataset_SUNN4.rds")

markers <- vroom("selected_markers/selected_markers_for_FigureS1.csv",
                 delim = ",")

## Help functions
my_UMAP <- function(gene_name = gene_name,
                    common_name = common_name,
                    celltype = celltype) {
    
    pa17 <- plot_cells(a17_cds,
                       genes = gene_name,
                       group_label_size = 0,
                       cell_size = 0.8,
                       scale_to_range = F) + 
        theme(strip.text = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              plot.title = element_blank(),
              legend.title=element_blank(),
              line = element_blank() ) +
        viridis::scale_color_viridis(option = "A") + 
        theme(legend.key.size = unit(1.0, 'cm')) 
    
    psunn4 <- plot_cells(sunn4_cds,
                         genes = gene_name,
                         group_label_size = 0,
                         cell_size = 0.8,
                         scale_to_range = F)  + 
        theme(strip.text = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              plot.title = element_blank(),
              legend.title=element_blank(),
              line = element_blank() ) +
        viridis::scale_color_viridis(option = "A") + 
        theme(legend.key.size = unit(1.0, 'cm')) 
    
    title <- ggdraw() + 
        draw_label(common_name,
            fontface = 'bold',
            size = 30,
            x = 0,
            hjust = 0) +
        theme( plot.margin = margin(0, 0, 0, 7) )
    
    p_combined <- cowplot::plot_grid(pa17,
                                     NULL,
                                     psunn4,
                                     labels = c("A17 (WT)","", "sunn-4"),
                                     label_size = 24,
                                     hjust = c(-5, 0, -2.5),
                                     nrow = 1,
                                     ncol = 3,
                                     rel_widths = c(1, 0.15, 1)) 
    
    p_combined2 <- cowplot::plot_grid(
        title, p_combined,
        ncol = 1,
        rel_heights = c(0.1, 1)
    ) 
    
    ggplot2::ggsave(filename = paste0("comparing_expression_datasets/",
                                      gene_name,
                                      "-",
                                      common_name,
                                      "-",
                                      celltype,
                                      ".svg"),
                    p_combined2,
                    height=15,
                    width=34,
                    units="cm",
                    bg = "#FFFFFF",
                    dpi = 300)
    
    p_combined2
}

my_UMAP2 <- function() {
    
    pa17 <- plot_cells(a17_cds,
                       group_label_size = 6,
                       cell_size = 0.5,
                       scale_to_range = F) + 
        theme(strip.text = element_blank(),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_blank(),
              legend.title=element_blank()) +
        theme(legend.key.size = unit(1.0, 'cm')) 
    
    psunn4 <- plot_cells(sunn4_cds,
                         group_label_size = 6,
                         cell_size = 0.5,
                         scale_to_range = F) + 
        theme(strip.text = element_blank(),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_blank(),
              legend.title=element_blank()) +
        theme(legend.key.size = unit(1.0, 'cm'))   
    
    p_combined <- cowplot::plot_grid(pa17,
                                     psunn4,
                                     labels = c("A17 (WT)", "sunn-4"),
                                     label_size = 24,
                                     hjust = c(-3, -3),
                                     nrow = 1,
                                     ncol = 2,
                                     rel_widths = c(1, 1)) 
    
    ggplot2::ggsave(
        filename = paste0("comparing_expression_datasets/A17_sunn-4_comparison.svg"),
        p_combined,
        height=12,
        width=32,
        units="cm",
        bg = "#FFFFFF",
        dpi = 300)
    
    p_combined
} 

system("mkdir -p comparing_expression_datasets/")

for ( i in 1:nrow( markers ) ) {
    
    ( combined_UMAP <- my_UMAP( gene_name = markers[i, 1],
                                common_name = markers[i, 2],
                                celltype = markers[i, 3]) )
    
}

my_UMAP2()