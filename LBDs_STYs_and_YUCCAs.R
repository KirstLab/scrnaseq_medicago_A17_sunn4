set.seed(1407)
suppressMessages( require(docopt) )
suppressMessages( require(monocle3) )
suppressMessages( require(vroom) )
suppressMessages( require(cowplot) )
suppressMessages( require(ggplot2) )
suppressMessages( require(dplyr) )

cds <- readRDS("rds_files/batched_integrated_clustered_complete_dataset.rds")

all_genes <- c(
    "MtrunA17Chr1g0207471", # MtLBD1
    "MtrunA17Chr4g0043421", # MtLBD4
    "MtrunA17Chr4g0030001", # MtLBD11
    "MtrunA17Chr7g0260971", # MtLBD16
    "MtrunA17Chr8g0351701", # MtLBD18
    
    "MtrunA17Chr8g0372461", # MtSTY2
    "MtrunA17Chr5g0404781", # MtSTY3 
    "MtrunA17Chr3g0082511", # MtSTY4
    "MtrunA17Chr3g0142171", # MtSTY5
    "MtrunA17Chr4g0035591", # MtSTY6
    "MtrunA17Chr5g0441921", # MtSTY7
    "MtrunA17Chr1g0155791", # MtSTY8
    "MtrunA17Chr8g0353111", # MtSTY9
    
    "MtrunA17Chr6g0485621", # MtYUC2
    "MtrunA17Chr1g0169111", # MtYUC3
    "MtrunA17Chr7g0262591", # MtYUC8
    "MtrunA17Chr1g0182991", # MtYUC9
    "MtrunA17Chr5g0448621", # MtNIN,
    "MtrunA17Chr8g0392301"  # MtCRE1
)

common_n <- c(
    "MtLBD1",
    "MtLBD4",
    "MtLBD11",
    "MtLBD16",
    "MtLBD18",
    
    "MtSTY2",
    "MtSTY3",
    "MtSTY4",
    "MtSTY5",
    "MtSTY6",
    "MtSTY7",
    "MtSTY8",
    "MtSTY9",
    
    "MtYUC2",
    "MtYUC3",
    "MtYUC8",
    "MtYUC9",
    "MtNIN",
    "MtCRE1"
)

all_genes_p <- plot_cells(cds,
                          genes = all_genes,
                          cell_size = 1,
                          scale_to_range = F,
                          label_leaves = T,
                          label_cell_groups = F,
                          group_label_size = 5)  +
    theme(strip.text = element_text(size = 16),
          text = element_text(size = 16),
          plot.title = element_text(size=14,
                                    face = "bold") ) +
    viridis::scale_color_viridis(option = "A") 

system( paste0( "mkdir -p LBDs_STYs_YUCCAs") )
ggplot2::ggsave(filename ="LBDs_STYs_YUCCAs/all_genes_UMAP_combined_timepoints_v2.svg",
                all_genes_p,
                height=48,
                width=80,
                units="cm",
                bg = "#FFFFFF")

ggplot2::ggsave(filename ="LBDs_STYs_YUCCAs/all_genes_UMAP_combined_timepoints_v2.png",
                all_genes_p,
                height=48,
                width=80,
                units="cm",
                bg = "#FFFFFF",
                dpi = 300)

for( i in 1:length(all_genes) ) {
    
    ( p_by_time <- plot_cells(cds,
                              genes = all_genes[i],
                              labels_per_group = F,
                              cell_size = 1,
                              scale_to_range = F) +
          facet_wrap(~timepoint, nrow = 1) + 
          theme(strip.text = element_text(size = 20),
                plot.title = element_text(size=20,
                                          face = "bold") ) +
          viridis::scale_color_viridis(option = "A") )
    
    ggplot2::ggsave( paste0(
        "LBDs_STYs_YUCCAs/",
        common_n[i],
        "__all_genes_UMAP_combined_timepoints_v2.svg"),
        p_by_time,
        height=12,
        width=60,
        units="cm",
        bg = "#FFFFFF")
    
    ggplot2::ggsave( paste0(
        "LBDs_STYs_YUCCAs/",
        common_n[i],
        "__all_genes_UMAP_combined_timepoints_v2.png"),
        p_by_time,
        height=48,
        width=80,
        units="cm",
        bg = "#FFFFFF",
        dpi = 300)
}
