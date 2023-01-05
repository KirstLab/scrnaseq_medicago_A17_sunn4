set.seed(1407)
suppressMessages( require(docopt) )
suppressMessages( require(monocle3) )
suppressMessages( require(vroom) )
suppressMessages( require(cowplot) )
suppressMessages( require(ggplot2) )
suppressMessages( require(dplyr) )

threads <- as.numeric(8)
umi_theshold <- as.numeric( 400 )

## A17 WT
cds_0h <- monocle3::load_cellranger_data("data/A17_sep_2022_0h_10k/",
                                         umi_cutoff = umi_theshold)
cds_24h <- monocle3::load_cellranger_data("data/A17_sep_2022_24h_10k/",
                                          umi_cutoff = umi_theshold)
cds_48h <- monocle3::load_cellranger_data("data/A17_sep_2022_48h_10k/",
                                          umi_cutoff = umi_theshold)
cds_96h <- monocle3::load_cellranger_data("data/A17_sep_2022_96h_10k/",
                                          umi_cutoff = umi_theshold)
## sunn-4
cds_sunn_0h <- monocle3::load_cellranger_data("data/Sunn_sep_2022_0h_10k/",
                                              umi_cutoff = umi_theshold)
cds_sunn_24h <- monocle3::load_cellranger_data("data/Sunn_sep_2022_24h_10k/",
                                               umi_cutoff = umi_theshold)
cds_sunn_48h <- monocle3::load_cellranger_data("data/Sunn_sep_2022_48h_10k/",
                                               umi_cutoff = umi_theshold)
cds_sunn_96h <- monocle3::load_cellranger_data("data/Sunn_sep_2022_96h_10k/",
                                               umi_cutoff = umi_theshold)

colData(cds_0h)$timepoint <- "0h"
colData(cds_24h)$timepoint <- "24h"
colData(cds_48h)$timepoint <- "48h"
colData(cds_96h)$timepoint <- "96h"
colData(cds_sunn_0h)$timepoint <- "0h"
colData(cds_sunn_24h)$timepoint <- "24h"
colData(cds_sunn_48h)$timepoint <- "48h"
colData(cds_sunn_96h)$timepoint <- "96h"

colData(cds_0h)$batch2 <- "A17_0h"
colData(cds_24h)$batch2 <- "A17_24h"
colData(cds_48h)$batch2 <- "A17_48h"
colData(cds_96h)$batch2 <- "A17_96h"
colData(cds_sunn_0h)$batch2 <- "Sunn4_0h"
colData(cds_sunn_24h)$batch2 <- "Sunn4_24h"
colData(cds_sunn_48h)$batch2 <- "Sunn4_48h"
colData(cds_sunn_96h)$batch2 <- "Sunn4_96h"

colData(cds_0h)$Group <- "A17-WT"
colData(cds_24h)$Group <- "A17-WT"
colData(cds_48h)$Group <- "A17-WT"
colData(cds_96h)$Group <- "A17-WT"
colData(cds_sunn_0h)$Group <- "Sunn-4"
colData(cds_sunn_24h)$Group <- "Sunn-4"
colData(cds_sunn_48h)$Group <- "Sunn-4"
colData(cds_sunn_96h)$Group <- "Sunn-4"

cds <- monocle3::combine_cds(list(cds_0h,
                                  cds_24h,
                                  cds_48h,
                                  cds_96h,
                                  cds_sunn_0h,
                                  cds_sunn_24h,
                                  cds_sunn_48h,
                                  cds_sunn_96h))

rm(cds_0h,
   cds_24h,
   cds_48h, 
   cds_96h,
   cds_sunn_0h,
   cds_sunn_24h,
   cds_sunn_48h,
   cds_sunn_96h)

################
## Clustering ##
################

cds <- monocle3::preprocess_cds(cds, num_dim = 100)
monocle3::plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds, verbose = T)

( p1 <- monocle3::plot_cells(cds,
                             color_cells_by="batch2",
                             label_cell_groups=FALSE) )

# removes possible batch effects
cds = align_cds(cds,
                num_dim = 100,
                alignment_group = "batch2", verbose = T)

cds <- reduce_dimension(cds, verbose = T)

( p2 <- plot_cells(cds,
                   color_cells_by="batch2",
                   label_cell_groups=FALSE) )

( p_time <- plot_cells(cds,
                   color_cells_by="timepoint",
                   label_cell_groups=FALSE) )

ggplot2::ggsave(filename = paste0("images/cells_distributed_by_timepoint_",
                                  umi_theshold, "_UMI.png"),
                p_time,
                height=10,
                width=15,
                units="cm",
                bg = "#FFFFFF")

p3 <- cowplot::plot_grid(p1, p2)

system("mkdir -p images/")
ggplot2::ggsave(filename = paste0("images/batch_effect_comparison_",
                                  umi_theshold, "_UMI.png"),
                p3,
                height=10,
                width=26,
                units="cm",
                bg = "#FFFFFF")

system("mkdir -p rds_files")
cds <- cluster_cells( cds, 
                              verbose = T,
                              random_seed = 1407
)

rowData(cds)$gene_short_name <- row.names(rowData(cds))
saveRDS( cds,
         paste0("rds_files/", "batched_integrated_clustered_complete_dataset.rds") )

( p5 <- plot_cells(cds, group_label_size = 11, cell_size = 0.75, alpha = 0.7) )

ggplot2::ggsave(filename = paste0("images/clustered_dataset_", umi_theshold, "_UMI.png"),
                p5,
                height=20,
                width=30,
                units="cm",
                bg = "#FFFFFF")

ggplot2::ggsave(filename = paste0("images/clustered_dataset_", umi_theshold, "_UMI.svg"),
                p5,
                height=20,
                width=30,
                units="cm",
                bg = "#FFFFFF")

( p8 <- plot_cells(cds,
                   group_label_size = 4,
                   cell_size = 0.5, 
                   alpha = 0.6) + 
        facet_wrap(~timepoint, nrow = 1) + 
        theme(strip.text = element_text(size = 20)) )

ggplot2::ggsave(filename = paste0("images/clusteres_by_time_",
                                  umi_theshold,
                                  "_UMI.png"),
                p8,
                height=10,
                width=50,
                units="cm",
                bg = "#FFFFFF")

( p9 <- plot_cells(cds,
                   group_label_size = 4,
                   cell_size = 0.5, 
                   alpha = 0.6) + 
        facet_wrap(~Group + timepoint, nrow = 2) + 
        theme(strip.text = element_text(size = 20)) )

ggplot2::ggsave(filename = paste0("images/clusteres_by_time_and_genotype",
                                  umi_theshold,
                                  "_UMI.png"),
                p9,
                height=22,
                width=50,
                units="cm",
                dpi = 300,
                bg = "#FFFFFF")

ggplot2::ggsave(filename = paste0("images/clusteres_by_time_and_genotype",
                                  umi_theshold,
                                  "_UMI.svg"),
                p9,
                height=22,
                width=50,
                units="cm",
                bg = "#FFFFFF")

################################################
## Identification of markers for each cluster ##
################################################

### Top 100
system("mkdir -p top_100_markers_per_cluster")
marker_test_res <- top_markers(cds, 
                               group_cells_by="cluster",
                               genes_to_test_per_group = 300,
                               marker_sig_test = T,
                               cores = threads,
                               verbose = T)

write.csv(marker_test_res, 
          "top_100_markers_per_cluster/top_300_markers_per_cluster.csv",
          row.names = F)

filtering_criteria = "specificity"

folder_name <- paste0("top_100_markers_per_cluster/", "by_", filtering_criteria )
system( paste0("mkdir -p ", folder_name) )
for (c in unique(marker_test_res$cell_group) ) {
    
    top_specific_markers <- marker_test_res %>%
        dplyr::filter(marker_test_q_value < 0.05) %>%
        dplyr::filter(fraction_expressing >= 0.20) %>%
        dplyr::filter(cell_group == c) %>%
        dplyr::arrange( desc(specificity) )
    
    top_specific_markers_top100 <- top_specific_markers[1:100, ]
    
    write.csv(top_specific_markers_top100, 
              paste0(folder_name,
                     "/list_of_top_100_for_cluster_",
                     c,
                     "_",filtering_criteria,
                     ".csv"),
              row.names = F)
    
    write.table(top_specific_markers_top100, 
                paste0(folder_name,
                       "/list_of_top_100_for_cluster_",
                       c,
                       "_",filtering_criteria,
                       ".tsv"),
                row.names = F, col.names = T, quote = F, sep = "\t")
}

### Top 1000
system("mkdir -p top_1000_markers_per_cluster")
marker_test_res <- top_markers(cds, 
                               group_cells_by="cluster",
                               genes_to_test_per_group = 1000,
                               marker_sig_test = T,
                               cores = threads,
                               verbose = T)

write.csv(marker_test_res, 
          "top_1000_markers_per_cluster/top_1000_markers_per_cluster.csv",
          row.names = F)

filtering_criteria = "specificity"

folder_name <- paste0("top_1000_markers_per_cluster/", "by_", filtering_criteria )
system( paste0("mkdir -p ", folder_name) )
for (c in unique(marker_test_res$cell_group) ) {
    
    top_specific_markers <- marker_test_res %>%
        dplyr::filter(marker_test_q_value < 0.05) %>%
        dplyr::filter(fraction_expressing >= 0.20) %>%
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
    
    write.table(top_specific_markers_top1000, 
                paste0(folder_name,
                       "/list_of_top_1000_for_cluster_",
                       c,
                       "_",filtering_criteria,
                       ".tsv"),
                row.names = F, col.names = T, quote = F, sep = "\t")
}