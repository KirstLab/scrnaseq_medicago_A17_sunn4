set.seed(1407)

outfile <- "logs/combining_all_samples_and_clustering.out" # File name of output log
#Check its existence
if ( file.exists(outfile) ) {
    #Delete file if it exists
    file.remove(outfile)
}

my_log <- file(outfile) 
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

suppressMessages( require(monocle3) )
suppressMessages( require(vroom) )
suppressMessages( require(cowplot) )
suppressMessages( require(ggplot2) )
suppressMessages( require(dplyr) )

threads <- as.numeric(8)
umi_theshold <- as.numeric(400)

## A17 WT
cds_0h <- readRDS("rds_files/cds_A17_0h_after_scDblFinder.rds")
cds_24h <- readRDS("rds_files/cds_A17_24h_after_scDblFinder.rds")
cds_48h <- readRDS("rds_files/cds_A17_48h_after_scDblFinder.rds")
cds_96h <- readRDS("rds_files/cds_A17_96h_after_scDblFinder.rds")
## sunn-4
cds_sunn_0h <- readRDS("rds_files/cds_sunn_0h_after_scDblFinder.rds")
cds_sunn_24h <- readRDS("rds_files/cds_sunn_24h_after_scDblFinder.rds")
cds_sunn_48h <- readRDS("rds_files/cds_sunn_48h_after_scDblFinder.rds")
cds_sunn_96h <- readRDS("rds_files/cds_sunn_96h_after_scDblFinder.rds")

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
cds <- align_cds(cds,
                 num_dim = 100,
                 alignment_group = "batch2",
                 verbose = T)

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

cds <- cluster_cells( cds,
                      verbose = T,
                      random_seed = 1407
)

rowData(cds)$gene_short_name <- row.names(rowData(cds))

( p5.1 <- plot_cells(cds,
                     group_label_size = 11,
                     cell_size = 0.75,
                     alpha = 0.7) )

print("Dimentions of the dataset BEFORE removing the doublets")
dim(colData(cds))
( pdoub <- plot_cells(cds,
                      color_cells_by="scDblFinder",
                      label_cell_groups=FALSE) + facet_wrap(~scDblFinder) )

ggplot2::ggsave(filename = paste0("images/clustered_dataset_", umi_theshold, "_UMI_singlet_vs_doublet.png"),
                pdoub,
                height=20,
                width=30,
                units="cm",
                bg = "#FFFFFF")

ggplot2::ggsave(filename = paste0("images/clustered_dataset_", umi_theshold, "_UMI_singlet_vs_doublet.svg"),
                pdoub,
                height=20,
                width=30,
                units="cm",
                bg = "#FFFFFF")


## Remove doublets
doublets <- cds@colData
doublets <- doublets[doublets$scDblFinder == "doublet", ]
cells_that_are_doublets <- unique(rownames(doublets))

cds <- cds[, !colnames(cds) %in% cells_that_are_doublets ]

( p5.2 <- plot_cells(cds,
                     group_label_size = 11,
                     cell_size = 0.75,
                     alpha = 0.7) )

print("Dimentions of the dataset AFTER removing the doublets")
dim(colData(cds))
plot_cells(cds,
           color_cells_by="scDblFinder",
           label_cell_groups=FALSE)

( p5 <- cowplot::plot_grid(p5.1, p5.2) )

saveRDS( cds,
         paste0("rds_files/", "batched_integrated_clustered_complete_dataset.rds") )

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
folder_name2 <- paste0("top_1000_markers_per_cluster/", "by_", filtering_criteria, "/top100" )
system( paste0("mkdir -p ", folder_name) )
system( paste0("mkdir -p ", folder_name2) )
for (c in unique(marker_test_res$cell_group) ) {
    
    top_specific_markers <- marker_test_res %>%
        dplyr::filter(marker_test_q_value < 0.05) %>%
        dplyr::filter(fraction_expressing >= 0.20) %>%
        dplyr::filter(cell_group == c) %>%
        dplyr::arrange( desc(specificity) )
    
    write.table(top_specific_markers, 
                paste0(folder_name,
                       "/list_of_top_1000_for_cluster_",
                       c,
                       "_",filtering_criteria,
                       ".tsv"),
                row.names = F, col.names = T, quote = F, sep = "\t")
    
    top_specific_markers_top100 <- top_specific_markers[1:100, ]
    
    write.table(top_specific_markers_top100, 
                paste0(folder_name2,
                       "/list_of_top_100_for_cluster_",
                       c,
                       "_",filtering_criteria,
                       ".tsv"),
                row.names = F, col.names = T, quote = F, sep = "\t")
}

closeAllConnections()