"Usage: step1_integrating_data_and_clustering_sunn4_only.R (--cores <cores>) (--UMI <umi>) (--out=<out>)
--cores=<cores> Number of cores [default:4].
--UMI <umi> Minimim UMI for a cell to be included [default:400].
--out=<out> Output prefix to the rds files [default:cds].
step1_integrating_data_and_clustering_sunn4_only.R -h | --help  show this message.
" -> doc

# retrieve the command-line arguments
suppressMessages( require(docopt) )
opts <- docopt(doc)

set.seed(1407)
suppressMessages( require(monocle3) )
suppressMessages( require(vroom) )
suppressMessages( require(cowplot) )
suppressMessages( require(ggplot2) )
suppressMessages( require(dplyr) )

threads <- as.numeric( opts$cores )
umi_theshold = as.numeric( opts$UMI )

cds_sunn_0h <- monocle3::load_cellranger_data("data/Sunn_sep_2022_0h_10k/",
                                              umi_cutoff = umi_theshold)
cds_sunn_24h <- monocle3::load_cellranger_data("data/Sunn_sep_2022_24h_10k/",
                                               umi_cutoff = umi_theshold)
cds_sunn_48h <- monocle3::load_cellranger_data("data/Sunn_sep_2022_48h_10k/",
                                               umi_cutoff = umi_theshold)
cds_sunn_96h <- monocle3::load_cellranger_data("data/Sunn_sep_2022_96h_10k/",
                                               umi_cutoff = umi_theshold)

colData(cds_sunn_0h)$timepoint <- "0h"
colData(cds_sunn_24h)$timepoint <- "24h"
colData(cds_sunn_48h)$timepoint <- "48h"
colData(cds_sunn_96h)$timepoint <- "96h"

colData(cds_sunn_0h)$batch2 <- "Sunn4_0h"
colData(cds_sunn_24h)$batch2 <- "Sunn4_24h"
colData(cds_sunn_48h)$batch2 <- "Sunn4_48h"
colData(cds_sunn_96h)$batch2 <- "Sunn4_96h"

colData(cds_sunn_0h)$Group <- "Sunn-4"
colData(cds_sunn_24h)$Group <- "Sunn-4"
colData(cds_sunn_48h)$Group <- "Sunn-4"
colData(cds_sunn_96h)$Group <- "Sunn-4"

cds <- monocle3::combine_cds(list(
    cds_sunn_0h,
    cds_sunn_24h,
    cds_sunn_48h,
    cds_sunn_96h))

rm(cds_sunn_0h,
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

cds_batched = align_cds(cds,
                        num_dim = 100,
                        alignment_group = "batch2")

cds_batched <- monocle3::preprocess_cds(cds_batched, num_dim = 100)
cds_batched = reduce_dimension(cds_batched)

( p2 <- plot_cells(cds_batched,
                   color_cells_by="batch2",
                   label_cell_groups=FALSE) )

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
cds_batched <- cluster_cells( cds_batched, 
                              verbose = T,
                              random_seed = 1407)
saveRDS( cds_batched,
         paste0("rds_files/",
                "batched_integrated_clustered_complete_dataset_SUNN4.rds") )

( p5 <- plot_cells(cds_batched,
                   group_label_size = 9, 
                   cell_size = 0.4,
                   alpha = 0.8) )
ggplot2::ggsave(filename = paste0("images/clustered_dataset_sunn4_",
                                  umi_theshold, "_UMI.png"),
                p5,
                height=10,
                width=15,
                units="cm",
                bg = "#FFFFFF")

ggplot2::ggsave(filename = paste0("images/clustered_dataset_sunn4_",
                                  umi_theshold, "_UMI.svg"),
                p5,
                height=10,
                width=15,
                units="cm",
                bg = "#FFFFFF")

( p8 <- plot_cells(cds_batched,
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

( p9 <- plot_cells(cds_batched,
                   group_label_size = 4,
                   cell_size = 0.5, 
                   alpha = 0.6) + 
        facet_wrap(~timepoint+Group, nrow = 2) + 
        theme(strip.text = element_text(size = 20)) )

ggplot2::ggsave(filename = paste0("images/clusteres_by_time_and_genotype",
                                  umi_theshold,
                                  "_UMI.png"),
                p9,
                height=22,
                width=50,
                units="cm",
                bg = "#FFFFFF")