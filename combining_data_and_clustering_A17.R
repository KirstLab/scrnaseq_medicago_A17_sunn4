"Usage: combining_data_and_clustering_A17.R (--cores <cores>) (--UMI <umi>) (--out=<out>)
--cores=<cores> Number of cores [default:4].
--UMI <umi> Minimim UMI for a cell to be included [default:400].
--out=<out> Output prefix to the rds files [default:cds].
combining_data_and_clustering_A17.R -h | --help  show this message.
" -> doc

set.seed(1407)

outfile <- "logs/combining_data_and_clustering_A17.out"
if ( file.exists(outfile) ) {
    file.remove(outfile)
}

my_log <- file(outfile)
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

suppressMessages( require(docopt) )
opts <- docopt(doc)

suppressMessages( require(monocle3) )
suppressMessages( require(vroom) )
suppressMessages( require(cowplot) )
suppressMessages( require(ggplot2) )
suppressMessages( require(dplyr) )

############################################################
## Loading the data and combining it in a single cds file ##
############################################################
threads <- as.numeric( opts$cores )
umi_theshold = as.numeric( opts$UMI )

## A17 WT
cds_0h <- readRDS("rds_files/A17_0h_after_scDblFinder.rds")
cds_24h <- readRDS("rds_files/A17_24h_after_scDblFinder.rds")
cds_48h <- readRDS("rds_files/A17_48h_after_scDblFinder.rds")
cds_96h <- readRDS("rds_files/A17_96h_after_scDblFinder.rds")

colData(cds_0h)$timepoint <- "0h"
colData(cds_24h)$timepoint <- "24h"
colData(cds_48h)$timepoint <- "48h"
colData(cds_96h)$timepoint <- "96h"

colData(cds_0h)$batch2 <- "A17_0h"
colData(cds_24h)$batch2 <- "A17_24h"
colData(cds_48h)$batch2 <- "A17_48h"
colData(cds_96h)$batch2 <- "A17_96h"

colData(cds_0h)$Group <- "A17-WT"
colData(cds_24h)$Group <- "A17-WT"
colData(cds_48h)$Group <- "A17-WT"
colData(cds_96h)$Group <- "A17-WT"

cds <- monocle3::combine_cds( list(cds_0h,
                                   cds_24h,
                                   cds_48h,
                                   cds_96h) )

rm(cds_0h,
   cds_24h,
   cds_48h, 
   cds_96h)

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

system("mkdir -p rds_files")
cds_batched <- cluster_cells( cds_batched, 
                              verbose = T,
                              random_seed = 1407
)
saveRDS( cds_batched,
         paste0("rds_files/",
                "batched_integrated_clustered_complete_dataset_A17.rds") )

( p5 <- plot_cells(cds_batched,
                   group_label_size = 9,
                   cell_size = 0.4,
                   alpha = 0.8) )

ggplot2::ggsave(filename = paste0("images/clustered_dataset_A17_",
                                  umi_theshold, "_UMI.png"),
                p5,
                height=10,
                width=15,
                units="cm",
                bg = "#FFFFFF")

ggplot2::ggsave(filename = paste0("images/clustered_dataset_A17_",
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

closeAllConnections()