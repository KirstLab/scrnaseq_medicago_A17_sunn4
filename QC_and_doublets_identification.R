set.seed(1407)

system("mkdir -p logs")
system("mkdir -p rds_files")

outfile <- "logs/QC_and_doublets_identification.out" # File name of output log
#Check its existence
if ( file.exists(outfile) ) {
    #Delete file if it exists
    file.remove(outfile)
}

my_log <- file(outfile) 
#sink(my_log, append = TRUE, type = "output")
#sink(my_log, append = TRUE, type = "message")

suppressMessages(require(scDblFinder))
suppressMessages(require(monocle3))
suppressMessages( require(vroom) )
suppressMessages( require(cowplot) )
suppressMessages( require(ggplot2) )
suppressMessages( require(dplyr) )

umi_theshold <- as.numeric(400)

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

########################
## Find doublet cells ##
########################
sample_list <- list("cds_A17_0h" = cds_0h,
                    "cds_A17_24h" = cds_24h,
                    "cds_A17_48h" = cds_48h,
                    "cds_A17_96h" = cds_96h,
                    "cds_sunn_0h" = cds_sunn_0h,
                    "cds_sunn_24h" = cds_sunn_24h,
                    "cds_sunn_48h" = cds_sunn_48h,
                    "cds_sunn_96h" = cds_sunn_96h)

for( s in 1:length(sample_list) ) {
    
    print(paste0( "Searching for doublets on sample: ", names(sample_list)[s] ))
    
    sce <- scDblFinder(
        sce = monocle3::exprs(sample_list[[s]]),
        threshold = T
    )
    
    cds_set <- sample_list[[s]]
    colData(cds_set)$scDblFinder <- sce$scDblFinder.class
    
    saveRDS(cds_set, 
            paste0( "rds_files/",
                    names(sample_list)[s], "_after_scDblFinder.rds") )
}

closeAllConnections()
