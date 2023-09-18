set.seed(1407)

outfile <- "logs/produce_gifs_from_umaps_video3_cortex_nodule.out" # File name of output log
#Check its existence
if ( file.exists(outfile) ) {
    #Delete file if it exists
    file.remove(outfile)
}

my_log <- file(outfile)
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

suppressMessages( require(vroom) )
suppressMessages( require(tidyverse) )
suppressMessages( require(monocle3) )
suppressMessages( require(animation) )

# Reads the RDS file with the clustered data
cds <- readRDS("RECLUSTERING/CORTEX_NOD_sunn4_all_clusters_but_CN5/rds_file_subset/medicago_integrated_subset_cortex_nodule.rds")

colData(cds)$Group <- ifelse(colData(cds)$Group == "A17-WT", "A17", "sunn-4")

## To filter the cells of each time point, we need to split the dataset before plotting. I am not sure if this is the easiest way, but I was not able to find an alternative.
cells_on_time <- cds@colData

#### 0h
cells_on_time <- cells_on_time[cells_on_time$timepoint == "0h", ]
cells_on_time_names <- unique(rownames(cells_on_time))

cds_subset <- cds[, colnames(cds) %in% cells_on_time_names ]

( p0 <- plot_cells(cds_subset,
                   #genes = gene_name,
                   group_label_size = 0,
                   cell_size = 1.1,
                   scale_to_range = F)+ 
        ggtitle("0h") + 
        xlim(-10, 10) +
        ylim(-6, 6) +
        theme(strip.text = element_text(size = 24),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_text(size = 40, 
                                        hjust = 0.5),
              legend.title=element_blank()) )

#### 24h
cells_on_time <- cds@colData
cells_on_time <- cells_on_time[cells_on_time$timepoint == "24h", ]
cells_on_time_names <- unique(rownames(cells_on_time))

cds_subset <- cds[, colnames(cds) %in% cells_on_time_names ]

( p24h <- plot_cells(cds_subset,
                     #genes = gene_name,
                     group_label_size = 0,
                     cell_size = 1.1,
                     scale_to_range = F) + 
        ggtitle("24h") +
        xlim(-6, 8) +
        xlim(-10, 10) +
        ylim(-6, 6) +
        theme(strip.text = element_text(size = 24),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_text(size = 40, 
                                        hjust = 0.5),
              legend.title=element_blank()) )

#### 48h
cells_on_time <- cds@colData
cells_on_time <- cells_on_time[cells_on_time$timepoint == "48h", ]
cells_on_time_names <- unique(rownames(cells_on_time))

cds_subset <- cds[, colnames(cds) %in% cells_on_time_names ]

( p48h <- plot_cells(cds_subset,
                     #genes = gene_name,
                     group_label_size = 0,
                     cell_size = 1.1,
                     scale_to_range = F) + 
        ggtitle("48h") + 
        xlim(-10, 10) +
        ylim(-6, 6) +
        theme(strip.text = element_text(size = 24),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_text(size = 40, 
                                        hjust = 0.5),
              legend.title=element_blank()) )

#### 96h
cells_on_time <- cds@colData
cells_on_time <- cells_on_time[cells_on_time$timepoint == "96h", ]
cells_on_time_names <- unique(rownames(cells_on_time))

cds_subset <- cds[, colnames(cds) %in% cells_on_time_names ]

( p96h <- plot_cells(cds_subset,
                     #genes = gene_name,
                     group_label_size = 0,
                     cell_size = 1.1,
                     scale_to_range = F) + 
        ggtitle("96h") + 
        xlim(-10, 10) +
        ylim(-6, 6) +
        theme(strip.text = element_text(size = 24),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              plot.title = element_text(size = 40, 
                                        hjust = 0.5),
              legend.title=element_blank()) )

## Here, we create a new folder to save the UMAP plots
folder_name = "medicago_GIFs/VIDEO3/"

system( paste0("mkdir -p ", folder_name) )

ggsave(plot = p0,
       height = 15,
       width = 20,
       units = "cm",
       filename = paste0(folder_name, "p0h.svg") )

ggsave(plot = p24h,
       height = 15,
       width = 20,
       units = "cm",
       filename = paste0(folder_name, "p24h.svg") )

ggsave(plot = p48h,
       height = 15,
       width = 20,
       units = "cm",
       filename = paste0(folder_name, "p48h.svg") )

ggsave(plot = p96h,
       height = 15,
       width = 20,
       units = "cm",
       filename = paste0(folder_name, "p96h.svg") )

## This package requires magick to be installed on mac (you can do it using homebrew)
require(magick)

# Reads the UMAPs from the folders we save them. WARNING: Make sure that only the plots you want are saved on that folder. all images will be included on the gif
imgs <- list.files(folder_name,
                   full.names = TRUE)

## list file names and read in
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 0.5 frames per second
img_animated <- image_animate(img_joined, fps = 0.5)

## view animated image
img_animated

## save to disk
image_write(image = img_animated,
            path = paste0(folder_name, "Pereira_el_al_2023_video3.GIF") )

closeAllConnections()
