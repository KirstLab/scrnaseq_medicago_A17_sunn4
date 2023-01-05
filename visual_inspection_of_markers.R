"Usage: visual_inspection_of_markers.R (--rds <rds>) (--markers <markers>) [--common=<common>] (--out <out>)
--rds    RDS file containing the clustered dataset.
--common=<common> Indicates if the second column in the markers file contains the common name of the gene (or any character string) to be used to label the plots [default: FALSE].
--markers   File containing the list of markers to be used, must be in the csv format.
--out <out> Output directory where the images will be saved.
visual_inspection_of_markers.R -h | --help  show this message.
" -> doc

suppressMessages( require(docopt) )
opts <- docopt(doc)

set.seed(1407)
suppressMessages( require(monocle3) )
suppressMessages( require(vroom) )
suppressMessages( require(cowplot) )
suppressMessages( require(ggplot2) )
suppressMessages( require(dplyr) )
suppressMessages( require(stringr) )

cds <- readRDS(opts$rds)

( pc <- plot_cells(cds,
                   cell_size = 0.45,
                   group_label_size = 4,
                   scale_to_range = F) )

ext <- tools::file_ext(opts$markers)

top_specific_marker_ids_df <- switch(ext,
                                     csv = vroom::vroom(opts$markers, delim = ","),
                                     tsv = vroom::vroom(opts$markers, delim = "\t"),
                                     validate("Invalid file; Please upload a .csv or .tsv file") )

top_specific_marker_ids_df <- top_specific_marker_ids_df[complete.cases(top_specific_marker_ids_df[, 1]), ]

colnames(top_specific_marker_ids_df)[1] <- "genes"

top_specific_marker_ids <- top_specific_marker_ids_df %>%
    dplyr::filter(genes %in% rownames(cds))

print( paste( "From the ", nrow(unique(top_specific_marker_ids_df[, 1] ) ),
              "markers, a total of ",
              nrow(unique(top_specific_marker_ids[, 1] ) ),
              "were found in the dataset."
) )

system( paste0("mkdir -p ", opts$out) )
for( m in 1:nrow(top_specific_marker_ids) ) {
    
    name_2_title <- if( opts$common == T ) {
        
        top_specific_marker_ids[m, 2]
        
    } else if( opts$common == F ) {
        
        top_specific_marker_ids[m, 1]
    }
    
    geneID <- top_specific_marker_ids[m, 1]
    
    ( p_by_time <- plot_cells(cds,
                              genes = geneID,
                              labels_per_group = F,
                              cell_size = 0.6,
                              scale_to_range = F) +
            facet_wrap(~Group+timepoint, nrow = 2) + 
            theme(legend.title = element_blank()) +
            theme(strip.text = element_text(size = 14),
                  plot.title = element_text(size=20,
                                            face = "bold") ) +
            viridis::scale_color_viridis(option = "A")  )
    
    ( p_by_sample <- plot_cells(cds,
                                genes = geneID,
                                labels_per_group = F,
                                cell_size = 0.6,
                                scale_to_range = F) +
            theme(legend.position="none") +
            facet_wrap(~Group, nrow = 2) +
            theme(strip.text = element_text(size = 14),
                  plot.title = element_text(size=20,
                                            face = "bold") ) +
            viridis::scale_color_viridis(option = "A")  )
    
    ( p_combined <- plot_cells(cds,
                               genes = geneID,
                               labels_per_group = F,
                               cell_size = 0.6,
                               scale_to_range = F) +
            theme(legend.position="none") +
            theme(strip.text = element_blank() ) +
            viridis::scale_color_viridis(option = "A")  )
    
    ( p_combined2 <- cowplot::plot_grid(pc, p_combined, nrow = 2) )
    
    ( p <- cowplot::plot_grid(p_combined2,
                              p_by_sample,
                              p_by_time,
                              nrow = 1,
                              ncol = 3,
                              rel_widths = c(1.1,1.1,4) ) )
    
    title <- ggdraw() + 
        draw_label(
            name_2_title,
            fontface = 'bold',
            size = 16,
            x = 0,
            hjust = 0
        ) +
        theme(
            plot.margin = margin(0, 0, 0, 7)
        )
    
    ( p_combined <- cowplot::plot_grid(
        title, p,
        ncol = 1,
        rel_heights = c(0.075, 1)
    ) ) 
    
    ggplot2::ggsave(paste0(opts$out,"/", name_2_title, ".png"),
                    p_combined,
                    height=20,
                    width=50,
                    units="cm",
                    dpi = 300,
                    bg = "#FFFFFF")
}
