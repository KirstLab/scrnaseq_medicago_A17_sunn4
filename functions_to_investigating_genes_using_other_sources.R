my_complexheatmap <- function(scaled_data_frame = scaled_data_frame,
                              arg_show_row_names = opts$show_gene_names,
                              arg_column_title = arg_column_title) {
    
    if ( opts$color_scheme == "MAGMA" ) {
        
        color_scheme <- scales::viridis_pal(option = "A")(11)
        
    } else if ( opts$color_scheme == "VIRIDIS" ) {
        
        color_scheme <- scales::viridis_pal(option = "D")(11)
        
    } else {
        
        print("The --color_scheme option must be 'VIRIDIS' or 'MAGMA'!")
        quit()
        
    }
    
    col_fun = circlize::colorRamp2( breaks = seq(0, 1, 0.1),
                                    colors = color_scheme)
    
    scaled_data_frame <- scaled_data_frame[ c(t_genes), ]
    
    ComplexHeatmap::Heatmap(scaled_data_frame,
                            border = TRUE,
                            rect_gp = grid::gpar(col = "white", lwd = 0.01),
                            name = "Expression",
                            cluster_rows = T,
                            cluster_columns = F,
                            col = col_fun,
                            show_row_names = arg_show_row_names,
                            show_column_names = T,
                            show_row_dend = F,
                            use_raster = F,
                            row_names_gp = grid::gpar(fontsize = 12),
                            column_title = arg_column_title)
    
}

LCM_func_v5 <- function(genes_names = genes_names,
                        plot_name = plot_name,
                        show_names = opts$show_gene_names) {
    
    LCM_sub <- LCMv5 %>%
        filter(gene %in% genes_names)
    
    print( paste ("From the", length(genes_names),
                  "input genes, a total of",
                  length(unique(LCM_sub$gene)),
                  "genes IDs had correspondents in the LCM data"))

    LCM_sub <- LCM_sub %>% 
        distinct() %>%
        column_to_rownames("gene")
    
    LCM_q_scaled <- t(dynutils::scale_quantile(t(LCM_sub)))
    
    LCM_0h <- LCM_q_scaled[, grepl(colnames(LCM_q_scaled), pattern = "0") ]
    LCM_12h <- LCM_q_scaled[, grepl(colnames(LCM_q_scaled), pattern = "12") ]
    LCM_24h <- LCM_q_scaled[, grepl(colnames(LCM_q_scaled), pattern = "24") ]
    LCM_48h <- LCM_q_scaled[, grepl(colnames(LCM_q_scaled), pattern = "48") ]
    LCM_72h <- LCM_q_scaled[, grepl(colnames(LCM_q_scaled), pattern = "72") ]
    
    ( p_0h <- my_complexheatmap(LCM_0h, arg_column_title = "0h") )
    ( p_12h <- my_complexheatmap(LCM_12h, arg_column_title = "12h") )
    ( p_24h <- my_complexheatmap(LCM_24h, arg_column_title = "24h") )
    ( p_48h <- my_complexheatmap(LCM_48h, arg_column_title = "48h") )
    ( p_72h <- my_complexheatmap(LCM_72h, arg_column_title = "72h") )
    
    (lcm_plot <- suppressWarnings( ComplexHeatmap::draw(p_0h + p_12h + p_24h + p_48h + p_72h,
                                                        ht_gap = unit(0.5, "cm") ) ) )
    
    if (show_names == TRUE) {
        
        plots_height <- 5 + ( 0.2 * nrow(p_0h) )
        
    } else if (show_names == FALSE) {
        
        plots_height <- 5 + ( 0.1 * nrow(p_0h) )
        
    }
    
    if(opts$image_format == "png") {
        
        png(  paste0(opts$out_images,
                     "/",
                     plot_name,
                     "_", "LCM_plot.png" ),
              width = (ncol(p_0h) * 1),
              height = plots_height,
              pointsize = 12,
              res = 300, units = "in")
        ComplexHeatmap::draw(lcm_plot)
        dev.off()
        
    } else if (opts$image_format == "svg") {
        
        svg(  paste0(opts$out_images,
                     "/",
                     plot_name,
                     "_", "LCM_plot.svg" ),
              width = (ncol(p_0h) * 1),
              height = plots_height,
              pointsize = 12 )
        ComplexHeatmap::draw(lcm_plot)
        dev.off()
        
    } else {
        print("--image_format must be 'png' or 'svg'!")
        quit()
    }
    
}

# Generates the heatmaps of the Roux_2014 dataset
atlast_v2_SRP028599 <- function(genes_names = genes_names,
                                plot_name = plot_name,
                                SRP = "SRP028599",
                                show_names = opts$show_gene_names) {
    
    atlas_v2_sub <- atlas_v2 %>%
        dplyr::filter(locus_tag %in% genes_names) %>%
        column_to_rownames("locus_tag")
    
    print( paste ("A total of", nrow(atlas_v2_sub), "genes IDs, from",
                  length(genes_names), "had correspondent gene IDs in the dataset"))
    
    # Select the data from Rox et al 2014, based on the project number
    SRP_complete <- atlas_v2_sub[, grepl(colnames(atlas_v2_sub), pattern = SRP)]
    SRP_scaled <- t(dynutils::scale_quantile( t(SRP_complete) ) )
    
    ( p_SRP <- my_complexheatmap(SRP_scaled,
                                 arg_column_title = SRP) )
    
    if (show_names == TRUE) {
        
        plots_height <- 5 + ( 0.2 * nrow(p_SRP) )
        
    } else if (show_names == FALSE) {
        
        plots_height <- 5 + ( 0.1 * nrow(p_SRP) )
        
    }
    
    png( paste0(opts$out_images, "/",
                plot_name, "_",
                SRP, "_Roux_2014_complete.png" ),
         width = (ncol(p_SRP) * 0.25),
         height = plots_height,
         pointsize = 12,
         res = 300, units = "in")
    ComplexHeatmap::draw(p_SRP)
    dev.off()
    
    SRP_complete_LCM <- SRP_complete[, grepl(colnames(SRP_complete), pattern = "Sm_RbmL_Nodule")]
    SRP_scaled <- t(dynutils::scale_quantile( t(SRP_complete_LCM) ) )
    
    ( p_SRP <- my_complexheatmap(SRP_scaled,
                                 arg_column_title = SRP) )
    
    if (show_names == TRUE) {
        
        plots_height <- 5 + ( 0.2 * nrow(p_SRP) )
        
    } else if (show_names == FALSE) {
        
        plots_height <- 5 + ( 0.1 * nrow(p_SRP) )
        
    }
    
    png( paste0(opts$out_images, "/",
                plot_name, "_",
                SRP, "_Roux_2014_LCM_of_nodules.png" ),
         width = (ncol(p_SRP) * 1),
         height = plots_height,
         pointsize = 12,
         res = 300, units = "in")
    ComplexHeatmap::draw(p_SRP)
    dev.off()
    
    SRP_complete_roots <- SRP_complete[, !grepl(colnames(SRP_complete), pattern = "Sm_RbmL_Nodule")]
    SRP_scaled <- t(dynutils::scale_quantile( t(SRP_complete_roots) ) )
    
    ( p_SRP <- my_complexheatmap(SRP_scaled,
                                 arg_column_title = SRP) )
    if (show_names == TRUE) {
        
        plots_height <- 5 + ( 0.2 * nrow(p_SRP) )
        
    } else if (show_names == FALSE) {
        
        plots_height <- 5 + ( 0.1 * nrow(p_SRP) )
        
    }
    png( paste0(opts$out_images, "/",
                plot_name, "_",
                SRP, "_Roux_2014_roots.png" ),
         width = (ncol(p_SRP) * 1),
         height = plots_height,
         pointsize = 12,
         res = 300, units = "in")
    ComplexHeatmap::draw(p_SRP)
    dev.off()
}

## Schiessl et al., 2019
atlast_v2_SRP212693 <- function(genes_names = genes_names,
                                plot_name = plot_name,
                                SRP = "SRP212693",
                                show_names = opts$show_gene_names) {
    
    atlas_v2_sub <- atlas_v2 %>%
        dplyr::filter(locus_tag %in% genes_names) %>%
        column_to_rownames("locus_tag")
    
    print( paste ("A total of", nrow(atlas_v2_sub), "genes IDs, from",
                  length(genes_names), "had correspondent gene IDs in the dataset"))
    
    SRP_complete <- atlas_v2_sub[, grepl( colnames( atlas_v2_sub ), pattern = SRP ) ]
    
    SRP_WT <- SRP_complete[, grepl( colnames( SRP_complete ), pattern = "WT" ) ]
    SRP_sub <- SRP_WT[, grepl( colnames( SRP_WT ), pattern = "Smeliloti" ) ]
    
    # Reorder to correct the timepoint
    SRP_sub <- SRP_sub[, c(1:4, 35:40, 51:55, 60:63, 5:8, 13:16, 17:21,
                           26:30, 31:34, 41:50, 56:59, 64:67, 9:12, 22:25) ] 
    
    SRP_scaled <- t( dynutils::scale_quantile( t(SRP_sub) ) )
    
    ( p_SRP <- my_complexheatmap(SRP_scaled,
                                 arg_column_title = "SRP212693 - Schiessl et al., 2019.") )
    
    if (show_names == TRUE) {
        
        plots_height <- 5 + ( 0.2 * nrow(p_SRP) )
        
    } else if (show_names == FALSE) {
        
        plots_height <- 5 + ( 0.1 * nrow(p_SRP) )
        
    }
    
    png( paste0( opts$out_images, "/",
                 plot_name, "_", 
                 SRP, "_Schiessl_2019_complete.png" ),
         width = (ncol(p_SRP) * 0.25),
         height = plots_height,
         pointsize = 12,
         res = 300, units = "in")
    ComplexHeatmap::draw(p_SRP)
    dev.off()
}