my_complexheatmap <- function(scaled_data_frame = scaled_data_frame,
                              arg_show_row_names = opts$show_gene_names,
                              arg_column_title = arg_column_title) {
    
    if ( opts$color_scheme == "MAGMA" ) {
        
        color_scheme <- scales::viridis_pal(option = "A")(11)
        
    } else if ( opts$color_scheme == "VIRIDIS" ) {
        
        color_scheme <- scales::viridis_pal(option = "plasma")(11)
        
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
                            show_heatmap_legend = T,
                            show_row_dend = F,
                            use_raster = F,
                            row_names_gp = grid::gpar(fontsize = 12),
                            column_title = arg_column_title)
    
}

LCM_func_v5 <- function(genes_names = genes_names,
                        plot_name = plot_name,
                        show_names = opts$show_gene_names ) {
    
    LCM_sub <- LCMv5 %>%
        filter(gene %in% genes_names)
    
    LCM_sub <- LCM_sub %>%
        rowwise() %>% 
        dplyr::mutate( EPI_0h = mean( c(JFED0, JFED0B, JFED0C) ),
                       EPI_12h = mean( c(JFED12, JFED12B, JFED12C) ),
                       EPI_24h = mean( c(JFED24, JFED24B, JFED24C) ),
                       EPI_48h = mean( c(JFED48, JFED48B, JFED48C) ),
                       EPI_72h = mean( c(JFED72, JFED72B) ),
                       
                       ICA_0h = mean( c(JFICA0, JFICA0B, JFICA0C) ),
                       ICA_12h = mean( c(JFICA12, JFICA12B, JFICA12C) ),
                       ICA_24h = mean( c(JFICA24, JFICA24B, JFICA24C) ),
                       ICA_48h = mean( c(JFICA48, JFICA48B, JFICA48C) ),
                       ICA_72h = mean( c(JFICA72, JFICA72B, JFICA72C) ),
                       
                       ICB_0h = mean( c(JFICB0, JFICB0B, JFICB0C) ),
                       ICB_12h = mean( c(JFICB12, JFICB12B, JFICB12C) ),
                       ICB_24h = mean( c(JFICB24, JFICB24B, JFICB24C) ),
                       ICB_48h = mean( c(JFICB48, JFICB48B, JFICB48C) ),
                       ICB_72h = mean( c(JFICB72, JFICB72B, JFICB72C) ),
                       
                       NOD_48h = mean( c(JFNOD48, JFNOD48B, JFNOD48C) ),
                       NOD_72h = mean( c(JFNOD72, JFNOD72B, JFNOD72C) ),
                       
                       OC_0h = mean( c(JFOC0, JFOC0B, JFOC0C) ),
                       OC_12h = mean( c(JFOC12, JFOC12B, JFOC12C) ),
                       OC_24h = mean( c(JFOC24, JFOC24B, JFOC24C) ),
                       OC_48h = mean( c(JFOC48, JFOC48B, JFOC48C) ),
                       OC_72h = mean( c(JFOC72, JFOC72B, JFOC72C) ),
                       
                       VS_0h = mean( c(JFVS0, JFVS0B, JFVS0C) ),
                       VS_12h = mean( c(JFVS12, JFVS12B, JFVS12C) ),
                       VS_24h = mean( c(JFVS24, JFVS24B, JFVS24C) ),
                       VS_48h = mean( c(JFVS48, JFVS48B, JFVS48C) ),
                       VS_72h = mean( c(JFVS72, JFVS72B, JFVS72C) ),
        ) %>%
        dplyr::select(gene, 
                      EPI_0h, EPI_12h, EPI_24h, EPI_48h, EPI_72h,
                      ICA_0h, ICA_12h, ICA_24h, ICA_48h, ICA_72h,
                      ICB_0h, ICB_12h, ICB_24h, ICB_48h, ICB_72h,
                      NOD_48h, NOD_72h,
                      OC_0h, OC_12h, OC_24h, OC_48h, OC_72h,
                      VS_0h, VS_12h, VS_24h, VS_48h, VS_72h) 
    
    print( paste ("From the", length(genes_names),
                  "input genes, a total of",
                  length(unique(LCM_sub$gene)),
                  "genes IDs had correspondents in the LCM data"))
    
    LCM_sub <- LCM_sub %>% 
        distinct() %>%
        column_to_rownames("gene")
    
    LCM_q_scaled <- t(dynutils::scale_quantile(t(LCM_sub)))
    
    LCM_0h <- LCM_q_scaled[, grepl(colnames(LCM_q_scaled), pattern = "0") ]
    LCM_0h <- LCM_0h[, c(1,4,2,3,5)]
    LCM_0h <- LCM_0h %>%
        as.data.frame() %>%
        dplyr::rename(EPI = EPI_0h,
                      ICA = ICA_0h,
                      ICB = ICB_0h,
                      OC = OC_0h,
                      VS = VS_0h)
    
    LCM_12h <- LCM_q_scaled[, grepl(colnames(LCM_q_scaled), pattern = "12") ]
    LCM_12h <- LCM_12h[, c(1,4,2,3,5)]
    LCM_12h <- LCM_12h %>%
        as.data.frame() %>%
        dplyr::rename(EPI = EPI_12h,
                      ICA = ICA_12h,
                      ICB = ICB_12h,
                      OC = OC_12h,
                      VS = VS_12h)
    
    LCM_24h <- LCM_q_scaled[, grepl(colnames(LCM_q_scaled), pattern = "24") ]
    LCM_24h <- LCM_24h[, c(1,4,2,3,5)]
    LCM_24h <- LCM_24h %>%   
        as.data.frame() %>%
        dplyr::rename(EPI = EPI_24h,
                      ICA = ICA_24h,
                      ICB = ICB_24h,
                      OC = OC_24h,
                      VS = VS_24h)
    
    LCM_48h <- LCM_q_scaled[, grepl(colnames(LCM_q_scaled), pattern = "48") ]
    LCM_48h <- LCM_48h[, c(1,5,2,3,6,4)]
    LCM_48h <- LCM_48h %>%
        as.data.frame() %>%
        dplyr::rename(EPI = EPI_48h,
                      ICA = ICA_48h,
                      ICB = ICB_48h,
                      NOD = NOD_48h,
                      OC = OC_48h,
                      VS = VS_48h)
    
    LCM_72h <- LCM_q_scaled[, grepl(colnames(LCM_q_scaled), pattern = "72") ]
    LCM_72h <- LCM_72h[, c(1,5,2,3,6,4)]
    LCM_72h <- LCM_72h %>%
        as.data.frame() %>%
        dplyr::rename(EPI = EPI_72h,
                      ICA = ICA_72h,
                      ICB = ICB_72h,
                      NOD = NOD_72h,
                      OC = OC_72h,
                      VS = VS_72h)
    
    ( p_0h <- my_complexheatmap(LCM_0h, arg_column_title = "0h") )
    ( p_12h <- my_complexheatmap(LCM_12h, arg_column_title = "12h") )
    ( p_24h <- my_complexheatmap(LCM_24h, arg_column_title = "24h") )
    ( p_48h <- my_complexheatmap(LCM_48h, arg_column_title = "48h") )
    ( p_72h <- my_complexheatmap(LCM_72h, arg_column_title = "72h") )
    
    (lcm_plot <- suppressWarnings( ComplexHeatmap::draw(p_0h + p_12h + p_24h + p_48h + p_72h,
                                                        ht_gap = unit(0.5, "cm") ) ) )
    
    if (show_names == TRUE) {
        
        plots_height <- 5 + ( 0.2 * nrow(LCM_0h) )
        
    } else if (show_names == FALSE) {
        
        plots_height <-  0.03 * nrow(LCM_0h) 
        
    }
    
    if(opts$image_format == "png") {
        
        png(  paste0(opts$out_images,
                     "/",
                     plot_name,
                     "_", "LCM_plot.png" ),
              width = (ncol(p_0h)/3)*5,
              height = plots_height,
              pointsize = 12,
              res = 300, units = "in")
        ComplexHeatmap::draw(lcm_plot)
        dev.off()
        
    } else if ( opts$image_format == "svg" ) {
        
        svg(  paste0(opts$out_images,
                     "/",
                     plot_name,
                     "_", "LCM_plot.svg" ),
              width = (ncol(p_0h)/3)*5,
              height = plots_height,
              pointsize = 12 )
        ComplexHeatmap::draw(lcm_plot)
        dev.off()
        
    } else {
        print("--image_format must be 'png' or 'svg'!")
        quit()
    }
    
}

# Generates the heatmaps of the top 100 expressed genes per cluster, comparing among clusters.
sc_cluster <- function(genes_names = genes_names,
                       plot_name = plot_name,
                       show_names = opts$show_gene_names,
                       cds_path = cds_path,
                       c_order = c_order) {
    
    cds <- readRDS(cds_path)
    
    cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                    cell_group=clusters(cds)[colnames(cds)])
    
    agg_mat <- aggregate_gene_expression(cds, NULL,
                                         cell_group_df,
                                         scale_agg_values = F)
    
    agg_mat2 <- agg_mat[rownames(agg_mat) %in% genes_names, ]
    
    agg_mat2 <- as.matrix(agg_mat2)
    
    agg_q_scaled <- t(dynutils::scale_quantile(t(agg_mat2)))
    agg_q_scaled <- agg_q_scaled[, c_order]
    
    ( p_agg <- my_complexheatmap( scaled_data_frame = agg_q_scaled,
                                  arg_column_title = "Clusters" ) )
    
    if (show_names == TRUE) {
        
        plots_height <- 5 + ( 0.2 * nrow(p_agg) )
        
    } else if (show_names == FALSE) {
        
        plots_height <-  0.03 * nrow(p_agg) 
        
    }
    
    if(opts$image_format == "png") {
        
        png(  paste0(opts$out_images,
                     "/",
                     plot_name,
                     "_", "HEATMAP_among_clusters.png" ),
              width = ncol(p_agg)/3,
              height = plots_height,
              pointsize = 12,
              res = 300, units = "in")
        ComplexHeatmap::draw(p_agg)
        dev.off()
        
    } else if (opts$image_format == "svg") {
        
        svg(  paste0(opts$out_images,
                     "/",
                     plot_name,
                     "_", "HEATMAP_among_clusters.svg" ),
              width = ncol(p_agg)/3,
              height = plots_height,
              pointsize = 12 )
        ComplexHeatmap::draw(p_agg)
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
    
    SRP_complete <- atlas_v2_sub[, grepl(colnames(atlas_v2_sub), pattern = SRP)]
    SRP_complete_LCM <- SRP_complete[, grepl(colnames(SRP_complete), pattern = "RbmL")]
    
    SRP_complete_LCM$genes <- rownames(SRP_complete_LCM)
    
    SRP_complete_LCM <- SRP_complete_LCM %>%
        rowwise() %>% 
        dplyr::mutate(
            Roots = mean( c(SRP028599.Mt_RbmL_Root_R1,
                            SRP028599.Mt_RbmL_Root_R2,
                            SRP028599.Mt_RbmL_Root_R3) ) ,
            
            'N10d' = mean( c(SRP028599.Mt_RbmL_Nodule_R1,
                             SRP028599.Mt_RbmL_Nodule_R2,
                             SRP028599.Mt_RbmL_Nodule_R3) ) ,
            
            FI = mean( c( SRP028599.Mt_Sm_RbmL_Nodule_ZI_R1,
                          SRP028599.Mt_Sm_RbmL_Nodule_ZI_R2,
                          SRP028599.Mt_Sm_RbmL_Nodule_ZI_R3 ) ) ,
            
            FIId = mean( c(SRP028599.Mt_Sm_RbmL_Nodule_ZIId_R1,
                           SRP028599.Mt_Sm_RbmL_Nodule_ZIId_R2,
                           SRP028599.Mt_Sm_RbmL_Nodule_ZIId_R3) ),
            
            FIIp = mean( c(SRP028599.Mt_Sm_RbmL_Nodule_ZIIp_R1,
                           SRP028599.Mt_Sm_RbmL_Nodule_ZIIp_R2,
                           SRP028599.Mt_Sm_RbmL_Nodule_ZIIp_R3) ),
            
            IZ = mean( c( SRP028599.Mt_Sm_RbmL_Nodule_IZ_R1,
                          SRP028599.Mt_Sm_RbmL_Nodule_IZ_R2,
                          SRP028599.Mt_Sm_RbmL_Nodule_IZ_R3) ),
            
            ZIII = mean( c(SRP028599.Mt_Sm_RbmL_Nodule_ZIII_R1,
                           SRP028599.Mt_Sm_RbmL_Nodule_ZIII_R2,
                           SRP028599.Mt_Sm_RbmL_Nodule_ZIII_R3) )
        )  %>%
        dplyr::select(genes, Roots, 'N10d', FI, FIId, FIIp, IZ, ZIII) %>%
        tibble::column_to_rownames("genes")
    
    SRP_scaled <- t(dynutils::scale_quantile( t(SRP_complete_LCM) ) )
    
    ( p_SRP <- my_complexheatmap(scaled_data_frame = SRP_scaled,
                                 arg_column_title = SRP) )
    
    if (show_names == TRUE) {
        
        plots_height <- 5 + ( 0.2 * nrow(p_SRP) )
        
    } else if (show_names == FALSE) {
        
        plots_height <-  0.03 * nrow(p_SRP) 
        
    }
    
    if(opts$image_format == "png") {
        
        png( paste0(opts$out_images, "/",
                    plot_name, "_",
                    SRP, "_Roux_2014_LCM_of_nodules.png" ),
             width = (ncol(p_SRP)/2),
             height = plots_height,
             pointsize = 12,
             res = 300, units = "in")
        ComplexHeatmap::draw(p_SRP)
        dev.off()
        
    } else if (opts$image_format == "svg") {
        
        svg(  paste0(opts$out_images, "/",
                     plot_name, "_",
                     SRP, "_Roux_2014_LCM_of_nodules.svg" ),
              width = (ncol(p_SRP)/2),
              height = plots_height,
              pointsize = 12 )
        ComplexHeatmap::draw(p_SRP)
        dev.off()
        
    } else {
        print("--image_format must be 'png' or 'svg'!")
        quit()
    }
    
}

## Schiessl et al., 2019
atlast_v2_SRP212693_opt3 <- function(genes_names = genes_names,
                                     plot_name = plot_name,
                                     SRP = "SRP212693",
                                     show_names = opts$show_gene_names) {
    
    atlas_v2_sub <- atlas_v2 %>%
        dplyr::filter(locus_tag %in% genes_names) %>%
        column_to_rownames("locus_tag")
    
    print( paste ("A total of", nrow(atlas_v2_sub), "genes IDs, from",
                  length(genes_names), "had correspondent gene IDs in the dataset"))
    
    # Select the data from Rox et al 2014, based on the project number
    SRP_complete <- atlas_v2_sub[, grepl( colnames( atlas_v2_sub ), pattern = SRP ) ]
    
    # Select the datasets that are related with the treatment of the WT.
    SRP_WT <- SRP_complete[, grepl( colnames( SRP_complete ),
                                    pattern = "WT" ) ]
    SRP_WT <- SRP_WT[, grepl( colnames( SRP_WT ),
                              pattern = "Smeliloti" ) ]
    
    SRP_WT$genes <- rownames(SRP_WT)
    SRP_WT <- SRP_WT %>%
        rowwise() %>% 
        dplyr::mutate( "0h" = mean( c(SRP212693.WT_0h_Smeliloti_R1,
                                      SRP212693.WT_0h_Smeliloti_R2,
                                      SRP212693.WT_0h_Smeliloti_R3,
                                      SRP212693.WT_0h_Smeliloti_R4) ),
                       
                       "2h" = mean( c(SRP212693.WT_2h_Smeliloti_R1,
                                      SRP212693.WT_2h_Smeliloti_R2,
                                      SRP212693.WT_2h_Smeliloti_R3,
                                      SRP212693.WT_2h_Smeliloti_R4,
                                      SRP212693.WT_2h_Smeliloti_R5,
                                      SRP212693.WT_2h_Smeliloti_R6
                       ) ),  
                       
                       "4h" = mean( c(SRP212693.WT_4h_Smeliloti_R1,
                                      SRP212693.WT_4h_Smeliloti_R2,
                                      SRP212693.WT_4h_Smeliloti_R3,
                                      SRP212693.WT_4h_Smeliloti_R4,
                                      SRP212693.WT_4h_Smeliloti_R5) ),
                       
                       "8h" = mean( c(SRP212693.WT_8h_Smeliloti_R1,
                                      SRP212693.WT_8h_Smeliloti_R2,
                                      SRP212693.WT_8h_Smeliloti_R3,
                                      SRP212693.WT_8h_Smeliloti_R4) ),
                       
                       "10h" = mean( c(SRP212693.WT_10h_Smeliloti_R1,
                                       SRP212693.WT_10h_Smeliloti_R2,
                                       SRP212693.WT_10h_Smeliloti_R3,
                                       SRP212693.WT_10h_Smeliloti_R4) ),
                       
                       "12h" = mean( c(SRP212693.WT_12h_Smeliloti_R1,
                                       SRP212693.WT_12h_Smeliloti_R2,
                                       SRP212693.WT_12h_Smeliloti_R3,
                                       SRP212693.WT_12h_Smeliloti_R4) ),
                       
                       "14h" = mean( c(SRP212693.WT_14h_Smeliloti_R1,
                                       SRP212693.WT_14h_Smeliloti_R2,
                                       SRP212693.WT_14h_Smeliloti_R3,
                                       SRP212693.WT_14h_Smeliloti_R4,
                                       SRP212693.WT_14h_Smeliloti_R5) ), 
                       
                       "16h" = mean( c(SRP212693.WT_16h_Smeliloti_R1,
                                       SRP212693.WT_16h_Smeliloti_R2,
                                       SRP212693.WT_16h_Smeliloti_R3,
                                       SRP212693.WT_16h_Smeliloti_R4,
                                       SRP212693.WT_16h_Smeliloti_R5) ), 
                       
                       "24h" =  mean( c(SRP212693.WT_24h_Smeliloti_R1,
                                        SRP212693.WT_24h_Smeliloti_R2,
                                        SRP212693.WT_24h_Smeliloti_R3,
                                        SRP212693.WT_24h_Smeliloti_R4) ),
                       
                       "36h" = mean( c(SRP212693.WT_36h_Smeliloti_R1,
                                       SRP212693.WT_36h_Smeliloti_R2,
                                       SRP212693.WT_36h_Smeliloti_R3,
                                       SRP212693.WT_36h_Smeliloti_R4) ),
                       
                       "48h" = mean( c(SRP212693.WT_48h_Smeliloti_R1,
                                       SRP212693.WT_48h_Smeliloti_R2,
                                       SRP212693.WT_48h_Smeliloti_R3,
                                       SRP212693.WT_48h_Smeliloti_R4,
                                       SRP212693.WT_48h_Smeliloti_R5,
                                       SRP212693.WT_48h_Smeliloti_R6 ) ),
                       
                       "72h" = mean( c(SRP212693.WT_72h_Smeliloti_R1,
                                       SRP212693.WT_72h_Smeliloti_R2,
                                       SRP212693.WT_72h_Smeliloti_R3,
                                       SRP212693.WT_72h_Smeliloti_R4) ),
                       
                       "96h" = mean( c(SRP212693.WT_96h_Smeliloti_R1,
                                       SRP212693.WT_96h_Smeliloti_R2,
                                       SRP212693.WT_96h_Smeliloti_R3,
                                       SRP212693.WT_96h_Smeliloti_R4) ), 
                       
                       "5d" = mean( c(SRP212693.WT_120h_Smeliloti_R1,
                                      SRP212693.WT_120h_Smeliloti_R2,
                                      SRP212693.WT_120h_Smeliloti_R3,
                                      SRP212693.WT_120h_Smeliloti_R4) ),
                       
                       "7d" = mean( c(SRP212693.WT_168h_Smeliloti_R1,
                                      SRP212693.WT_168h_Smeliloti_R2,
                                      SRP212693.WT_168h_Smeliloti_R3,
                                      SRP212693.WT_168h_Smeliloti_R4) ) ) %>% 
        select( c("genes", "0h", "2h", "4h", "8h", "10h", "12h", "14h", "24h", "36h", "48h", "72h", "5d", "7d") ) %>% 
        tibble::column_to_rownames("genes")
    
    SRP_scaled <- t(dynutils::scale_quantile( t(SRP_WT) ) )
    
    ( p_SRP <- my_complexheatmap(SRP_scaled,
                                 arg_column_title = "SRP212693 - Schiessl et al., 2019.") )
    
    if (show_names == TRUE) {
        
        plots_height <- 5 + ( 0.2 * nrow(p_SRP) )
        
    } else if (show_names == FALSE) {
        
        plots_height <-  0.03 * nrow(p_SRP) 
        
    }
    
    if(opts$image_format == "png") {
        
        png( paste0( opts$out_images, "/",
                     plot_name, "_", 
                     SRP, "_Schiessl_2019_complete_opt3.png" ),
             width = ncol(p_SRP)/2,
             height = plots_height,
             pointsize = 12,
             res = 300, units = "in")
        ComplexHeatmap::draw(p_SRP)
        dev.off()
        
    } else if (opts$image_format == "svg") {
        
        svg(  paste0(opts$out_images, "/",
                     plot_name, "_", 
                     SRP, "_Schiessl_2019_complete_opt3.svg" ),
              width = ncol(p_SRP)/2,
              height = plots_height,
              pointsize = 12 )
        ComplexHeatmap::draw(p_SRP)
        dev.off()
        
    } else {
        print("--image_format must be 'png' or 'svg'!")
        quit()
    }
    
    
}