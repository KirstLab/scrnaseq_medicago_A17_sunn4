"Usage: number_of_cells_per_cluster_over_time.R (--rds <rds>) (--prefix <prefix>)
--rds    RDS file containing the clustered dataset.
--prefix   Prefix to identify the output files
number_of_cells_per_cluster_over_time.R -h | --help  show this message.
" -> doc

set.seed(1407)

outfile <- "logs/number_of_cells_per_cluster_over_time.out" # File name of output log
#Check its existence
if ( file.exists(outfile) ) {
    #Delete file if it exists
    file.remove(outfile)
}

my_log <- file(outfile) 
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

suppressMessages( require(docopt) )
suppressMessages( require(tidyverse) )
suppressMessages( require(monocle3) )

opts <- docopt(doc)
cds <- readRDS(opts$rds)
colData( cds )$cluster <- monocle3::clusters( cds )

cds_meta <- as.data.frame( colData( cds ) )

counts_per_cluster <- cds_meta %>% 
    dplyr::count(cluster, Group, timepoint)

write_excel_csv(counts_per_cluster,
                file = paste0( opts$prefix, "/",
                              "raw_number_of_cells_per_cluster_by_genotype_and_timepoint.csv") )

counts_per_group <- cds_meta %>% 
    dplyr::count(Group, timepoint)

write_excel_csv(counts_per_cluster, 
                file = paste0( opts$prefix, "/",
                              "raw_number_of_cells_per_genotype_by_timepoint.csv") )

p <- ggplot(data=counts_per_cluster,
            aes(x=timepoint, y=n, group=Group)) +
    geom_line( aes( color = Group ) )  + 
    geom_point(aes(color = Group)) +
    facet_wrap(~cluster, scales = "free") +
    theme_bw() + 
    theme(strip.text = element_text(size = 15),
          axis.text = element_text(size = 15),
          axis.text.x=element_text(angle = 0, hjust = 0.5, vjust = 0),
          plot.title = element_blank(),
          legend.title=element_blank(),
          legend.text = element_text(size = 17) ) +
    scale_color_manual( values=c("#fc9003","#031cfc") )

ggplot2::ggsave( 
    paste0("images/", 
           opts$prefix,"/",
           "_raw_number_of_cells_per_clusters.png"),
    p,
    height=40,
    width=60,
    units="cm",
    dpi = 300,
    bg = "#FFFFFF")

ggplot2::ggsave(    
    paste0("images/", 
           opts$prefix,"/",
           "_raw_number_of_cells_per_clusters.svg"),
    p,
    height=40,
    width=60,
    units="cm",
    dpi = 300,
    bg = "#FFFFFF")

p <- ggplot(data=counts_per_cluster,
            aes(x=timepoint, y=n, group=Group)) +
    geom_line( aes( color = Group ) )  + 
    geom_point(aes(color = Group)) +
    facet_wrap(~cluster, scales = "free") +
    theme_bw() + 
    ylim(0, max(counts_per_cluster$n) + 5) +
    theme(strip.text = element_text(size = 15),
          axis.text = element_text(size = 15),
          axis.text.x=element_text(angle = 0, hjust = 0.5, vjust = 0),
          plot.title = element_blank(),
          legend.title=element_blank(),
          legend.text = element_text(size = 17) ) +
    scale_color_manual( values=c("#fc9003","#031cfc") )

ggplot2::ggsave(    
    paste0("images/", 
           opts$prefix,"/",
           "_raw_number_of_cells_per_clusters_fixed_scale.png"),
    p,
    height=40,
    width=60,
    units="cm",
    dpi = 300,
    bg = "#FFFFFF")

ggplot2::ggsave(
    paste0("images/", 
           opts$prefix,"/",
           "_raw_number_of_cells_per_clusters_fixed_scale.svg"),
    p,
    height=40,
    width=60,
    units="cm",
    dpi = 300,
    bg = "#FFFFFF")

counts_per_cluster_wide <- cds_meta %>% 
    dplyr::count(cluster, Group, timepoint) %>% 
    tidyr::pivot_wider(names_from = c(Group, timepoint), 
                       values_from = n) 

write.table(counts_per_cluster_wide,
            paste0(opts$prefix,"/",
                   "_raw_cells_per_cluster_integrated_29_C_cds_clustered.csv"),
            col.names = T,
            row.names = F,
            quote = F,
            sep = ",")

## After normalizing the number of cells 
### For each cluster and genotype, the counts on each time is divided by the total of cells of that genotype in the timepoint. Then the normalized value id multiplied by 1000 to create the counts per thousand metric.
WT <- counts_per_cluster%>%
    dplyr::filter(Group == "A17-WT")
WT_G <- counts_per_group%>%
    dplyr::filter(Group == "A17-WT")

sunn4 <- counts_per_cluster%>%
    dplyr::filter(Group == "Sunn-4")
sunn4_G <- counts_per_group%>%
    dplyr::filter(Group == "Sunn-4")

WT$CT <- ifelse(
    WT$timepoint == "0h",  round( ( (WT$n / WT_G[WT_G$timepoint == "0h", "n"] ) * 1000), 2),
    ifelse( WT$timepoint == "24h", round( ( (WT$n / WT_G[WT_G$timepoint == "24h", "n"] ) * 1000), 2),
            ifelse( WT$timepoint == "48h",  round( ( (WT$n / WT_G[WT_G$timepoint == "48h", "n"] ) * 1000), 2),
                    ifelse( WT$timepoint == "96h", round( ( (WT$n / WT_G[WT_G$timepoint == "96h", "n"] ) * 1000), 2),
                            "NA"
                    ) 
            )
    )
)

sunn4$CT <- ifelse(
    sunn4$timepoint == "0h",  round( ( (sunn4$n / sunn4_G[sunn4_G$timepoint == "0h", "n"] ) * 1000), 2),
    ifelse( sunn4$timepoint == "24h", round( ( (sunn4$n / sunn4_G[sunn4_G$timepoint == "24h", "n"] ) * 1000), 2),
            ifelse( sunn4$timepoint == "48h",  round( ( (sunn4$n / sunn4_G[sunn4_G$timepoint == "48h", "n"] ) * 1000), 2),
                    ifelse( sunn4$timepoint == "96h", round( ( (sunn4$n / sunn4_G[sunn4_G$timepoint == "96h", "n"] ) * 1000), 2),
                            "NA"
                    ) 
            )
    )
)

counts_per_cluster <- rbind(sunn4, WT) %>%
    dplyr::arrange(cluster, Group, timepoint)
counts_per_cluster$CT <- as.numeric(counts_per_cluster$CT)

write_excel_csv(counts_per_cluster,
                file = paste0( opts$prefix, "/",
                              "number_of_cells_per_cluster_by_genotype_and_timepoint_including_CT.csv") )

write_excel_csv(counts_per_cluster,
                file = paste0( opts$prefix, "/",
                              "raw_number_of_cells_per_cluster_by_genotype_and_timepoint.csv") )

p <- ggplot(data = counts_per_cluster,
            aes(x=timepoint, y=CT, group=Group)) +
    geom_line( aes( color = Group ) )  + 
    geom_point(aes(color = Group)) +
    facet_wrap(~cluster, scales = "free") +
    theme_bw() + 
    ylim(0, max(counts_per_cluster$CT) + 5) +
    theme(strip.text = element_text(size = 15),
          #axis.title = element_text(size = 16),
          axis.text = element_text(size = 15),
          axis.text.x=element_text(angle = 0, hjust = 0.5, vjust = 0),
          plot.title = element_blank(),
          legend.title=element_blank(),
          legend.text = element_text(size = 17) ) +
    scale_color_manual( values=c("#fc9003","#031cfc") )

p

ggplot2::ggsave(
    paste0("images/",
           opts$prefix, "/",
           "_CT_number_of_cells_per_clusters.png"),
    p,
    height=40,
    width=60,
    units="cm",
    dpi = 300,
    bg = "#FFFFFF")

ggplot2::ggsave(
    paste0("images/",
           opts$prefix, "/",
           "_CT_number_of_cells_per_clusters.svg"),
    p,
    height=40,
    width=60,
    units="cm",
    dpi = 300,
    bg = "#FFFFFF")

closeAllConnections()