suppressMessages( require(monocle3) )
suppressMessages( library(corrplot) )
suppressMessages( library(ggplot2) )

## Correlation with combined cluster
cds <- readRDS( "rds_files/batched_integrated_clustered_complete_dataset.rds")

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=clusters(cds)[colnames(cds)])

agg_mat <- aggregate_gene_expression(cds, NULL, cell_group_df, scale_agg_values = F)
agg_mat <- as.matrix(agg_mat)
agg_mat <- agg_mat[rowSums(agg_mat) > 0, ]

M = cor(agg_mat)
corrplot.mixed(M,
               lower = 'number',
               upper = 'pie',
               order = 'hclust',
               number.digits = 2,
               number.cex = 0.8,
               tl.pos = "lt",
               tl.col = "black",
               tl.offset=0.5,
               tl.srt = 0)

png("images/correlation_among_clusters_whole_dataset.png",
    width = 40, 
    height = 20,
    res = 300,
    units = "cm")
corrplot.mixed(M,
               lower = 'number',
               upper = 'pie',
               order = 'hclust',
               number.digits = 2,
               number.cex = 0.8,
               tl.pos = "lt",
               tl.col = "black",
               tl.offset=0.5,
               tl.srt = 0)
dev.off()

write.csv(M, "correlation_among_clusters_whole_dataset.csv")

## Correlation by cluster by and timepoints
colData(cds)$clusters <- as.numeric(clusters(cds))
colData(cds)$c_time <- paste0( "C", colData(cds)$clusters,"-", colData(cds)$timepoint )

cell_group_df <- tibble::tibble( cell=row.names( colData(cds) ), 
                                 cell_group= colData(cds)$c_time )

agg_mat <- aggregate_gene_expression(cds, NULL, cell_group_df,scale_agg_values = F)
agg_mat <- as.matrix(agg_mat)
agg_mat <- agg_mat[rowSums(agg_mat) > 0, ]

M = cor(agg_mat)

COL2( diverging = c("RdBu",
                    "BrBG",
                    "PiYG",
                    "PRGn",
                    "PuOr",
                    "RdYlBu"),
      n = 200)

corrplot(M,
         method = 'square',
         order = 'hclust',
         number.digits = 2,
         col = COL2('RdYlBu', 10),
         tl.col = 'black')

png("images/correlation_among_clusters_timepoint_whole_dataset.png",
    width = 60, 
    height = 50,
    res = 300,units = "cm")
corrplot(M,
         method = 'square',
         order = 'hclust',
         number.digits = 2,
         col = COL2('RdYlBu', 10),
         tl.col = 'black')
dev.off()

write.csv(M, "correlation_among_clusters_timepoint_whole_dataset.csv")
