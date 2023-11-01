set.seed(1407)
##########################################
### Trajectory Differential Expression ###
##########################################

outfile <- "logs/TradeSeq.out"
if ( file.exists(outfile) ) {
    file.remove(outfile)
}

my_log <- file(outfile)
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

require(tradeSeq)
require(slingshot)
require(tidyverse)

my_write_csv <- function( my_obj = my_obj,
                          name = name,
                          col.names = T,
                          row.names = F,
                          quote = F,
                          sep = "\t") {
    
    write.table(my_obj,
                name,
                col.names = T,
                row.names = F,
                quote = F,
                sep = "\t")
}

my_ggsave <- function( name = name,
                       plot = plot,
                       height = 12,
                       width = 16) {
    
    ggplot2::ggsave(
        filename = name,
        plot = plot,
        height = height,
        width = width,
        units = "in",
        bg = "#FFFFFF", 
        dpi = 300)
    
}

annot_names <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv")
annot_summary <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv")
colnames(annot_names)[1] <- "gene_id"
colnames(annot_summary)[1] <- "gene_id"

annot_summary <- merge( annot_summary,
                        annot_names,
                        by = "gene_id",
                        all.x = T)

annot_summary <- annot_summary[, c(1, 14, 2:13)]

## Using TradeSeq
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 6

sling_res <- readRDS("RECLUSTERING/CORTEX_NOD/rds_file_subset/cortex_nodules_slingshot_trajectory.rds")
cds_subset <- readRDS("RECLUSTERING/CORTEX_NOD/rds_file_subset/clusters_cortex_nodule_selected_for_trajectory.rds")

plot_cells(cds_subset,
           cell_size = 0.75,
           group_label_size = 8) 

sds2 <- SlingshotDataSet(sling_res)

counts_mt <- as.matrix(assay(cds_subset))
counts_mt <- counts_mt[rowSums(counts_mt) > 0, ]

icMat <- evaluateK(counts = counts_mt,
                   sds = sds2,
                   verbose = F)

sce_fitted <- tradeSeq::fitGAM(counts = counts_mt,
                               sds = sds2,
                               nknots = 6,
                               verbose = T,
                               parallel = TRUE,
                               BPPARAM = BPPARAM)

system( "mkdir -p rds_files")
saveRDS(sce_fitted,
        file = "rds_files/TRADEseq_cortex_nodules.rds")

## Association test to check if each genes is DEG within the trajectories of each lineage.
foldername4 <- "tradeseq_markers_cortex_nod/"
system( paste0("mkdir -p ", foldername4) )

feat_importances <- tradeSeq::associationTest(sce_fitted,
                                              lineages = TRUE)

# Adjust the p-value detected in each lineage to correct for multiple test (by using FDR).
feat_importances$fdr_1 <- stats::p.adjust(feat_importances$pvalue_1,
                                          method = "fdr",
                                          n = length(feat_importances$pvalue_1))

feat_importances$genes <- rownames(feat_importances)

DEGs_within_traj_of_lineage_1 <- feat_importances %>%
    dplyr::filter(fdr_1 < 0.001) %>%
    dplyr::select(waldStat_1, df_1, pvalue_1, fdr_1, genes) %>%
    dplyr::rename(gene_id = genes) %>%
    dplyr::arrange( desc( waldStat_1 ) )

DEGs_within_traj_of_lineage_1 <- merge(DEGs_within_traj_of_lineage_1,
                                       annot_summary, 
                                       by = "gene_id",
                                       all.x = T )

write.table( DEGs_within_traj_of_lineage_1,
              paste0(foldername4, "/tradeseq_DEGs_lineage1.tsv"), col.names = T, row.names = F, sep = "\t" )

## Select among the known RNS 
RNS_genes <- vroom::vroom("selected_markers/RNS_genes_curated_list.tsv")

DEG_lineage_1 <- DEGs_within_traj_of_lineage_1 %>%
    dplyr::filter(gene_id %in% RNS_genes$gene_id)

### Heatmap ###
sling_pt_heat <- slingshot::slingPseudotime(sling_res) %>% 
    as.data.frame()

#### Lineage1 - DEGs
Lin1 <- data.frame(gene_id =  DEG_lineage_1$gene_id )

Lin1 <- merge( Lin1,
               annot_summary, 
               by = "gene_id",
               all.x = T )

pst.ord <- order(sling_pt_heat$Lineage1, na.last = NA)
heatdata <- counts_mt
heatdata <- heatdata[rownames(heatdata) %in% Lin1$gene_id, pst.ord]

cells_order_in_traject <- tibble(cells = colnames(heatdata))

# Replaces the gene ID by its acronym if it is available.
for( i in 1:nrow(heatdata) ) {
    
    gene <- rownames(heatdata)[i]
    
    rownames(heatdata)[i] <- ifelse(
        
        rownames(heatdata)[i] %in% Lin1$gene_id,
        unique(Lin1[Lin1$gene_id %in% rownames(heatdata)[i], "acronym"]),
        Lin1$gene_id)
    
    rownames(heatdata)[i] <- ifelse(
        
        is.na(rownames(heatdata)[i]), 
        gene,
        rownames(heatdata)[i] )
    
}

heatdata_log <- log1p(heatdata)
heatdata_scaled <- dynutils::scale_quantile(heatdata_log) 

p_heatmap <- pheatmap::pheatmap(heatdata_scaled,
                                show_rownames = T,
                                show_colnames = F,
                                cluster_rows = T,
                                cluster_cols = FALSE,
                                scale="none",
                                color = viridis::viridis(n = 20, option = "D"),
                                clustering_method="ward.D2",
                                fontsize=15)

foldername5 <- "tradeseq_markers_cortex_nod/images/"

system( paste0("mkdir -p ", foldername5) )
my_ggsave(
    paste0(foldername5, "TRADESEQ_RNS_genes_lineage1_trajectory.svg"), 
    p_heatmap,
    height= nrow(heatdata_scaled)*0.2,
    width=8)

my_ggsave(
    paste0(foldername5, "TRADESEQ_RNS_genes_lineage1_trajectory.png"), 
    p_heatmap,
    height= nrow(heatdata_scaled)*0.2,
    width=8)

## Binned heatmap
N_GROUPS = 50
cells_per_group <- floor( ncol(heatdata) / (N_GROUPS) ) 
extra_cells <- ncol(heatdata) - (N_GROUPS * cells_per_group)

colnames_heatmap <- c( 
    
    rep( 1:extra_cells,
         each =  cells_per_group + 1),
    
    rep( (extra_cells + 1) : N_GROUPS, each = cells_per_group )
    
)

print(table(colnames_heatmap))

cells_per_bin <- table(colnames_heatmap)

colnames(heatdata) <- colnames_heatmap

averaged_heatdata <- data.frame(
    row.names = make.names( rownames(heatdata), unique = T ) )
for(c in unique(colnames_heatmap) ) {
    
    sub <- heatdata[, colnames(heatdata) == c]
    sub2 <- as.data.frame(rowMeans(sub))
    
    averaged_heatdata <- cbind(averaged_heatdata, sub2)
}

averaged_heatdata <- as.matrix(averaged_heatdata)
colnames(averaged_heatdata) <- paste("bin", 1:N_GROUPS, sep = "-")
averaged_heatdata_log <- log1p(averaged_heatdata)
averaged_heatdata_scaled <- t(dynutils::scale_quantile(t(averaged_heatdata_log)))

bined_heatmap <- pheatmap::pheatmap(averaged_heatdata_scaled,
                                    show_rownames = T,
                                    show_colnames = F,
                                    cluster_rows = T,
                                    cluster_cols = FALSE,
                                    scale="none",
                                    color = viridis::viridis(n = 20, option = "D"),
                                    clustering_method="ward.D2",
                                    fontsize=15)

my_ggsave(paste0(foldername5, "TRADESEQ_RNS_genes_lineage1_trajectory_50BINs.svg"), 
          bined_heatmap,
          height= nrow(heatdata_scaled)*0.2,
          width=8)

my_ggsave(paste0(foldername5, "TRADESEQ_RNS_genes_lineage1_trajectory_50BINs.png"), 
          bined_heatmap,
          height= nrow(heatdata_scaled)*0.2,
          width=8)