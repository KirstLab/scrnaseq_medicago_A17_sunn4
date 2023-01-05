set.seed(1407)

suppressMessages( require( monocle3 ) )
suppressMessages( require( Seurat ) )
suppressMessages( require( dynwrap ) )
suppressMessages( require( tidyverse ) )
suppressMessages( require( dynplot ) )
suppressMessages( require( assertthat ) )
suppressMessages( require( grDevices ) )
suppressMessages( require( slingshot ) )
suppressMessages( require( RColorBrewer ) )
suppressMessages( require( ComplexHeatmap ) )
suppressMessages( require( SingleCellExperiment ) )
suppressMessages( require( dynfeature ) )
suppressMessages( require( ggthemes ) )
suppressMessages( require( tradeSeq ) )
suppressMessages( require(viridis) )
suppressMessages( require(dynutils) )

source("dynplot/R/my_plot_heatmap.R")
source("dynplot/R/add_density_coloring.R")
source("dynplot/R/dummy_proofing.R")
source("dynplot/R/expect_ggplot.R")            
source("dynplot/R/is_colour_vector.R")
source("dynplot/R/linearise_cells.R")
source("dynplot/R/milestone_palette.R")
source("dynplot/R/mix_colors.R")
source("dynplot/R/my_plot_heatmap.R")
source("dynplot/R/optimize_order.R")
source("dynplot/R/package.R")
source("dynplot/R/plot_edge_flips.R")
source("dynplot/R/plot_linearised_comparison.R")
source("dynplot/R/plot_strip.R")              
source("dynplot/R/project_waypoints.R")  
source("dynplot/R/add_milestone_coloring.R")
source("dynplot/R/theme_clean.R")
source("dynplot/R/add_cell_coloring.R")
source("dynplot/R/plot_topology.R")     
source("dynplot/R/plot_onedim.R")
source("dynplot/R/plot_graph.R") 
source("dynplot/R/plot_dendro.R")
source("dynplot/R/plot_dimred.R")

cds <- readRDS("RECLUSTERING/CORTEX_NOD_sunn4/rds_file_subset/medicago_integrated_subset_cortex_nodule.rds")

plot_cells(cds,
           label_cell_groups=T,
           graph_label_size=1.5,
           cell_size = 1,
           group_label_size = 7) + 
    facet_wrap(~ timepoint, nrow = 1)

plot_cells(cds,
           label_cell_groups=T,
           graph_label_size=1.5,
           cell_size = 1,
           group_label_size = 7) + 
    facet_wrap(~Group + timepoint, nrow = 2)

### Removing the subcluster that won't be part of the trajectory
to_remove <- c(5, 6)

colData(cds)$cluster <- monocle3::clusters(cds)
cds_sub_meta <- as.data.frame( colData(cds) )
cds_sub_meta <- cds_sub_meta[!cds_sub_meta$cluster %in% to_remove, ]

to_remove_cells <- as.character(rownames(cds_sub_meta))
cds <- cds[, colnames(cds) %in% to_remove_cells ]

plot_cells(cds,
           label_cell_groups=T,
           graph_label_size=1.5,
           cell_size = 1,
           group_label_size = 7) + 
    facet_wrap(~ timepoint, nrow = 1)

############################################################################
### Using dynverse to identify the importance of genes in the trajectory ###
############################################################################

## Preparing the data to the format accepted by dynverse. First we draw the trajectory using the docker image of slingshot, then we find the importance of each gene to explain the trajectory.
colData(cds)$clusters = monocle3::clusters(cds)

COUNTS <- as.matrix(assay(cds))
COUNTS <- COUNTS[rowSums(COUNTS) > 0, ]

DATA_counts <- Seurat::NormalizeData(COUNTS)

expressed_genes <- rownames(DATA_counts)

object_expression <- Matrix::t(as(DATA_counts, 'sparseMatrix'))
object_counts <-  Matrix::t(as(COUNTS, 'sparseMatrix'))

dataset_inf <- dynwrap::wrap_expression(
    
    counts = object_counts,
    expression =  object_expression
    
)

end_cells <- NULL

sc_meta <- as.data.frame(colData(cds))
sc_meta <- sc_meta %>%
    mutate(cells = rownames(.))

sc_meta_cluster <- data.frame(cell_id = sc_meta$cells,
                              group_id = sc_meta$clusters)

dimred <- SingleCellExperiment::reducedDims(cds)[["PCA"]]

STARTING_CLUSTER = 2
start_cells <- sc_meta[sc_meta$clusters == STARTING_CLUSTER, ]
start_cells <- as.character(start_cells$cells)

# fetch newest version of the slingshot within dynverse 
## Docker is required here!
method_id <- paste0("dynverse/ti_", "slingshot", ":latest")
methods_selected <- create_ti_method_container(method_id)

dataset <- add_prior_information(
    
    dataset_inf,
    start_id = start_cells,
    groups_id = sc_meta_cluster,
    dimred = dimred
)

sds <- infer_trajectory( dataset,
                         methods_selected(),
                         verbose = T,
                         give_priors = c(
                             "start_id",
                             "groups_id",
                             "dimred"
                         ) )

foldername1 <- "RECLUSTERING/CORTEX_NOD_sunn4_cortex_meristem_trajectory/"
system( paste0("mkdir -p ", foldername1, "/rds_file_subset") )
saveRDS(sds, 
        paste0(foldername1,
               "rds_file_subset/sunn4_trajectory_cortex_nodule.rds")
)

( p <- dynplot::plot_dimred(sds,
                            grouping = sc_meta_cluster,
                            color_density = "grouping") +
        ggtitle("Cell grouping") )

foldername2 <- "RECLUSTERING/CORTEX_NOD_sunn4_cortex_meristem_trajectory/dynverse_markers_cortex_nod/images/"
system( paste0( "mkdir -p ", foldername2) )
ggplot2::ggsave(
    paste0(foldername2, "pericycle_trajectory_dynverse_grouping.svg"),
    p,
    height=12,
    width=16,
    units="in",
    bg = "#FFFFFF", 
    dpi = 300)

( p1 <- dynplot::plot_dendro(sds,
                             grouping = sc_meta_cluster) +
        ggtitle("Trajectory") )

ggplot2::ggsave(
    paste0(foldername2,
           "pericycle_trajectory_dynverse_dendro.svg"),
    p1,
    height=12,
    width=16,
    units="in",
    bg = "#FFFFFF", 
    dpi = 300)

foldername3 <- "RECLUSTERING/CORTEX_NOD_sunn4_cortex_meristem_trajectory/dynverse_markers_cortex_nod/"
system( paste0( "mkdir -p ", foldername3 ) )

## Reads the files from the genome to add the annotation to the genes lists generated below.
annot_names <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv")
annot_summary <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv")
colnames(annot_names)[1] <- "gene_id"
colnames(annot_summary)[1] <- "gene_id"

annot_summary <- merge( annot_summary,
                        annot_names,
                        by = "gene_id",
                        all.x = T)

annot_summary <- annot_summary[, c(1, 14, 2:13, 15:16)]

## Overall importance
feat_importances <- dynfeature::calculate_overall_feature_importance( 
    trajectory = sds,
    expression_source = dataset)

colnames(feat_importances)[1] <- "gene_id"
feat_importances_annot <- merge( feat_importances,
                                 annot_summary, 
                                 by = "gene_id",
                                 all.x = T )

feat_importances_annot <- feat_importances_annot %>%
    dplyr::arrange( desc(importance) )

feat_importances_annot$rank <- seq( 1 : nrow(feat_importances_annot) )

write.table(feat_importances_annot,
            paste0(foldername3, "cortex_nodule_trajectory_dynverse_whole_trajectory.tsv"),
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t")

feat_importances_annot_sub <- feat_importances_annot[c(1:100), ]
( top_100_heatmap <- my_plot_heatmap(
    sds,
    expression_source = dataset,
    features_oi = as.character(feat_importances_annot_sub$gene_id),
    scale = T) )

ggplot2::ggsave(
    paste0(foldername2, "dynfeatures_top_100_whole_trajectory.png" ),
    top_100_heatmap,
    height= nrow(feat_importances_annot_sub)*0.2,
    width=20,
    units="in",
    bg = "#FFFFFF", 
    dpi = 300)

## extracts features that are specifically up-regulated or down-regulated in a specific branch
feat_importances_branch <- dynfeature::calculate_branch_feature_importance( 
    trajectory = sds,
    expression_source = dataset)

feat_importances_branch_sub <- feat_importances_branch %>%
    dplyr::group_by(to) %>% 
    dplyr::slice_max(order_by = importance, n = 100)

colnames(feat_importances_branch_sub)[1] <- "gene_id"
feat_importances_branch_sub <- merge( feat_importances_branch_sub,
                                      annot_summary, 
                                      by = "gene_id",
                                      all.x = T )

feat_importances_branch_sub <- feat_importances_branch_sub %>%
    dplyr::arrange( desc(importance) )

feat_importances_branch_sub$rank <- seq( 1 : nrow(feat_importances_branch_sub) )

write.table(feat_importances_branch_sub,
            paste0(foldername3, "cortex_nodule_trajectory_dynverse_top_100_each_branch.tsv"),
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t")

######################
### Using TradeSeq ###
######################
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 7

sce <- SingleCellExperiment(assays = List(counts = COUNTS ) )

reducedDims(sce) <- SimpleList(PCA = SingleCellExperiment::reducedDims(cds)[["PCA"]])
colData(sce)$clusters <- monocle3::clusters(cds)

sds_trade <- slingshot::slingshot(
    sce,
    clusterLabels = "clusters",
    start.clus = STARTING_CLUSTER,
    reducedDim = 'PCA'
)

slingshot::slingLineages(sds_trade)

sds2 <- SlingshotDataSet(sds_trade)
counts <- as.matrix(assays(sds_trade)$counts)

# To evaluated the number of knots to be used on the modeling phase
icMat <- evaluateK(counts = counts, 
                   sds = sds2,
                   k = 3:10, 
                   nGenes = 200,
                   verbose = T)

sce_fitted <- fitGAM(counts = counts,
                     sds = sds2,
                     nknots = 6,
                     verbose = T,
                     parallel = TRUE,
                     BPPARAM = BPPARAM)

saveRDS(sce_fitted,
        file = paste0("RECLUSTERING/CORTEX_NOD_sunn4_cortex_meristem_trajectory/rds_file_subset/TRADEseq_cortex_nodules.rds")
)

feat_importances <- tradeSeq::associationTest(sce_fitted)

feat_importances$fdr <- stats::p.adjust(feat_importances$pvalue,
                                        method = "fdr",
                                        n = length(feat_importances$pvalue))

feat_importances$genes <- rownames(feat_importances)
feat_importances2 <- feat_importances[ ,c(ncol(feat_importances), 1:( ncol(feat_importances) - 1)) ]

foldername4 <- "RECLUSTERING/CORTEX_NOD_sunn4_cortex_meristem_trajectory/tradeseq_markers_cortex_nod/"
system( paste0("mkdir -p ", foldername4) )
write.table( feat_importances2,
             paste0(foldername4, "TRADEseq_cortex_nodules.tsv"),
             col.names = T,
             row.names = F,
             sep = "\t",
             quote = F
)

### Crossing the results of dynfeatures and tradseq ###
tradseq <- vroom::vroom(
    paste0(foldername4, "TRADEseq_cortex_nodules.tsv"),
    delim = "\t")

dynverse <- vroom::vroom(
    paste0(foldername3,
           "cortex_nodule_trajectory_dynverse_whole_trajectory.tsv"),
    delim = "\t", na = "NA")

tradseq <- tradseq %>%
    dplyr::rename(gene_id = genes) %>%
    dplyr::filter(fdr < 0.01)

combined <- merge(tradseq, dynverse, by = "gene_id")
combined <- combined %>%
    dplyr::arrange( desc( waldStat ) )

write.table(combined,
            paste0(foldername4, 
                   "tradeseq_significant_genes_with_dynverse_ranks.tsv"),
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t"
)

combined2 <- combined[1:100, ]
write.table(combined2,
            paste0(foldername4, 
                   "tradeseq_significant_genes_with_dynverse_ranks_seleceted_to_heatmap.tsv"),
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t"
)

( top_100_tradeseq <- my_plot_heatmap(
    sds,
    expression_source = dataset,
    features_oi = as.character(combined2$gene_id),
    scale = T) )

ggplot2::ggsave(
    paste0(foldername2, "TRADESEQ_top_100_whole_trajectory.png" ),
    top_100_tradeseq,
    height= nrow(feat_importances_annot_sub)*0.2,
    width=20,
    units="in",
    bg = "#FFFFFF", 
    dpi = 300)

### Heatmap ###

## Gather the cells order within the trajectory
pst.ord <- order(sds_trade$slingPseudotime_1, na.last = NA)
heatdata <- assays(sds_trade)$counts
heatdata <- heatdata[rownames(heatdata) %in% combined2$gene_id, pst.ord]

cells_order_in_traject <- tibble(cells = colnames(heatdata))

# Replaces the gene ID by its acronym if it is available.
for( i in 1:nrow(heatdata) ) {
    
    gene <- rownames(heatdata)[i]
    
    rownames(heatdata)[i] <- ifelse(
        
        rownames(heatdata)[i] %in% combined2$gene_id,
        unique(combined2[combined2$gene_id %in% rownames(heatdata)[i], "acronym"]),
        combined2$gene_id)
    
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
                                color = viridis(n = 20, option = "D"),
                                clustering_method="ward.D2",
                                fontsize=15)

foldername5 <- "RECLUSTERING/CORTEX_NOD_sunn4_cortex_meristem_trajectory/tradeseq_markers_cortex_nod/images/"
system( paste0("mkdir -p ", foldername5) )
ggplot2::ggsave(
    paste0(foldername5, "TRADESEQ_top100_lineage1_trajectory.svg"), 
    p_heatmap,
    height= nrow(heatdata_scaled)*0.2,
    width=20,
    units="in",
    bg = "#FFFFFF", 
    dpi = 300)

foldername5 <- "RECLUSTERING/CORTEX_NOD_sunn4_cortex_meristem_trajectory/tradeseq_markers_cortex_nod/images/"
system( paste0("mkdir -p ", foldername5) )
ggplot2::ggsave(
    paste0(foldername5, "TRADESEQ_top100_lineage1_trajectory.png"), 
    p_heatmap,
    height= nrow(heatdata_scaled)*0.2,
    width=20,
    units="in",
    bg = "#FFFFFF", 
    dpi = 300)

## Divides the trajectory in bins to facilitate the visualization
N_GROUPS = 50
colnames_heatmap <- c( 
    rep( 1:(N_GROUPS-1),
         each =  floor( ncol(heatdata) / (N_GROUPS) ) ),
    
    rep(N_GROUPS, ncol(heatdata) - ( (N_GROUPS-1) * floor(ncol(heatdata)/(N_GROUPS) ) ) )
)

print(table(colnames_heatmap))

cells_per_bin <- table(colnames_heatmap)

colnames(heatdata) <- colnames_heatmap

averaged_heatdata <- data.frame(row.names = rownames(heatdata))
for(c in unique(colnames_heatmap) ) {
    
    sub <- heatdata[, colnames(heatdata) == c]
    sub2 <- as.data.frame(rowMeans(sub))
    
    averaged_heatdata <- cbind(averaged_heatdata, sub2)
}

averaged_heatdata <- as.matrix(averaged_heatdata)
averaged_heatdata_log <- log1p(averaged_heatdata)
averaged_heatdata_scaled <- t(dynutils::scale_quantile(t(averaged_heatdata_log)))

bined_heatmap <- pheatmap::pheatmap(averaged_heatdata_scaled,
                                    show_rownames = T,
                                    show_colnames = F,
                                    cluster_rows = T,
                                    cluster_cols = FALSE,
                                    scale="none",
                                    color = viridis(n = 20, option = "D"),
                                    clustering_method="ward.D2",
                                    fontsize=15)

ggplot2::ggsave(paste0(foldername5, "TRADESEQ_top100_lineage1_trajectory_50BINs.svg"), 
                bined_heatmap,
                height= nrow(heatdata_scaled)*0.2,
                width=16,
                units="in",
                bg = "#FFFFFF", 
                dpi = 300)

ggplot2::ggsave(paste0(foldername5, "TRADESEQ_top100_lineage1_trajectory_50BINs.png"),
                bined_heatmap,
                height= nrow(heatdata_scaled)*0.2,
                width=16,
                units="in",
                bg = "#FFFFFF", 
                dpi = 300)

## Expression table of each bin in the trajectory but for all DEGs

pst.ord <- order(sds_trade$slingPseudotime_1, na.last = NA)
heatdata <- assays(sds_trade)$counts
heatdata <- heatdata[rownames(heatdata) %in% combined$gene_id, pst.ord]

cells_order_in_traject <- tibble(cells = colnames(heatdata))
for( i in 1:nrow(heatdata) ) {
    
    gene <- rownames(heatdata)[i]
    
    rownames(heatdata)[i] <- ifelse(
        
        rownames(heatdata)[i] %in% combined$gene_id,
        unique(combined[combined$gene_id %in% rownames(heatdata)[i], "acronym"]),
        combined$gene_id)
    
    rownames(heatdata)[i] <- ifelse(
        
        is.na(rownames(heatdata)[i]), 
        gene,
        rownames(heatdata)[i] )
    
}

N_GROUPS = 50
colnames_heatmap <- c( 
    rep( 1:(N_GROUPS-1),
         each =  floor( ncol(heatdata) / (N_GROUPS) ) ),
    
    rep(N_GROUPS, ncol(heatdata) - ( (N_GROUPS-1) * floor(ncol(heatdata)/(N_GROUPS) ) ) )
)

print( table(colnames_heatmap) )

cells_per_bin <- table(colnames_heatmap)

colnames(heatdata) <- colnames_heatmap

heatdata <- as.data.frame(heatdata)

averaged_heatdata <- data.frame(row.names = rownames(heatdata) )
for(c in unique(colnames_heatmap) ) {
    
    sub <- heatdata[, colnames(heatdata) == c]
    sub2 <- as.data.frame(rowMeans(sub))
    
    averaged_heatdata <- cbind(averaged_heatdata, sub2)
}

averaged_heatdata <- as.matrix(averaged_heatdata)
averaged_heatdata_log <- log1p(averaged_heatdata)
averaged_heatdata_scaled <- t(dynutils::scale_quantile(t(averaged_heatdata_log)))

averaged_heatdata_scaled <- as.data.frame(averaged_heatdata_scaled)
colnames(averaged_heatdata_scaled) <- paste0("bin_", seq(1, 50,1))

averaged_heatdata_scaled$gene <- rownames(averaged_heatdata_scaled)
averaged_heatdata_scaled <- averaged_heatdata_scaled[, c(51, 1:50)]

write.table(averaged_heatdata_scaled,
            paste0(foldername4, "expression_of_DEGs_within_trajectory_cortex_nodule_meristem.tst"),
            col.names = T, row.names = F, sep = "\t", quote = F)
