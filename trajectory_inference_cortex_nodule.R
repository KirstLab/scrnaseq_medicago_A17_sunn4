set.seed(1407)

outfile <- "logs/trajectory_inference_cortex_nodule.out" # File name of output log
#Check its existence
if ( file.exists(outfile) ) {
    #Delete file if it exists
    file.remove(outfile)
}

my_log <- file(outfile)
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

require(ggplot2)         # plot utilities 
require(monocle3)
require(tidyverse)
require(paletteer)

## Help functions ##
theme_umap <- function(base.size = 14) {
    ggplot2::theme_classic(base_size = base.size) + 
        ggplot2::theme(axis.ticks = ggplot2::element_blank(), 
                       axis.text = ggplot2::element_blank(), 
                       plot.subtitle = ggplot2::element_text(face = "italic", size = 11), 
                       plot.caption = ggplot2::element_text(face = "italic", size = 11))
}
guide_umap <- function(key.size = 4) {
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = key.size, alpha = 1)))
}

## Help functions
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

palette_cluster <- paletteer::paletteer_d("ggsci::default_jama")
palette_celltype <- paletteer::paletteer_d("ggsci::category20_d3")
palette_heatmap <- paletteer::paletteer_d("wesanderson::Zissou1")

# Reads the RDS file containing the cluster for the trajectory. This is the clustered data as output by Seurat.
cds <- readRDS("RECLUSTERING/CORTEX_NOD_sunn4_all_clusters_but_CN5/rds_file_subset/medicago_integrated_subset_cortex_nodule.rds")

# ### Removing the subcluster that won't be part of the trajectory
to_remove <- c(3)

colData(cds)$cluster <- monocle3::clusters(cds)
cds_sub_meta <- as.data.frame( colData(cds) )
cds_sub_meta <- cds_sub_meta[!cds_sub_meta$cluster %in% to_remove, ]

to_remove_cells <- as.character(rownames(cds_sub_meta))
cds <- cds[, colnames(cds) %in% to_remove_cells ]

( pp <- monocle3::plot_cells(cds,
                             label_cell_groups=T,
                             graph_label_size=1.5,
                             cell_size = 1,
                             group_label_size = 7) )

# Reclustering
cds <- clear_cds_slots(cds)

cds <- monocle3::preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds)

cds <- cluster_cells(cds,
                     verbose = T, 
                     resolution = 1e-03, 
                     random_seed = 1407)

monocle3::plot_cells(cds,
                     label_cell_groups=T,
                     graph_label_size=1.5,
                     cell_size = 1,
                     group_label_size = 7)

system("mkdir -p RECLUSTERING/CORTEX_NOD_sunn4_all_clusters_but_CN5/tradeseq_markers_cortex_nod/images/")

ggsave(filename = "RECLUSTERING/CORTEX_NOD_sunn4_all_clusters_but_CN5/tradeseq_markers_cortex_nod/images/cortex_nod_without_cn5_UMAP.svg",
       pp, width = 12, height = 8, bg = "white")

UMAP <- SingleCellExperiment::reducedDims(cds)[["UMAP"]]

colData(cds)$cluster <- monocle3::clusters(cds)

### trajectory inference with slingshot
sling_res <- slingshot::slingshot( UMAP, 
                                   clusterLabels = colData(cds)$cluster, 
                                   start.clus = "5",  # select starting point
                                   approx_points = 1000,
                                   extend = "n")

slingshot::slingLineages(sling_res)

sling_pt <- slingshot::slingPseudotime(sling_res) %>% 
    as.data.frame() 

# Check how many columns (lineages) exist and generate the colnames
col_n <- paste("PT", seq(1:ncol(sling_pt)), sep = "")
sling_pt <- sling_pt %>% 
    magrittr::set_colnames(col_n) 

# For each lineage, generate a new UMAP showing the pseudotime.
# Cell ordering based on slingshot
for (c in 1:ncol(sling_pt)) {
    
    pt_lin <- paste0("PT", c)
    
    ( ps_plot <- UMAP %>% 
            as.data.frame() %>% 
            magrittr::set_colnames(c("UMAP_1", "UMAP_2")) %>% 
            mutate(PT = sling_pt[[pt_lin]]) %>% 
            ggplot(aes(x = UMAP_1, y = UMAP_2, color = PT)) + 
            geom_point(size = 1, alpha = 0.75) + 
            labs(x = "UMAP 1", 
                 y = "UMAP 2", 
                 color = paste("Pseudotime lineage",c )) + 
            scale_color_gradientn(colors = palette_heatmap, 
                                  labels = scales::label_number(accuracy = 0.1),
                                  na.value="white") + 
            theme_umap() )
    
    ggsave(filename = paste("RECLUSTERING/CORTEX_NOD_sunn4_all_clusters_but_CN5/tradeseq_markers_cortex_nod/images/Pseudotime_lineage", c, "UMAP.png"),
           ps_plot, width = 12, height = 8, dpi = 300, bg = "white")
    
    ggsave(filename = paste("RECLUSTERING/CORTEX_NOD_sunn4_all_clusters_but_CN5/tradeseq_markers_cortex_nod/images/Pseudotime_lineage", c, "UMAP.svg"),
           ps_plot, width = 12, height = 8, dpi = 300, bg = "white")
    
}

## Reads the files from the genome to add the annotation to the genes lists generated below.
annot_names <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv")
annot_summary <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv")
colnames(annot_names)[1] <- "gene_id"
colnames(annot_summary)[1] <- "gene_id"

annot_summary <- merge( annot_summary,
                        annot_names,
                        by = "gene_id",
                        all.x = T)

annot_summary <- annot_summary[, c(1, 14, 2:13)]

# ########################################
### Trajectory Differential Expression ###
##########################################

require(tradeSeq)
require(slingshot)

## Using TradeSeq
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 5

sds2 <- SlingshotDataSet(sling_res)

counts_mt <- as.matrix(assay(cds))
counts_mt <- counts_mt[rowSums(counts_mt) > 0, ]

# To evaluated the number of knots to be used on the modeling phase
icMat <- evaluateK(counts = counts_mt,
                   sds = sds2,
                   verbose = F)

sce_fitted <- tradeSeq::fitGAM(counts = counts_mt,
                               sds = sds2,
                               nknots = 6,
                               verbose = T,
                               parallel = TRUE,
                               BPPARAM = BPPARAM)

saveRDS(sce_fitted,
        file = paste0("RECLUSTERING/CORTEX_NOD_sunn4_all_clusters_but_CN5/rds_file_subset/TRADEseq_cortex_nodules.rds")
)

#sce_fitted <- readRDS("RECLUSTERING/CORTEX_NOD_sunn4_all_clusters_but_CN5/rds_file_subset/TRADEseq_cortex_nodules.rds")

## Association test to check if each genes is DEG within the trajectories of each lineage.
foldername4 <- "RECLUSTERING/CORTEX_NOD_sunn4_all_clusters_but_CN5/tradeseq_markers_cortex_nod/"
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

my_write_csv( DEGs_within_traj_of_lineage_1,
              paste0(foldername4, "/tradeseq_DEGs_lineage1_c214356.csv") )

## Select among the known RNS 
RNS_genes <- vroom::vroom("selected_markers/RNS_genes_curated_list.tsv")

DEG_lineage_1_exc_RNS <- DEGs_within_traj_of_lineage_1 %>%
    dplyr::filter(gene_id %in% RNS_genes$gene_id)

### Heatmap ###

## Gather the cells order within the trajectory
sling_pt_heat <- slingshot::slingPseudotime(sling_res) %>% 
    as.data.frame() 

#### Lineage1 - exclusive DEGs
Lin1 <- data.frame(gene_id =  DEG_lineage_1_exc_RNS$gene_id )

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

foldername5 <- "RECLUSTERING/CORTEX_NOD_sunn4_all_clusters_but_CN5/tradeseq_markers_cortex_nod/images/"

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

## Divides the trajectory in bins to facilitate the visualization
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
closeAllConnections()