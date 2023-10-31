set.seed(1407)

outfile <- "hdWGCNA.out"

if ( file.exists(outfile) ) {
  file.remove(outfile)
}

my_log <- file(outfile) 
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

require(Seurat)
require(tidyverse)
require(cowplot)
require(patchwork)
require(WGCNA)
require(hdWGCNA)

theme_set(theme_cowplot())

my_ggsave <- function( name = name,
                       plot = plot,
                       height = 12,
                       width = 14) {
  
  ggplot2::ggsave(
    filename = name,
    plot = plot,
    height = height,
    width = width,
    units = "in",
    bg = "#FFFFFF",
    dpi = 300,
    limitsize = F)
  
}

ncores = 24
INPUT = "data/mt_whole_dataset/seurat_formated_whole_dataset.rds"
sample_name = "mt_whole_dataset_pearson"
group_of_interest = "whole_dataset"

use_harmony = NULL
harmony_group = "timepoint"

metacells_grouping = c("cell_type", harmony_group)

n_of_hub_genes = 100
my_reduction = "pca"

enableWGCNAThreads(nThreads = ncores)

dir_name <- paste0("mkdir -p ", sample_name)
system( dir_name)

seurat_obj <- readRDS(INPUT)

p1 <- DimPlot(seurat_obj, group.by='seurat_clusters', label=TRUE) +
  umap_theme() + ggtitle('') + NoLegend()

my_ggsave(name = paste0(sample_name, "/UMAP_plot_of_input_dataset.png"),
          plot = p1)

my_ggsave(name = paste0(sample_name, "/UMAP_plot_of_input_dataset.svg"),
          plot = p1)

seurat_obj <- Seurat::FindVariableFeatures(object = seurat_obj,
                                           assay = "RNA")
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

if ( is.null(seurat_obj@meta.data$cell_type) ) {
  
  seurat_obj@meta.data$cell_type <- seurat_obj@meta.data$seurat_clusters
  
}

seurat_obj@meta.data$cell_type <- "whole_dataset"

p2 <- DimPlot(seurat_obj, group.by='cell_type', label=TRUE) +
  umap_theme() + ggtitle('') + NoLegend()

my_ggsave(name = paste0(sample_name, "/UMAP_plot_of_cell_types.png"),
          plot = p2)

my_ggsave(name = paste0(sample_name, "/UMAP_plot_of_cell_types.svg"),
          plot = p2)

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.001,
  wgcna_name = sample_name
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = metacells_grouping,
  reduction = my_reduction,
  k = 25,
  max_shared = 10,
  ident.group = 'cell_type'
)

seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj2 <- NormalizeMetacells(seurat_obj)
seurat_obj2 <- ScaleMetacells(seurat_obj2, features=VariableFeatures(seurat_obj2))
seurat_obj2 <- RunPCAMetacells(seurat_obj2, features=VariableFeatures(seurat_obj2))
seurat_obj2 <- RunHarmonyMetacells(seurat_obj2, group.by.vars='timepoint')
seurat_obj2 <- RunUMAPMetacells(seurat_obj2, reduction='harmony', dims=1:15)

p3 <- DimPlotMetacells(seurat_obj2, group.by='cell_type',label = T) + 
  umap_theme() +
  ggtitle("Cell Type") +
  theme(legend.position = "none")
p4 <- DimPlotMetacells(seurat_obj2, group.by='cell_type') + 
  umap_theme() +
  ggtitle("Cell Type") 
p5 <- p3|p4

my_ggsave(name = paste0(sample_name, "/UMAP_plot_of_metacells.png"),
          plot = p5)

my_ggsave(name = paste0(sample_name, "/UMAP_plot_of_metacells.svg"),
          plot = p5)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = group_of_interest,
  group.by='cell_type',
  assay = 'RNA',
  slot = 'data'
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed'
)

plot_list <- PlotSoftPowers(seurat_obj)

p6 <- wrap_plots(plot_list, ncol=2)
my_ggsave(name = paste0(sample_name, "/soft_power_decision.png"),
          plot = p6)

my_ggsave(name = paste0(sample_name, "/soft_power_decision.svg"),
          plot = p6)

seurat_obj <- ConstructNetwork(
  seurat_obj, 
  setDatExpr=FALSE,
  tom_name = sample_name,
  overwrite_tom = TRUE
)

png(filename = paste0(sample_name, "/hdWGCNA_Dendrogram.png"),
    height = 10, width = 20,units = "cm",res = 300, bg = "white")
PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
dev.off()

svg(filename = paste0(sample_name, "/hdWGCNA_Dendrogram.svg"),
    height = 10, width = 20, bg = "white")
PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
dev.off()

PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')

seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars=harmony_group
)

hMEs <- GetMEs(seurat_obj)
write.csv(hMEs,
          file = paste0(sample_name, "/harmonized_module_eigengenes.csv"),
          row.names = TRUE)

MEs <- GetMEs(seurat_obj, harmonized=FALSE)
write.csv(MEs, file = paste0(sample_name, "/module_eigengenes.csv"),
          row.names = TRUE)

seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type',
  group_name = group_of_interest,
  corFnc = 'cor',
  corOptions = "use = 'p', method = 'pearson'")

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = paste0(group_of_interest, "-M")
)

n_whole_datasets <- length( unique(seurat_obj@misc$mt_whole_dataset_pearson$wgcna_modules$module) )
p7 <- PlotKMEs(seurat_obj, ncol = round(n_whole_datasets/3) )

my_ggsave(name = paste0(sample_name, "/genes_ranked_by_kME_per_module.png"),
          plot = p7, width = round(n_whole_datasets/3)*5, height = 15)

my_ggsave(name = paste0(sample_name, "/genes_ranked_by_kME_per_module.svg"),
          plot = p7, width = round(n_whole_datasets/3)*5, height = 15)

modules <- GetModules(seurat_obj)
write.csv(modules, 
          file = paste0(sample_name, "/module_assignment_table.csv"),
          row.names = T)

hub_df <- GetHubGenes(seurat_obj, n_hubs = n_of_hub_genes)
write_csv(hub_df, file = paste0(sample_name, "/hub_genes_of_modules.csv") )

seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = n_of_hub_genes,
  method='UCell')

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  reduction = "UMAP",
  features='hMEs',
  order=TRUE,
  ucell = T,
  raster = F
)

p8 <- wrap_plots(plot_list, ncol=round(n_whole_datasets/3)) 
my_ggsave(name = paste0(sample_name, "/UMAP_hMEs.png"),
          plot = p8, width = round(n_whole_datasets/3)*3, 
          height = 15)

my_ggsave(name = paste0(sample_name, "/UMAP_hMEs.svg"),
          plot = p8, width = round(n_whole_datasets/3)*3, 
          height = 15)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  reduction = "UMAP",
  features='scores',
  order='shuffle',
  ucell = TRUE,
  raster = F
)

p9 <- wrap_plots(plot_list, ncol=round(n_whole_datasets/3)) 

my_ggsave(name = paste0(sample_name, "/UMAP_hub_scores_provided_by_UCell.png"),
          plot = p9, width = round(n_whole_datasets/3)*3, height = 15)

my_ggsave(name = paste0(sample_name, "/UMAP_hub_scores_provided_by_UCell.svg"),
          plot = p9, width = round(n_whole_datasets/3)*3, height = 15)

png(filename = paste0(sample_name, "/module_correlagram.png"),
    height = round(n_whole_datasets/3)*3,
    width = round(n_whole_datasets/3)*3,
    units = "cm",res = 300, bg = "white")
ModuleCorrelogram(seurat_obj)
dev.off()

svg(filename = paste0(sample_name, "/module_correlagram.svg"),
    height = round(n_whole_datasets/3)*3,
    width = round(n_whole_datasets/3)*3, bg = "white")
ModuleCorrelogram(seurat_obj)
dev.off()

MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs)
mods <- mods[mods != 'grey']

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

p12 <- DotPlot(seurat_obj, features=mods, group.by = 'cell_type') +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

my_ggsave(name = paste0(sample_name, "/dot_plot_of_module_expression.png"),
          plot = p12, width = round(n_whole_datasets/3)*3, height = 15)
my_ggsave(name = paste0(sample_name, "/dot_plot_of_module_expression.svg"),
          plot = p12, width = round(n_whole_datasets/3)*3, height = 15)

for (i in 1:length(mods) ) {
  
  md_name <- paste0("whole_dataset-M", i)
  
  md_plot <- VlnPlot(
    seurat_obj,
    features = md_name,
    pt.size = 0,
    group.by = 'cell_type') +
    geom_boxplot(width=.25, fill='white')+
    xlab('') + ylab('hME') + NoLegend()
  
  my_ggsave(name = paste0(sample_name, 
                          "/Violin_plot_", md_name, ".png"),
            plot = md_plot, height = 8, width = 16)
  my_ggsave(name = paste0(sample_name, 
                          "/Violin_plot_", md_name, ".svg"),
            plot = md_plot, height = 8, width = 16)
  
}

require(igraph)
ModuleNetworkPlot(seurat_obj, outdir = sample_name)

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 25,
  n_neighbors=15,
  min_dist=0.1
)

umap_df <- GetModuleUMAP(seurat_obj)

p13 <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color,
    size=umap_df$kME*2
  ) +
  umap_theme()

my_ggsave(name = paste0(sample_name,
                        "/UMAP_plot_25hubgenes_15_neighbors.png"),
          plot = p13, width = 12, height = 12)
my_ggsave(name = paste0(sample_name,
                        "/UMAP_plot_25hubgenes_15_neighbors.svg"),
          plot = p13, width = 12, height = 12)

p14 <- ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1,
  label_hubs=2,
  keep_grey_edges=FALSE,
)

my_ggsave(name = paste0(sample_name,
                        "/UMAP_plot_25hubgenes_15_neighbors_with_labeled_top_3.png"),
          plot = p14, width = 12, height = 12)
my_ggsave(name = paste0(sample_name,
                        "/UMAP_plot_25hubgenes_15_neighbors_with_labeled_top_3.svg"),
          plot = p14, width = 12, height = 12)

md_comparisons <- as_tibble( expand.grid( c(0,24, 48, 96),
                                          c(0,24, 48, 96) ) ) %>%
  dplyr::filter(!Var1 == Var2) %>%
  dplyr::filter(!Var2 < Var1) %>%
  dplyr::filter(!Var2 == 0) %>%
  dplyr::mutate(Var1 = paste0(Var1, "h")) %>%
  dplyr::mutate(Var2 = paste0(Var2, "h"))

DMEs_all <- data.frame()
for (c in 1:nrow(md_comparisons) ) {
  
  g1 <- md_comparisons[c, 1] %>% pull()
  g2 <- md_comparisons[c, 2] %>% pull()
  
  group1 <- seurat_obj@meta.data %>% 
    subset(cell_type == group_of_interest & timepoint == g1) %>%
    rownames
  
  group2 <- seurat_obj@meta.data %>% 
    subset(cell_type == group_of_interest & timepoint == g2) %>%
    rownames
  
  DMEs <- FindDMEs(
    seurat_obj,
    barcodes1 = group1,
    barcodes2 = group2,
    test.use='wilcox',
    wgcna_name=NULL, 
    verbose = T ) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::mutate(comparison = paste0(g1, "_vs_", g2) ) %>%
    dplyr::mutate(p_val_adj = round(p_val_adj, 5) )
  
  DMEs_all <- rbind(DMEs_all, DMEs)
  
}

write_csv(DMEs_all, file = paste0(sample_name, "/Sig_Diff_Modules_between_time_points.csv") )

seurat_obj <- Seurat::DietSeurat(seurat_obj,counts = T,data = T,scale.data = F)

saveRDS(seurat_obj,
        file = paste0(sample_name, "/processed_seurat_obg.rds") )

treatment <- as.data.frame(seurat_obj@meta.data$timepoint)
colnames(treatment)[1] <- "timepoint"
treatment$treatment <- ifelse(treatment$timepoint == "0h",
                              "Control",
                              "Treated")

seurat_obj@meta.data$treatment <- treatment$treatment

seurat_obj$timepoint <- as.factor(seurat_obj$timepoint)
seurat_obj$batch2 <- as.factor(seurat_obj$batch2)
seurat_obj$Group <- as.factor(seurat_obj$Group)
seurat_obj$treatment <- as.factor(seurat_obj$treatment)

cur_traits <- c("timepoint",
                "Group",
                "treatment")

seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = cur_traits,
  group.by='cell_type'
)

mt_cor <- GetModuleTraitCorrelation(seurat_obj)
CORR <- as.data.frame(mt_cor$cor$whole_dataset)
FDR <- as.data.frame(mt_cor$fdr$whole_dataset)

cor_whole_dataset <- rbind(CORR, FDR)
write.csv(cor_whole_dataset,
          file = "correlations_module_timepoint_and_treatment.csv", 
          row.names = T)

png(filename = paste0(sample_name, "/correlations_module_timepoint_and_treatment.png"),
    height = 50, width = 20,units = "cm",res = 300, bg = "white")

PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = 'yellow',
  mid_color = 'black',
  low_color = 'purple',
  plot_max = 0.2,
  combine=TRUE
)

dev.off()

svg(filename = paste0(sample_name, "/correlations_module_timepoint_and_treatment.svg"),
    height = 50, width = 20, bg = "white")

PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = 'yellow',
  mid_color = 'black',
  low_color = 'purple',
  plot_max = 0.2,
  combine=TRUE
)

dev.off()

net <- GetTOM(seurat_obj)
saveRDS(net, "mt_whole_dataset_pearson/network_TOM_file.rds")

closeAllConnections()