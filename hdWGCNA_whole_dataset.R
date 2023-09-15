set.seed(1407)

outfile <- "hdWGCNA.out" # File name of output log
#Check its existence
if ( file.exists(outfile) ) {
  #Delete file if it exists
  file.remove(outfile)
}

my_log <- file(outfile) 
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")
# single-cell analysis package
require(Seurat)

# plotting and data science packages
require(tidyverse)
require(cowplot)
require(patchwork)

# co-expression network analysis packages:
require(WGCNA)
require(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

## Help functions
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

#############################
## Setting main parameters ##
#############################

## Number of cores
ncores = 24

## RDS file containing the Seurat (clustered) dataset
INPUT = "data/mt_whole_dataset/seurat_formated_whole_dataset.rds"

## Name of the output folder
sample_name = "mt_whole_dataset_pearson"

## What cell type (or cluster) is of interest
group_of_interest = "whole_dataset"

## what variable to use for batch correction (in the modules)
## >> Create a variable here telling with you should normalize or not
use_harmony = NULL
harmony_group = "timepoint"

## What variables to consider when building the metacells.
metacells_grouping = c("cell_type", harmony_group)

## Number of hub genes to be reported
n_of_hub_genes = 100

## What dimension reduction to use during the analysis (this if for the network construction, not data representation. Plots will always use UMAP).
my_reduction = "pca"

# optionally enable multithreading
enableWGCNAThreads(nThreads = ncores)

## Creates dir to save all outputs
dir_name <- paste0("mkdir -p ", sample_name)
system( dir_name)

# load the snRNA-seq dataset
seurat_obj <- readRDS(INPUT)

p1 <- DimPlot(seurat_obj, group.by='seurat_clusters', label=TRUE) +
  umap_theme() + ggtitle('') + NoLegend()

my_ggsave(name = paste0(sample_name, "/UMAP_plot_of_input_dataset.png"),
          plot = p1)

my_ggsave(name = paste0(sample_name, "/UMAP_plot_of_input_dataset.svg"),
          plot = p1)

seurat_obj <- Seurat::FindVariableFeatures(object = seurat_obj,
                                           assay = "RNA")
# Scaling the object
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

if ( is.null(seurat_obj@meta.data$cell_type) ) {
  
  seurat_obj@meta.data$cell_type <- seurat_obj@meta.data$seurat_clusters
  
}

seurat_obj@meta.data$cell_type <- "whole_dataset"

## dimplot 
p2 <- DimPlot(seurat_obj, group.by='cell_type', label=TRUE) +
  umap_theme() + ggtitle('') + NoLegend()

my_ggsave(name = paste0(sample_name, "/UMAP_plot_of_cell_types.png"),
          plot = p2)

my_ggsave(name = paste0(sample_name, "/UMAP_plot_of_cell_types.svg"),
          plot = p2)

## Set up Seurat object for WGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.001, # fraction of cells that a gene needs to be expressed in order to be included in the network.
  wgcna_name = sample_name # the name of the hdWGCNA experiment
)

# Construct metacells
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = metacells_grouping, # specify the columns in seurat_obj@meta.data to group by
  reduction = my_reduction, # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cell_type' # set the Idents of the metacell seurat object
)

seurat_obj <- NormalizeMetacells(seurat_obj)

###########################################
## Processing the Metacell Seurat Object ##
###########################################

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

####################################
## Co-expression network analysis ##
####################################

## Set up the expression matrix
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = group_of_interest, # the name of the group of interest in the group.by column
  group.by='cell_type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups#
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Select soft-power threshold
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
p6 <- wrap_plots(plot_list, ncol=2)
my_ggsave(name = paste0(sample_name, "/soft_power_decision.png"),
          plot = p6)

my_ggsave(name = paste0(sample_name, "/soft_power_decision.svg"),
          plot = p6)

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  setDatExpr=FALSE,
  tom_name = sample_name, # name of the topological overlap matrix written to disk
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

##########################################
### Module Eigengenes and Connectivity ###
##########################################

# Compute harmonized module eigengenes
# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars=harmony_group ###>>> We may remove this normalization, if it is affecting the biological variation that is expected
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)
write.csv(hMEs,
          file = paste0(sample_name, "/harmonized_module_eigengenes.csv"),
          row.names = TRUE)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
write.csv(MEs, file = paste0(sample_name, "/module_eigengenes.csv"),
          row.names = TRUE)

###################################
### Compute module connectivity ###

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type',
  group_name = group_of_interest,
  corFnc = 'cor',
  corOptions = "use = 'p', method = 'pearson'")

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = paste0(group_of_interest, "-M")
)

# plot genes ranked by kME for each module
n_whole_datasets <- length( unique(seurat_obj@misc$mt_whole_dataset_pearson$wgcna_modules$module) )
p7 <- PlotKMEs(seurat_obj, ncol = round(n_whole_datasets/3) )

my_ggsave(name = paste0(sample_name, "/genes_ranked_by_kME_per_module.png"),
          plot = p7, width = round(n_whole_datasets/3)*5, height = 15)

my_ggsave(name = paste0(sample_name, "/genes_ranked_by_kME_per_module.svg"),
          plot = p7, width = round(n_whole_datasets/3)*5, height = 15)

## Getting the module assignment table
# get the module assignment table:
modules <- GetModules(seurat_obj)
write.csv(modules, 
          file = paste0(sample_name, "/module_assignment_table.csv"),
          row.names = T)

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = n_of_hub_genes)
write_csv(hub_df, file = paste0(sample_name, "/hub_genes_of_modules.csv") )

seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = n_of_hub_genes,
  method='UCell')



###########################
### Basic Visualization ###
###########################

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  reduction = "UMAP",
  features='hMEs', # plot the hMEs
  order=TRUE, # order so the points with highest hMEs are on top
  ucell = T, # depending on Seurat vs UCell for gene scoring
  raster = F
)

# stitch together with patchwork
p8 <- wrap_plots(plot_list, ncol=round(n_whole_datasets/3)) 
my_ggsave(name = paste0(sample_name, "/UMAP_hMEs.png"),
          plot = p8, width = round(n_whole_datasets/3)*3, 
          height = 15)

my_ggsave(name = paste0(sample_name, "/UMAP_hMEs.svg"),
          plot = p8, width = round(n_whole_datasets/3)*3, 
          height = 15)

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  reduction = "UMAP",
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE, # depending on Seurat vs UCell for gene scoring
  raster = F
)

# stitch together with patchwork
p9 <- wrap_plots(plot_list, ncol=round(n_whole_datasets/3)) 

my_ggsave(name = paste0(sample_name, "/UMAP_hub_scores_provided_by_UCell.png"),
          plot = p9, width = round(n_whole_datasets/3)*3, height = 15)

my_ggsave(name = paste0(sample_name, "/UMAP_hub_scores_provided_by_UCell.svg"),
          plot = p9, width = round(n_whole_datasets/3)*3, height = 15)

# plot module correlagram
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

#################################
### Seurat plotting functions ###
# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs)
mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p12 <- DotPlot(seurat_obj, features=mods, group.by = 'cell_type') +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
my_ggsave(name = paste0(sample_name, "/dot_plot_of_module_expression.png"),
          plot = p12, width = round(n_whole_datasets/3)*3, height = 15)
my_ggsave(name = paste0(sample_name, "/dot_plot_of_module_expression.svg"),
          plot = p12, width = round(n_whole_datasets/3)*3, height = 15)

#Plot violin oplots per module
for (i in 1:length(mods) ) {
  
  md_name <- paste0("whole_dataset-M", i)
  
  md_plot <- VlnPlot(
    seurat_obj,
    features = md_name,
    pt.size = 0,
    group.by = 'cell_type') +
    # add box-and-whisker plots on top:
    geom_boxplot(width=.25, fill='white')+
    # change axis labels and remove legend:
    xlab('') + ylab('hME') + NoLegend()
  
  my_ggsave(name = paste0(sample_name, 
                          "/Violin_plot_", md_name, ".png"),
            plot = md_plot, height = 8, width = 16)
  my_ggsave(name = paste0(sample_name, 
                          "/Violin_plot_", md_name, ".svg"),
            plot = md_plot, height = 8, width = 16)
  
}

#############################
### Network Visualization ###
require(igraph)

ModuleNetworkPlot(seurat_obj, outdir = sample_name)

## Applying UMAP to co-expression networks
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 25, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
p13 <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
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
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE,
)

my_ggsave(name = paste0(sample_name,
                        "/UMAP_plot_25hubgenes_15_neighbors_with_labeled_top_3.png"),
          plot = p14, width = 12, height = 12)
my_ggsave(name = paste0(sample_name,
                        "/UMAP_plot_25hubgenes_15_neighbors_with_labeled_top_3.svg"),
          plot = p14, width = 12, height = 12)

####################################################
### Differential module eigengene (DME) analysis ###

## Comparing all possibilities for all time-points
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
    dplyr::mutate(p_val_adj = round(p_val_adj, 5) )# %>%
  #tibble::rownames_to_column("Module")
  
  DMEs_all <- rbind(DMEs_all, DMEs)
  
}

write_csv(DMEs_all, file = paste0(sample_name, "/Sig_Diff_Modules_between_time_points.csv") )

## This removed the scaled dataset from the seurat object, to reduce it size before writting it to the disk.
seurat_obj <- Seurat::DietSeurat(seurat_obj,counts = T,data = T,scale.data = F)

saveRDS(seurat_obj,
        file = paste0(sample_name, "/processed_seurat_obg.rds") )

#########################
### Trait correlation ###

treatment <- as.data.frame(seurat_obj@meta.data$timepoint)
colnames(treatment)[1] <- "timepoint"
treatment$treatment <- ifelse(treatment$timepoint == "0h",
                              "Control",
                              "Treated")

seurat_obj@meta.data$treatment <- treatment$treatment

# Covert the character to factors
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

# get the mt-correlation results
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

## save the network 
net <- GetTOM(seurat_obj)
saveRDS(net, "mt_whole_dataset_pearson/network_TOM_file.rds")

closeAllConnections()