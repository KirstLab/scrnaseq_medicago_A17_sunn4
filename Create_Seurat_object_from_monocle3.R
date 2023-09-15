set.seed(1407)

outfile <- "logs/Create_Seurat_object_from_monocle3.out" # File name of output log
#Check its existence
if ( file.exists(outfile) ) {
    #Delete file if it exists
    file.remove(outfile)
}

my_log <- file(outfile) 

sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

require(Seurat)
require(monocle3)
require(vroom)
require(tidyverse)

my.cds <- readRDS( "data/mt_whole_dataset/batched_integrated_clustered_complete_dataset.rds" )
colData(my.cds)$cell_type <- monocle3::clusters(my.cds)

plot_cells(my.cds,
           label_cell_groups = T,
           graph_label_size = 1.5,
           cell_size = 1,
           group_label_size = 7)

## transfer the counts from the monocle object
rna_counts <- monocle3::exprs(my.cds)

## Transfer the normalized data
rna_data <- monocle3::normalized_counts(my.cds)

## Create a new Seurat Object
adata <- CreateSeuratObject(counts = rna_counts,
                            assay = "RNA",
                            project = "scMedicago")

## Merge the data (normalized) dataset to the object
rna_counts_assay <- CreateAssayObject(counts = rna_counts)
rna_counts_assay@key <- "RNA_"
adata@assays$RNA <- rna_counts_assay

adata@assays$RNA@data <- rna_data

# Add the cell embeddings for UMAP
mat_emb <- reducedDims(my.cds)$UMAP

UMA_assay <- CreateDimReducObject(embeddings = mat_emb, key = 'UMAP_' )
UMA_assay@key <- "UMAP_"
adata@reductions$UMAP <- UMA_assay

## Adds embeding from PCA
# Add the cell embeddings for UMAP
mat_emb <- reducedDims(my.cds)$PCA

pca_assay <- CreateDimReducObject(embeddings = mat_emb, key = 'pca_' )
pca_assay@key <- "pca_"
adata@reductions$pca <- pca_assay

## Add clustering and other meta information
colData( my.cds )$seurat_clusters <- colData( my.cds )$cell_type
adata <- AddMetaData(adata, as.data.frame( colData( my.cds ) ) )

Seurat::DimPlot(adata, reduction = 'UMAP',
                group.by = "seurat_clusters",
                pt.size = 1 )

## Test to see if the expression is matching the expected
###Counts
Seurat::FeaturePlot(adata,
                    features = "MtrunA17Chr1g0197491",
                    pt.size = 1,
                    slot = "counts")

### Normalized
Seurat::FeaturePlot(adata,
                    features = "MtrunA17Chr1g0197491",
                    pt.size = 1,
                    slot = "data")

saveRDS(adata, "data/mt_whole_dataset/seurat_formated_whole_dataset.rds")

closeAllConnections()