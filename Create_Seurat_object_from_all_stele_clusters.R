set.seed(1407)

require(Seurat)
require(monocle3)
require(vroom)
require(tidyverse)

gene_names <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv")

monocle_to_seurat <- function( my.cds ) {
    
    colData(my.cds)$cell_type <- monocle3::clusters(my.cds)
    
    rna_counts <- monocle3::exprs(my.cds)
    rna_data <- monocle3::normalized_counts(my.cds)
    
    adata <- CreateSeuratObject(counts = rna_counts,
                                assay = "RNA",
                                project = "scMedicago")
    
    rna_counts_assay <- CreateAssayObject(counts = rna_counts)
    rna_counts_assay@key <- "RNA_"
    adata@assays$RNA <- rna_counts_assay
    adata@assays$RNA@data <- rna_data
    
    mat_emb <- reducedDims(my.cds)$UMAP
    
    UMA_assay <- CreateDimReducObject(embeddings = mat_emb, key = 'UMAP_' )
    UMA_assay@key <- "UMAP_"
    adata@reductions$UMAP <- UMA_assay
    
    mat_emb <- reducedDims(my.cds)$PCA
    
    pca_assay <- CreateDimReducObject(embeddings = mat_emb, key = 'pca_' )
    pca_assay@key <- "pca_"
    adata@reductions$pca <- pca_assay
    
    colData( my.cds )$seurat_clusters <- colData( my.cds )$cell_type
    adata <- AddMetaData(adata, as.data.frame( colData( my.cds ) ) )
    adata
}

cds <- readRDS( "RECLUSTERING/PERICYCLE_v2/rds_file_subset/medicago_integrated_subset_pericycle.rds" )
cds_seurat <- monocle_to_seurat(cds)

saveRDS(cds_seurat, 
        "RECLUSTERING/PERICYCLE_v2/rds_file_subset/seurat_formated_pericycle.rds")