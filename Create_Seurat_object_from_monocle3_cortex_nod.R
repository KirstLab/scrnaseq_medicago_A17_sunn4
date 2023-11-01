set.seed(1407)

outfile <- "logs/Create_Seurat_object_from_monocle3_cortex_nod.out" # File name of output log
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

gene_names <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv")

monocle_to_seurat <- function( my.cds ) {
    
    colData(my.cds)$cell_type <- monocle3::clusters(my.cds)
    
    ## transfer the counts from the monocle object
    rna_counts <- monocle3::exprs(my.cds)

    ## Transfer the normalized data
    rna_data <- monocle3::normalized_counts(my.cds)
    
    for(i in 1:nrow(rna_data)) {
        
        rownames(rna_data)[1] <- ifelse( rownames(rna_data)[1] %in% gene_names$locus_tag,
                                           gene_names$acronym,
                                           gene_names$locus_tag)
    }
    
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
    adata
}

# read data
cds <- readRDS( "RECLUSTERING/CORTEX_NOD/rds_file_subset/clusters_cortex_nodule_selected_for_trajectory.rds" )

cds_seurat <- monocle_to_seurat(cds)

saveRDS(cds_seurat,
        "RECLUSTERING/CORTEX_NOD/rds_file_subset/seurat_formated_selected_clusters_cortex_nodule.rds")

# read data
cds <- readRDS( "RECLUSTERING/CORTEX_NOD/rds_file_subset/cortex_nod_selected_clusters_RECLUSTERED.rds" )

cds_seurat <- monocle_to_seurat(cds)

saveRDS(cds_seurat,
        "RECLUSTERING/CORTEX_NOD/rds_file_subset/seurat_formated_cortex_nod_selected_clusters_RECLUSTERED.rds")

closeAllConnections()