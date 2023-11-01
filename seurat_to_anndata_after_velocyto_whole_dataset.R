# Read the rds file
seurat_obj <- readRDS("data/seurat_formated_batched_integrated_clustered_complete_dataset.rds")

# save metadata table:
seurat_obj$barcode <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$UMAP@cell.embeddings[colnames(seurat_obj),1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$UMAP@cell.embeddings[colnames(seurat_obj),2]

system("mkdir -p ./tmp_whole")
write.csv(seurat_obj@meta.data, file='./tmp_whole/metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file='./tmp_whole/counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat_obj@reductions$pca@cell.embeddings, file='./tmp_whole/pca.csv', quote=F, row.names=F)

# write gene names
write.table(
    data.frame('gene'=rownames(seurat_obj)),file='./tmp_whole/gene_names.csv',
    quote=F,row.names=F,col.names=F
)