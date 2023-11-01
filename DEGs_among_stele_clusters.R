set.seed(1407)

require(Seurat)
require(tidyverse)
seurat_obj <- readRDS("RECLUSTERING/PERICYCLE_v2/rds_file_subset/seurat_formated_pericycle.rds")

seurat_obj <- NormalizeData(seurat_obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

DimPlot(seurat_obj, reduction = "UMAP", group.by = "seurat_clusters")

Idents(seurat_obj) <- seurat_obj@meta.data$seurat_clusters

freq_cells <- table( as.data.frame( Idents(seurat_obj) ) )

comp_df <- FindMarkers(seurat_obj,
                       ident.1 = 16,
                       ident.2 = c(5, 7, 13, 14, 16, 17,18, 19),
                       logfc.threshold = 1,
                       min.pct = 0.05,
                       min.cells.group = 20) %>% 
    dplyr::mutate(comparison = paste0("C16", " vs ", "stele") ) %>%
    dplyr::filter(p_val_adj < 0.05)

comp_df$genes <- rownames(comp_df)

## Adds annotation
annot_names <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv")
annot_summary <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv")
colnames(annot_names)[1] <- "gene_id"
colnames(annot_summary)[1] <- "gene_id"

annot_summary <- merge( annot_summary,
                        annot_names,
                        by = "gene_id",
                        all.x = T)
annot_summary <- annot_summary[, c(1, 14, 2:13, 15:16)]

all_degs <- merge( comp_df,
                   annot_summary, by.x = "genes",
                   by.y = "gene_id",
                   all.x = T ) %>%
    arrange(comparison, desc(avg_log2FC))

write.table(all_degs,
            "all_DEGs_C16_vs_Stele_using_Seurat.tsv",
            col.names = T,
            row.names = F,
            quote = T,
            sep = ",")
