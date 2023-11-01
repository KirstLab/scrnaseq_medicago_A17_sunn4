set.seed(1407)

outfile <- "logs/DEGs_across_timepoints_whole_dataset.out" # File name of output log
#Check its existence
if ( file.exists(outfile) ) {
    #Delete file if it exists
    file.remove(outfile)
}

my_log <- file(outfile) 
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

require(Seurat)
require(tidyverse)
seurat_obj <- readRDS("rds_files/seurat_formated_whole_dataset.rds")

seurat_obj <- NormalizeData(seurat_obj,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

seurat_obj@meta.data$c_time <- paste0("C",
                                      seurat_obj@meta.data$seurat_clusters, "_",
                                      seurat_obj@meta.data$timepoint)

DimPlot(seurat_obj, reduction = "UMAP", group.by = "seurat_clusters")
DimPlot(seurat_obj, reduction = "UMAP", group.by = "timepoint")
DimPlot(seurat_obj, reduction = "UMAP", group.by = "c_time")

Idents(seurat_obj) <- seurat_obj@meta.data$c_time
all_ctime <- unique(seurat_obj@meta.data$c_time)

freq_cells <- as.data.frame( table(seurat_obj@meta.data$c_time) )

# help function
my_wilcox <- function(g1 = g1, g2 = g2) { 
    
    n1 <- freq_cells %>%
        filter(Var1 == g1) %>%
        select(Freq) %>% as.numeric()
    n2 <- freq_cells %>%
        filter(Var1 == g2) %>%
        select(Freq) %>% as.numeric()
    
    if ( is.na(n1) | is.na(n2) ) { 
    } else {
        
        if (n1 < 20 | n2 < 20) { 
            
            print( paste0( "One of both groups in the comparison '",
                           g1, " vs ", g2,
                           "' have less than 20 cells"))
            
            empty_df <- data.frame(matrix(NA, nrow = 1, ncol = 7))
            
            colnames(empty_df) <- c("p_val",
                                    "avg_log2FC",
                                    "pct.1",
                                    "pct.2",
                                    "p_val_adj",
                                    "comparison",
                                    "genes")
            
            empty_df["comparison"] <- paste0(g1, " vs ", g2)
            empty_df
        } else {
            
            comp_df <- FindMarkers(seurat_obj,
                                   ident.1 = g1,
                                   ident.2 = g2,
                                   logfc.threshold = 1,
                                   min.pct = 0.05,
                                   min.cells.group = 20) %>% 
                dplyr::mutate(comparison = paste0(g1, " vs ", g2) ) %>%
                dplyr::mutate(comparison = paste0(
                    g1, " vs ",
                    g2 ) ) %>%
                dplyr::filter(p_val_adj < 0.05)
            
            comp_df$genes <- rownames(comp_df)
            
            comp_df
        } 
    }
}
all_degs <- data.frame()
for (c in paste0("C", 1:24) ) {
    
    print(c)
    
    dge_0_vs_24 <- my_wilcox(g1 = paste0(c, "_", "24h"),
                             g2 = paste0(c, "_", "0h") )
    
    dge_24_vs_48 <- my_wilcox(g1 = paste0(c, "_", "48h"),
                              g2 = paste0(c, "_", "24h") )
    
    dge_0_vs_48 <- my_wilcox(g1 = paste0(c, "_", "48h"),
                             g2 = paste0(c, "_", "0h") )
    
    dge_48_vs_96 <- my_wilcox(g1 = paste0(c, "_", "96h"),
                              g2 = paste0(c, "_", "48h") )
    
    dge_24_vs_96 <- my_wilcox(g1 = paste0(c, "_", "96h"),
                              g2 = paste0(c, "_", "24h") )
    
    dge_0_vs_96 <- my_wilcox(g1 = paste0(c, "_", "96h"),
                             g2 = paste0(c, "_", "0h") )
    
    
    combined <- rbind(dge_0_vs_24, dge_24_vs_48, dge_0_vs_48,
                      dge_48_vs_96,dge_24_vs_96, dge_0_vs_96)
    
    all_degs <- rbind(all_degs, combined)
    
}

all_degs <- all_degs %>%
    filter(!is.na(all_degs$genes))

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

all_degs <- merge( all_degs,
                   annot_summary, by.x = "genes",
                   by.y = "gene_id",
                   all.x = T ) %>%
    arrange(comparison, desc(avg_log2FC))

write.table(all_degs,
            "all_DEGs_whole_dataset_among_timepoints_using_Seurat.tsv",
            col.names = T,
            row.names = F,
            quote = T,
            sep = ",")

