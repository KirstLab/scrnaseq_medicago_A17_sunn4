require(hdWGCNA)
require(WGCNA)
require(tidyverse)

net <- readRDS("hdWGCNA_whole_dataset_pearson/network_TOM_file.rds")

gene_name <- "MtrunA17Chr5g0448621"
gene_common_name <- "MtNIN"
N_genes = 30

net_sub <- as.data.frame(net[gene_name,])
colnames(net_sub) <-  gene_common_name

net_formt <- net_sub %>% 
    as.data.frame() %>%
    rownames_to_column( "Target" ) %>%
    pivot_longer(!Target, names_to = "Source", values_to = "kME") %>%
    dplyr::filter( !is.na(kME) ) %>%
    dplyr::arrange( desc( kME ) )

annot_names <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv") %>%
    dplyr::rename(Target = locus_tag)
annot_summary <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv") %>%
    dplyr::rename(Target = LOCUS_TAG)

annot_summary <- merge( annot_summary,
                        annot_names,
                        by = "Target",
                        all.x = T)

net_formt_ann <- merge(net_formt, annot_summary, by = "Target", all.x = T)
net_formt_ann <- net_formt_ann %>%
    dplyr::arrange(desc(kME)) %>%
    slice_head(n = N_genes)

write_csv(net_formt_ann, 'MtNIN_test_whole_network.csv')

closeAllConnections()