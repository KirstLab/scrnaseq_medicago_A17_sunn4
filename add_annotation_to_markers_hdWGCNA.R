set.seed(1407)

require(vroom)
require(tidyverse)

annot_names <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv") %>%
    dplyr::rename(gene_name = locus_tag)
annot_summary <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv") %>%
    dplyr::rename(gene_name = LOCUS_TAG)

annot_summary <- merge( annot_summary,
                        annot_names,
                        by = "gene_name",
                        all.x = T)

annot_summary2 <- annot_summary %>%
    dplyr::rename(gene_id = gene_name) %>%
    dplyr::rename(gene_name = acronym) %>%
    dplyr::rename(acronym = gene_id)

sub <- vroom("hdWGCNA_whole_dataset_pearson/module_assignment_table.csv") %>%
    dplyr::filter(!module == "grey") %>%
    dplyr::select(-"...1")

for (m in 1:length( unique(sub$module) ) ) {
    
    var_n <- paste0("M", m,"$")
    keep_c <- c("gene_name", "module", "color")
    
    sub2 <- sub %>%
        dplyr::filter( grepl(var_n, module) ) %>%
        dplyr::select_if( grepl( var_n, names(.) ) | names(.) %in% keep_c )
    
    sub2a <- sub2 %>%
        dplyr::filter( grepl("MtrunA17", gene_name) ) %>% 
        merge( annot_summary, 
               by = "gene_name",
               all.x = T )
    
    sub2b <- sub2 %>%
        dplyr::filter( !grepl("MtrunA17", gene_name) ) %>%
        merge( annot_summary2, 
               by = "gene_name",
               all.x = T )
    
    sub2 <- rbind(sub2a, sub2b) %>%
        dplyr::arrange( desc( pick( starts_with( "kME" ) ) ) )
    
    system("mkdir -p hdWGCNA_whole_dataset_pearson/gene_annotation_of_modules/")
    write_csv(sub2,
              paste0("hdWGCNA_whole_dataset_pearson/gene_annotation_of_modules/module_annotation_M", m, ".csv") )
}