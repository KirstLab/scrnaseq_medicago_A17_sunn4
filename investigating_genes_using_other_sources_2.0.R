"Usage: investigating_genes_using_other_sources.R (--out_images <out_imag>) [--markers_file <markers_file>] [--markers_dir <markers_dir>] [--color_scheme=<color_scheme>] [--plot_name <p_name>] [--show_gene_names=<gene_names>] [--image_format=<img_f>]
--out_images <out_imag> Output directory where the images will be saved.
--markers_file <markers_file>   File containing the list of markers to be used, must be in the csv format and have the gene IDs in the first column. Use EITHER \"--markers_file\" OR \"--markers_dir\" [default: NULL].
--markers_dir <markers_dir> Folder containing the files that should be processed. Each file must be in the csv format and have the gene IDs in the first column. Use EITHER \"--markers_file\" OR \"--markers_dir\" [default: NULL].
--color_scheme=<color_scheme>   Color scheme to use in the heatmaps. Must be 'VIRIDIS' or 'MAGMA' [default: MAGMA].
--plot_name=<p_name> Name to be used to identify the plots. Only works with  \"--markers_file\" and will be ignored if using  \"--markers_dir\" [default: Plot].
--show_gene_names=<gene_names>  Indicates if the gene names is to be show in the rows of the heatmaps or not [default: TRUE].
--image_format=<img_f>  Format to export the image, can be svg or png [default: png].
investigating_genes_using_other_sources.R -h | --help  show this message.
" -> doc

suppressMessages( require(docopt) )
opts <- docopt(doc)

suppressMessages( require(vroom) )
suppressMessages( require(tidyverse) )
suppressMessages( require(circlize) )
suppressMessages( require(ComplexHeatmap) )
suppressMessages( require(monocle3) )

source("functions_to_investigating_genes_using_other_sources_2.0.R")

LCMv5 <- suppressMessages( vroom("data_from_other_sources/LCM/Parsed_JFLCM_v5FPKM.csv", delim = ",") )
colnames(LCMv5)[1] <- "gene"

atlas_v2 <- suppressMessages( vroom::vroom("data_from_other_sources/20211023_MtExpressV2-Dataset/data/gene/mtexpress_v2.log2_tmm.matrix.gene.tsv") )

if ( opts$markers_file == "NULL" && opts$markers_dir == "NULL" ) {
    
    print("Please add a file containing the markers file ('--markers_file') OR a directory contaning one or more markers files ('--markers_dir')")
    quit()
    
} else if ( opts$markers_file == "NULL" || opts$markers_dir == "NULL" ) {

} else {
    
    print("Please add a file containing the markers file ('--markers_file') OR a directory contaning one or more markers files ('--markers_dir'), but not both.")
    quit()
    
}

system( paste0("mkdir -p ", opts$out_images) )
if ( !opts$markers_dir == "NULL" ) {
    
    markers_files <- list.files(opts$markers_dir, full.names = T)
    
    for ( m in 1:length(markers_files) ) {
        print(m)
        ext <- tools::file_ext(markers_files[m])
        
        markers <- switch(ext,
                          csv = vroom::vroom(markers_files[m], delim = ","),
                          tsv = vroom::vroom(markers_files[m], delim = "\t"),
                          validate("Invalid file; Please upload a .csv or .tsv file") )
        
        markers <- as.data.frame( markers[, 1] )
        
        t_genes <- unique( markers[complete.cases(markers), 1] )
        
        markers_files_sub <- strsplit(markers_files[m], split = "/")
        markers_files_sub <- markers_files_sub[[1]][length(markers_files_sub[[1]])]
        markers_files_sub <- strsplit(markers_files_sub, split = "\\.")[[1]][1]
        
        print( paste0( "Comparing gene list on the file '", markers_files_sub, "' with the LCM dataset.") )
        
        LCM_func_v5(genes_names = t_genes,
                    plot_name = markers_files_sub)
        
        c_order = c(
            1, 3, 14, 22, # Epidermis / Root hair
            12, # Root hair
            15, 29, # Lateral root
            26, # Xylem
            10, 23, # Stele
            18, # Cell division vasculature
            8, # Pericycle
            7, 24, # Endodermis
            2, 4, 9, 11, 13, 16, # Cortex
            6, # Nodule
            5, # NA
            17, #NA
            19, 20, 21, #NA 
            25, # NA
            27, 28 # NA
        )
        
        cds_path = "rds_files/batched_integrated_clustered_complete_dataset.rds"
        
        sc_cluster( genes_names = t_genes,
                    plot_name = markers_files_sub,
                    cds_path = cds_path,
                    c_order = c_order)
        
        print( paste0( "Comparing gene list on the file '", markers_files_sub, "' with the SRP028599 dataset.") )
        atlast_v2_SRP028599(genes_names = t_genes,
                            plot_name = markers_files_sub)
        
        print( paste0( "Comparing gene list on the file '", markers_files_sub, "' with the SRP212693 dataset.") )
        
        atlast_v2_SRP212693_opt3(genes_names = t_genes,
                                 plot_name = markers_files_sub)
        
    }
}