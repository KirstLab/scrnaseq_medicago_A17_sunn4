"Usage: add_annotation_to_markers.R (--dir <dir>) (--log_name <log>)
--dir   Directory containing the files (CSV format) to be annotated. Each file must have the gene IDs on the first column.
--log_name  Name for the log file.
add_annotation_to_markers.R -h | --help  show this message.
" -> doc

set.seed(1407)

suppressMessages( require(docopt) )
suppressMessages( require(vroom) )
opts <- docopt(doc)

outfile <- paste0("logs/", opts$log)
if ( file.exists(outfile) ) {
    file.remove(outfile)
}

my_log <- file(outfile) 
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

folder_containing_markers <- paste0(
    opts$dir,
    "/")

files <- list.files(folder_containing_markers, pattern = "*.tsv")

annot_names <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv")
annot_summary <- vroom::vroom("MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv")
colnames(annot_names)[1] <- "gene_id"
colnames(annot_summary)[1] <- "gene_id"

annot_summary <- merge( annot_summary,
                        annot_names,
                        by = "gene_id",
                        all.x = T)
annot_summary <- annot_summary[, c(1, 14, 2:13, 15:16)]

for( i in 1:length( files ) ) {
    
    sub <- vroom( file = paste0(folder_containing_markers, files[i] ) )
    sub <- sub[!is.na(sub$gene_id), ]
    
    sub2 <- merge( sub,
                   annot_summary, 
                   by = "gene_id",
                   all.x = T )
    
    file_name <- strsplit(x = files[i], split = "\\.")[[1]][1]
    write.table(sub2,
                paste0(folder_containing_markers,
                       "including_annot_", file_name, ".tsv"),
                col.names = T,
                row.names = F,
                quote = F,
                sep = "\t")
    
}
closeAllConnections()