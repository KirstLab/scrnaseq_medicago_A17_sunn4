set.seed(1407)
system("mkdir -p logs")
system("mkdir -p images/QC")

outfile <- "logs/removing_empty_cells.out"
if ( file.exists(outfile) ) {
    file.remove(outfile)
}

my_log <- file(outfile) 
sink(my_log, append = TRUE, type = "output")
sink(my_log, append = TRUE, type = "message")

suppressMessages( require(docopt) )
suppressMessages(require(tidyverse) )
suppressMessages(require(DropletUtils) )
suppressMessages(require(scuttle) )
suppressMessages(require(scater) )

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 7

samples <- c("A17_0h",
             "A17_24h",
             "A17_48h",
             "A17_96h",
             "Sunn_0h",
             "Sunn_24h",
             "Sunn_48h",
             "Sunn_96h")

fdr <- 0.05

for ( s in samples) {
    
    sce <- DropletUtils::read10xCounts(
        paste0("data/", s, "/raw_feature_bc_matrix/"),
        col.names = T)
    mol.info <- read10xMolInfo(paste0("data/", s, "/molecule_info.h5") )
    
    br.out <- DropletUtils::barcodeRanks(sce)
    
    png( paste0("images/QC/knee_plot_", s, ".png"),
         width = 10, height = 6, units = "in", res = 300)
    
    par(mar = c(1, 1, 1, 1))
    plot(br.out$rank, br.out$total, log = "xy", xlab = "Rank", ylab = "Total")
    o <- order(br.out$rank)
    lines(br.out$rank[o], br.out$fitted[o], col = "red")
    
    abline(h = metadata(br.out)$knee, col = "dodgerblue", lty = 2)
    abline(h = metadata(br.out)$inflection, col = "forestgreen", lty = 2)
    legend("bottomleft",
           lty = 2, col = c("dodgerblue", "forestgreen"),
           legend = c("knee", "inflection")
    )
    
    dev.off()
    
    empty <- DropletUtils::emptyDrops(sce,
                                      BPPARAM = BPPARAM,
                                      niters = 50000)
    
    is.cell <- empty$FDR <= fdr
    (empty_and_limited <- table(Limited = empty$Limited,
                                Significant = is.cell))
    
    if (empty_and_limited[2, 1] > 0) {
        print("Warning: there are cells that will benefit of more iterations. Increasing the number of iterations to 100000")
        
        empty <- DropletUtils::emptyDrops(sce, niters = 100000)
        is.cell <- empty$FDR <= opts$FDR
        empty_and_limited <- table(Limited = empty$Limited, Significant = is.cell)
        empty_and_limited
        
        if (empty_and_limited[2, 1] > 0) {
            print("warning: there are still cells that would benefit of more iterations after 100000 iterations!!")
        }
    }
    
    png(paste0("images/QC/probability_plot_", s, ".png"),
        width = 5, height = 5, units = "in", res = 300)
    plot(empty$Total, -empty$LogProb,
         col = ifelse(is.cell, "blue", "red"),
         xlab = "Total UMI count", ylab = "-Log Probability"
    )
    dev.off()
    
    sce <- sce[, which(empty$FDR <= fdr)]
    sce_counts <- counts(sce)
    rownames(sce_counts) <- gsub(pattern = "MtrunA17_",
                                 replacement = "MtrunA17",
                                 x = rownames(sce_counts))
    
    print(paste("Dimension of the dataset:", s))
    print( dim(sce_counts) )
    folder_name <- paste0("data/filtered_",
                     s,"/outs/filtered_feature_bc_matrix")
    system( paste0("mkdir -p ", "data/filtered_",
                                      s,"/outs/filtered_feature_bc_matrix") )
    DropletUtils::write10xCounts(x = sce_counts,
                                 path = folder_name,
                                 version = "3", overwrite = T)
    
}

closeAllConnections()