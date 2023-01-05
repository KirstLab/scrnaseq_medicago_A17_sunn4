## Packages from CRAN
install.packages("assertthat")
install.packages("circlize")
install.packages("cowplot")
install.packages("docopt")
install.packages("dplyr")
install.packages("dynfeature")
install.packages("dynplot")
install.packages("dynutils")
install.packages("dynwrap")
install.packages("ggplot2")
install.packages("ggthemes")
install.packages("RColorBrewer")
install.packages("scales")
install.packages("Seurat")
install.packages("stringr")
install.packages("tidyverse")
install.packages("viridis")
install.packages("vroom")

## Packages from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

## Required for monocle3
BiocManager::install( c(
    'BiocGenerics',
    'batchelor',
    'DelayedArray',
    'DelayedMatrixStats',
    'limma',
    'Matrix.utils',
    'S4Vectors',
    'SingleCellExperiment',
    'SummarizedExperiment',
    "slingshot",
    "tradeSeq"
))

## Packages from devtools
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github('dynverse/dyno')