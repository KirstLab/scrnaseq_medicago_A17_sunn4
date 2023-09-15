<p align="center">
  <!-- <a href="https://github.com/othneildrew/Best-README-Template">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a> -->

  <!-- <h2 align="center">Asc-Seurat</h2> -->

  <p align="center">
    <h3 align="center"> Scripts to reproduce the scRNA-seq analysis of <i> Medicago truncatula </i> roots</h3>
  </p>
</p>


<!-- ABOUT THE PROJECT -->
## About the _M. truncatula_ scRNA-seq analysis

This repository contains the source code and instructions to reproduce the results of the manuscript "The Single-Cell Transcriptome Program of Nodule Development Cellular Lineages in _Medicago truncatula_" published by Pereira et al., 2023 (add doi once the publication is released).

## Installing software and list of necessary files to execute the analysis.

### Installing the necessary software.

All analyses were conducted using R v.4.1.2 in macOS 13.1. A list of necessary R packages, and the commands to install them, is present in the script _installing_necessary_packages.R_.

To run the script and perform the packages installation, open a terminal and execute the script below:

```sh
Rscript installing_necessary_packages.R
```

## Necessary inputs for the analysis

1. Single-cell RNA-seq data for all samples, which needs to be organized as shown below. To download the dataset, visit (add the link to the datasets once the paper is published.)

```sh
> tree data/

data
├── A17_sep_2022_0h_10k
│   └── outs
│       └── filtered_feature_bc_matrix
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
├── A17_sep_2022_24h_10k
│   └── outs
│       └── filtered_feature_bc_matrix
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
├── A17_sep_2022_48h_10k
│   └── outs
│       └── filtered_feature_bc_matrix
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
├── A17_sep_2022_96h_10k
│   └── outs
│       └── filtered_feature_bc_matrix
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
├── Sunn_sep_2022_0h_10k
│   └── outs
│       └── filtered_feature_bc_matrix
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
├── Sunn_sep_2022_24h_10k
│   └── outs
│       └── filtered_feature_bc_matrix
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
├── Sunn_sep_2022_48h_10k
│   └── outs
│       └── filtered_feature_bc_matrix
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
└── Sunn_sep_2022_96h_10k
    └── outs
        └── filtered_feature_bc_matrix
            ├── barcodes.tsv.gz
            ├── features.tsv.gz
            └── matrix.mtx.gz
```

## Analyse of the genotypes independently

First, combine all timepoints of each genotype and perform the clustering. Then, visualize the expression of marker genes.

```sh
Rscript combining_data_and_clustering_A17.R --cores 8 --UMI 400 --out Medicago_ONLY_A17_integ_cds_clustered.rds
```

```sh
Rscript combining_data_and_clustering_sunn4.R --cores 8 --UMI 400 --out Medicago_ONLY_SUNN4_integ_cds_clustered.rds
```

```sh
Rscript comparing_expression_in_separeted_datasets.R
```

## Combined analysis of all datasets

### Generating the combined dataset and clustering

The first step of the analysis is to perform quality control of the dataset, including removal of empty cells and identification of doublets. Next, the datasets are combined and the clustering is performed.

The script listed below will perform those steps and save the clustered datasets as an rds file.

```sh
Rscript QC_and_doublets_identification.R
Rscript combining_all_samples_and_clustering.R
```

During the execution, a couple of folders will be created to store the datasets and the log files. The folder _rds_files_ will contain the rds files with the clustered datasets. The folder _logs_ will contain the log files generated during the execution of the script.

After the clustering, the top 1000 most specific genes (_de novo_ gene markers) of each cluster are investigated and exported. When less than 1000 genes pass the thresholds, those will be reported instead. Those files are stored at the folder _top_1000_markers_per_cluster/by_specificity_.

### Number of cells per cluster and timepoint -- whole dataset

The script below calculates the number of cells per cluster and time point. The values are shown in cells per thousand, to account for variations in the number of cells due to sampling bias.

```sh
Rscript number_of_cells_per_cluster_over_time.R --rds rds_files/batched_integrated_clustered_complete_dataset.rds --prefix number_of_cells_per_cluster/whole_dataset
```

## Cluster annotation and cell type identification

### Adding gene annotation to the top 1000 markers of each cluster -- whole dataset.
The script below adds the gene annotation to the top-specific genes. It requires the presence of the files: "MtrunA17r5.0-ANR-EGN-r1.9.gene_names.tsv" and "MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv" in the directory. Both were originally obtained at https://medicago.toulouse.inra.fr/MtrunA17r5.0-ANR/ .

```sh
Rscript add_annotation_to_markers.R --dir top_1000_markers_per_cluster/by_specificity/
Rscript add_annotation_to_markers.R --dir top_1000_markers_per_cluster/by_specificity/top100
```

### Visualization of the expression of the top 100 markers of each cluster -- whole dataset.

The command below generates the expression profile of the most specific genes of each cluster in the LCM and other reference datasets. The first command below shows the expression, including the different replicates. The second shows only the averages when more than one replicate is present.

To use this functions, it is necessary to download the external datasets from FigShare (add the link to Figshare).

```sh
Rscript investigating_genes_using_other_sources.R --out_images images/annotation_using_other_sources_images --markers_dir top_1000_markers_per_cluster/by_specificity/top100 --color_scheme VIRIDIS --show_gene_names FALSE --image_format png
```

```sh
Rscript investigating_genes_using_other_sources_2.0.R --out_images images/annotation_using_other_sources_images --markers_dir top_1000_markers_per_cluster/by_specificity/top100 --color_scheme VIRIDIS --show_gene_names FALSE --image_format png
```

Next, we can generate the expression profile of each of the top specific genes per cluster, represented as a UMAP plot per gene.

```sh
mkdir -p images/expression_profiles_top100/
i=1
while [ $i -le 29 ]
do

Rscript visual_inspection_of_markers.R --rds rds_files/batched_integrated_clustered_complete_dataset.rds --markers top_1000_markers_per_cluster/by_specificity/top100/list_of_top_100_for_cluster_${i}_specificity.tsv --out images/expression_profiles_top100/top_100_for_cluster_${i} --common=FALSE

i=$(($i+1))
done
```

### Visualizing the expression profile of key RNS genes reported from other publications in the whole dataset.

```sh
Rscript visual_inspection_of_markers.R --rds rds_files/batched_integrated_clustered_complete_dataset.rds --markers selected_markers/figure_2_roy_et_al_2020.csv --out images/annotation_using_markers/figure_2_roy_et_al_2020 --common=TRUE
```

```sh
Rscript visual_inspection_of_markers.R --rds rds_files/batched_integrated_clustered_complete_dataset.rds --markers selected_markers/markers_Schiessl_2019_table_s2.csv --out images/annotation_using_markers/markers_Schiessl_2019_table_s2 --common=TRUE
```

```sh
Rscript visual_inspection_of_markers.R --rds rds_files/batched_integrated_clustered_complete_dataset.rds --markers selected_markers/RNS_genes_curated_list.tsv --out images/annotation_using_markers/RNS_genes_curated_list --common=TRUE
```

## high-dimention weighted gene co-expression network analysis (hdWGCNA)

Note, this step demans a large amoung of RAM memory. We recommend using a cluster with at least 150GB of RAM.

WGCNA is a powerful tool to identify modules of co-expressed genes. The output of this script is a list of modules, and the genes that are part of each module. 

```sh
Rscript Create_Seurat_object_from_monocle3.R
Rscript hdWGCNA.R
```

Next, we can subset the generated network to filter the genes with the highest connectivity (stronger co-expression) to any gene of interest. The script below will generate a list of top N genes with the highest connectivity to the gene of interest. 

```sh
Rscript subseting_the_network.R
```

Nothe that the output list can also be used as input to CytoScape to visualization and to generate a network plot.

## RNA velocity

```sh

```

<!-- ### Expression profile of key markers using heatmap -- whole dataset

Generates the expression of genes listed by Schiess et al. 2019 and Roy et al.,2020.
```sh
Rscript generates_expression_of_key_genes_as_heatmap.R
``` -->

## Reclustering analysis of the pericycle cells

The scripts below will extract the cells from the pericycle cluster, then recluster them using a higher resolution.

```sh
Rscript reclustering_pericycle.R 
```

```sh
Rscript add_annotation_to_markers.R --dir RECLUSTERING/PERICYCLE/top_1000_markers_per_cluster/by_specificity
```


## Reclustering analysis of the root hair and IT

```sh
Rscript reclustering_root_hair_IT.R
```

```sh
Rscript add_annotation_to_markers.R --dir RECLUSTERING/Epidermis_roothair/top_1000_markers_per_cluster/by_specificity
```

## Reclustering analysis of cortex and nodule cells (sunn-4)

```sh
Rscript open reclustering_cortex_nodule.R
```

```sh
Rscript add_annotation_to_markers.R --dir RECLUSTERING/CORTEX_NOD_sunn4/top_1000_markers_per_cluster/by_specificity
```

```sh
Rscript investigating_genes_using_other_sources.R --out_images RECLUSTERING/CORTEX_NOD_sunn4/annotation_using_other_sources_images --markers_dir RECLUSTERING/CORTEX_NOD_sunn4/top_1000_markers_per_cluster/by_specificity/ --color_scheme VIRIDIS --show_gene_names FALSE --image_format png
```

```sh
Rscript visual_inspection_of_markers.R --rds RECLUSTERING/CORTEX_NOD_sunn4/rds_file_subset/medicago_integrated_subset_cortex_nodule.rds  --markers selected_markers/RNS_genes_curated_list.tsv --out RECLUSTERING/CORTEX_NOD_sunn4/annotation_using_markers/RNS_genes_curated_list --common=TRUE
```

```sh
Rscript visual_inspection_of_markers.R --rds RECLUSTERING/CORTEX_NOD_sunn4/rds_file_subset/medicago_integrated_subset_cortex_nodule.rds --markers selected_markers/figure_2_roy_et_al_2020.csv --out RECLUSTERING/CORTEX_NOD_sunn4/annotation_using_markers/figure_2_roy_et_al_2020 --common=TRUE
```

```sh
Rscript visual_inspection_of_markers.R --rds RECLUSTERING/CORTEX_NOD_sunn4/rds_file_subset/medicago_integrated_subset_cortex_nodule.rds --markers selected_markers/markers_Schiessl_2019_table_s2.csv --out RECLUSTERING/CORTEX_NOD_sunn4/annotation_using_markers/markers_Schiessl_2019_table_s2 --common=TRUE
```

```sh
Rscript visual_inspection_of_markers.R --rds RECLUSTERING/CORTEX_NOD_sunn4/rds_file_subset/medicago_integrated_subset_cortex_nodule.rds --markers selected_markers/GusGenes.csv --out RECLUSTERING/CORTEX_NOD_sunn4/annotation_using_markers/GusGenes --common=FALSE
```


```sh
mkdir -p RECLUSTERING/CORTEX_NOD_sunn4/expression_profiles_top1000/
i=1
while [ $i -le 8 ]
do

Rscript visual_inspection_of_markers.R --rds RECLUSTERING/CORTEX_NOD_sunn4/rds_file_subset/medicago_integrated_subset_cortex_nodule.rds --markers RECLUSTERING/CORTEX_NOD_sunn4/top_1000_markers_per_cluster/by_specificity/list_of_top_1000_for_cluster_${i}_specificity.csv --out RECLUSTERING/CORTEX_NOD_sunn4/expression_profiles_top1000/top_1000_for_cluster_${i} --common=FALSE

Rscript visual_inspection_of_markers.R --rds rds_files/batched_integrated_clustered_complete_dataset.rds --markers RECLUSTERING/CORTEX_NOD_sunn4/top_1000_markers_per_cluster/by_specificity/list_of_top_1000_for_cluster_${i}_specificity.csv --out RECLUSTERING/CORTEX_NOD_sunn4/expression_profiles_top1000/top_1000_for_cluster_${i}_in_the_whole_dataset --common=FALSE

i=$(($i+1))
done
```

### Trajectory inference on the reclustered dataset

```sh
Rscript trajectory_inference_cortex_nodule.R
```

### Number of cells per cluster and timepoint -- cortex/nodule

The script below calculates the number of cells per cluster and time point. The values are shown in cells per thousand, to account for variations in the number of cells due to sampling bias.

```sh
Rscript number_of_cells_per_cluster_over_time.R --rds RECLUSTERING/CORTEX_NOD_sunn4/rds_file_subset/medicago_integrated_subset_cortex_nodule.rds --prefix cortex_nodule
```

## Additional analysis

### Generates heatmap for Figure 1.

```sh
Rscript generates_the_heatmap_containing_gene_markers_for_figure1.R
```

### Expression profiles of LBDs, STYs, and YUCCAS

```sh
Rscript LBDs_STYs_and_YUCCAs.R
```

### Generates UMAP plot coloring cells by the cell type

```sh
Rscript umap_colored_by_cell_type.R
```

### t-test number of nodules per plant _MtSTY4_

```sh
Rscript t_test_number_nodules_per_plants_RNAi_STY4.R
```
