import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
import scvelo as scv
from tqdm import tqdm as tqdm

sample_name = "whole_sample_dynamic_model"

# load sparse matrix:
X = io.mmread("./tmp_whole/counts.mtx")

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("./tmp_whole/metadata.csv")
cell_meta['seurat_clusters'] = cell_meta['seurat_clusters'].astype("category")

# load gene names:
with open("./tmp_whole/gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("./tmp_whole/pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color=['seurat_clusters'], frameon=False, save="_"+sample_name+"_plot_to_compare_to_original.svg", show=False)
sc.pl.umap(adata, color=['seurat_clusters'], frameon=False, save="_"+sample_name+"_plot_to_compare_to_original.png", show=False)

# Loads loom file (velocyto)
loom_obj = scv.read("data/velocyto/merged_samples.loom", cache=True, cleanup=True, validate=False)

# merge matrices into the original adata object
adata = scv.utils.merge(adata, loom_obj)
del loom_obj

# plot umap to check
sc.pl.umap(adata, color='seurat_clusters', legend_loc='on data', title='', frameon=False, save="_"+sample_name+"_plot_to_compare_to_original_after_merging_loom.svg", show=False)
sc.pl.umap(adata, color='seurat_clusters', legend_loc='on data', title='', frameon=False, save="_"+sample_name+"_plot_to_compare_to_original_after_merging_loom.png", show=False)

# proportion of spliced/unspliced reads
scv.pl.proportions(adata, groupby='seurat_clusters', save="_"+sample_name+"_splicing.svg", show=False, figsize = (6,6))
scv.pl.proportions(adata, groupby='seurat_clusters', save="_"+sample_name+"_splicing.png", show=False, figsize = (6,6))

# pre-process
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

# compute velocity
scv.tl.recover_dynamics(adata, n_jobs = 48)
scv.tl.velocity(adata, mode='dynamical', n_jobs = 48)
scv.tl.velocity_graph(adata, n_jobs = 48)

scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save="_"+sample_name+'embedding.svg')
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save="_"+sample_name+'embedding.png')

scv.pl.velocity_embedding_grid(adata, basis='umap', color='seurat_clusters', save="_"+sample_name+'embedding_grid.svg', title='')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='seurat_clusters', save="_"+sample_name+'embedding_grid.png', title='')

scv.pl.velocity_embedding_stream(adata, basis='umap', color='seurat_clusters', save="_"+sample_name+'embedding_stream.svg', title='')
scv.pl.velocity_embedding_stream(adata, basis='umap', color='seurat_clusters', save="_"+sample_name+'embedding_stream.png', title='')

scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)
df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()

pd.DataFrame.to_csv(df, "_"+sample_name+"rank_of_dynamic_genes_according_to_velocity.csv")
           
scv.pl.scatter(adata, df['6'][:5], ylabel='1', frameon=False, color='seurat_clusters', size=10, linewidth=1.5, save="_"+sample_name+'top5_ranked_scatter.svg')
scv.pl.scatter(adata, df['6'][:5], ylabel='1', frameon=False, color='seurat_clusters', size=10, linewidth=1.5, save="_"+sample_name+'top5_ranked_scatter.png')

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save="_"+sample_name+'velocity_confidence.svg')
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save="_"+sample_name+'velocity_confidence.png')

scv.pl.velocity_graph(adata, threshold=.1, color='seurat_clusters', save="_"+sample_name+'velocity_graph.svg')
scv.pl.velocity_graph(adata, threshold=.1, color='seurat_clusters', save="_"+sample_name+'velocity_graph.png')

x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=70)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax, save="_"+sample_name+'velocity_graph_transitions_scatter.svg')
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax, save="_"+sample_name+'velocity_graph_transitions_scatter.png')

scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save="_"+sample_name+'velocity_pseudotime.svg')
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save="_"+sample_name+'velocity_pseudotime.png')