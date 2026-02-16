from gpca import *
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import squidpy as sq
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances as pair
from sklearn.metrics import adjusted_rand_score as ari_score


def graphPCA_run(counts, locations):
    # counts = gene expression matrix in cell x gene format (otherwise transpose matrix first)
    # locations = spatial coordinates matrix
    adata = ad.AnnData(X = counts)
    adata.obsm['spatial'] = locations
    print("anndata file created")
    sc.experimental.pp.normalize_pearson_residuals(adata)
    print("normalisation completed.")
    sc.pp.scale(adata)
    print("scaling completed.")
    
    Z,W = Run_GPCA(adata, location = locations, n_components = 50, method = 'knn', platform = 'ST', _lambda = 0.5, save_reconstruction = True)
    print("graphPCA run completed.")
    adata.uns["GraphPCA"] = Z
    adata.uns["GraphPCA_loadings"] = W
    print("graphPCA output stored in adata.")
    
    return adata





# def graphPCA_run(sample, geneNames = "hvgs", sampleList = None, hvgs = 2000, platform = "Visium", nn = None):
    # sample = path to h5ad file, including filename
    # sampleList = list of h5ad filepaths, including file names
    # geneNames = whether to use hvgs or all genes 
    # hvgs = number of hvgs to use, default is 2000
    # platform = to use Visium or ST - this automatically defines the n_neighbours
    # nn = defines n_neighbours, by default automatically chosen for a platform
    # 
    # adata = ad.AnnData(X = counts)
    # adata.obsm['spatial'] = locations
    # print("anndata file created")
    # adata = sc.read_h5ad(sample)
    # adata.var_names_make_unique()
    # sc.pp.filter_genes(adata, min_cells = 20)
    # sc.experimental.pp.normalize_pearson_residuals(adata)
    # print("normalisation completed.")
    # sc.pp.scale(adata)
    # print("scaling completed.")
    # 
    # if geneNames == "hvgs":
    #     geneList = sc.pp.highly_variable_genes(adata, flavor = 'seurat_v3', n_top_genes = hvgs)
    #     adata = adata[:, geneList]
    # adata = adata[:, 0:5000]
    
    # x_array = adata.obs["array_row"].tolist()
    # y_array = adata.obs["array_col"].tolist()
    # location = adata.obsm["X_spatial"].astype(np.float32)

    # Z,W = Run_GPCA(adata, location = location, n_components = 50, method = "knn", platform = platform, _lambda = 0.5, n_neighbors = nn, save_reconstruction = True)
    # print("graphPCA run completed.")
    # adata.obsm["GraphPCA"] = Z
    # adata.obsm["GraphPCA_loadings"] = W
    # print("graphPCA output stored in adata.")
    # 
    # return adata
