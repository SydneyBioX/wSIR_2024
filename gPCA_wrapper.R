graphPCA_run <- function(cntPath, locPath){
    require(reticulate)
    require(here)
    
    # cntPath = path to counts.Rds - must be a matrix of gene x cell format
    # locPath = path to coords.Rds - matrix of spatial coordinates
    
    # import python packages
    py_run_string("import scanpy as sc")
    py_run_string("import anndata as ad")
    py_run_string("import numpy as np")
    py_run_string("import pandas as pd")
    py_run_string("import squidpy as sq")
    py_run_string("from sklearn.cluster import KMeans")
    py_run_string("from sklearn.metrics import pairwise_distances as pair")
    py_run_string("from sklearn.metrics import adjusted_rand_score as ari_score")
    py_run_string("import pyreadr")
    
    # read count and locations Rds files
    py_run_string(paste0("counts = pyreadr.read_r('", cntPath, "')[None]"))
    py_run_string(paste0("locations = pyreadr.read_r('", locPath, "')[None]"))
    print("Counts and location files read.")
    
    # convert count and locations variables to matrices (numpy array in python),
    # transpose count matrix for anndata object, and 
    # create anndata object from them 
    py_run_string(paste0("adata = ad.AnnData(X = counts.to_numpy().transpose())"))
    py_run_string(paste0("adata.obsm['spatial'] = locations.to_numpy()"))
    print("anndata file created")
    
    # normalise data as suggested for graphPCA
    py_run_string(paste0("sc.experimental.pp.normalize_pearson_residuals(adata)"))
    print("normalisation completed.")
    
    # scale the data
    py_run_string(paste0("sc.pp.scale(adata)"))
    print("scaling completed.")
    
    # run graphPCA to generate directions and loadings
    py_run_string(paste0("Z, W = Run_GPCA(adata, location = locations, n_components = 50, method = 'knn', platform = 'ST', _lambda = 0.5, save_reconstruction = True)"))
    print("graphPCA run completed.")
    
    # export directions and loadings to R environment from python environment
    Z = py_to_r(py_eval("Z"))
    W = py_to_r(py_eval("W"))
    
    # return directions and loadings
    return(list("Z" = Z, "W" = W))
}