# This file contains wrapper functions for spatial dimension reduction methods

library(tidyverse)
library(scran)
library(MouseGastrulationData)
library(energy)
library(reticulate)
library(glue)
library(doBy)
library(vctrs)
library(pryr)
library(nnet) # multinom function
library(DescTools) # PseudoR2 function
library(Seurat)

# PRECAST

precast_wrapper <- function(exprs, coords, K = 7, maxIter = 30) {
  seurat_obj <- CreateSeuratObject(counts = t(exprs), project = "SpatialProject")
  row.names(coords) <- rownames(exprs)
  colnames(coords) <- c("row", "col")
  seurat_obj <- AddMetaData(seurat_obj, metadata = coords)
  
  preobj <- CreatePRECASTObject(seuList = list(seurat_obj), 
                                selectGenesMethod = "HVGs", 
                                gene.number = ncol(exprs),
                                verbose = FALSE) 
  PRECASTObj <- AddAdjList(preobj, platform = "ST")
  PRECASTObj <- AddParSetting(PRECASTObj, 
                              Sigma_equal = FALSE, 
                              coreNum = 1, 
                              maxIter = maxIter, 
                              verbose = FALSE)
  
  PRECASTObj <- PRECAST(PRECASTObj, K = K)
  
  seuInt <- PRECASTObj@seulist[[1]]
  seuInt@meta.data$cluster <- factor(unlist(PRECASTObj@resList[[1]]$cluster))
  seuInt@meta.data$batch <- 1
  seuInt <- Add_embed(PRECASTObj@resList[[1]]$hZ[[1]], seuInt, embed_name = "PRECAST")
  posList <- lapply(PRECASTObj@seulist, function(x) cbind(x$row, x$col))
  seuInt <- Add_embed(posList[[1]], seuInt, embed_name = "position")
  Idents(seuInt) <- factor(seuInt@meta.data$cluster)
  return(seuInt@reductions$PRECAST@cell.embeddings)
}

# DRSC

create_adjacency <- function(coords, pct0 = 0.999) {
  dmat <- as.matrix(dist(coords))
  adj_mat <- dmat
  thres_val <- quantile(dmat, 1-pct0)
  adj_mat[dmat < thres_val] <- 1
  adj_mat[dmat >= thres_val] <- 0
  adj_mat <- as(adj_mat, "dgCMatrix")
  return(adj_mat)
}

drsc_wrapper <- function(exprs, coords, pct0 = 0.999, K = 7) {
  adj <- create_adjacency(coords = coords, pct0 = pct0)
  drsc_obj <- DR.SC_fit(X = exprs, Adj_sp = adj, K = K)
  lowdim <- drsc_obj$Objdrsc[[1]]$hZ
  return(lowdim)
}

# spatialPCA

spca_wrapper <- function(counts = NULL, # 1) These need to be genes x cells (standard counts slot sce)
                         # 2) These need to have colnames (cell names)
                         coords = NULL,
                         ncores = 10) {
  spca_res <- CreateSpatialPCAObject(counts = counts, # create SPCA object
                                     location = coords,
                                     numCores_spark = ncores)
  spca_res <- SpatialPCA_buildKernel(spca_res, # build kernel (from spca tutorial)
                                     kerneltype="gaussian", 
                                     bandwidthtype="SJ",
                                     bandwidth.set.by.user=NULL,
                                     sparseKernel = TRUE,
                                     sparseKernel_ncore = ncores)
  spca_res <- SpatialPCA_EstimateLoading(spca_res, # estimate loadings (spca tutorial)
                                         fast=TRUE,
                                         SpatialPCnum=20) 
  spca_res <- SpatialPCA_SpatialPCs(spca_res, fast=TRUE) # calculate spatialPCs (spca tutorial)
  spca_projected <- t(spca_res@SpatialPCs) # extract projected matrix
  kept_cells <- spca_res@normalized_expr %>% colnames() # spca removes some cells - check which are left
  kept_indices <- colnames(counts) %in% kept_cells # find which of original cells are kept
  coords_kept <- coords[kept_indices,] # subset original coords to just get the coords of cells spca kept
  return(list(projected = spca_projected, coords = coords_kept)) # return projected and (subsetted) coords
}

# STAMP

## need to set up your python

use_python("/Library/Frameworks/Python.framework/Versions/3.8/bin/python3")

## wrapper itself

test_file <- "stamp_wrapper.py"
source_python(test_file)

stamp_wrapper_complete <- function(X, coords) {
  out <- stamp_wrapper_fn(counts = as.matrix(X), 
                          coordinates = coords)
  stamp_cells_taken <- as.numeric(out[[2]]) + 1 # +1 since python 0-indexed
  stamp_coords_taken <- coords[stamp_cells_taken,]
  return(list(embedding = out[[1]],
              coords = stamp_coords_taken))
}