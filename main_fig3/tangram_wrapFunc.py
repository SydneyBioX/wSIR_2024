#!/usr/bin/env python

import os, sys
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import torch
import tangram as tg
import anndata as ad

# sys.path.append('/dski/nobackup/pratibha//pythonVENVs/TANGRAM_analysis/.venv/bin')
# requires python 3.8

def tangram_run(spe_counts, spe_locs, sce_counts):
    '''
    function to predict spatial coordinates for single cell data.
    Returns a data frame of predicted coordinates.
    '''
    
    # spe_counts, spe_locs, and sce_counts must be matrices
    
    ad_sp = ad.AnnData(X = spe_counts)
    ad_sp.obsm['spatial'] = spe_locs
    print("spatial anndata file created.")
    ad_sc = ad.AnnData(X = sce_counts)
    print("single cell anndata file created.")
    
    # use all genes to train data
    tg.pp_adatas(ad_sc, ad_sp, genes = None)
    print("Used all features for training.")
    # checking if genes in the training_genes fields 
    # of both annDatas follow the same order.
    # no value returned if assertion is met.
    assert ad_sc.uns['training_genes'] == ad_sp.uns['training_genes']
    
    # aligning sc data to reference spatial data
    ad_map = tg.map_cells_to_space(
        adata_sc = ad_sc, 
        adata_sp = ad_sp,
        num_epochs = 500, 
        mode = 'cells',
        # device = 'cuda:0', # for running using gpu
        device = 'cpu', # for running on cpu
        # density_prior='rna_count_based' # for sequencing-based dataset
        density_prior = 'uniform')
    print("Tangram mapped query cells to reference.")
    
    # predicting coordinates for sc data
    mapping_mat = ad_map.X # probability matrix
    sp_coords = ad_sp.obsm['spatial'] # reference coordinates
    pred_coords = [] # query coordinates to be predicted
    for i in range(mapping_mat.shape[0]):
        cur_mapping_probs = mapping_mat[i]
        cur_ind = np.argmax(cur_mapping_probs)
        cur_pred_coords = np.take(sp_coords, cur_ind, axis = 0)
        pred_coords.append(cur_pred_coords)
    pred_coords = np.asarray(pred_coords)
    print("Query coordinates predicted.")
    print("Tangram run finished.")
    return pred_coords
