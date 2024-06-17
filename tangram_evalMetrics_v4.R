## Tangram - predicting spatial coordinates for single cell data
tangram_run <- function(sp_filename, sc_filename){
  # requires Python 3.8
  require(reticulate)
  require(here)
  require(reshape)
  require(dplyr)
  require(tidyr)
  
  py_run_string("import numpy as np")
  py_run_string("import scanpy as sc")
  py_run_string("import tangram as tg")
  
  py_run_string(paste0("ad_sp = sc.read_h5ad('", sp_filename, "')"))
  py_run_string(paste0("ad_sc = sc.read_h5ad('", sc_filename, "')"))
  py_run_string("tg.pp_adatas(ad_sc, ad_sp, genes = None)")
  py_run_string("ad_map = tg.map_cells_to_space(adata_sc = ad_sc, adata_sp = ad_sp, num_epochs = 1000, device = 'cuda:0', density_prior = 'uniform')")
  
  runs = py_eval("ad_map.X.shape[0]")
  py_run_string("mapping_mat = ad_map.X")
  py_run_string("sp_coords = np.column_stack((ad_sp.obs['X'].values, ad_sp.obs['Y'].values))")
  py_run_string("pred_coords = []")
  for (i in 1:runs){
    py_run_string(paste0("cur_mapping_probs = mapping_mat[",i-1,"]"))
    py_run_string("cur_ind = np.argmax(cur_mapping_probs)")
    py_run_string("cur_pred_coords = np.take(sp_coords, cur_ind, axis = 0)")
    py_run_string("pred_coords.append(cur_pred_coords)")
  }
  py_run_string("pred_coords = np.asarray(pred_coords)")
  pred_coords = py_to_r(py_eval("pred_coords"))
  return(pred_coords)
}



eval_metrics <- function(true_coords, pred_coords, celltypes, neighbor_num, 
                         neighHits = TRUE, JSD = TRUE, spearCoef = TRUE){
  require(BiocNeighbors)
  require(reshape)
  require(dplyr)
  require(tidyr)
  
  lim = nrow(pred_coords)
  true_neigh_data = findKNN(true_coords, k = lim-1, BNPARAM = KmknnParam())
  true_dist = cbind(rep(0, times = lim), true_neigh_data$distance)
  true_ind = cbind(c(1:lim), true_neigh_data$index)
  pred_neigh_data <- findKNN(pred_coords, k = lim-1, BNPARAM = KmknnParam())
  pred_dist = cbind(rep(0, times = lim), pred_neigh_data$distance)
  pred_ind = cbind(c(1:lim), pred_neigh_data$index)
  
  
  if (neighHits == TRUE) {
    # function to calculate the number of similar neighbourhood
    # cells observed in ground truth vs predicted coordinates.
    neighbor_hits_all = data.frame(matrix(nrow = 0, ncol = 3))
    colnames(neighbor_hits_all) <- c("CellNum", "NN", "NumHits")
    avg_neighbor_hits = data.frame()
    for (k in seq(20, 220, by = 20)){
      neighbor_hits = list()
      for (i in c(1:lim)) {
        cur_true_neighbor_ind = as.character(unique(true_ind[i, 2:k]))
        cur_pred_neighbor_ind = as.character(unique(pred_ind[i, 2:k]))
        overlap_neighbors = length(intersect(cur_true_neighbor_ind, cur_pred_neighbor_ind))
        neighbor_hits[i] = overlap_neighbors
      }
      neighbor_hits_all = rbind(neighbor_hits_all, 
                                cbind("CellNum" = c(1:lim), 
                                      "NN" = rep(paste0('NN_', k), times = lim), 
                                      "NumHits" = unlist(neighbor_hits)))
      avg_neighbor_hits = rbind(avg_neighbor_hits, cbind("NN_select" = paste0('NN_', k),
                                                         "NN" = k,
                                                         "AvgHits" = mean(unlist(neighbor_hits)),
                                                         "prop_AvgHits" = mean(unlist(neighbor_hits)) / k))
    }
    neighbor_hits_all$CellNum = as.numeric(neighbor_hits_all$CellNum)
    neighbor_hits_all = cast(neighbor_hits_all, CellNum~NN) %>%
      select(avg_neighbor_hits$NN_select)
    avg_neighbor_hits$NN_select = factor(avg_neighbor_hits$NN_select, levels = avg_neighbor_hits$NN_select)
  } else {
    neighbor_hits_all = NULL
    avg_neighbor_hits = NULL
  }
  
  
  if(JSD == TRUE) {
    # function to calculate Jensen-Shannon distance to check
    # how well the local cell-type heterogeneity is preserved.
    celltype_set = unique(celltypes)
    celltype_df = data.frame("celltypes" = celltype_set)
    jsd_vect = c()
    for (i in c(1:lim)) {
      true_counts = celltypes[as.vector(true_ind[i, 2:neighbor_num])] %>%
        table() %>%
        prop.table() %>%
        as.data.frame()
      colnames(true_counts) = c("celltypes", "prop")
      true_data = merge(celltype_df, true_counts, by = "celltypes", row.names = NULL, all = TRUE) %>%
        mutate(prop = replace_na(prop, 0))
      pred_counts = celltypes[as.vector(pred_ind[i, 2:neighbor_num])] %>% 
        table() %>%
        prop.table() %>%
        as.data.frame()
      colnames(pred_counts) = c("celltypes", "prop")
      pred_data = merge(celltype_df, pred_counts, by = "celltypes", all = TRUE) %>%
        mutate(prop = replace_na(prop, 0))
      input_data = rbind(pred_data$prop, true_data$prop)
      jsd_val = sqrt(jensen_difference(pred_data$prop, true_data$prop, 
                                       testNA = TRUE, unit = "log2"))
      jsd_vect = append(jsd_vect, jsd_val)
    }
  } else {
    jsd_vect = NULL
  }
  
  
  if (spearCoef == TRUE) {
    # function to calculate Spearman's rank correlation coefficient.
    # to check whether global cell-cell relationships can be captured.
    true_rev_ind = apply(true_ind, 1, function(x){
      y = data.frame(index = c(1:length(x)), cell_index = x) %>%
        arrange(cell_index)
      return(y$index)
    })
    tmpT_ind = melt(true_rev_ind)[, c(1,3)]
    true_dist_un = matrix(true_dist[as.matrix(tmpT_ind)],
                          nrow = nrow(true_dist),
                          ncol = ncol(true_dist)) 
    pred_rev_ind = apply(pred_ind, 1, function(x){
      y = data.frame(index = c(1:length(x)), cell_index = x) %>%
        arrange(cell_index)
      return(y$index)
    })
    tmpP_ind = melt(pred_rev_ind)[, c(1,3)]
    pred_dist_un = matrix(pred_dist[as.matrix(tmpP_ind)],
                          nrow = nrow(pred_dist), 
                          ncol = ncol(pred_dist))
    spearman_df = data.frame(matrix(nrow = 0, ncol = 2))
    colnames(spearman_df) = c("Correlation", "P_value")
    for (i in c(1:lim)) {
      corr = cor.test(pred_dist_un[i,], true_dist_un[i,], method = 'spearman')
      spearman_df = rbind(spearman_df, cbind("Correlation" = corr$estimate, 
                                             "P_value" = corr$p.value))
    }
  } else {
    spearman_df = NULL
  }
  
  return (list("all_neighHits" = neighbor_hits_all, 
               "avg_neighHits" = avg_neighbor_hits,
               "JSdist" = jsd_vect,
               "spearmanCoef" = spearman_df))
}
