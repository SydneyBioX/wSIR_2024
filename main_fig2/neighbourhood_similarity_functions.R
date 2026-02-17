
# show jaccard values for neighbourhood similarity for varying neighbourhood sizes
# computes jaccard index for all cells rather than average over all cells
# jaccard_compute is better for single neighbourhoodd size value

# means input: if you want a line plot of the mean of each DR method at each neighbourhood size, or a boxplot over all cells
# verbose: print progress update as you move to each neighbourhood size

jaccard_from_lowdim <- function(coords_test, 
                                low_dim_obj, 
                                n_vals, 
                                means = TRUE, 
                                verbose = TRUE,
                                sd = TRUE) {
  dr_methods <- names(low_dim_obj)
  
  dist_true = dist(coords_test) # create true distance matrix
  
  # initialise results df, fill neighbourhood_size and dr_method columns (different cases whether you want just means or all cells)
  
  if (means) {
    jc_df <- matrix(NA, nrow = length(n_vals)*length(dr_methods), ncol = 4) %>% as.data.frame()
    colnames(jc_df) <- c("jaccard", "neighbourhood_size", "dr_method", "sd")
    jc_df$neighbourhood_size <- as.factor(vec_rep_each(n_vals, length(dr_methods))) # need package "vctrs" for this
    jc_df$dr_method <- rep(dr_methods, length(n_vals))
  } else {
    jc_df <- matrix(NA, nrow = nrow(coords_test)*length(n_vals)*length(dr_methods), ncol = 3) %>% as.data.frame()
    colnames(jc_df) <- c("jaccard", "neighbourhood_size", "dr_method")
    jc_df$neighbourhood_size <- as.factor(vec_rep_each(n_vals, nrow(coords_test)*length(dr_methods))) # need package "vctrs" for this
    jc_df$dr_method <- rep(vec_rep_each(dr_methods, nrow(coords_test)), length(n_vals))
  }
  
  # create low-dim distance matrices
  
  if ("sir" %in% dr_methods) {
    dist_dr_sir <- dist(low_dim_obj$sir[[2]]) %>% as.matrix()
  }
  if ("wsir" %in% dr_methods) {
    dist_dr_wsir <- dist(low_dim_obj$wsir[[2]]) %>% as.matrix()
  }
  if ("pca" %in% dr_methods) {
    dist_dr_pca <- dist(low_dim_obj$pca[[2]]) %>% as.matrix()
  }
  if ("lda" %in% dr_methods) {
    dist_dr_lda <- dist(low_dim_obj$lda[[2]]) %>% as.matrix()
  }
  if ("pls" %in% dr_methods) {
    dist_dr_pls <- dist(low_dim_obj$pls[[2]]) %>% as.matrix()
  }
  
  i = 0
  for (n in n_vals) {
    
    if (verbose) {
      print(paste("Starting neighbourhood size:", n))
    }
    
    true_neighbours <- find_neighbours(dist_true = dist_true, n = n) # compute n true neighbours
    
    jc_vals <- c()
    jc_means <- c()
    jc_sds <- c()
    if ("sir" %in% dr_methods) {
      jc_sir <- jaccard_from_dmat(neighbours_true = true_neighbours, dist_dr = dist_dr_sir)
      jc_vals <- jc_vals %>% append(jc_sir)
      jc_sir_mean <- mean(jc_sir)
      jc_means <- jc_means %>% append(jc_sir_mean)
      jc_sir_sd <- sd(jc_sir)
      jc_sds <- jc_sds %>% append(jc_sir_sd)
      
    }
    if ("wsir" %in% dr_methods) {
      jc_wsir <- jaccard_from_dmat(neighbours_true = true_neighbours, dist_dr = dist_dr_wsir)
      jc_vals <- jc_vals %>% append(jc_wsir)
      jc_wsir_mean <- mean(jc_wsir)
      jc_means <- jc_means %>% append(jc_wsir_mean)
      jc_wsir_sd <- sd(jc_wsir)
      jc_sds <- jc_sds %>% append(jc_wsir_sd)
    }
    if ("pca" %in% dr_methods) {
      jc_pca <- jaccard_from_dmat(neighbours_true = true_neighbours, dist_dr = dist_dr_pca)
      jc_vals <- jc_vals %>% append(jc_pca)
      jc_pca_mean <- mean(jc_pca)
      jc_means <- jc_means %>% append(jc_pca_mean)
      jc_pca_sd <- sd(jc_pca)
      jc_sds <- jc_sds %>% append(jc_pca_sd)
    }
    if ("lda" %in% dr_methods) {
      jc_lda <- jaccard_from_dmat(neighbours_true = true_neighbours, dist_dr = dist_dr_lda)
      jc_vals <- jc_vals %>% append(jc_lda)
      jc_lda_mean <- mean(jc_lda)
      jc_means <- jc_means %>% append(jc_lda_mean)
      jc_lda_sd <- sd(jc_lda)
      jc_sds <- jc_sds %>% append(jc_lda_sd)
    }
    if ("pls" %in% dr_methods) {
      jc_pls <- jaccard_from_dmat(neighbours_true = true_neighbours, dist_dr = dist_dr_pls)
      jc_vals <- jc_vals %>% append(jc_pls)
      jc_pls_mean <- mean(jc_pls)
      jc_means <- jc_means %>% append(jc_pls_mean)
      jc_pls_sd <- sd(jc_pls)
      jc_sds <- jc_sds %>% append(jc_pls_sd)
    }
    
    if (means) {
      jc_df$jaccard[(i*length(dr_methods)+1):((i+1)*length(dr_methods))] <- jc_means
      jc_df$sd[(i*length(dr_methods)+1):((i+1)*length(dr_methods))] <- jc_sds
    } else {
      jc_df$jaccard[(i*nrow(coords_test)*length(dr_methods)+1):((i+1)*nrow(coords_test)*length(dr_methods))] <- jc_vals 
    }
    
    i = i + 1
  }
  
  # make plot
  
  if (means) {
    plot <- ggplot(aes(x = neighbourhood_size, y = jaccard, colour = dr_method, group = dr_method), data = jc_df) + geom_line()
  } else {
    plot <- ggplot(aes(x = neighbourhood_size, y = jaccard, fill = dr_method), data = jc_df) + geom_boxplot()
  }
  
  return(list(plot = plot, 
              dataframe = jc_df))
}

# for a single neighbourhood size

jaccard_compute <- function(low_dim_obj,
                            coords_test,
                            n_size = 50) {
  
  dr_methods <- names(low_dim_obj)
  dist_true = dist(coords_test) # create true distance matrix
  
  true_neighbours <- find_neighbours(dist_true = dist_true, n = n_size) # compute n true neighbours
  
  jc_vals <- list()
  if ("sir" %in% dr_methods) {
    dist_dr_sir <- dist(low_dim_obj$sir[[2]]) %>% as.matrix()
    jc_sir <- jaccard_from_dmat(neighbours_true = true_neighbours, dist_dr = dist_dr_sir)
    jc_vals$sir <- mean(jc_sir)
    
  }
  if ("wsir" %in% dr_methods) {
    dist_dr_wsir <- dist(low_dim_obj$wsir[[2]]) %>% as.matrix()
    jc_wsir <- jaccard_from_dmat(neighbours_true = true_neighbours, dist_dr = dist_dr_wsir)
    jc_vals$wsir <- mean(jc_wsir)
  }
  if ("pca" %in% dr_methods) {
    dist_dr_pca <- dist(low_dim_obj$pca[[2]]) %>% as.matrix()
    jc_pca <- jaccard_from_dmat(neighbours_true = true_neighbours, dist_dr = dist_dr_pca)
    jc_vals$pca <- mean(jc_pca)
  }
  if ("lda" %in% dr_methods) {
    dist_dr_lda <- dist(low_dim_obj$lda[[2]]) %>% as.matrix()
    jc_lda <- jaccard_from_dmat(neighbours_true = true_neighbours, dist_dr = dist_dr_lda)
    jc_vals$lda <- mean(jc_lda)
  }
  if ("pls" %in% dr_methods) {
    dist_dr_pls <- dist(low_dim_obj$pls[[2]]) %>% as.matrix()
    jc_pls <- jaccard_from_dmat(neighbours_true = true_neighbours, dist_dr = dist_dr_pls)
    jc_vals$pls <- mean(jc_pls)
  }
  
  return(jc_vals)
}

# jaccard between two sets

jaccard <- function(a,b) {
  int <- length(intersect(a,b))
  union <- length(a) + length(b) - int
  return(int/union)
}

# jaccard from distance matrix and true neighbours set

jaccard_from_dmat <- function(neighbours_true, dist_dr) {
  dist_dr <- dist_dr %>% as.matrix()
  n_obs <- length(neighbours_true)
  n <- length(neighbours_true[[1]])
  
  jc_dr <- rep(0, n_obs)
  
  for (i in 1:n_obs) {
    neighbours_current <- neighbours_true[[i]]
    neighbours_dr <- which.minn(dist_dr[i,], n)[-1]
    jc_dr[i] <- jaccard(neighbours_dr, neighbours_current)
  }
  return(jc_dr)
}

# find n neighbours from true distance matrix

find_neighbours <- function(dist_true, n) {
  dist_true <- dist_true %>% as.matrix()
  n = n+1 # if want 10 neighbours, need 11 closest as that includes itself
  neighbours_true <- list()
  for (i in 1:dim(dist_true)[1]) {
    neighbours_true[[i]] <- which.minn(dist_true[i,], n)[-1]
  }
  return(neighbours_true)
}


