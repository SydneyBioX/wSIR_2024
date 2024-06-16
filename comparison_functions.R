# ABOUT

### This .R file contains the relevant files to facilitate comparison between various methods and datasets.
### It is made to save space in rmd files and prevent repeated comparison code chunks

# required packages

source("wsir_only.R") # this .R file contains the most updated WSIR functions (including QR scaling method)

library(tidyverse) # for pipes
library(MASS) # for lda function in lda_compare
library(mgcv) # for gam function in gam_compute
library(doBy) # for which.minn(x,n) function in knn_compute
library(pls) # for plsr function
library(Rfast) # for distCor metric

# functions for easy output of low-dim results from DR functions

lda_compare <- function(exprs_train, 
                        exprs_test, 
                        coords_train, 
                        slices = 3, 
                        varThreshold = 1) {
  tiles <- tile_allocator(coords = coords_train, slices = slices)
  sliceName <- "coordinate"
  labels <- tiles[,sliceName,drop=FALSE]
  
  exprs_train <- exprs_train %>% as.data.frame()
  exprs_train$tile <- labels$coordinate
  
  lda_mod <- lda(tile~., data = exprs_train)
  
  pred <- predict(lda_mod, as.data.frame(exprs_test))
  
  lda_proj <- pred$x
  
  return(list(predict(lda_mod)$x, lda_proj)) # 1st projected train, 2nd projected test
}

pls_compare <- function(exprs_train,
                        exprs_test,
                        coords_train,
                        coords_test) {
  
  pls_mod <- plsr(as.matrix(coords_train) ~ as.matrix(exprs_train))
  pred_test <- predict(pls_mod, newdata = as.matrix(exprs_test))
  pred_train <- predict(pls_mod, newdata = as.matrix(exprs_train))
  
  ncomps <- dim(pred_train)[3]
  pls_mses_train <- rep(0, ncomps)
  
  for (i in 1:ncomps) { # make a prediction for every number of components
    mse_train <- sum((rowSums((pred_train[,,i] - coords_train)^2))^0.5)
    pls_mses_train[i] <- mse_train
  }
  
  best_comp <- which.min(pls_mses_train)
  
  pred_coords <- pred_test[,,best_comp]
  
  rmse <- sum((rowSums((pred_coords - coords_test)^2))^0.5)
  
  return(list(pred_coords, rmse)) # predicted test coords
}

# exactly the same structure as WSIR() function
pca_obj <- function(X, varThreshold = 0.95) {
  
  # make pca object
  pca <- prcomp(X) 
  
  # find number of directions based on varThreshold
  cumvar_pca <- cumsum(pca$sdev^2)/sum(pca$sdev^2)
  d_pca <- which(cumvar_pca > varThreshold)[1]
  
  # compute scores, directions
  scores <- pca$x[,1:d_pca]
  directions <- pca$rotation[,1:d_pca]
  colnames(directions) <- NULL
  
  return(list(scores = scores,
              directions = directions,
              estd <- d_pca))
}

# extract_low_dim: find low-dim representations of train/test data according to dr methods

extract_low_dim <- function(exprs_train,
                          exprs_test,
                          coords_train,
                          directions = 40,
                          varThreshold = 0.95,
                          slices = 5,
                          SIR = FALSE,
                          WSIR = FALSE,
                          LDA = FALSE,
                          PCA = FALSE,
                          PLS = FALSE,
                          alpha = 4) {
  output <- list()
  
  if (SIR + WSIR + LDA + PCA + PLS == 0){
    return("error: need to set at least one method TRUE")
  }
  
  if (SIR) {
    sir <- WSIR(X = exprs_train,
                coords = coords_train,
                slices = slices,
                alpha = 0,
                varThreshold = varThreshold,
                maxDirections = directions)
    sir_project <- project_WSIR(wsir = sir, newdata = as.matrix(exprs_test))
    output$sir <- list(sir$scores, sir_project)
  }
  
  if (WSIR) {
    wsir <- WSIR(X = exprs_train,
                 coords = coords_train,
                 slices = slices,
                 varThreshold = varThreshold,
                 maxDirections = directions,
                 alpha = alpha)
    wsir_project <- project_WSIR(wsir = wsir, newdata = as.matrix(exprs_test))
    output$wsir <- list(wsir$scores, wsir_project)
  }
  
  if (LDA) {
    lda <- lda_compare(exprs_train = exprs_train,
                       exprs_test = exprs_test,
                       coords_train = coords_train,
                       slices = slices)
    output$lda <- lda
  }
  
  if (PCA) {
    pca <- pca_obj(X = exprs_train, varThreshold = varThreshold)
    pca_train <- pca$scores
    pca_test <- project_WSIR(wsir = pca, newdata = exprs_test)
    output$pca <- list(pca_train, pca_test)
  }
  
  if (PLS) { # this gives predicted coordinates of dimension nrow * 2
    pls_mod <- plsr(as.matrix(coords_train) ~ as.matrix(exprs_train))
    pred_train_all <- predict(pls_mod, newdata = as.matrix(exprs_train))
    pred_test_all <- predict(pls_mod, newdata = as.matrix(exprs_test))
    
    ncomps <- dim(pred_train_all)[3]
    pls_mses_train <- rep(0, ncomps)
    
    for (i in 1:ncomps) { # make a prediction for every number of components
      pred_current <- pred_train_all[,,i]
      mse_train <- sum((rowSums((pred_current - coords_train)^2))^0.5)
      pls_mses_train[i] <- mse_train
    }
    
    best_comp <- which.min(pls_mses_train)
    pred_test_best <- pred_test_all[,,best_comp]
    pred_train_best <- pred_train_all[,,best_comp]
    output$pls <- list(pred_train_best, pred_test_best)
  }
  
  return(output)
}

# code to compute performance metrics

gam_compute = function(low_dim_obj, 
                       coords_train, 
                       coords_test,
                       train_mse = FALSE) { # low_dim_obj needs to be just for one method, e.g low_dim_obj = low_dims$wsir
  low_train <- low_dim_obj[[1]]
  low_test <- low_dim_obj[[2]]
  df <- data.frame(low_train, coords_train)
  formula_x <- as.formula(paste("x~", paste(names(df)[1:(ncol(df)-2)], collapse="+")))
  formula_y <- as.formula(paste("y~", paste(names(df)[1:(ncol(df)-2)], collapse="+")))
  gam_mod <- mgcv::gam(formula = list(formula_x, 
                                       formula_y),
                        family = mvn(d = 2), 
                        data = df)
  low_test <- as.data.frame(low_test)
  colnames(low_test) <- colnames(df)[1:(ncol(df)-2)]
  mse = sum((rowSums((predict(gam_mod, newdata = low_test) - coords_test)^2))^0.5)
  
  if (train_mse) { # in case we want to see training mse. If not then return test mse only. Default is test only
    mse_train <- sum((rowSums((predict(gam_mod) - coords_train)^2))^0.5)
    return(mse_train)
  } else {
    return(mse)
  }
}

knn_positions <- function(low_train,
                          low_test,
                          coords_train,
                          k = 5) {
  n_train <- nrow(low_train)
  n_test <- nrow(low_test)
  
  predicted_positions <- matrix(NA, nrow = n_test, ncol = 2)
  for (i in 1:n_test) { # for each test observation, find its closest k train observations in low-dim space, average location between
    dists <- rep(0, n_train)
    for (j in 1:n_train) {
      dists[j] = sum((low_test[i,] - low_train[j,])^2)
    }
    close_trains <- which.minn(dists, n = k)
    close_positions <- coords_train[close_trains,]
    avg_position <- colMeans(close_positions)
    predicted_positions[i,] <- avg_position
  }
  return(predicted_positions)
}

knn_compute = function(low_dim_obj, # as in from extract_low_dim$wsir for example
                       coords_train,
                       coords_test,
                       k = 5) {
  
  predicted_positions <- knn_positions(low_train = low_dim_obj[[1]],
                                       low_test = low_dim_obj[[2]],
                                       coords_train = coords_train,
                                       k = k)
  
  mse = sum((rowSums((predicted_positions - coords_test)^2))^0.5)
  return(mse)
}

subset_lower.tri = function(m) {
  mm = m[lower.tri(m, diag = FALSE)]
  return(c(mm))
}

corDist = function(x,y) {cor(subset_lower.tri(dist(x)),
                             subset_lower.tri(dist(y)),
                             method = "spearman",
                             use = "pairwise.complete") 
} # this function is correlation of distances: cells close in Z are close in Y?
distCor = function(x,y) {Rfast::dcor(x = x,
                                     y = y)$dcor
} # this function is distance correlation:

# compute performance metric for each method selected

metric_gam <- function(SIR = FALSE,
                       WSIR = FALSE,
                       LDA = FALSE,
                       PCA = FALSE,
                       PLS = FALSE,
                       low_dims,
                       coords_train,
                       coords_test) {
  mses <- list()
  
  if (SIR) {
    MSE_sir <- gam_compute(low_dim_obj = low_dims$sir, coords_train = coords_train, coords_test = coords_test)
    mses$sir <- MSE_sir
  }
  
  if (WSIR) {
    MSE_wsir <- gam_compute(low_dim_obj = low_dims$wsir, coords_train = coords_train, coords_test = coords_test)
    mses$wsir <- MSE_wsir
  }
  
  if (LDA) {
    MSE_lda <- gam_compute(low_dim_obj = low_dims$lda, coords_train = coords_train, coords_test = coords_test)
    mses$lda <- MSE_lda
  }
  
  if (PCA) {
    MSE_pca <- gam_compute(low_dim_obj = low_dims$pca, coords_train = coords_train, coords_test = coords_test)
    mses$pca <- MSE_pca
  }
  if (PLS) {
    MSE_pls <- sum((rowSums((low_dims$pls[[2]] - coords_test)^2))^0.5)
    mses$pls <- MSE_pls
  }
  return(mses)
}

metric_knn <- function(SIR = FALSE,
                       WSIR = FALSE,
                       LDA = FALSE,
                       PCA = FALSE,
                       PLS = FALSE,
                       low_dims,
                       coords_train,
                       coords_test,
                       k = 1) {
  mses = list()
  
  if (SIR) {
    MSE_sir <- knn_compute(low_dim_obj = low_dims$sir, coords_train = coords_train, coords_test = coords_test, k = k)
    mses$sir <- MSE_sir
  }
  
  if (WSIR) {
    MSE_wsir <- knn_compute(low_dim_obj = low_dims$wsir, coords_train = coords_train, coords_test = coords_test, k = k)
    mses$wsir <- MSE_wsir
  }
  
  if (LDA) {
    MSE_lda <- knn_compute(low_dim_obj = low_dims$lda, coords_train = coords_train, coords_test = coords_test, k = k)
    mses$lda <- MSE_lda
  }
  
  if (PCA) {
    MSE_pca <- knn_compute(low_dim_obj = low_dims$pca, coords_train = coords_train, coords_test = coords_test, k = k)
    mses$pca <- MSE_pca
  }
  if (PLS) {
    MSE_pls <- sum((rowSums((low_dims$pls[[2]] - coords_test)^2))^0.5)
    mses$pls <- MSE_pls
  }
  return(mses)
}

metric_distCor <- function(SIR = FALSE,
                           WSIR = FALSE,
                           LDA = FALSE,
                           PCA = FALSE,
                           PLS = FALSE,
                           low_dims,
                           coords_test) {
  distCors <- list()
  if (SIR) {
    distCors$sir <- distCor(x = coords_test, y = low_dims$sir[[2]])
  }
  if (WSIR) {
    distCors$wsir <- distCor(x = coords_test, y = low_dims$wsir[[2]])
  }
  if (LDA) {
    distCors$lda <- distCor(x = coords_test, y = low_dims$lda[[2]])
  }
  if (PCA) {
    distCors$pca <- distCor(x = coords_test, y = low_dims$pca[[2]])
  }
  if (PLS) {
    distCors$pls <- distCor(x = coords_test, y = low_dims$pls[[2]])
  }
  return(distCors)
}

metric_corDist <- function(SIR = FALSE,
                           WSIR = FALSE,
                           LDA = FALSE,
                           PCA = TRUE,
                           PLS = FALSE,
                           low_dims,
                           coords_test) {
  corDists <- list()
  if (SIR) {
    corDists$sir <- corDist(x = coords_test, y = low_dims$sir[[2]])
  }
  if (WSIR) {
    corDists$wsir <- corDist(x = coords_test, y = low_dims$wsir[[2]])
  }
  if (LDA) {
    corDists$lda <- corDist(x = coords_test, y = low_dims$lda[[2]])
  }
  if (PCA) {
    corDists$pca <- corDist(x = coords_test, y = low_dims$pca[[2]])
  }
  if (PLS) {
    corDists$pls <- corDist(x = coords_test, y = low_dims$pls[[2]])
  }
  return(corDists)
}

# neighbourhood similarity function and helper functions

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

jaccard <- function(a,b) {
  int <- length(intersect(a,b))
  union <- length(a) + length(b) - int
  return(int/union)
}

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

find_neighbours <- function(dist_true, n) {
  dist_true <- dist_true %>% as.matrix()
  n = n+1 # if want 10 neighbours, need 11 closest as that includes itself
  neighbours_true <- list()
  for (i in 1:dim(dist_true)[1]) {
    neighbours_true[[i]] <- which.minn(dist_true[i,], n)[-1]
  }
  return(neighbours_true)
}

# put it all together in a function (including train/test split, various metrics)

analysis_all <- function(SIR = FALSE,
                         WSIR = FALSE,
                         LDA = FALSE,
                         PCA = FALSE,
                         PLS = FALSE,
                         exprs,
                         coords,
                         slices = 4,
                         directions = 50,
                         varThreshold = 0.95,
                         GAM = FALSE,
                         distCor = FALSE,
                         corDist = FALSE,
                         neighbourhood_similarity = FALSE,
                         KNN = FALSE,
                         k = 1,
                         nrep = 3,
                         plot_show = TRUE,
                         alpha = 4,
                         n_size = 50) {
  
  # initialise results vectors
  
  gam_mses_sir <- rep(0, nrep)
  gam_mses_wsir <- rep(0, nrep)
  gam_mses_lda <- rep(0, nrep)
  gam_mses_pca <- rep(0, nrep)
  gam_mses_pls <- rep(0, nrep)
  
  knn_mses_sir <- rep(0, nrep)
  knn_mses_wsir <- rep(0, nrep)
  knn_mses_lda <- rep(0, nrep)
  knn_mses_pca <- rep(0, nrep)
  knn_mses_pls <- rep(0, nrep)
  
  distCors_sir <- rep(0, nrep)
  distCors_wsir <- rep(0, nrep)
  distCors_lda <- rep(0, nrep)
  distCors_pca <- rep(0, nrep)
  distCors_pls <- rep(0, nrep)
  
  corDists_sir <- rep(0, nrep)
  corDists_wsir <- rep(0, nrep)
  corDists_lda <- rep(0, nrep)
  corDists_pca <- rep(0, nrep)
  corDists_pls <- rep(0, nrep)
  
  jaccard_sir <- rep(0, nrep)
  jaccard_wsir <- rep(0, nrep)
  jaccard_lda <- rep(0, nrep)
  jaccard_pca <- rep(0, nrep)
  jaccard_pls <- rep(0, nrep)
  
  # perform analysis nrep times
  for(i in 1:nrep) {
    print(i)
    
    # split data into train / test
    keep <- sample(c(TRUE, FALSE), nrow(exprs), replace = TRUE)
    exprs_train <- exprs[keep,]
    coords_train <- coords[keep,]
    
    exprs_test <- exprs[!keep,]
    coords_test <- coords[!keep,]
    
    # find low-dimensional representations
    low_dims <- extract_low_dim(exprs_train = exprs_train,
                                exprs_test = exprs_test,
                                coords_train = coords_train,
                                slices = slices,
                                SIR = SIR, 
                                WSIR = WSIR, 
                                LDA = LDA,
                                PCA = PCA,
                                PLS = PLS,
                                varThreshold = varThreshold,
                                directions = directions,
                                alpha = alpha)
    
    # calculate results from selected metric, save in relevant vectors
    ## GAM
    if (GAM) {
      gam_mses <- metric_gam(SIR = SIR,
                         WSIR = WSIR,
                         LDA = LDA,
                         PCA = PCA,
                         PLS = PLS,
                         coords_test = coords_test,
                         coords_train = coords_train,
                         low_dims = low_dims)
      if (SIR) {
        gam_mses_sir[i] <- gam_mses$sir
      }
      if (WSIR) {
        gam_mses_wsir[i] <- gam_mses$wsir
      }
      if (LDA) {
        gam_mses_lda[i] <- gam_mses$lda
      }
      if (PCA) {
        gam_mses_pca[i] <- gam_mses$pca
      }
      if (PLS) {
        gam_mses_pls[i] <- gam_mses$pls
      }
    }
    ## KNN
    if (KNN) {
      knn_mses <- metric_knn(SIR = SIR,
                             WSIR = WSIR,
                             LDA = LDA,
                             PCA = PCA,
                             PLS = PLS,
                             k = k,
                             coords_test = coords_test,
                             coords_train = coords_train,
                             low_dims = low_dims)
      if (SIR) {
        knn_mses_sir[i] <- knn_mses$sir
      }
      if (WSIR) {
        knn_mses_wsir[i] <- knn_mses$wsir
      }
      if (LDA) {
        knn_mses_lda[i] <- knn_mses$lda
      }
      if (PCA) {
        knn_mses_pca[i] <- knn_mses$pca
      }
      if (PLS) {
        knn_mses_pls[i] <- knn_mses$pls
      }
    }
    ## distCor
    if (distCor) {
      distCors <- metric_distCor(SIR = SIR,
                                 WSIR = WSIR,
                                 LDA = LDA,
                                 PCA = PCA,
                                 PLS = PLS,
                                 coords_test = coords_test,
                                 low_dims = low_dims)
      if (SIR) {
        distCors_sir[i] <- distCors$sir
      }
      if (WSIR) {
        distCors_wsir[i] <- distCors$wsir
      }
      if (LDA) {
        distCors_lda[i] <- distCors$lda
      }
      if (PCA) {
        distCors_pca[i] <- distCors$pca
      }
      if (PLS) {
        distCors_pls[i] <- distCors$pls
      }
    }
    ## corDist
    if (corDist) {
      corDists <- metric_corDist(SIR = SIR,
                                 WSIR = WSIR,
                                 LDA = LDA,
                                 PCA = PCA,
                                 PLS = PLS,
                                 coords_test = coords_test,
                                 low_dims = low_dims)
      if (SIR) {
        corDists_sir[i] <- corDists$sir
      }
      if (WSIR) {
        corDists_wsir[i] <- corDists$wsir
      }
      if (LDA) {
        corDists_lda[i] <- corDists$lda
      }
      if (PCA) {
        corDists_pca[i] <- corDists$pca
      }
      if (PLS) {
        corDists_pls[i] <- corDists$pls
      }
    }
    if (neighbourhood_similarity) {
      n_sims <- jaccard_compute(low_dim_obj = low_dims,
                                coords_test = coords_test,
                                n_size = n_size)
      if (SIR) {
        jaccard_sir[i] <- n_sims$sir
      }
      if (WSIR) {
        jaccard_wsir[i] <- n_sims$wsir
      }
      if (LDA) {
        jaccard_lda[i] <- n_sims$lda
      }
      if (PCA) {
        jaccard_pca[i] <- n_sims$pca
      }
      if (PLS) {
        jaccard_pls[i] <- n_sims$pls
      }
    }
  }
  
  # collect results into a dataframe
  results_df <- matrix(NA, nrow = nrep*5, ncol = 1) %>% as.data.frame() # nrep*5 since there are 5 methods
  colnames(results_df) <- c("method")
  results_df$method <- c(rep("SIR", nrep), rep("WSIR", nrep), rep("LDA", nrep), rep("PCA", nrep), rep("PLS", nrep))
  if (GAM) {
    results_df$mse_GAM <- c(gam_mses_sir, gam_mses_wsir, gam_mses_lda, gam_mses_pca, gam_mses_pls)
  }
  if (distCor) {
    results_df$distCor <- c(distCors_sir, distCors_wsir, distCors_lda, distCors_pca, distCors_pls)
  }
  if (corDist) {
    results_df$corDist <- c(corDists_sir, corDists_wsir, corDists_lda, corDists_pca, corDists_pls)
  }
  if (KNN) {
    results_df$mse_KNN <- c(knn_mses_sir, knn_mses_wsir, knn_mses_lda, knn_mses_pca, knn_mses_pls)
  }
  if (neighbourhood_similarity) {
    results_df$neighbourhood_similarity <- c(jaccard_sir, jaccard_wsir, jaccard_lda, jaccard_pca, jaccard_pls)
  }
  methods_used <- c("SIR"[SIR], "WSIR"[WSIR], "LDA"[LDA], "PCA"[PCA], "PLS"[PLS])
  results_df <- results_df[(results_df$method %in% methods_used),]
  
  # prepare outputs 
  output = list()
  output$dataframe <- results_df
  if (plot_show == TRUE) {
    if (GAM) {
      plot_gam = ggplot(aes(x = method, y = mse_GAM), data = results_df) +
        geom_point() +
        theme_classic() +
        ggtitle("MSE from GAM model from DR methods with nrep iterations")
      output$figure_gam = plot_gam
    }
    if (distCor) {
      plot_distCor = ggplot(aes(x = method, y = distCor), data = results_df) +
        geom_point() +
        theme_classic() +
        ggtitle("distCor between DR methods and test set coordinates with nrep iterations")
      output$figure_distCor = plot_distCor
    }
    if (corDist) {
      plot_corDist = ggplot(aes(x = method, y = corDist), data = results_df) +
        geom_point() +
        theme_classic() +
        ggtitle("corDist between DR methods and test set coordinates with nrep iterations")
      output$figure_corDist = plot_corDist
    }
    if (KNN) {
      plot_knn = ggplot(aes(x = method, y = mse_KNN), data = results_df) +
        geom_point() +
        theme_classic() +
        ggtitle("MSE from KNN location predictor from DR methods with nrep iterations")
      output$figure_knn = plot_knn
    }
    if (neighbourhood_similarity) {
      plot_nsim = ggplot(aes(x = method, y = neighbourhood_similarity), data = results_df) +
        geom_point() +
        theme_classic() +
        ggtitle("jaccard index for neighbourhood similarity between DR methods with nrep iterations")
      output$figure_nsim = plot_nsim
    }
  }
  
  return(output)
}
















