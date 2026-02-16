
library(wSIR)
library(pls) # for plsr function
#library(Rfast) # for distCor metric
library(tidyverse) # for pipes
library(MASS) # for lda function in lda_compare
library(vctrs) # for vec_rep_each function

tile_allocator <- function(coords, slices = 3) {
  sliced <- lapply(coords, function(x) as.integer(cut(rank(x), slices)))
  allocation <- as.factor(do.call(paste, c(sliced, sep = ", ")))
  
  newcoords <- cbind(coords, coordinate = allocation)
  
  return(newcoords)
}

lda_compare <- function(exprs_train, 
                        exprs_test, 
                        coords_train, 
                        slices = 3, 
                        varThreshold = 1) {
  tiles <- tile_allocator(coords = coords_train, slices = slices) # replace this
  sliceName <- "coordinate"
  labels <- tiles[,sliceName,drop=FALSE]
  
  exprs_train <- exprs_train %>% as.matrix() %>% as.data.frame()
  exprs_train$tile <- labels$coordinate
  
  lda_mod <- lda(tile~., data = exprs_train)
  
  pred <- predict(lda_mod, as.data.frame(as.matrix(exprs_test)))
  
  lda_proj <- pred$x
  
  return(list(predict(lda_mod)$x, lda_proj)) # 1st projected train, 2nd projected test
}

pls_compare <- function(exprs_train,
                        exprs_test,
                        coords_train,
                        coords_test) {
  
  pls_mod <- plsr(as.matrix(coords_train) ~ as.matrix(exprs_train))
  pred_test_scores_all <- predict(pls_mod, newdata = as.matrix(exprs_test), type = "scores")
  pred_train_coords <- predict(pls_mod, newdata = as.matrix(exprs_train))
  
  ncomps <- dim(pred_train_coords)[3]
  pls_mses_train <- rep(0, ncomps)
  
  for (i in 1:ncomps) { # make a prediction for every number of components
    mse_train <- sum((rowSums((pred_train_coords[,,i] - coords_train)^2))^0.5)
    pls_mses_train[i] <- mse_train
  }
  
  best_comp <- which.min(pls_mses_train)
  
  pred_test_scores_selected <- pred_test_scores_all[,1:best_comp]
  
  return(pred_test_scores_selected) # predicted test coords
}

pca_obj <- function(X, varThreshold = 0.95, maxDirections = 40) {
  
  # make pca object
  pca <- prcomp(X) 
  
  # find number of directions based on varThreshold
  cumvar_pca <- cumsum(pca$sdev^2)/sum(pca$sdev^2)
  d_pca <- which(cumvar_pca >= varThreshold)[1]
  
  d_pca <- min(maxDirections, d_pca)
  
  # compute scores, directions
  scores <- pca$x[,1:d_pca]
  directions <- pca$rotation[,1:d_pca]
  colnames(directions) <- NULL
  
  return(list(scores = scores,
              directions = directions,
              estd <- d_pca))
}

graphPCA_obj <- function(X, coords, k = 6, lambda = 0.5, n_components = 50){
  require(FNN)
  require(Matrix)
  
  knn_result <- get.knn(coords, k = k)
  # Create the adjacency matrix (sparse format)
  adj_matrix <- sparseMatrix(
    i = rep(1:nrow(coords), each = k),
    j = as.vector(knn_result$nn.index),
    x = 1
  )
  # Create the degree matrix (diagonal matrix)
  degree_matrix <- Diagonal(x = rowSums(adj_matrix))
  # Compute the Laplacian matrix
  laplacian_matrix <- degree_matrix - adj_matrix
  # Convert to full matrix if needed
  n <- nrow(X)
  graphL <- as.matrix(laplacian_matrix)
  G <- Diagonal(n) + lambda*graphL
  Ginv <- Rfast::spdinv(as.matrix(G))
  Expr <- scale(X)
  C <- t(Expr) %*% Ginv %*% Expr
  eigenC <- eigen(C, symmetric = T)
  out_GPCA <- list()
  out_GPCA$loadings <- eigenC$vectors[, 1:n_components]
  out_GPCA$scores <- Ginv %*% Expr %*% out_GPCA$loadings
  return(out_GPCA)
}

create_low_dim <- function(count_train, 
                           exprs_train,
                           exprs_test,
                           coords_train,
                           count_test,
                           coords_test,
                           maxDirections = 40,
                           varThreshold = 0.95,
                           slices = 5,
                           optim_params = TRUE,
                           verbose = TRUE,
                           alpha_vals = c(0,2,4,8),
                           slice_vals = c(5,7,10,15),
                           SIR = FALSE,
                           WSIR = FALSE,
                           LDA = FALSE,
                           PCA = FALSE,
                           PLS = FALSE,
                           GPCA = FALSE,
                           alpha = 4) {
  output <- list()
  
  if (WSIR) {
    if (verbose) {
      print("Starting wSIR")
    }
    if (optim_params) {
      wsir_optim_obj <- exploreWSIRParams(X = exprs_train,
                                          coords = coords_train,
                                          alpha_vals = alpha_vals,
                                          slice_vals = slice_vals,
                                          varThreshold = varThreshold,
                                          maxDirections = maxDirections)
      alpha <- wsir_optim_obj$best_alpha
      slices <- wsir_optim_obj$best_slices
    }
    
    wsir <- wSIR:::wSIRSpecifiedParams(X = exprs_train,
                       coords = coords_train,
                       alpha = alpha,
                       slices = slices,
                       varThreshold = varThreshold,
                       maxDirections = maxDirections)
    wsir_project <- projectWSIR(wsir = wsir, newdata = as.matrix(exprs_test))
    output$wsir <- list(wsir$scores, wsir_project)
  }
  
  if (SIR) {
    if (verbose) {
      print("Starting SIR")
    }
    sir <- wSIR:::wSIRSpecifiedParams(X = exprs_train,
                coords = coords_train,
                slices = slices,
                alpha = 0,
                varThreshold = varThreshold,
                maxDirections = maxDirections)
    sir_project <- projectWSIR(wsir = sir, newdata = as.matrix(exprs_test))
    output$sir <- list(sir$scores, sir_project)
  }
  
  if (LDA) {
    if (verbose) {
      print("Starting LDA")
    }
    lda <- lda_compare(exprs_train = exprs_train,
                       exprs_test = exprs_test,
                       coords_train = coords_train,
                       slices = slices)
    output$lda <- lda
  }
  
  if (GPCA) {
    if (verbose) {
      print("Starting GPCA")
    }
    out <- graphPCA_run(count_train, as.matrix(coords_train))  
    gpca_loadings <- out$uns$GraphPCA_loadings
    # for the normalization version on the test set
    out_test <- graphPCA_run(count_test, as.matrix(coords_test))
    gpca_test <- out_test$X %*% gpca_loadings
    output$gpca <- list(out$uns$GraphPCA, gpca_test)
  }
    
  if (PCA) {
    if (verbose) {
        print("Starting PCA")
      }
      pca <- pca_obj(X = exprs_train, varThreshold = varThreshold, maxDirections = maxDirections)
      pca_train <- pca$scores
      pca_test <- projectWSIR(wsir = pca, newdata = as.matrix(exprs_test))
      output$pca <- list(pca_train, pca_test)
  }
  
  if (PLS) { # this gives predicted coordinates of dimension nrow * 2
    if (verbose) {
      print("Starting PLS")
    }
    pls_mod <- plsr(as.matrix(coords_train) ~ as.matrix(exprs_train))
    pred_train_coords_all <- predict(pls_mod, newdata = as.matrix(exprs_train))
    pred_test_scores_all <- predict(pls_mod, newdata = as.matrix(exprs_test), type = "scores")
    pred_train_scores_all <- predict(pls_mod, newdata = as.matrix(exprs_train), type = "scores")
    
    pls_mses_train <- rep(0, maxDirections)
    
    for (i in 1:maxDirections) { # make a prediction for every number of components
      pred_current <- pred_train_coords_all[,,i]
      mse_train <- sum((rowSums((pred_current - coords_train)^2))^0.5)
      pls_mses_train[i] <- mse_train
    }
    
    best_comp <- which.min(pls_mses_train)
    pred_test_best <- pred_test_scores_all[,1:best_comp]
    pred_train_best <- pred_train_scores_all[,1:best_comp]
    output$pls <- list(pred_train_best, pred_test_best)
  }
  
  return(list(output, alpha, slices))
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
    distCors$sir <- distCor(x = coords_test, y = as.matrix(low_dims$sir[[2]]))
  }
  if (WSIR) {
    distCors$wsir <- distCor(x = coords_test, y = as.matrix(low_dims$wsir[[2]]))
  }
  if (LDA) {
    distCors$lda <- distCor(x = coords_test, y = as.matrix(low_dims$lda[[2]]))
  }
  if (PCA) {
    distCors$pca <- distCor(x = coords_test, y = as.matrix(low_dims$pca[[2]]))
  }
  if (PLS) {
    distCors$pls <- distCor(x = coords_test, y = as.matrix(low_dims$pls[[2]]))
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

corDist = function(x,y) {cor(subsetLowerTri(dist(x)),
                             subsetLowerTri(dist(y)),
                             method = "spearman",
                             use = "pairwise.complete") 
} # this function is correlation of distances: cells close in Z are close in Y?
distCor = function(x,y) {Rfast::dcor(x = x,
                                     y = y)$dcor
} # this function is distance correlation:

subsetLowerTri = function(m) {
  mm = m[lower.tri(m, diag = FALSE)]
  return(c(mm))
}

metric_bc_distCor <- function(low_dims, coords_test) {
  lapply(low_dims, function(onemethod){
    bc_distCor(x = coords_test, y = onemethod[[2]])
  })
}  
  

bc_distCor <- function(x, y) {
  Rfast::dcor.ttest(x, y)[["BC dcor"]]
}
