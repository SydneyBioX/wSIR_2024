# This .R file contains functions for analysis of WSIR outputs.

# They take in a WSIR output to discover and explore spatial genes (genes with high loading).

library(tidyverse)
library(vctrs)
library(umap)

# visualise_wsir: show each cell coloured by its value for WSIR1 / WSIR2 / ... 

visualise_wsir <- function(coords, WSIR, dirs = ncol(WSIR$scores), mincol = "blue", maxcol = "red") {
  dirs <- min(dirs, ncol(WSIR$scores)) # make sure it is a valid value
  
  # initiliase empty long df
  vis_df_long <- matrix(NA, nrow = dirs*nrow(coords), ncol = 4) %>% as.data.frame() 
  colnames(vis_df_long) <- c("x", "y", "value", "WSIR_direction")
  
  # fill columns of long df with relevant WSIR1/2/... values
  vis_df_long$x <- rep(coords$x, dirs)
  vis_df_long$y <- rep(coords$y, dirs)
  vis_df_long$value <- as.vector(WSIR$scores[,1:dirs])
  vis_df_long$WSIR_direction <- as.factor(vec_rep_each(c(1:dirs), nrow(coords)))
  
  # produce plot
  plot <- ggplot(aes(x = x, y = y, color = value), data = vis_df_long) + 
    geom_point() +
    theme_classic() +
    facet_wrap(~WSIR_direction, scales = "fixed") + # one panel per WSIR direction
    ggtitle("Cells at true positions coloured by WSIR values") +
    scale_color_gradient(low = mincol, high = maxcol)
  return(plot)
}

# Show genes with high/low loading in WSIR1 in barplot

top_genes <- function(WSIR, highest = 10) {
  
  # make df with directions from WSIR object
  wsir_dirs_df <- WSIR$directions %>% as.data.frame()
  wsir_dirs_df$gene <- rownames(WSIR$directions)
  
  # retain in df only genes with high/low loading
  top_n <- doBy::which.maxn(abs(wsir_dirs_df$V1), highest)
  wsir_dirs_top_n <- wsir_dirs_df[top_n,]
  
  top_genes_ls <- wsir_dirs_df$V1[top_n]
  names(top_genes_ls) <- wsir_dirs_df$gene[top_n]
  top_genes_ls <- sort(top_genes_ls, decreasing = TRUE)
  
  # produce plot
  loadings_plot <- ggplot(aes(x = reorder(gene, -V1, sum), y = V1), data = wsir_dirs_top_n) +
    geom_col(width = 0.6) +
    theme_minimal() +
    labs(x = "Gene", y = "LowDim1") +
    geom_hline(yintercept = 0) +
    ggtitle(paste("Top", highest, "genes with highest/lowest loading in WSIR1"))
  
  return(list(plot = loadings_plot,
              genes = top_genes_ls))
}

# function to generate umap plots on the WSIRS (or SIRS or PCs etc)

vis_umap <- function(exprs, # same as the X input to WSIR function
                     WSIR, # output from WSIR function
                     highest_genes, # output from top_genes function
                     n_genes = length(names(highest_genes$genes))) { # integer for number of genes you want to show
  n_genes <- min(n_genes, length(names(highest_genes$genes))) # make sure it is a valid number of genes
  gene_names <- names(highest_genes$genes)[1:n_genes]
  
  umap_obj <- umap(WSIR$scores) # create umap object
  gene_inds <- match(gene_names, colnames(exprs)) # identify the relevant columns in the gene expression matrix
  
  # create and fill umap_df
  umap_df <- matrix(NA, nrow = n_genes*nrow(exprs), ncol = 4) %>% as.data.frame()
  colnames(umap_df) <- c("UMAP1", "UMAP2", "gene", "expression")
  
  umap_df$UMAP1 <- rep(umap_obj$layout[,1], n_genes)
  umap_df$UMAP2 <- rep(umap_obj$layout[,2], n_genes)
  umap_df$gene <- vec_rep_each(gene_names, nrow(exprs))
  umap_df$expression <- exprs[, gene_inds] %>% as.vector()
  
  # create plot
  plot = ggplot(data = umap_df, aes(x = UMAP1, y = UMAP2, colour = expression)) + 
    geom_point() + 
    facet_wrap(~gene)
  return(plot)
}

