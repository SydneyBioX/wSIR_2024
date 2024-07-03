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
