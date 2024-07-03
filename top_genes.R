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
