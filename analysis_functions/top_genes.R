#' top_genes
#'
#' @description
#' A function to find and visualise the genes with the highest (in absolute value) loading in WSIR1. These genes contribute the most to the first low-dimensional direction. 
#'
#' @param WSIR wsir object as output of wSIR function. If you want to analyse a different DR method, ensure the slot named 'directions' contains the loadings as a matrix and gene names as rownames (e.g could be PC directions)
#' @param highest integer for how many of the top genes you would like to see. Recommend no more than 20 for ease of visualisation. Default is 10.
#'
#' @return plot 
#'
#' @examples
#' # need example?
#'
#' @export

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
