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
