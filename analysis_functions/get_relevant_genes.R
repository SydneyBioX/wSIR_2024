# get_relevant_genes

#' Identify the relevant genes from wSIR
#'
#' @description
#' This function calculates standardized z-scores for gene expression values, converts them to corresponding p-values, identifies genes with p-values below a specified threshold, and returns their names as relevant genes.
#'
#' @param out output from running wSIR function
#' @param p.threshold threshold value for unadjusted p-value
#'
#' @return character vector of relevant features
#'
#' @examples
#' #hello()
#'
#' @export
get_relevant_genes = function(out,
                              p.threshold = 0.01) {

  # out is the output from wSIR function
  # p.threshold is the (unadjusted) P-value threshold

  out_z = apply(out[[2]], 2, function(x) (x - mean(x)) / sd(x))
  p_z = pnorm(abs(out_z), lower.tail = FALSE)
  svgs = names(which(apply(p_z, 1, function(x) any(x < p.threshold))))

  return(svgs)
}
