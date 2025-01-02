#' Calculate transcripts per million for each isoform
#'
#' @param isoforms description
#' @param min_tpm description
#'
#' \deqn{
#'  TPM = \frac{q_i/l_i}{\sum q_j/l_j} 10^6
#' }
calc_tpm <- function(isoforms) {
  expr <- dplyr::distinct(isoforms, read_seq, .keep_all=TRUE)
  total <- sum(expr$read_num / nchar(expr$read_seq))
  tpm <- isoforms$read_num / nchar(isoforms$read_seq) / total * 1e6
  return(tpm)
}


#' Calculate tissue specificity index of an expression profile.
#'
#' \code{calculate_tsi} returns the tau value of a gene expression profile
#'   across tissues. The index \eqn{\tau} is defined as:
#'   \deqn{
#'    \tau=\frac{\sum_{i}^N (1-x_i)}/N-1
#'   }
#'
#' @param exp Expression profile of a gene
#' @return tau
calc_one_tsi <- function(exp){
  max_exp <- max(exp)
  tsi <- sum(1 - exp/max_exp) / (length(exp) - 1)
  return(tsi)
}
