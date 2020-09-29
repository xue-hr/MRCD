#' Generate matrix M_s
#'
#' Generate matrix M_s in formula 2.9 in
#' "The Asymptotic Variance Matrix of the Sample Correlation Matrix".
#'
#' @param n dimension.
#'
#' @return M_s M_s matrix.
#' @export
generate_M_s <- function(n)
{
  # For n, generate matrix M_s in formula 2.9
  K = matrix(0,n^2,n^2)
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      ind_row = (i-1)*n + j
      ind_col = (j-1)*n + i
      K[ind_row,ind_col] = K[ind_row,ind_col] + 1
    }
  }
  M_s = 1/2 * (diag(n^2) + K)
  return(M_s)
}
