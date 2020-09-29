#' Generate matrix M_d
#'
#' Generate matrix M_d in formula 2.13 in
#' "The Asymptotic Variance Matrix of the Sample Correlation Matrix".
#'
#' @param n dimension.
#'
#' @return M_d M_d matrix.
#' @export
generate_M_d <- function(n)
{
  # For n, generate matrix M_d in formula 2.13
  K = matrix(0,n^2,n^2)
  for(i in 1:n)
  {
    ind_rc = (i-1)*n+i
    K[ind_rc,ind_rc] = K[ind_rc,ind_rc] + 1
  }

  return(K)
}
