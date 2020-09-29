#' Estimate Variance
#'
#' Use bisection to search root, which is the MLE of sigma^2.
#'
#' @param v vector of length p.
#' @param Vg covariance matrix.
#'
#' @return Estiamte of sigma^2.
#' @export
estimate_sigma2 <- function(v,Vg)
{
  eigen_Vg = eigen(Vg)
  Q = t(eigen_Vg$vectors)
  eigen_values = eigen_Vg$values

  ###Now: Vg = t(Q)%*%diag(eigen_values)%*%Q

  x = Q%*%v
  left_end = max(min(c(x^2 - eigen_values)),max(-eigen_values)+1e-14)
  right_end = max(c(x^2 - eigen_values))
  f = function(r) {
    sum(-(x^2 / (eigen_values+r)^2)) + sum ( 1 / (eigen_values+r))
  }

  if(f(left_end)>0 | is.na(f(left_end)))
  {
    root = left_end
  } else {
    while(f(right_end)<0)
    {
      right_end = right_end+1
    }
    root_find = pracma::bisect(f, left_end, right_end,maxiter = 500)
    root = root_find$root
  }

  return(max(0,root))

}
