#' Calculate Asymptotic Variance Matrix
#'
#' Calculate the asymptotic covariance matrix of sample correlations of n
#' SNPs with X, Use Theorem 2 in
#' "The Asymptotic Variance Matrix of the Sample Correlation Matrix".
#'
#' @param SNP n by p matrix, genotype data for n individuals of p SNPs.
#' @param rho (n+1) by (n+1) matrix, sample correlation matrix of n SNPs and X.
#'
#' @return n by n matrix, which is the asymptotic covariance matrix of
#' sample correlations of n SNPs with X.
#' @export
calculate_asymptotic_variance <- function(SNP,rho)
{
  n = ncol(rho)-1
  SNP = scale(SNP)
  SNP = SNP * sqrt(nrow(SNP)) / sqrt(nrow(SNP)-1)
  M_s = generate_M_s(n + 1)
  M_d = generate_M_d(n + 1)
  ###
  B = ( kronecker( diag(n+1) , rho) )
  C = eigenMapMatMult(M_s, B)
  C = eigenMapMatMult(C,M_d)
  M_1 = diag((n+1)^2) - C

  V = generate_V(SNP,rho)

  C = eigenMapMatMult(M_1,V)
  B = t(M_1)
  asymp_cov = eigenMapMatMult(C,B)

  target_ind1 = n * (n+1) + 1
  target_ind2 = (n+1)^2-1

  return( asymp_cov[(target_ind1:target_ind2) , (target_ind1:target_ind2)] )
}
