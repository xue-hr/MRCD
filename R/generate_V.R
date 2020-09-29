#' Generate matrix V
#'
#' Generate matrix V in formula 3.6 in
#' "The Asymptotic Variance Matrix of the Sample Correlation Matrix".
#'
#' @param SNP n by p matrix, genotype data for n individuals of p SNPs. And it is
#' standardized with column means 0 and variances 1.
#' @param rho (n+1) by (n+1) matrix, sample correlation matrix of n SNPs and X.
#'
#' @return V V matrix.
#' @export
generate_V <- function(SNP,rho)
{
  n = ncol(SNP)
  Sigma = rho[1:n,1:n]
  inv_Sigma = solve(Sigma,tol = 0)
  rho_X = as.matrix(rho[1:n,(n+1)])
  alpha =  inv_Sigma %*% rho_X
  e2 = 1 - t(rho_X) %*% inv_Sigma %*% rho_X

  if(nrow(alpha)>1)
  {
    SNP_alpha = SNP %*% diag(as.numeric(alpha))
  } else{
    SNP_alpha = SNP * as.numeric(alpha)
  }
  S = rowSums(SNP_alpha)
  ###

  sampleszie = length(S)
  sim_X = S + rnorm(sampleszie,0,sqrt(e2))
  sim_X = scale(sim_X)*sqrt(sampleszie) / sqrt(sampleszie-1)
  M_SX = cbind(SNP,sim_X)

  BigM = NULL
  for(i in 1:(n+1))
  {
    BigM = cbind(BigM,M_SX[,i]*M_SX)
  }
  ###
  B = t(BigM)
  C = eigenMapMatMult(B,BigM)
  V = C / sampleszie
  ###
  #V = t(BigM)%*%BigM / sampleszie
  ###

  vSigma = as.matrix(c(rho))
  V = V - vSigma%*%t(vSigma)

  return(V)
}
