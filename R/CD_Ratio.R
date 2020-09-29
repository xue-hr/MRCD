#' Perform CD-Ratio Method with Multiple SNPs
#'
#' Perform CD-Ratio method for two traits T1 and T2,
#' assuming all SNPs are valid instruments.
#'
#' @param pruned a list, has "loci_bed" and "sig_part".
#' "loci_bed" is the reference panel, n by p matrix, genotype of
#' n individuals for p SNPs;
#' sig_part is a data.frame for the summary statistics, p by 12,
#' columns are "chr,pos,rsid,A1,A2,beta_T1,se_T1,N_T1,beta_T2,se_T2,N_T2,loci".
#' @param V_T1 asymptotic covariance matrix for correlations of p SNPs with T1,
#' default is NULL.
#' @param V_T2 asymptotic covariance matrix for correlations of p SNPs with T2,
#' default is NULL.
#'
#' @return a list with four elements "T1toT2", "T2toT1", "Q_T1toT2" and
#' "Q_T2toT1". "T1toT2" is a vector of 2 elements: estimated K and se for
#' direction T1 to T2.
#' "Q_T1toT2" is the test statistic for Goodness-of-Fit Tests for direction
#' T1 to T2. "T2toT1" and "Q_T2toT1" are similar.
#' @export
CD_Ratio <- function(pruned,V_T1 = NULL,V_T2 = NULL)
{
  SNP = pruned$loci_bed
  sig_part = pruned$sig_part
  p = ncol(SNP)
  #
  N_T1 = sig_part[,8]
  N_T2 = sig_part[,11]
  T1_T = sig_part[,6] / sig_part[,7]
  T1_r = T1_T / sqrt(N_T1 - 2 + T1_T^2)
  T2_T = sig_part[,9] / sig_part[,10]
  T2_r = T2_T / sqrt(N_T2 - 2 + T2_T^2)
  #
  rho_T1 = matrix(0, ncol = (p+1), nrow = (p+1))
  rho_T1[1:p,1:p] = cor(SNP)
  rho_T1[1:p,p+1] = T1_r
  rho_T1[p+1,1:p] = T1_r
  rho_T1[p+1,p+1] = 1
  #
  rho_T2 = matrix(0, ncol = (p+1), nrow = (p+1))
  rho_T2[1:p,1:p] = cor(SNP)
  rho_T2[1:p,p+1] = T2_r
  rho_T2[p+1,1:p] = T2_r
  rho_T2[p+1,p+1] = 1
  #
  if(is.null(V_T1))
  {
    V_T1 = calculate_asymptotic_variance(SNP,rho_T1)
  }
  if(is.null(V_T2))
  {
    V_T2 = calculate_asymptotic_variance(SNP,rho_T2)
  }
  #

  #T1 to T2
  jacobian = cbind( diag(1/T1_r) ,
                    -diag(T2_r/T1_r^2) )
  combined_V = rbind( cbind(V_T2,matrix(0,ncol = p, nrow = p)) / mean(N_T2),
                      cbind(matrix(0,ncol = p, nrow = p),V_T1) / mean(N_T1)
  )
  V = jacobian %*% combined_V %*% t(jacobian)
  inv_V = solve(V, tol = 0)
  est_vec = T2_r/T1_r
  gls_est = sum(inv_V %*% est_vec) / sum(inv_V)
  gls_var = 1 / sum(inv_V)

  T1toT2 = c(gls_est,sqrt(gls_var))
  Q_T1toT2 = (est_vec - gls_est)%*%inv_V%*%(est_vec - gls_est)

  #T2 to T1
  jacobian = cbind( diag(1/T2_r) ,
                    -diag(T1_r/T2_r^2) )
  combined_V = rbind( cbind(V_T1,matrix(0,ncol = p, nrow = p)) / mean(N_T1),
                      cbind(matrix(0,ncol = p, nrow = p),V_T2) / mean(N_T2)
  )
  V = jacobian %*% combined_V %*% t(jacobian)
  inv_V = solve(V, tol = 0)
  est_vec = T1_r/T2_r
  gls_est = sum(inv_V %*% est_vec) / sum(inv_V)
  gls_var = 1 / sum(inv_V)

  T2toT1 = c(gls_est,sqrt(gls_var))
  Q_T2toT1 = (est_vec - gls_est)%*%inv_V%*%(est_vec - gls_est)

  names(T1toT2) = c("K","se(K)")
  names(T2toT1) = c("K","se(K)")

  return(list(T1toT2 = T1toT2,
              T2toT1 = T2toT1,
              Q_T1toT2 = Q_T1toT2,
              Q_T2toT1 = Q_T2toT1))
}
