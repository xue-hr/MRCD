#' Perform CD-Egger Method with Multiple Independent SNPs
#'
#' Perform CD-Egger method for two traits T1 and T2.
#'
#' @param sig_part is a data.frame for the summary statistics, p by 12,
#' columns are: chr, pos, rsid, A1, A2, beta_T1, se_T1, N_T1, beta_T2, se_T2, N_T2, loci.
#' @param num_iteration number of iteration, default is 20.
#' @param V_T1 asymptotic covariance matrix for correlations of p SNPs with T1,
#' default is NULL.
#' @param V_T2 asymptotic covariance matrix for correlations of p SNPs with T2,
#' default is NULL.
#'
#' @return a list with four elements "T1toT2", "T2toT1", "Q_T1toT2" and
#' "Q_T2toT1". "T1toT2" is a vector of 4 elements: estimated K and se,
#' estimated intercept b0 and se, for
#' direction T1 to T2.
#' "Q_T1toT2" is the test statistic for Goodness-of-Fit Tests for direction
#' T1 to T2. "T2toT1" and "Q_T2toT1" are similar.
#' @export
CD_Egger_Independent <- function(sig_part,num_iteration = 20,V_T1 = NULL,V_T2 = NULL)
{
  p = nrow(sig_part)
  SNP = matrix(rnorm(p*p*100),p*100)
  SNP = scale(SNP)

  SIGMA_SNP = cov(SNP)
  SIGMA_INV = solve(SIGMA_SNP)
  v_coef = as.numeric(SIGMA_SNP%*%rep(1,p))
  SIGMA_SQUARE = SIGMA_SNP%*%SIGMA_SNP

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

  ### T1 to T2
  var_pleio = 0
  for(i in 1:num_iteration)
  {
    cov_mat = V_T2/mean(N_T2) + var_pleio*SIGMA_SQUARE
    M_tem = rbind(v_coef,T1_r)

    I_InfoMat = M_tem%*%solve(cov_mat,tol = 0)%*%t(M_tem)

    b0_K = solve(I_InfoMat,tol = 0)%*%M_tem%*%solve(cov_mat,tol = 0)%*%T2_r

    vv =
      as.numeric(SIGMA_INV%*%(T2_r - b0_K[1]*v_coef - b0_K[2]*T1_r))
    var_pleio = estimate_sigma2(v = vv,
                                Vg = SIGMA_INV%*%
                                  (V_T2/mean(N_T2))%*%SIGMA_INV)
  }
  vec_tem = T2_r - b0_K[1]*v_coef - b0_K[2]*T1_r
  Q_T1toT2 = vec_tem%*%solve(cov_mat,tol = 0)%*%vec_tem
  T1toT2 = c(b0_K,sqrt(diag(solve(I_InfoMat,tol = 0))))

  ### T2 to T1
  var_pleio = 0
  for(i in 1:num_iteration)
  {
    cov_mat = V_T1/mean(N_T1) + var_pleio*SIGMA_SQUARE
    M_tem = rbind(v_coef,T2_r)

    I_InfoMat = M_tem%*%solve(cov_mat,tol = 0)%*%t(M_tem)

    b0_K = solve(I_InfoMat,tol = 0)%*%M_tem%*%solve(cov_mat,tol = 0)%*%T1_r

    vv =
      as.numeric(SIGMA_INV%*%(T1_r - b0_K[1]*v_coef - b0_K[2]*T2_r))
    var_pleio = estimate_sigma2(v = vv,
                                Vg = SIGMA_INV%*%
                                  (V_T1/mean(N_T1))%*%SIGMA_INV)
  }
  vec_tem = T1_r - b0_K[1]*v_coef - b0_K[2]*T2_r
  Q_T2toT1 = vec_tem%*%solve(cov_mat,tol = 0)%*%vec_tem
  T2toT1 = c(b0_K,sqrt(diag(solve(I_InfoMat,tol = 0))))

  names(T1toT2) = c("b0","K","se(b0)","se(K)")
  names(T2toT1) = c("b0","K","se(b0)","se(K)")

  return(list(T1toT2 = T1toT2,
              T2toT1 = T2toT1,
              Q_T1toT2 = Q_T1toT2,
              Q_T2toT1 = Q_T2toT1))
}
