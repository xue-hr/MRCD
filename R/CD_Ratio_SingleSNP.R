#' Perform CD-Ratio Method For with A Single SNP
#'
#' Perform CD-Ratio method for two traits T1 and T2,
#' with one SNP.
#'
#' @param pruned a list, has "loci_bed" and "sig_part".
#' "loci_bed" is the reference panel, n by 1 matrix, genotype of
#' n individuals for 1 SNP;
#' sig_part is a data.frame for the summary statistics, 1 by 12,
#' columns are: chr, pos, rsid, A1, A2, beta_T1, se_T1, N_T1, beta_T2, se_T2, N_T2, loci.
#'
#' @return a list with two elements "T1toT2", "T2toT1".
#' "T1toT2" is a vector of 2 elements: estimated K and se for
#' direction T1 to T2.
#' "T2toT1" is similar.
#' @export
CD_Ratio_SingleSNP <- function(pruned)
{
  SNP = pruned$loci_bed
  SNP = scale(SNP)
  sig_part = pruned$sig_part
  p = ncol(SNP)
  #
  N_T1 = sig_part[,8]
  N_T2 = sig_part[,11]
  T1_T = sig_part[,6] / sig_part[,7]
  T1_r = T1_T / sqrt(N_T1 - 2 + T1_T^2)
  T2_T = sig_part[,9] / sig_part[,10]
  T2_r = T2_T / sqrt(N_T2 - 2 + T2_T^2)
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
  V_T1 = calculate_asymptotic_variance(SNP,rho_T1)
  V_T2 = calculate_asymptotic_variance(SNP,rho_T2)
  #

  ###T1 to T2
  est_vec = T2_r/T1_r
  var_vec = est_vec^2 * ( V_T2/N_T2/T2_r^2 +
                            V_T1/N_T1/T1_r^2 )
  T1toT2 = c(est_vec,sqrt(var_vec))

  ###T2 to T1
  est_vec = T1_r/T2_r
  var_vec = est_vec^2 * ( V_T1/N_T1/T1_r^2 +
                            V_T2/N_T2/T2_r^2 )
  T2toT1 = c(est_vec,sqrt(var_vec))

  names(T1toT2) = c("K","se(K)")
  names(T2toT1) = c("K","se(K)")

  return(list(T1toT2 = T1toT2,
              T2toT1 = T2toT1
  ))

}
