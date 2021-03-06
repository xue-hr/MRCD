#' Perform 3 CD Methods with Multiple SNPs
#'
#' Perform CD-Ratio, CD-Egger and CD-GLS method for two traits T1 and T2,
#' only calculate asymptotic covariance matrices of sample correlations once to
#' save time.
#'
#' @param pruned a list, has "loci_bed" and "sig_part".
#' "loci_bed" is the reference panel, n by p matrix, genotype of
#' n individuals for p SNPs;
#' sig_part is a data.frame for the summary statistics, p by 12,
#' columns are: chr, pos, rsid, A1, A2, beta_T1, se_T1, N_T1, beta_T2, se_T2, N_T2, loci.
#' @param num_iteration number of iteration, default is 20.
#'
#' @return a list with 3 elements "CD_Ratio_result", "CD_Egger_result",
#'  "CD_GLS_result".
#' @export
CD_3_methods<-function(pruned,num_iteration = 20)
{
  ### perform all 3 CD methods: CD_Ratio, CD_Egger, CD_GLS
  ### "pruned" is a list, has "loci_bed" and "sig_part"
  ### loci_bed is the reference panel, n by p
  ### sig_part is the summary statistics, p by 12,
  ### columns are "chr,pos,rsid,A1,A2,beta_T1,se_T1,N_T1,beta_T2,se_T2,N_T2,loci"

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
  V_T1 = calculate_asymptotic_variance(SNP,rho_T1)
  V_T2 = calculate_asymptotic_variance(SNP,rho_T2)
  #
  CD_Ratio_result =
    CD_Ratio(pruned,V_T1 = V_T1,V_T2 = V_T2)

  CD_Egger_result =
    CD_Egger(pruned,num_iteration,V_T1 = V_T1,V_T2 = V_T2)

  CD_GLS_result = CD_GLS(pruned,num_iteration,V_T1 = V_T1,V_T2 = V_T2)

  return(list(CD_Ratio_result = CD_Ratio_result,
              CD_Egger_result = CD_Egger_result,
              CD_GLS_result = CD_GLS_result
  ))
}
