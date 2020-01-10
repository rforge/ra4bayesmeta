
# compute parameter m of the SGC(m) prior based on target RLMC value
m_inf_sgc <- function(rlmc, alpha=0.5){
  # inputs:
  # rlmc: target RLMC value
  # alpha: quantile of the SGC distribution to use as reference threshold U_ref
  # output:
  # parameter m of the SGC distribution with C=sigma_ref^{-2}
  
  log.quot <- log(rlmc) - log(1-rlmc)
  # return(1- log(alpha)/log(1+exp(log.quot)))
  log.res <- log(-log(alpha)) - log(log(1+exp(log.quot)))
  return(1 + exp(log.res))
}