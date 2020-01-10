
M_inf_sigc <- function(rlmc, df, alpha=0.5){
  # inputs:
  # rlmc: target RLMC value
  # df: data frame in bayesmeta format
  # alpha: quantile of the SGC distribution to use as reference threshold U_ref
  # output:
  # parameter m of the SGC distribution with C=sigma_ref^{-2}
  sigma.ref <- sigma_ref(df)
  log.quot <- -6*log(sigma.ref) + 2*log(1-rlmc) -2*log(rlmc) 
  # return(1- log(alpha)/log(1+exp(log.quot)))
  log.res <- log(-log(1-alpha)) - log(log(1+exp(log.quot)))
  return(1 + exp(log.res))
}