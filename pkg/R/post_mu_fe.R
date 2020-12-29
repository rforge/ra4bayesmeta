
####---- Marginal posterior for mu in the fixed-effects model (tau=0 benchmark) ----####

# implementation: analytical tau=0 benchmark with the prior and posterior predictive distributions
# model y_i\sim N(mu, sigma_i^2), sigma_i fixed, known, i=1,...,k, prior mu\sim N(nu, gamma^2)

post_mu_fe <- function(df, mu.mean = 0, mu.sd = 4){
  # input:
  # df: data frame in bayesmeta format containing y and sigma
  # mu.mean: mean of the normal prior for mu
  # mu.sd: standard deviation of the normal prior for mu
  
  y_i <- df$y
  sigma_i <- df$sigma
  
  # mean of the normal posterior of mu
  numerator <- sum(y_i/sigma_i^2)+mu.mean/mu.sd^2
  denominator <- sum(1/sigma_i^2)+1/mu.sd^2
  mean <- numerator/denominator
  
  # sd of the normal posterior of mu
  sd <- sqrt(1/(sum(1/sigma_i^2)+ 1/mu.sd^2))

  return(list(mean=mean, sd=sd))
}
