# compute fits for the HN0 and J benchmarks and the actual priors specified
fit_models_RA <- function (df, tau.prior = list(), scale.hn0 = 1/500,
                           mu.mean = 0,  mu.sd = 4, interval.type = "central") {
  
  # compute the posterior for the Jefferys prior
  fit.j <- bayesmeta(y = df[, "y"], sigma = df[, "sigma"], 
                     mu.prior.mean = mu.mean, mu.prior.sd = mu.sd, tau.prior = "Jeffreys",
                     interval.type= interval.type)
  
  # compute the posterior for the HN0 = HN(0, 0.01) prior
  # tau_sd_AA0 <- 1/500
  # tau_sd_AA0 <- 1/100
  fit.hn0 <- bayesmeta(y = df[, "y"], sigma = df[, "sigma"], 
                       mu.prior.mean = mu.mean, mu.prior.sd = mu.sd, 
                       tau.prior = function (x) dhalfnormal(x, scale = scale.hn0),
                       interval.type= interval.type)
  
  if (length(tau.prior) > 0) {
    n.act <- length(tau.prior)
    fit.actual <- list()
    fit.bms <- list(fit.hn0, fit.j)
    for (i in 1:n.act) {
      fit.actual[[i]] <- bayesmeta(y = df[, "y"], 
                                   sigma = df[, "sigma"], mu.prior.mean = mu.mean, 
                                   mu.prior.sd = mu.sd, tau.prior = tau.prior[[i]],
                                   interval.type= interval.type)
      # interval.type="central")
    }
    fits <- c(fit.bms, fit.actual)
    names.fit.actual <- rep(NA, times = n.act)
    for (j in 1:n.act) {
      names.fit.actual[j] <- paste0("fit.actual_", 
                                    j, sep = "")
    }
    names(fits) <- c("fit.hn0", "fit.j", names.fit.actual)
    return(fits)
  }
  else fits <- list(fit.hn0 = fit.hn0, fit.j = fit.j)
  return(fits)
}