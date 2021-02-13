# Hellinger distance between the marginal posteriors of 2 bayesmeta fits
# input: fit0, fit1 are bayesmeta fits
# parameter: determines which marginal posteriors to use;
#            possible values are "mu", "tau", "theta_new" 
#            and "theta" (for random effects, specify the number of the study in the individual argument)
# individual: used only if parameter = "theta", specifies the number of the study in the data
# method: method for computation, either "integral" for numerical computation of integrals
#         or "moment" for a moment-based approxmate analytical formula
#         the moment-based method relies on the normal approximation of the densities
H_fits <- function(fit1, fit2, parameter="mu", individual = NA,
                   method = "integral"){
  
  # number of studies
  k <- fit1$k
  if(parameter == "theta" & is.na(individual))
    stop(paste0("individual should be an integer between 1 and ", k, ".",
                sep = ""))
  # if individual is larger than k, give an error message
  if(parameter == "theta" & individual > k)
    stop(paste0("individual is too large. There are only ", k, " studies in the data set.",
                sep = ""))
  
  if(method == "integral"){
    # lower integration limit is -Inf except for tau (will be reset in that case)
    lower <- -Inf
    if(parameter == "mu"){
      dens1 <- function(x) fit1$dposterior(mu = x)
      dens2 <- function(x) fit2$dposterior(mu = x)
    }
    if(parameter == "tau"){
      # lower integration limit
      lower <- 0
      dens1 <- function(x) fit1$dposterior(tau = x)
      dens2 <- function(x) fit2$dposterior(tau = x)
    }
    if(parameter == "theta_new"){
      dens1 <- function(x) fit1$dposterior(mu = x, predict = TRUE)
      dens2 <- function(x) fit2$dposterior(mu = x, predict = TRUE)
    }
    if(parameter == "theta"){
      dens1 <- function(x) fit1$dposterior(theta = x, individual = individual)
      dens2 <- function(x) fit2$dposterior(theta = x, individual = individual)
    }
    h <- H(dens1, dens2, lower = lower)
  }
  if(method == "moment"){
    if(parameter == "mu"){
      mean1 <- fit1$summary["mean", parameter]
      sd1 <- fit1$summary["sd", parameter]
      mean2 <- fit2$summary["mean", parameter]
      sd2 <- fit2$summary["sd", parameter]
    }
    if(parameter == "tau"){
      # computation of mean and sd for log_tau
      log_tau_estimates <- function(fit, ilim = 700){
        mean_post_log_tau <- integrate(function(x) {x * fit$dposterior(tau = exp(x)) * exp(x)}, lower = -ilim, upper = ilim)$value
        sd_post_log_tau <- sqrt(integrate(function(x) {x^2 * fit$dposterior(tau = exp(x)) * exp(x)}, lower = -ilim, upper = ilim)$value - mean_post_log_tau^2)
        return(c(mean_post_log_tau, sd_post_log_tau))
      }
      
      estimates1 <- log_tau_estimates(fit1)
      estimates2 <- log_tau_estimates(fit2)
      mean1 <- estimates1[1]
      sd1 <- estimates1[2]
      mean2 <- estimates2[1]
      sd2 <- estimates2[2]
    }
    if(parameter == "theta_new"){
      mean1 <- fit1$summary["mean", "theta"]
      sd1 <- fit1$summary["sd", "theta"]
      mean2 <- fit2$summary["mean", "theta"]
      sd2 <- fit2$summary["sd", "theta"]
    }
    if(parameter == "theta"){
        mean1 <- fit1$theta["mean", individual]
        sd1 <- fit1$theta["sd", individual]
        mean2 <- fit2$theta["mean", individual]
        sd2 <- fit2$theta["sd", individual]
    }
    h <- H_normal(mean1 = mean1, sd1 = sd1, mean2 = mean2, sd2 = sd2)
  }
  return(h)
}