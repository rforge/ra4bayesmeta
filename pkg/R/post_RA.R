# compute post-RA table using the benchmarks HN0 and J for the parameters mu, tau, theta_new
# (the random effects theta_i are included if show.re=TRUE)

post_RA <- function(df, tau.prior = list(function(x) dhalfnormal(x, scale = 1)), 
                    show.re = FALSE, estimate = "median", ci.method = "central", 
                    scale.hn0 = 1/500, mu.mean = 0, mu.sd = 4) {
  
  fits <- fit_models_RA(df = df, tau.prior = tau.prior, 
                        scale.hn0 = scale.hn0,
                        mu.mean = mu.mean, mu.sd = mu.sd)
  
  fits.bm <- fits[c("fit.hn0", "fit.j")]
  n.act <- length(tau.prior) + 2
  fits.actual <- fits[c(1, 3:n.act, 2)]
  k <- fits.actual[[1]]$k
  n.bm <- length(fits.bm)
  if(show.re == FALSE)
    table <- matrix(nrow = 3 * n.act + 1, ncol = 4 + n.bm)
  else
    table <- matrix(nrow = (3 + k) * n.act + 1, ncol = 4 + n.bm)
  
  # H distances betw. po_HN0 and po_J benchmarks to compute signed informativeness
  H_dist_post_bms_mu <- H(function(x) fits.bm[["fit.hn0"]]$dposterior(mu = x), 
                          function(x) fits.bm[["fit.j"]]$dposterior(mu = x))
  H_dist_post_bms_tau <- H(function(x) fits.bm[["fit.hn0"]]$dposterior(tau = x), 
                           function(x) fits.bm[["fit.j"]]$dposterior(tau = x), lower=0)
  H_dist_post_bms_theta_new <- H(function(x) fits.bm[["fit.hn0"]]$dposterior(mu = x, predict = TRUE), 
                                 function(x) fits.bm[["fit.j"]]$dposterior(mu = x, predict = TRUE))
  
  # compute mean and sd for the FE benchmark
  mean.fe <- post_mu_fe(df = df,  mu.mean = mu.mean, mu.sd = mu.sd)$"mean"
  sd.fe <- post_mu_fe(df = df,  mu.mean = mu.mean,  mu.sd = mu.sd)$"sd"
  
  # add FE benchmark in first row
  # median = mean for normal distribution
  table[1,1] <- mean.fe
  table[1,2] <- qnorm(p=0.025, mean = mean.fe, sd = sd.fe)
  table[1,3] <- qnorm(p=0.975, mean = mean.fe, sd = sd.fe)
  
  j <- 1
  table[1, 4+j] <- H(function(x) dnorm(x, mean = mean.fe, sd = sd.fe), 
                     function(x) fits.bm[[j]]$dposterior(mu = x))
  j <- 2
  table[1, 4+j] <- sign(table[1, 5] - H_dist_post_bms_mu) * H(function(x) dnorm(x, mean = mean.fe, sd = sd.fe), 
                                                              function(x) fits.bm[[j]]$dposterior(mu = x))
  
  for (i in 1:n.act) {
    m <- i + 1
    table[m,1] <- fits.actual[[i]]$summary[estimate, "mu"]
    table[m,2] <- fits.actual[[i]]$post.interval(mu.level=0.95, method=ci.method)[1]
    table[m,3] <- fits.actual[[i]]$post.interval(mu.level=0.95, method=ci.method)[2]
    
    table[m, 5] <- H(function(x) fits.actual[[i]]$dposterior(mu = x), 
                     function(x) fits.bm[["fit.hn0"]]$dposterior(mu = x))
    table[m, 6] <- sign(table[m, 5] - H_dist_post_bms_mu) * H(function(x) fits.actual[[i]]$dposterior(mu = x), 
                                                              function(x) fits.bm[["fit.j"]]$dposterior(mu = x))
    
    m <- n.act + i + 1  
    table[m,1] <- fits.actual[[i]]$summary[estimate, "tau"]
    table[m,2] <- fits.actual[[i]]$post.interval(tau.level=0.95, method=ci.method)[1]
    table[m,3] <- fits.actual[[i]]$post.interval(tau.level=0.95, method=ci.method)[2]
    
    table[m, 5] <- H(function(x) fits.actual[[i]]$dposterior(tau = x), 
                     function(x) fits.bm[["fit.hn0"]]$dposterior(tau = x), 
                     lower = 0)
    table[m, 6] <- sign(table[m, 5] - H_dist_post_bms_tau) * H(function(x) fits.actual[[i]]$dposterior(tau = x), 
                                                               function(x) fits.bm[["fit.j"]]$dposterior(tau = x), 
                                                               lower = 0)
    if(show.re == FALSE){
      m <- 2*n.act + i + 1
      table[m,1] <- fits.actual[[i]]$summary[estimate, "theta"]
      table[m,2] <- fits.actual[[i]]$post.interval(theta.level=0.95, method=ci.method, predict=TRUE)[1]
      table[m,3] <- fits.actual[[i]]$post.interval(theta.level=0.95, method=ci.method, predict=TRUE)[2]
      
      table[m, 5] <- H(function(x) fits.actual[[i]]$dposterior(mu = x, predict = TRUE), 
                       function(x) fits.bm[["fit.hn0"]]$dposterior(mu = x))
      table[m, 6] <- sign(table[m, 5] - H_dist_post_bms_theta_new) * H(function(x) fits.actual[[i]]$dposterior(mu = x, predict = TRUE), 
                                                                function(x) fits.bm[["fit.j"]]$dposterior(mu = x, predict = TRUE))
    } else {
    for(j in 1:k){
      m <- 2*n.act + 1 + (j-1)*n.act + i
      table[m,1] <- fits[[i]]$qposterior(theta.p=0.5, individual=j)
      table[m,2] <- fits[[i]]$post.interval(theta.level=0.95, individual=j, method= ci.method)[1]
      table[m,3] <- fits[[i]]$post.interval(theta.level=0.95, individual=j, method= ci.method)[2]
      
      H_dist_post_bms_theta_j <- H(function(x) fits.bm[["fit.hn0"]]$dposterior(theta = x, 
                                                                               individual = j), 
                                   function(x) fits.bm[["fit.j"]]$dposterior(theta = x, individual = j))
      
      table[m, 5] <- H(function(x) fits.actual[[i]]$dposterior(theta = x, 
                                                               individual = j), 
                       function(x) fits.bm[["fit.hn0"]]$dposterior(theta = x, individual = j))
      table[m, 6] <- sign(table[m, 5] - H_dist_post_bms_theta_j) * H(function(x) fits.actual[[i]]$dposterior(theta = x, 
                                                                                                             individual = j), 
                                                                     function(x) fits.bm[["fit.j"]]$dposterior(theta = x, individual = j))
      
    }
    
    m <- (k+2)*n.act + 1 + i
    table[m,1] <- fits.actual[[i]]$summary[estimate, "theta"]
    table[m,2] <- fits.actual[[i]]$post.interval(theta.level=0.95, method=ci.method, predict=TRUE)[1]
    table[m,3] <- fits.actual[[i]]$post.interval(theta.level=0.95, method=ci.method, predict=TRUE)[2]
    
    table[m, 5] <- H(function(x) fits.actual[[i]]$dposterior(mu = x, predict = TRUE), 
                     function(x) fits.bm[["fit.hn0"]]$dposterior(mu = x, predict = TRUE))
    table[m, 6] <- sign(table[m, 5] - H_dist_post_bms_theta_new) * H(function(x) fits.actual[[i]]$dposterior(mu = x, predict = TRUE), 
                                                                   function(x) fits.bm[["fit.j"]]$dposterior(mu = x, predict = TRUE))
    }
  }
  
  # add length of credible intervals
  table[ ,4] <- table[ ,3] - table[ ,2]
  
  # set H distances under the same benchmarks to exactly 0
  pos.j <- n.act
  pos.hn0 <- 1
  if(show.re == FALSE)
    for(i in 1:3){
      table[n.act*(i-1) + pos.j + 1, 6] <- 0
      table[n.act*(i-1) + pos.hn0 + 1, 5] <- 0
    }
  else {
  for(i in 1:(k+3)){
    table[n.act*(i-1) + pos.j + 1, 6] <- 0
    table[n.act*(i-1) + pos.hn0 + 1, 5] <- 0
    }
  }
  
  colnames(table) <- c("estimate", "CrI_low", "CrI_up", "length_CrI",
                       "H(po_HN0, po_act)", "signed_inf")
  
  rownames.mu <- rep(NA, times = n.act + 1)
  rownames.tau <- rep(NA, times = n.act)
  rownames.theta.new <- rep(NA, times = n.act)
  rownames.theta.i <- rep(NA, times = n.act * k)
  
  rownames.mu[1] <- "mu, FE"
  rownames.mu[2] <- "mu, HN0"
  rownames.mu[n.act + 1] <- "mu, J"
  
  rownames.tau[1] <- "tau, HN0"
  rownames.tau[n.act] <- "tau, J"
  rownames.theta.new[1] <- "theta_new, HN0"
  rownames.theta.new[n.act] <- "theta_new, J"
  
  for(j in 1:k){
    rownames.theta.i[(j-1)*(n.act)+1] <- paste0("theta_", j, ", HN0", sep = "")
    rownames.theta.i[j*n.act] <- paste0("theta_", j, ", J", sep = "")
  }
  
  for (i in 1:(n.act-2)) {
    rownames.mu[i+2] <- paste0("mu, pri_act_", i, sep = "")
    rownames.tau[i+1] <- paste0("tau, pri_act_", i, sep = "")
    rownames.theta.new[i+1] <- paste0("theta_new, pri_act_", 
                                      i, sep = "")
    for (j in 1:k) {
      ii <- (j - 1) * n.act
      rownames.theta.i[ii + i + 1] <- paste0("theta_", j, ", pri_act_", 
                                             i, sep = "")
    }
  }
  if(show.re == FALSE)
    rownames <- c(rownames.mu, rownames.tau,
                  rownames.theta.new)
  else
    rownames <- c(rownames.mu, rownames.tau, rownames.theta.i, 
                  rownames.theta.new)
  
  row.names(table) <- rownames
  return(table)
}


