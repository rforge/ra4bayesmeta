
post_RA_fits <- function(fits.actual, fits.bm, 
                         H.dist.method = "integral"){
  # inputs:
  # fits.actual: a list of bayesmeta fits under actual priors
  # fits.bm: a list of bayesmeta fits under benchmark priors, 
  # with "most extreme" benchmarks first and last
  # output:
  # matrix with (k+3) rows and length(fits.bm)+1 columns giving
  # the Hellinger distances betw. marginal posteriors for different parameters
  # under the actual prior vs. one benchmark prior,
  # and the Hellinger distance betw. the 2 most extreme benchmarks
  
  k <- fits.actual[[1]]$k
  n.bm <- length(fits.bm)
  n.act <- length(fits.actual)
  n.row <- (3+k)*n.act
  # p <- length(fits.bm)+1
  tab.H <- matrix(nrow=n.row, ncol=n.bm)
  
  # loop over parameters
  # for(m in 1:(k+3)){
  #    ii <- (m-1)*n.act
  for(i in 1:n.act){
    # fit.actual <- fits.actual[[i]]
    # ii <- (i -1)*(k+3)
    for(j in 1:n.bm){
      ii <- 0
      tab.H[ii+i,j] <- H_fits(fits.actual[[i]], fits.bm[[j]],
                              parameter = "mu", method = H.dist.method)
      
      ii <- 1*n.act
      tab.H[ii+i,j] <- H_fits(fits.actual[[i]], fits.bm[[j]],
                              parameter = "tau", method = H.dist.method)
      for(l in 2:(k+1)){
        ii <- l*n.act
        tab.H[ii+i,j] <- H_fits(fits.actual[[i]], fits.bm[[j]],
                                parameter = "theta", individual=l-1,
                                method = H.dist.method)
      }
      ii <- (k+2)*n.act
      tab.H[ii+i,j] <- H_fits(fits.actual[[i]], fits.bm[[j]],
                              parameter = "theta_new", method = H.dist.method)
    }
  }
  
  colnames <- rep(NA, times = n.bm)
  for(j in 1:n.bm){
    colnames[j] <- paste0("H(po_{bm_",j,"}, po_act)", sep="")
  }
  
  rownames.mu <- rep(NA, times = n.act)
  rownames.tau <- rep(NA, times = n.act)
  rownames.theta.new <- rep(NA, times = n.act)
  rownames.theta.i <- rep(NA, times = n.act*k)
  for(j in 1:n.act){
    rownames.mu[j] <- paste0("mu, pri_act_",j, sep="")
    rownames.tau[j] <- paste0("tau, pri_act_",j, sep="")
    rownames.theta.new[j] <- paste0("theta_new, pri_act_",j, sep="")
    for(i in 1:k){
      ii <- (i-1)*n.act
      rownames.theta.i[ii+j] <- paste0("theta_",i, ", pri_act_",j, sep="")
    } 
  }
  rownames <- c(rownames.mu, rownames.tau, 
                rownames.theta.i, rownames.theta.new)
  
  colnames(tab.H) <- colnames
  rownames(tab.H) <- rownames
  
  return(tab.H)
}