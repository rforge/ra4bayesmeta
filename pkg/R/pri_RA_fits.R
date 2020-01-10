
# computes a table of Hellinger distances betw. the actual priors for tau
# in the specified fits fits.actual and the benchmark priors for tau in fits.bm
pri_RA_fits <- function(fits.actual, fits.bm){
  # inputs:
  # fits.actual: list of bayesmeta fits under actual priors for tau of interest
  # fits.bm: list of bayesmeta fits under benchmark priors for tau (usually 4 fits)
  # output:
  # matrix with one row per actual prior giving the Hellinger distances
  # between the actual priors and the benchmark priors
  
  n.act <- length(fits.actual)
  n.bm <- length(fits.bm)
  tab.H <- matrix(nrow=n.act, ncol=n.bm)
  
  for(i in 1:n.act){
    for(j in 1:n.bm){
      tab.H[i,j] <- H(function(x) fits.actual[[i]]$dprior(tau=x), 
                      function(x) fits.bm[[j]]$dprior(tau=x), lower=0)
    }
  }
  
  colnames <- rep(NA, times = n.bm)
  for(j in 1:n.bm){
    colnames[j] <- paste0("H(pri_bm_",j,", pri_act)", sep="")
  }
  
  rownames <- rep(NA, times = n.act)
  for(j in 1:n.act){
    rownames[j] <- paste0("pri_act_",j, sep="")
  }
  
  colnames(tab.H) <- colnames
  rownames(tab.H) <- rownames
  
  return(tab.H)
}