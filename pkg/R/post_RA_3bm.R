
post_RA_3bm <- function(df, tau.prior=list(function(x) dhalfnormal(x, scale=1)),
                        H.dist.method = "integral",
                        m_inf=NA, M_inf=NA, rlmc0=0.0001, rlmc1=0.9999,
                        mu.mean=0, mu.sd=4){
  
  if(is.na(m_inf))
    m_inf <- m_inf_sgc(rlmc=rlmc0)
  if(is.na(M_inf))
    M_inf <- M_inf_sigc(rlmc=rlmc1, df=df)
  
  
  thres <- 5*10^6
  
  if(m_inf > thres)
    warning(paste0("m_inf=", round(m_inf,0), 
                   ">5e+06. This may cause numerical problems in the bayesmeta() function.", sep=""))
  if(M_inf > thres)
    warning(paste0("M_inf=", round(M_inf,0), 
                   ">5e+06. This may cause numerical problems in the bayesmeta() function.", sep=""))
  
  
  res.fit <- fit_models_RA_5bm(df=df, tau.prior=tau.prior, compute.J.bm=FALSE,
                               m_J=NA, M_J=NA, upper.J=3, digits.J=1,
                               m_inf=m_inf, M_inf=M_inf, rlmc0=rlmc0, rlmc1=rlmc1,
                               mu.mean=mu.mean, mu.sd=mu.sd)
  fits <- res.fit[[1]]
  par <- res.fit[[2]][c(1,4,5)]
  
  fits.bm <- fits[c("fit.SGC.m_inf", "fit.j", "fit.SIGC.M_inf")]
  n.act <- length(tau.prior)
  fits.actual <- fits[c(6:(n.act+5))]
  
 
  k <- fits.actual[[1]]$k
  n.bm <- length(fits.bm)
  
  # p <- length(fits.bm)+1
  table <- matrix(nrow=(3+k)*n.act, ncol=n.bm)
  
  # loop over parameters
  # for(m in 1:(k+3)){
  #    ii <- (m-1)*n.act
    for(i in 1:n.act){
    # fit.actual <- fits.actual[[i]]
    # ii <- (i -1)*(k+3)
      for(j in 1:n.bm){
        ii <- 0
        table[ii+i,j] <- H_fits(fits.actual[[i]], fits.bm[[j]],
                                parameter = "mu", method = H.dist.method)
        
        ii <- 1*n.act
        table[ii+i,j] <- H_fits(fits.actual[[i]], fits.bm[[j]],
                                parameter = "tau", method = H.dist.method)
        
        for(l in 2:(k+1)){
          ii <- l*n.act
          table[ii+i,j] <- H_fits(fits.actual[[i]], fits.bm[[j]],
                                  parameter = "theta", individual=l-1,
                                  method = H.dist.method)
        }
        ii <- (k+2)*n.act
        table[ii+i,j] <- H_fits(fits.actual[[i]], fits.bm[[j]],
                                parameter = "theta_new", method = H.dist.method)
     }
    }
  
  colnames(table) <- c("H(po_{m_inf}, po_act)", "H(po_J, po_act)",
                       "H(po_{M_inf}, po_act)")
  
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
  row.names(table) <- rownames
  
  return(list(table, par))
}