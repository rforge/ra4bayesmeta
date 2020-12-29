
# TODO: test this function!
  
# fit models for the reference analysis
# compute fits for the reference analysis using bayesmeta
fit_models_RA_5bm <- function(df, tau.prior=list(), compute.J.bm=TRUE,
                          m_J=NA, M_J=NA, upper.J=3, digits.J=2,
                          m_inf=NA, M_inf=NA, rlmc0=0.0001, rlmc1=0.9999,
                          mu.mean=0, mu.sd=4){
  # inputs:
  # df: data frame in bayesmeta format
  # M_J, m_J, M_inf, m_inf: parameter values m for SICG (M) and SGC (m) benchmark priors
  #     SIGC(M_J) and SGC(m_J) are benchmarks for Jeffreys prior;
  #     SGC(m_inf) has RLMC close to 0; SIGC(M_inf) has RLMC close to 1
  # tau.prior: list of densities of the actual priors for tau (up to 10)
  # RLMC0, RLMC1: RLMC target values for SGC(m_inf) or SIGC(M_inf) benchmarks
  # mu.mean, mu.sd: parameter values of the normal prior for mu
  # output:
  # a list of 6 (or 5 of no actual tau prior specified) bayesmeta fits, 
  # in the following order: list(fit.SGC.m_inf, fit.SIGC.M_J, fit.SGC.m_J, fit.SIGC.M_inf, 
  #                              fit.j, fit.main)
  
  if(is.na(m_inf))
    m_inf <- m_inf_sgc(rlmc=rlmc0)
  if(is.na(M_inf))
    M_inf <- M_inf_sigc(rlmc=rlmc1, df=df)
  
  if(compute.J.bm==TRUE){
    if(is.na(m_J))
      m_J <- m_j_sgc(df=df, upper=upper.J, digits=digits.J, mu.mean=mu.mean, mu.sd=mu.sd)
    if(is.na(M_J))
      M_J <- M_j_sigc(df=df, upper=upper.J, digits=digits.J, mu.mean=mu.mean, mu.sd=mu.sd)
  }
  
  thres <- 5*10^6
  
  if(m_inf > thres)
    warning(paste0("m_inf=", round(m_inf,0), 
                   ">5e+06. This may cause numerical problems in the bayesmeta() function.", sep=""))
  if(M_inf > thres)
    warning(paste0("M_inf=", round(M_inf,0), 
                   ">5e+06. This may cause numerical problems in the bayesmeta() function.", sep=""))
  if(!is.na(m_J) && m_J > thres)
    warning(paste0("m_J=", round(m_J,0), 
                   ">5e+06. This may cause numerical problems in the bayesmeta() function.", sep=""))
  if(!is.na(M_J) && M_J > thres)
    warning(paste0("M_J=", round(M_J,0), 
                   ">5e+06. This may cause numerical problems in the bayesmeta() function.", sep=""))
  
  C <- sigma_ref(df)^{-2}
  
  par <- c(m_inf=m_inf, M_J=M_J, m_J=m_J, M_inf=M_inf, C=C)
  
  
  fit.j <- bayesmeta(y=df[,"y"], sigma=df[,"sigma"],
                     mu.prior.mean=mu.mean, mu.prior.sd=mu.sd,
                     tau.prior="Jeffreys")
  
  fit.SGC.m_inf <- bayesmeta(y=df[,"y"], sigma=df[,"sigma"],
                             mu.prior.mean=mu.mean, mu.prior.sd=mu.sd,
                             tau.prior=function(x){dsgc(x, m=m_inf, C=C)})
                             # tau.prior=function(t){dsgc(t, m=m_inf, C=C)})
  
  if(compute.J.bm==TRUE){
    fit.SGC.m_J <- bayesmeta(y=df[,"y"], sigma=df[,"sigma"],
                             mu.prior.mean=mu.mean, mu.prior.sd=mu.sd,
                             tau.prior=function(x){dsgc(x, m=m_J, C=C)})
                              # tau.prior=function(t){dsgc(t, m=m_J, C=C)})
  }else{
    fit.SGC.m_J <- NA
  }
  
  fit.SIGC.M_inf <- bayesmeta(y=df[,"y"], sigma=df[,"sigma"],
                              mu.prior.mean=mu.mean, mu.prior.sd=mu.sd,
                              tau.prior=function(x){dsigc(x, M=M_inf, C=C)})
                              # tau.prior=function(t){dsigc(t, M=M_inf, C=C)})
  
  if(compute.J.bm==TRUE){
    fit.SIGC.M_J <- bayesmeta(y=df[,"y"], sigma=df[,"sigma"],
                              mu.prior.mean=mu.mean, mu.prior.sd=mu.sd,
                              tau.prior=function(x){dsigc(x, M=M_J, C=C)})
                              # tau.prior=function(t){dsigc(t, M=M_J, C=C)})
  }else{
    fit.SIGC.M_J <- NA
  }
  
  if(length(tau.prior)>0){
    n.act <- length(tau.prior)
    fit.actual <- list()
    
     fit.bms <- list(fit.SGC.m_inf, fit.SIGC.M_J, fit.SGC.m_J, 
                     fit.SIGC.M_inf, fit.j)
    for(i in 1:n.act){
      fit.actual[[i]] <- bayesmeta(y=df[,"y"], sigma=df[,"sigma"],
                            mu.prior.mean=mu.mean, mu.prior.sd=mu.sd,
                            tau.prior=tau.prior[[i]])
      # res <- c(res, fit.actual)
    }
    fits <- c(fit.bms, fit.actual)
    names.fit.actual <- rep(NA, times = n.act)
    for(j in 1:n.act){
      names.fit.actual[j] <- paste0("fit.actual_",j, sep="")
    }
    # names.fit.actual <- c("fit.actual1", "fit.actual2", "fit.actual3", "fit.actual4", "fit.actual5",
    #                       "fit.actual6", "fit.actual7", "fit.actual8", "fit.actual9", "fit.actual10")[1:n.act]
    # 
    
    names(fits) <- c("fit.SGC.m_inf", "fit.SIGC.M_J", "fit.SGC.m_J", 
                    "fit.SIGC.M_inf", "fit.j", names.fit.actual)
    res <- list(fits=fits, par=par)
    return(res)
  }
  else
    fits <- list(fit.SGC.m_inf=fit.SGC.m_inf, fit.SIGC.M_J=fit.SIGC.M_J, fit.SGC.m_J=fit.SGC.m_J, 
                fit.SIGC.M_inf=fit.SIGC.M_inf, fit.j=fit.j)
    res <- list(fits=fits, par=par)
    return(res)
}