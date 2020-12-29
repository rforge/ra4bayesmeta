
# computes a table of Hellinger distances betw. the actual priors for tau
# in the specified fits fits.actual and the benchmark priors for tau in fits.bm
pri_RA_5bm <- function(df, tau.prior=list(function(x) dhalfnormal(x, scale=1)),
                       m_J=NA, M_J=NA, upper.J=3, digits.J=2,
                       m_inf=NA, M_inf=NA, rlmc0=0.0001, rlmc1=0.9999,
                       mu.mean=0, mu.sd=4){
  # inputs:
  
  # output:
  # matrix with one row per actual prior giving the Hellinger distances
  # between the actual priors and the benchmark priors
  
  # rlmc0 <- 0.0001
  # rlmc1 <- 0.9999
  
   if(is.na(m_inf))
    m_inf <- m_inf_sgc(rlmc=rlmc0)
   if(is.na(M_inf))
    M_inf <- M_inf_sigc(rlmc=rlmc1, df=df)
  
    # upper.J <- 3
    # digits.J <- 2
   if(is.na(m_J))
    m_J <- m_j_sgc(df=df, upper=upper.J, digits=digits.J, mu.mean=mu.mean, mu.sd=mu.sd)
   if(is.na(M_J))
    M_J <- M_j_sigc(df=df, upper=upper.J, digits=digits.J, mu.mean=mu.mean, mu.sd=mu.sd)
   
   # thres <- 5*10^6
   # 
   # if(m_inf > thres)
   #   warning(paste0("m_inf=", round(m_inf,0), 
   #                  ">5e+06. This may cause numerical problems in the bayesmeta() function.", sep=""))
   # if(M_inf > thres)
   #   warning(paste0("M_inf=", round(M_inf,0), 
   #                  ">5e+06. This may cause numerical problems in the bayesmeta() function.", sep=""))
   # if(m_J > thres)
   #   warning(paste0("m_J=", round(m_J,0), 
   #                  ">5e+06. This may cause numerical problems in the bayesmeta() function.", sep=""))
   # if(M_J > thres)
   #   warning(paste0("M_J=", round(M_J,0), 
   #                  ">5e+06. This may cause numerical problems in the bayesmeta() function.", sep=""))
    
    C <- sigma_ref(df)^{-2}
    
    par <- c(m_inf=m_inf, M_J=M_J, m_J=m_J, M_inf=M_inf, C=C)
  
  n.act <- length(tau.prior)
  table <- matrix(nrow=n.act, ncol=4)
  
  for(i in 1:length(tau.prior)){
    # for(j in 1:length(fits.b){
    if(is.function(tau.prior[[i]])){
      table[i,1] <- H( tau.prior[[i]], 
                      function(x) dsgc(x, m=m_inf, C=C), lower=0)
      table[i,2] <- H( tau.prior[[i]], 
                       function(x) dsigc(x, M=M_J, C=C), lower=0)
      table[i,3] <- H( tau.prior[[i]], 
                       function(x) dsgc(x, m=m_J, C=C), lower=0)
      table[i,4] <- H( tau.prior[[i]], 
                       function(x) dsigc(x, M=M_inf, C=C), lower=0)
    }
    if(is.character(tau.prior[[i]])){
      fit <- bayesmeta(y=df[,"y"], sigma=df[,"sigma"],
                       mu.prior.mean=mu.mean, mu.prior.sd=mu.sd,
                       tau.prior=tau.prior[[i]])
      table[i,1] <- H( function(x) fit$tau.prior(x), 
                       function(x) dsgc(x, m=m_inf, C=C), lower=0)
      table[i,2] <- H( function(x) fit$tau.prior(x), 
                       function(x) dsigc(x, M=M_J, C=C), lower=0)
      table[i,3] <- H( function(x) fit$tau.prior(x), 
                       function(x) dsgc(x, m=m_J, C=C), lower=0)
      table[i,4] <- H( function(x) fit$tau.prior(x), 
                       function(x) dsigc(x, M=M_inf, C=C), lower=0)
    }
  }

  colnames(table) <- c("H(SGC(m_inf), pri_act)", "H(SIGC(M_J), pri_act)",
                       "H(SGC(m_J), pri_act)", "H(SIGC(M_inf), pri_act)")
  
  rownames <- rep(NA, times = n.act)
  for(j in 1:n.act){
    rownames[j] <- paste0("pri_act_",j, sep="")
  }
  rownames(table) <- rownames
  # rownames(table) <- c("pri_act_1", "pri_act_2", "pri_act_3", "pri_act_4", "pri_act_5",
  #                      "pri_act_6", "pri_act_7", "pri_act_8", "pri_act_9", "pri_act_10",
  #                      "pri_act_11", "pri_act_12", "pri_act_13", "pri_act_14", "pri_act_15",
  #                      "pri_act_16", "pri_act_17", "pri_act_18", "pri_act_19", "pri_act_20")[1:n.act]
  
  
  return(list(table, par))
}