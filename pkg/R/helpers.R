
# functions to compute Hellinger distances between posterior benchmarks

# J and m_J benchmarks
H_dist_J_mJ <- function(df, mJ, 
                     mu.mean=0, mu.sd=4){
  # inputs:
  # df: data frame in bayesmeta format
  # mJ: value of m for the SGC(m) prior
  # mu.mean, mu.sd: parameter values of the normal prior for mu
  # output:
  # Hellinger distance betw. the marg. posteriors for tau
  # under Jeffreys prior and prior SGC(mJ)
 
  
  C <- sigma_ref(df=df)^{-2}
  
  fit.j <- bayesmeta(y=df[,"y"], sigma=df[,"sigma"],
                     mu.prior.mean=mu.mean, mu.prior.sd=mu.sd,
                     tau.prior="Jeffreys")
  
  fit.SGC.mJ <- bayesmeta(y=df[,"y"], sigma=df[,"sigma"],
                           mu.prior.mean=mu.mean, mu.prior.sd=mu.sd,
                           tau.prior=function(t){dsgc(t, m=mJ, C=C)})
  
  
  h <- H(function(x) fit.SGC.mJ$dposterior(tau=x), 
                   function(x) fit.j$dposterior(tau=x), lower=0)
  
  
  return(h)
}


# J and M_J benchmarks
H_dist_J_MJ <- function(df, MJ, 
                        mu.mean=0, mu.sd=4){
  # inputs:
  # df: data frame in bayesmeta format
  # mJ: value of m for the SIGC(m) prior
  # mu.mean, mu.sd: parameter values of the normal prior for mu
  # output:
  # Hellinger distance betw. the marg. posteriors for tau
  # under Jeffreys prior and prior SIGC(MJ)
 
  
  
  C <- sigma_ref(df=df)^{-2}
  
  fit.j <- bayesmeta(y=df[,"y"], sigma=df[,"sigma"],
                     mu.prior.mean=mu.mean, mu.prior.sd=mu.sd,
                     tau.prior="Jeffreys")
  
  fit.SIGC.MJ <- bayesmeta(y=df[,"y"], sigma=df[,"sigma"],
                            mu.prior.mean=mu.mean, mu.prior.sd=mu.sd,
                            tau.prior=function(t){dsigc(t, M=MJ, C=C)})
  
  h <- H(function(x) fit.j$dposterior(tau=x), 
                    function(x) fit.SIGC.MJ$dposterior(tau=x), lower=0)
  
  return(h)
}