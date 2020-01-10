
# TODO: ok to have no lower border for m as input?
  
M_j_sigc <- function(df, upper=3, digits=2,
                    mu.mean=0, mu.sd=4){
  
  if(!(digits==1 | digits==2 | digits==3)){
    stop("digits must be equal to 1, 2 or 3")
  }
  
  # if(digits==1 | digits==2){
    xx <- seq(from=1.1, to=upper, by=0.1)
    Hres <- rep(NA, times=length(xx))
    
    for(i in seq_along(xx)){
      Hres[i] <- H_dist_J_MJ(df=df, MJ=xx[i], mu.mean=mu.mean, mu.sd=mu.sd)
    } 
    
    index.min <- which(Hres==min(Hres))
    m_J <- xx[index.min]
  # }
  
  if(digits==2 | digits==3){
    m_J_global <- m_J
    if(m_J > 1.1)
      xx <- seq(from=m_J_global-0.05, to=m_J_global+0.05, by=0.01)
    else{
      xx <- seq(from=m_J_global-0.09, to=m_J_global+0.05, by=0.01)
    }
    Hres <- rep(NA, times=length(xx))
    
    for(i in seq_along(xx)){
      Hres[i] <- H_dist_J_MJ(df=df, MJ=xx[i], mu.mean=mu.mean, mu.sd=mu.sd)
    } 
    
    index.min <- which(Hres==min(Hres))
    m_J <- xx[index.min]
  }
  
  if(digits==3){
    m_J_global <- m_J
    xx <- seq(from=m_J_global-0.005, to=m_J_global+0.005, by=0.001)
    Hres <- rep(NA, times=length(xx))
    
    for(i in seq_along(xx)){
      Hres[i] <- H_dist_J_MJ(df=df, MJ=xx[i], mu.mean=mu.mean, mu.sd=mu.sd)
    } 
    
    index.min <- which(Hres==min(Hres))
    m_J <- xx[index.min]
  }
  
  return(m_J)
}



