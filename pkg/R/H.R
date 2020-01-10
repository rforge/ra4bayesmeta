
# Hellinger distance
H <- function(dens1, dens2, lower=-Inf){
  # dens1 and dens2: densities as functions of x
  # lower: lower limit of integration (usually -Inf or 0)
  
  # compute the Bhattacharyya coefficient (BC)
  integrand <- function(x) exp(1/2*(log(dens1(x))+log(dens2(x))))
  BC <- integrate(integrand, lower = lower, upper = Inf)$value
  
  if(BC>1)
    BC <- 1
  
  return(sqrt(1-BC))
}