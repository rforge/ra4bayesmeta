
####----  normal calibration function for Hellinger distance ---####

# computes either the shift(s) between two unit variance normal distribution
# with the given Hellinger distance(s) (output = "shift", see Roos et al. (2015) for details)
# or the corresponding area(s) of overlap between these two normal distributions (output = "ao")

cal_h_dist <- function (h, output = "shift"){
  
  shift <- sqrt(-8 * log(1 - h^2))
  
  if(output == "shift")
    return(shift)
  if(output == "ao"){
    ao <- pnorm(shift/2, mean=shift, sd=1) + (1-pnorm(shift/2, mean=0, sd=1))
    return(ao)
  }
}