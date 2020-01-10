
####----  normal calibration function for Hellinger distance: N(0,1) vs N(muh,1) for muh in (0,7) ----####

cal_h_dist <-function(h){
  # Calibration of the Hellinger distance with respect to a N(0,1) distribution
  # see Roos et al. (2015) for deteils
  return(sqrt(-8*log(1-h^2)))
}