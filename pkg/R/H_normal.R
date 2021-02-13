# Hellinger distance based on mean and sd estimates of two normal distributions

H_normal <- function(mean1, sd1, mean2, sd2){

# the sd-part of the BC under normal approximation
# this part quantifies the spread modification
BC_normal_sd_part <- sqrt((2*sd1*sd2)/(sd1^2+sd2^2))

# the mean-part of the BC under normal approximation
# the mean-part is adjusted for standard deviations (it corresponds to the Mahalanobis distance)
# this part quantifies the location modification
BC_normal_mean_part <- exp(-((mean1-mean2)^2/(4*(sd1^2+sd2^2))))

# computation of the total BC under normality assumption
BC_normal_total <- BC_normal_sd_part * BC_normal_mean_part

  ## computation of the total Hellinger distance H under normality assumption 
  return(sqrt(1 - BC_normal_total))
}