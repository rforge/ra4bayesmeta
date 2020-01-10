
# computation of the reference standard deviation as suggested in equation (7) by Sorbye and Rue (scaling)
# sigma_ref by geometric mean of sigmai
sigma_ref<-function(df){
  # input:
  # df: data frame with one column "sigma" containing the standard deviations sigmai in each study
  # output:
  # reference standard deviation as suggested in equation (7) by Sorbye and Rue (scaling)
  return(exp(mean(log(df$sigma))))
}