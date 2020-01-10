
# density function of the SGC(m,C) prior
dsgc <- function(x, m, C) {
  res <- rep(-Inf, length(x))
  res[x > 0] <- log(2)+log(C)+log(m-1)+log(x[x > 0])-m*log(1+C*x[x > 0]^2)
  result <- exp(res)
  return(result)
  # return(exp(log(2)+log(C)+log(m-1)+log(x)-m*log(1+C*x^2)))
}