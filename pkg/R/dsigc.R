
# density function of the SIGC(M,C) prior
dsigc <- function(x, M, C) {
  res <- rep(-Inf, length(x))
  res[x > 0] <- log(4)+log(C)+log(M-1)-5*log(x[x > 0])-M*log(1+C*x[x > 0]^{-4})
  result <- exp(res)
  return(result)
}