
# plot prior for tau and marginal posteriors for mu, tau, theta_new
plot_RA_fits <- function(fits.actual, fits.bm, type="pri.tau", xlim, ylim,
                         legend=FALSE, pos.legend="topright", legend.tau.prior=c(), bty="o",
                         col.actual=c("red","lightpink3","darkgreen","green", 
                                      "violetred")[1:length(fits.actual)], 
                         col.bm=c("cyan","black","blue","darkgray",
                                  "dodgerblue")[1:length(fits.bm)], 
                         lty.actual=rep(2, times=length(col.actual)), 
                         lty.bm=rep(1, times=length(col.bm)),
                         lwd.actual=rep(2, times=length(col.actual)) , 
                         lwd.bm=rep(2, times=length(col.bm))){
  # inputs:
  # fits: list of 6, 7, 8 or 9 (or 5 if no additional actual prior for tau is specified) bayesmeta fits
  # order: list(fit.SGC.m_inf, fit.SIGC.M_J, fit.SGC.m_J, fit.SIGC.M_inf, 
  #             fit.j, fit.main), where fit.j is the fit under Jeffreys' prior and 
  #             fit.main the fit under the actual prior for tau
  #         and 
  # type: specifies the parameter of interest and 
  #       if the prior or marg. posterior should be plotted
  # availbable options: "pri.tau", "post.mu", "post.tau", "post.theta.new"
  # legend: if TRUE, a legend will be added to the plot
  # pos.legend: legend position in the plot
  # legend.tau.prior: legend entry for the actual prior for tau in fit.main (may be NA)
  # xlim, ylim: determines the range of x-values/y-values shown
  # output:
  # plot of the prior densities for tau (of type=="pri.tau") or marg. posterior densities of mu/tau/theta_new
  # under different priors for tau (including 4 benchmark priors and Jeffreys' prior)
  
  fits <- c(fits.bm, fits.actual)
  n <- length(fits)
  y <- list()
  if(type=="pri.tau"){
    x <- c(seq(from=0, to=0.2, length=1000), seq(from=0.2, to=xlim[2], length=5000))
    xlab <- expression("heterogeneity "*tau)
    ylab <- "prior density"
    
    for(i in 1:n){
      y[[i]] <- fits[[i]]$tau.prior(x)
    }
    
    # if(is.na(xlim[1]) || is.na(xlim[2])){
    #   xlim <- c(0,5)
    # }
  }
  
  if(type=="post.mu"){
    x <- seq(from=xlim[1], to=xlim[2], length=2000)
    xlab <- expression("effect "*mu)
    ylab <- "posterior density"
    
    for(i in 1:n){
      y[[i]] <- fits[[i]]$dposterior(mu=x)
    }
    
    # if(is.na(xlim[1]) || is.na(xlim[2])){
    #   xlim <- c(-3,3)
    # }
  }
  
  if(type=="post.tau"){
    x <- c(seq(from=0, to=0.2, length=4000), seq(from=0.2, to=xlim[2], length=5000))
    xlab <- expression("heterogeneity "*tau)
    ylab <- "posterior density"
    
    for(i in 1:n){
      y[[i]] <- fits[[i]]$dposterior(tau=x)
    }     
    
    # if(is.na(xlim[1]) || is.na(xlim[2])){
    #   xlim <- c(0,5)
    # }
  }
  
  if(type=="post.theta.new"){
    x <- seq(from=xlim[1], to=xlim[2], length=2000)
    xlab <- expression("effect "*theta[new])
    ylab <- "posterior density"
    
    for(i in 1:n){
      y[[i]] <- fits[[i]]$dposterior(mu=x, predict=TRUE)
    }
    
    # if(is.na(xlim[1]) || is.na(xlim[2])){
    #   xlim <- c(-3,3)
    # }
  }
  
  n.bm <- length(fits.bm)
  n.act <- length(fits.actual)
  # lwd <- 2
  
  ind.bm <- 1:n.bm
  ind.act <- (n.bm+1):(n.bm+n.act)
  ind <- c(ind.bm, ind.act)
  
  plot(x, y[[2]], type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
  
  # plot fits under the actual priors
  for(i in 1:n.act){
    lines(x, y[[n.bm+i]] , type="l", lty=lty.actual[i], lwd=lwd.actual[i], col=col.actual[i])
  }
  # plot fits under the benchmark priors
  for(i in 1:n.bm){
    lines(x, y[[i]] , type="l", lty=lty.bm[i], lwd=lwd.bm[i], col=col.bm[i])
  }
  
  
  # for(i in 1:n){
  #   lines(x, y[[i]] , type="l", lty=mylty2[i], lwd=lwd, col=mycol2[i])
  # }
  
  # add sigma_i values to plot of prior for tau
  if(type=="pri.tau"){
    # extract the number of studies in the data set
    k <- fits[[1]]$k
    points(x=fits[[1]]$sigma, y=rep(0, times=k), pch=3, cex=0.5)
  }
  
  cex.legend <- 1.0
  
  mycol <- c(col.bm[1:n.bm], col.actual[1:n.act])
  mylty <- c(lty.bm[1:n.bm], lty.actual[1:n.act])
  mylwd <- c(lwd.bm[1:n.bm], lwd.actual[1:n.act])
  
  if(legend==TRUE){
    if(length(legend.tau.prior)==0)
      warning("legend entries need to be specified")
    
    
    if(length(legend.tau.prior)>0){
      # ind.act <- seq(from=n, to=6, by=-1)
      
      ind.leg <- c(ind.act, ind.bm)
      legend(pos.legend, legend.tau.prior,
             col=mycol[ind.leg], lwd=mylwd[ind.leg], lty=mylty[ind.leg], bg="white", bty=bty,
             cex=cex.legend)
    }
  }
}