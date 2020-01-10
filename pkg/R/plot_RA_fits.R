
# TODO: add generic legend in df vers. of this func. only?

# plot prior for tau and marginal posteriors for mu, tau, theta_new
plot_RA_fits <- function(fits.actual, fits.bm, type="pri.tau", xlim, ylim,
                    legend=FALSE, pos.legend="topright", legend.tau.prior=c(), bty="o"){
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
    x <- c(seq(from=0, to=0.2, length=1000), seq(from=0.2, to=5, length=5000))
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
    x <- seq(from=-10, to=10, length=2000)
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
    x <- c(seq(from=0, to=0.2, length=4000), seq(from=0.2, to=5, length=5000))
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
    x <- seq(from=-10, to=10, length=2000)
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
  lwd <- 2
  
  ind.act <- 6:(6+n.act-1)
  ind.bm <- 1:n.bm
  ind <- c(ind.bm, ind.act)
  mycol <- c(5,1,"blue","darkgray","orange","red", "lightpink3", "darkgreen", "green", "violetred")
  # mycol <- c(5,1,"blue","darkgray","orange","red", "lightpink3", "magenta4", "green")
  # mycol <- c(5,1,"blue","darkgray","orange","red", "lightpink3", "magenta4", "darkgoldenrod")
  # mycol <- c(5,1,"blue","darkgray","orange","red", "firebrick", "green")
  mylty <- c(rep(1, times=5), rep(2, times=5))
  
  mycol2 <- c(5,1,"blue","darkgray","orange","red", "lightpink3", "darkgreen", "green", "violetred")[ind]
  # mycol <- c(5,1,"blue","darkgray","orange","red", "lightpink3", "magenta4", "green")
  # mycol <- c(5,1,"blue","darkgray","orange","red", "lightpink3", "magenta4", "darkgoldenrod")
  # mycol <- c(5,1,"blue","darkgray","orange","red", "firebrick", "green")
  mylty2 <- c(rep(1, times=5), rep(2, times=5))[ind]
  # if(is.na(ylim[1]) || is.na(ylim[2]))
  #   plot(x, y[[2]], type="n", xlab=xlab, ylab=ylab, xlim=xlim)
  # else
  plot(x, y[[2]], type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
  
  # plot fits under the benchmark priors
  for(i in 1:n){
    lines(x, y[[i]] , type="l", lty=mylty2[i], lwd=lwd, col=mycol2[i])
  }
  # add fit under Jeffreys prior except for prior plot
  # if(type!="pri.tau"){
  #   lines(x, y[[5]] , type="l", lty=mylty[5], lwd=lwd, col=mycol[5])
  # }
  # add fits under specified prior(s) for tau if available
  # if(n.act>0){
  #   for(i in (n.bm+1):n)
  #     lines(x, y[[i]], type="l", lty=mylty[i], lwd=lwd, col=mycol[i])
  # }
  # if(length(y)>=6){
  #   lines(x, y[[6]], type="l", lty=mylty[6], lwd=lwd, col=mycol[6])
  # }
  # if(length(y)>=7){
  #   lines(x, y[[7]], type="l", lty=mylty[7], lwd=lwd, col=mycol[7])
  # }
  # if(length(y)>=8){
  #   lines(x, y[[8]], type="l", lty=mylty[8], lwd=lwd, col=mycol[8])
  # }
  # if(length(y)>=9){
  #   lines(x, y[[9]], type="l", lty=mylty[9], lwd=lwd, col=mycol[9])
  # }
  
  # add sigma_i values to plot of prior for tau
  if(type=="pri.tau"){
    # extract the number of studies in the data set
    k <- fits[[1]]$k
    points(x=fits[[1]]$sigma, y=rep(0, times=k), pch=3, cex=0.5)
  }
  
  cex.legend <- 1.0
  # cex.legend <- 0.8
  # add a legend if legend==TRUE
  # if(legend==TRUE && type=="pri.tau"){
  #   if(length(legend.tau.prior)==0)
  #     legend(pos.legend, c(expression(SCG(m[infinity])), expression(SIGC(M[J])), 
  #                          expression(SGC(m[J])), expression(SICG(M[infinity]))),
  #            col=mycol[1:4], lwd=2, lty=mylty[1:4], bg="white")
  #   
  #   if(length(legend.tau.prior)>0 && n>5){
  #     # ind.act <- seq(from=n, to=6, by=-1)
  #     ind.act <- 6:n
  #     ind <- c(ind.act, 1:4)
  #     
  #     legend(pos.legend, c(legend.tau.prior,  expression(SCG(m[infinity])), expression(SIGC(M[J])), 
  #                          expression(SGC(m[J])), expression(SICG(M[infinity]))),
  #            col=mycol[ind], lwd=2, lty=mylty[ind], bg="white",
  #            cex=cex.legend)
      
      # if(n==6)
      #   legend(pos.legend, c(legend.tau.prior,  expression(SCG(m[infinity])), expression(SIGC(M[J])), 
      #                        expression(SGC(m[J])), expression(SICG(M[infinity]))),
      #          col=mycol[c(6,1:4)], lwd=2, lty=mylty[c(6,1:4)], bg="white",
      #          cex=cex.legend)
      # if(n==7)
      #   legend(pos.legend, c(legend.tau.prior,  expression(SCG(m[infinity])), expression(SIGC(M[J])), 
      #                        expression(SGC(m[J])), expression(SICG(M[infinity]))),
      #          col=mycol[c(7,6,1:4)], lwd=2, lty=mylty[c(7,6,1:4)], bg="white",
      #          cex=cex.legend)
      # if(n==8)
      #   legend(pos.legend, c(legend.tau.prior,  expression(SCG(m[infinity])), expression(SIGC(M[J])), 
      #                        expression(SGC(m[J])), expression(SICG(M[infinity]))),
      #          col=mycol[c(8,7,6,1:4)], lwd=2, lty=mylty[c(8,7,6,1:4)], bg="white", 
      #          cex=cex.legend)
      # if(n==9)
      #   legend(pos.legend, c(legend.tau.prior,  expression(SCG(m[infinity])), expression(SIGC(M[J])), 
      #                        expression(SGC(m[J])), expression(SICG(M[infinity]))),
      #          col=mycol[c(9,8,7,6,1:4)], lwd=2, lty=mylty[c(9,8,7,6,1:4)], bg="white", 
      #          cex=cex.legend)
    
  
  if(legend==TRUE){
    if(length(legend.tau.prior)==0)
      warning("legend entries need to be specified")
      
    
    if(length(legend.tau.prior)>0){
      # ind.act <- seq(from=n, to=6, by=-1)
      
      ind.leg <- c(ind.act, ind.bm)
      legend(pos.legend, legend.tau.prior,
             col=mycol[ind.leg], lwd=2, lty=mylty[ind.leg], bg="white", bty=bty,
             cex=cex.legend)
    }
  }
      # if(n==6)
      #   legend(pos.legend, c(legend.tau.prior, expression(SCG(m[infinity])), expression(SIGC(M[J])),
      #                        "Jeffreys", expression(SGC(m[J])), expression(SICG(M[infinity]))),
      #          col=mycol[c(6,1:2,5,3:4)], lwd=2, lty=mylty[c(6,1:2,5,3:4)], bg="white",
      #          cex=cex.legend)
      # if(n==7)
      #   legend(pos.legend, c(legend.tau.prior, expression(SCG(m[infinity])), expression(SIGC(M[J])),
      #                        "Jeffreys", expression(SGC(m[J])), expression(SICG(M[infinity]))),
      #          col=mycol[c(7,6,1:2,5,3:4)], lwd=2, lty=mylty[c(7,6,1:2,5,3:4)], bg="white",
      #          cex=cex.legend)
      # if(n==8)
      #   legend(pos.legend, c(legend.tau.prior, expression(SCG(m[infinity])), expression(SIGC(M[J])),
      #                        "Jeffreys", expression(SGC(m[J])), expression(SICG(M[infinity]))),
      #          col=mycol[c(8,7,6,1:2,5,3:4)], lwd=2, lty=mylty[c(8,7,6,1:2,5,3:4)], bg="white",
      #          cex=cex.legend)
      # if(n==9)
      #   legend(pos.legend, c(legend.tau.prior, expression(SCG(m[infinity])), expression(SIGC(M[J])),
      #                        "Jeffreys", expression(SGC(m[J])), expression(SICG(M[infinity]))),
      #          col=mycol[c(9,8,7,6,1:2,5,3:4)], lwd=2, lty=mylty[c(9,8,7,6,1:2,5,3:4)], bg="white",
      #          cex=cex.legend)
}