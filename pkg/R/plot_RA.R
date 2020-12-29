
plot_RA <- function(df, tau.prior = list(),
                    type = "pri.tau", improper.prior = NULL, show.sigma.i = FALSE,
                    xlim, ylim, legend = FALSE, 
                    pos.legend = "topright", legend.tau.prior = c(), 
                    xlab = NULL, bty = "o", 
                    scale.hn0 = 1/500, mu.mean = 0, mu.sd = 4) {
  
  # default xlab depends on the parameter to show specified in argument "type"
  if(is.null(xlab)){
    xlab <- switch(type, 
                   pri.tau = expression("heterogeneity "*tau),
                   post.mu = expression("effect "*mu),
                   post.tau = expression("heterogeneity "*tau),
                   post.theta.new = expression("effect "*theta[new])
    )
  }

  fits <- fit_models_RA(df = df, tau.prior = tau.prior, scale.hn0 = scale.hn0,
                        mu.mean = mu.mean, mu.sd = mu.sd)
  n <- length(fits)
  y <- list()
  if (type == "pri.tau") {
    x <- c(seq(from = 0, to = 0.2, length = 1000), seq(from = 0.2, 
                                                       to = xlim[2], length = 5000))
    xlab <- expression("heterogeneity " * tau)
    ylab <- "prior density"
    for (i in 1:n) {
      y[[i]] <- fits[[i]]$tau.prior(x)
    }
  }
  if (type == "post.mu") {
    x <- seq(from = xlim[1], to = xlim[2], length = 2000)
    # xlab <- expression("effect " * mu)
    ylab <- "posterior density"
    for (i in 1:n) {
      y[[i]] <- fits[[i]]$dposterior(mu = x)
    }
    mean.fe <- post_mu_fe(df = df,  mu.mean = mu.mean, mu.sd = mu.sd)$"mean"
    sd.fe <- post_mu_fe(df = df,  mu.mean = mu.mean, mu.sd = mu.sd)$"sd"
    # y[[n+1]] <- function(x) dnorm(x, mean = mean.fe, sd = sd.fe)
  }
  if (type == "post.tau") {
    x <- c(seq(from = 0, to = 0.2, length = 4000), seq(from = 0.2, 
                                                       to = xlim[2], length = 5000))
    # xlab <- expression("heterogeneity " * tau)
    ylab <- "posterior density"
    for (i in 1:n) {
      y[[i]] <- fits[[i]]$dposterior(tau = x)
    }
  }
  if (type == "post.theta.new") {
    x <- seq(from = xlim[1], to = xlim[2], length = 2000)
    # xlab <- expression("effect " * theta[new])
    ylab <- "posterior density"
    for (i in 1:n) {
      y[[i]] <- fits[[i]]$dposterior(mu = x, predict = TRUE)
    }
  }
  lwd <- 2
  mycol <- c(5, "blue",
             # 1, 
             # "blue", "darkgray", 
             "orange", "red", "lightpink3", 
             "darkgreen", "green", "violetred")
  mycol.fe <- 1
  mylty <- c(rep(1, times = 2), 2, 2, 2, 2)
  mylty.fe <- 4
  plot(x, y[[2]], type = "n", xlab = xlab, ylab = ylab, xlim = xlim, 
       ylim = ylim)
  
  
  
  i <- 1
  lines(x, y[[i]], type = "l", lty = mylty[i], lwd = lwd, 
        col = mycol[i])
  
  if(type == "post.tau" | type == "post.theta.new")  {
    # if (type != "pri.tau") {
    lines(x, y[[2]], type = "l", lty = mylty[2], lwd = lwd, 
          col = mycol[2])
    if (n > 2) {
      for (i in 3:n) 
        lines(x, y[[i]], type = "l", lty = mylty[i], 
              lwd = lwd, col = mycol[i])
    }
  }
  if(type == "post.mu") {
    lines(x, y[[2]], type = "l", lty = mylty[2], lwd = lwd, 
          col = mycol[2])
    if (n > 2) {
      for (i in 3:n) 
        lines(x, y[[i]], type = "l", lty = mylty[i], 
              lwd = lwd, col = mycol[i])
    }
    # add the FE benchmark
    lines(x, dnorm(x, mean = mean.fe, sd = sd.fe), type = "l", lty = mylty.fe, lwd = lwd,
          # lty = mylty[2], lwd = lwd, 
          col = mycol.fe)
  }
  if (type == "pri.tau") {
    if (n > 2) {
      for (i in 3:n) 
        if(!((i-2) %in% improper.prior))
          lines(x, y[[i]], type = "l", lty = mylty[i], 
                lwd = lwd, col = mycol[i])
    }
    if(show.sigma.i == TRUE){
      k <- fits[[1]]$k
      points(x = fits[[1]]$sigma, y = rep(0, times = k), pch = 3, 
             cex = 0.5)
    }
  }
  cex.legend <- 1
  if (legend == TRUE && type == "pri.tau") {
    if (length(legend.tau.prior) == 0) 
      legend(pos.legend, c("HN0"), col = mycol[1], 
             lwd = 2, lty = mylty[1], bg = "white", bty = bty, 
             cex = cex.legend)
    if (length(legend.tau.prior) > 0 && n > 2) {
      ind.act <- 3:n
      ind <- c(1, ind.act)
      legend(pos.legend, c("HN0", legend.tau.prior), col = mycol[ind], 
             lwd = 2, lty = mylty[ind], bg = "white", bty = bty, 
             cex = cex.legend)
    }
  }
  if (legend == TRUE && type == "post.mu") {
    if (length(legend.tau.prior) == 0) 
      legend(pos.legend, c("FE", "HN0", "J"), col = c(mycol.fe, mycol[1:2]), 
             lwd = 2, lty = c(mylty.fe, mylty[1:2]), bg = "white", bty = bty, 
             cex = cex.legend)
    if (length(legend.tau.prior) > 0 && n > 2) {
      ind.act <- 3:n
      ind <- c(1, ind.act, 2)
      legend(pos.legend, c("FE", "HN0", legend.tau.prior, "J"), col = c(mycol.fe, mycol[ind]), 
             lwd = 2, lty = c(mylty.fe, mylty[ind]), bg = "white", bty = bty, 
             cex = cex.legend)
    }
  }
  if (legend == TRUE && (type == "post.tau" | type == "post.theta.new")) {
    if (length(legend.tau.prior) == 0) 
      legend(pos.legend, c("HN0", "J"), col = mycol[1:2], 
             lwd = 2, lty = mylty[1:2], bg = "white", bty = bty, 
             cex = cex.legend)
    if (length(legend.tau.prior) > 0 && n > 2) {
      ind.act <- 3:n
      ind <- c(1, ind.act, 2)
      legend(pos.legend, c("HN0", legend.tau.prior, "J"), col = mycol[ind], 
             lwd = 2, lty = mylty[ind], bg = "white", bty = bty, 
             cex = cex.legend)
    }
  }
}