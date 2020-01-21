\name{plot_RA_fits}
\alias{plot_RA_fits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reference analysis plot based on bayesmeta fits:
Plot heterogeneity benchmark priors and the
corresponding marginal posteriors
%%  ~~function to do ... ~~
}
\description{
Plots the benchmark heterogeneity priors inducing the specified model fits \code{fits.bm}
and the corresponding marginal benchmark posteriors for different parameters in the 
NNHM. 
Also displays the actual heterogeneity priors inducing the specified actual model fits \code{fits.actual}
and the corresponding marginal posteriors.
All bayesmeta fits should be based on the same data set.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_RA_fits(fits.actual, fits.bm, type="pri.tau", xlim, ylim,
             legend=FALSE, pos.legend="topright", 
             legend.tau.prior=c(), bty="o",
             col.actual=c("red", "lightpink3", "darkgreen", "green", 
                          "violetred"), 
             col.bm=c("cyan","black","blue","darkgray","dogerblue"), 
             lty.actual=rep(2, times=5), lty.bm=rep(1, times=5),
             lwd.actual=rep(2, times=5) , lwd.bm=rep(2, times=5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{fits.actual}{
    a list of model fits (max. 5 fits) of class bayesmeta, computed with
    the \code{bayesmeta} function in the package bayesmeta.
}
\item{fits.bm}{
    a list of model fits (max. 5 fits) of class bayesmeta, computed with
    the \code{bayesmeta} function in the package bayesmeta. 
    To be used as benchmarks.
}
\item{type}{
    specifies if heterogeneity priors or marginal posterior densities for
    a given parameter should be plotted. 
    Options are \code{"pri.tau"} for heterogeneity priors, \code{"post.mu"} for the marg. posterior of the main effect mu, \code{"post.tau"} for the marg. posterior of the heterogeneity tau and
    \code{"post.theta.new"} for the marg. posterior of the predicted effect theta_new in a new study.
}
\item{xlim}{
    a vector of two real numbers. Limits of the x-axis. (First number >= 0 for densities for tau.)
}
\item{ylim}{
    a vector of two real non-negative numbers. Limits of the y-axis.
}
\item{legend}{
    logical. Specifies if a legend should be added to the plot. Defaults to \code{FALSE}.
}
\item{pos.legend}{
    a character string specifing the position of the legend in the plot. 
    Options are all the keywords accepted by
    the \code{legend} function in the graphics package, i.e. 
    \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, \code{"topright"}, \code{"right"} and \code{"center"}.
}
\item{legend.tau.prior}{
   a vector of character strings or expressions with one
   entry for each fit in \code{fits.actual} and \code{fits.bm}
   (i.e. a vector of length length(\code{fits.actual})+length(\code{fits.bm})).
   Determines the legend entries corresponding to
   the actual und benchmark fits specified,
   in the same order as in the list c(\code{fits.actual}, \code{fits.bm}).
}
\item{bty}{
    the type of box to be drawn around the legend. The allowed values are "o" (the default) and "n".
}
\item{col.actual}{
    vector of color specifications for the lines representing the actual fits in
    \code{fits.actual}, in the same order as in that list.
    Vector of minimum length \code{length(fits.actual)}.
}
\item{col.bm}{
    vector of color specifications for the lines representing the benchmark fits in
    \code{fits.bm}, in the same order as in that list.
    Vector of minimum length \code{length(fits.bm)}.
}
\item{lty.actual}{
    vector of line type specifications for the lines representing the actual fits in
    \code{fits.actual}, in the same order as in that list.
    Vector of minimum length \code{length(fits.actual)}.
}
\item{lty.bm}{
    vector of line type specifications for the lines representing the benchmark fits in
    \code{fits.bm}, in the same order as in that list.
    Vector of minimum length \code{length(fits.bm)}.
}  
  \item{lwd.actual}{
    numeric vector specifying the width of the lines representing the actual fits in
    \code{fits.actual}, in the same order as in that list.
    Vector of minimum length \code{length(fits.actual)}.
}
\item{lwd.bm}{
    numeric vector specifying the width of the lines representing the benchmark fits in
    \code{fits.bm}, in the same order as in that list.
    Vector of minimum length \code{length(fits.bm)}.
}
%\item{lty}{
%    line types of the lines representing the different fits specified in
%    \code{c(fits.actual, fits.bm)}, in the same order as in that list.
%    Vector of length \code{length(fits.actual) + length(fits.bm)}.
%}
%\item{lwd}{
%    line widths of the lines representing the different fits specified in
%    \code{c(fits.actual, fits.bm)}, in the same order as in that list.
%    Vector of length \code{length(fits.actual) + length(fits.bm)}.
%}
}
\details{
Suggestions for prior and posterior benchmarks are provided
in Ott et al. (2020, Sections 2.5-2.7).
These benchmarks can be
computed using the function \code{fit_models_RA}.
}
\references{
Ott, M., Plummer, M., Roos, M. How vague is vague? How informative is
informative? Reference analysis for
Bayesian meta-analysis. Manuscript submitted to \emph{Statistics in Medicine}. 2020.
}
\seealso{
\code{\link[bayesmeta]{bayesmeta}} in the package bayesmeta,
\code{\link{fit_models_RA}}, \code{\link{plot_RA}} 
}
\examples{
# for aurigular acupuncture (AA) data set with two
# actual half-normal and half-Cauchy heterogeneity priors
data(aa)
# compute the model fits
fits <- fit_models_RA(df=aa, tau.prior=list(function(t)dhalfnormal(t, scale=1),
                                 function(t)dhalfcauchy(t, scale=1)))[[1]]
# benchmark fits under SGC(m_inf), SIGC(M_J), SGC(m_J) and SIGC(M_inf) priors
fits.bm.pri <- fits[1:4]
# benchmark fits under SGC(m_inf), Jeffreys' and SIGC(M_inf) priors
fits.bm.post <- fits[c(1,5,4)]
fits.actual <- fits[c(6,7)]

# prior densities
plot_RA_fits(fits.actual=fits.actual, fits.bm=fits.bm.pri,
             type="pri.tau", xlim=c(0,2), ylim=c(0,3), 
             legend=TRUE,
             legend.tau.prior=c("HN(1)", "HC(1)", 
              expression(SCG(m[infinity])), expression(SIGC(M[J])), 
              expression(SGC(m[J])), 
              expression(SICG(M[infinity]))))
          
# marginal posterior for the effect mu
plot_RA_fits(fits.actual=fits.actual, fits.bm=fits.bm.post, 
             type="post.mu", xlim=c(-1.5,2), ylim=c(0,3),
             legend=TRUE,
             legend.tau.prior=c("HN(1)", "HC(1)",
              expression(SCG(m[infinity])), "Jeffreys", 
              expression(SICG(M[infinity]))))        
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.