\name{plot_RA}
\alias{plot_RA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reference analysis plot based on a data frame using 2 benchmarks:
Plot heterogeneity benchmark priors and the
corresponding marginal posteriors
%%  ~~function to do ... ~~
}
\description{
Plots the proper half-normal benchmark heterogeneity prior proposed in Ott et al. (2021)
and the corresponding marginal benchmark posteriors for different parameters in the 
NNHM. Also displays the marginal Jeffreys reference posterior,
which is induced by an improper reference prior.
Displays the specified actual heterogeneity priors in \code{tau.prior}
and the corresponding marginal posteriors as well.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_RA(df, tau.prior = list(),
        type = "pri.tau", improper.prior = NULL, show.sigma.i = FALSE,
        xlim, ylim, legend = FALSE, 
        pos.legend = "topright", legend.tau.prior = c(), 
        xlab = NULL, bty = "o", 
        scale.hn0 = 1/500, mu.mean = 0, mu.sd = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{df}{
    data frame with one column "y" containing the (transformed) effect estimates for the individual studies and
    one column "sigma" containing the standard errors of these estimates.
%%     ~~Describe \code{x} here~~
}
\item{tau.prior}{
    list (maximum length 5) of prior specifications, which are \emph{either} functions returning the probability densities of the actual priors of interest for the heterogeneity parameter tau \emph{or} character strings specifying priors implemented in the \code{bayesmeta} function. See the documentation of the argument \code{tau.prior} of the \code{bayesmeta} function for details.
}
\item{type}{
    specifies if heterogeneity priors or marginal posterior densities for
    a given parameter should be plotted. 
    Options are \code{"pri.tau"} for heterogeneity priors, \code{"post.mu"} for the marg. posterior of the main effect mu, \code{"post.tau"} for the marg. posterior of the heterogeneity tau and
    \code{"post.theta.new"} for the marg. posterior of the predicted effect theta_new in a new study.
}
\item{improper.prior}{
    vector of integers specifying the positions of the improper priors in the list \code{tau.prior}
    or NULL (the default). Determines which prior densities should not be plotted. Applies only to \code{type = "pri.tau"}.
}
\item{show.sigma.i}{
    logical. Specifies if the standard errors of the estimates for the individual studies (the entries of the column "sigma" in the data frame \code{df}) should be added to the density plot as dots. Defaults to FALSE. Applies only to \code{type = "pri.tau"}.
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
   a vector of character strings or expressions of the same length
   as the list \code{tau.prior}. Determines the legend entries corresponding to
   the actual heterogeneity priors specified in \code{tau.prior},
   in the same order as the priors in that list.
}
\item{xlab}{
    x-axis label. If unspecified or set to NULL (the default), a default label specifying
    the paramter of interest is used.
}
\item{bty}{
    the type of box to be drawn around the legend. The allowed values are "o" (the default) and "n".
}
\item{scale.hn0}{
  scale parameter of the half-normal benchmark prior (usually a small number to 
  obtain a conservative benchmark which has most of its probability mass close to 0).
  Defaults to 0.002.
}
\item{mu.mean}{
  mean of the normal prior for the effect mu.
}
\item{mu.sd}{
  standard deviation of the normal prior for the effect mu.
}
}
\details{
The two heterogeneity benchmark priors used are
introduced in Ott et al. (2021, Section 3.4)
and they are denoted by 
HN0 and J.
Note that the latter Jeffreys reference prior J
is improper, but leads to a proper posterior if
there are at least two studies in the meta-analysis data set.

The default values for \code{mu.mean} and \code{mu.sd}
are suitable for effects mu on the log odds (ratio) scale.

If you wish to adapt the colors, line types or line widths of the plotted densities,
we suggest to compute the models fits first by using the function \code{fit_models_RA}
and then to apply the more flexible function \code{plot_RA_fits}
operating on model fits of class \code{bayesmeta}.
%%  ~~ If necessary, more details than the description above ~~
}
\references{
Ott, M., Plummer, M., Roos, M. How vague is vague? How informative is
informative? Reference analysis for
Bayesian meta-analysis. Manuscript revised for \emph{Statistics in Medicine}. 2021.
}
%\section{Warning }{
%This function may take several minutes to run.
%}
\seealso{
\code{\link{plot_RA_fits}},
\code{\link{plot_RA_5bm}},
\code{\link{fit_models_RA}}
}
\examples{
# for aurigular acupuncture (AA) data set with two
# actual half-normal and half-Cauchy heterogeneity priors
data(aa)
# prior densities
plot_RA(df=aa, tau.prior=list(function(t)dhalfnormal(t, scale=1),
                              function(t)dhalfcauchy(t, scale=1)),
        type="pri.tau", xlim=c(0,2), ylim=c(0,3), legend=TRUE,
        legend.tau.prior=c("HN(1)", "HC(1)"))
          
# marginal posterior for the effect mu
plot_RA(df=aa, tau.prior=list(function(t)dhalfnormal(t, scale=1),
                                  function(t)dhalfcauchy(t, scale=1)),
        type="post.mu", xlim=c(-1.5,2), ylim=c(0,3), legend=TRUE,
        legend.tau.prior=c("HN(1)", "HC(1)"))         
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.