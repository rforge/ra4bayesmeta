\name{plot_RA}
\alias{plot_RA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reference analysis plot based on a data frame:
Plot heterogeneity benchmark priors and the
corresponding marginal posteriors
%%  ~~function to do ... ~~
}
\description{
Plots the four proper benchmark heterogeneity priors proposed in Ott et al. (2020)
and the corresponding marginal benchmark posteriors for different parameters in the 
NNHM. Also displays the marginal Jeffreys' reference posterior,
which is induced by an improper reference prior.
Displays the specified actual heterogenitiy priors in \code{tau.prior}
and the corresponding marginal posteriors as well.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_RA(df, tau.prior=list(), type="pri.tau", xlim, ylim,
        m_J=NA, M_J=NA, upper.J=3, digits.J=2,
        m_inf=NA, M_inf=NA, rlmc0=0.0001, rlmc1=0.9999,
        legend=FALSE, pos.legend="topright", legend.tau.prior=c(), 
        bty="o", mu.mean=0, mu.sd=4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{df}{
    data frame with one column "y" containing the (transformed) effect estimates for the individual studies and
    one column "sigma" containing the standard errors of these estimates.
%%     ~~Describe \code{x} here~~
}
\item{tau.prior}{
    list (maximum length 5) of prior specifications, which are \emph{either} functions returning the probability densities of the actual priors of interest for the heterogeneity parameter tau \emph{or} character strings specifying priors implemented in the \code{bayesmeta}() function. See the documentation of the argument \code{tau.prior} of the \code{bayesmeta}() function for details.
}
%\item{tau.prior}{
%    list of functions returning the probability densities of the actual priors of interest for the heterogeneity parameter tau 
%    (max. 5 functions).
%}
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
\item{m_J}{
    parameter value \eqn{m=m_J} of the SGC(\eqn{m}) prior,
    which induces a marignal posterior for the heterogeneity standard deviation tau   close to Jeffreys' reference posterior (wrt the Hellinger distance).
    If set to \code{NA} (the default), this parameter is computed 
    using the function \code{m_j_sgc} and some other parameters.
}
\item{M_J}{
    parameter value \eqn{M=M_J} of the SIGC(\eqn{M}) prior,
    which induces a marignal posterior for the heterogeneity standard deviation tau   close to Jeffreys' reference posterior (wrt the Hellinger distance).
    If set to \code{NA} (the default), this parameter is computed 
    using the function \code{M_j_sigc} and some other parameters.
}
\item{upper.J}{
    upper bound for the parameters \code{M_J} and \code{m_J}. Real number in \eqn{(1,\infty)}. Is required only if \code{M_J=NA} or \code{m_J=NA}.
}
\item{digits.J}{
    specifies the desired precision of the parameter values \eqn{M_J} and \eqn{m_J}, i.e. to how many digits these two values
   should be determined. Possible values are 1,2,3. Defaults to 2.
   Is required only if \code{M_J=NA} or \code{m_J=NA}.
}
\item{m_inf}{
    parameter value \eqn{m=m_{inf}} of the SGC(\eqn{m}) prior,
such that the median relative latent model complexity (RLMC) is close to 0.
  If set to \code{NA} (the default), this parameter is computed 
    using the function \code{m_inf_sgc}, such that the median RLMC is
    approximately equal to \code{rlmc0}.
}
\item{M_inf}{
    parameter value \eqn{M=M_{inf}} of the SIGC(\eqn{M}) prior,
such that the median relative latent model complexity (RLMC) is close to 1.
  If set to \code{NA} (the default), this parameter is computed 
    using the function \code{M_inf_sigc}, such that the median RLMC is
    approximately equal to \code{rlmc1}.
}
\item{rlmc0}{
    RLMC target value for the SGC(\eqn{m_{inf}}) benchmark prior (typically close to 0).
    Is required only if \code{m_inf=NA}.
}
\item{rlmc1}{
    RLMC target value for the SIGC(\eqn{M_{inf}}) benchmark prior (typically close to 1).
    Is required only if \code{M_inf=NA}.
}
\item{legend}{
    logical. Specifies if a legend should be added to the plot. Defaults to \code{FALSE}.
}
\item{pos.legend}{
    a character string specifing the position of the legend in the plot. 
    Options are all the keywords accepted by
    the \code{legend}() function in the graphics package, i.e. 
    \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, \code{"topright"}, \code{"right"} and \code{"center"}.
}
\item{legend.tau.prior}{
   a vector of character strings or expressions of the same length
   as the list \code{tau.prior}. Determines the legend entries corresponding to
   the actual heterogeneity priors specified in \code{tau.prior},
   in the same order as the priors in that list.
}
\item{bty}{
    the type of box to be drawn around the legend. The allowed values are "o" (the default) and "n".
}
\item{mu.mean}{
  mean of the normal prior for the effect mu.
}
\item{mu.sd}{
  standard deviation of the normal prior for the effect mu.
}
}
\details{
The five heterogeneity benchmark priors used are
introduced in Ott et al. (2020, Sections 2.3.2, 2.5 and 2.7)
and they are denoted by 
SGC(\code{m_inf}), SIGC(\code{M_J}), SGC(\code{m_J}), SIGC(\code{M_inf})
and Jeffreys.
Note that the latter Jeffreys' reference prior
is improper, but leads to a proper posterior if
there are at least two studies in the meta-analysis data set.

The default values for \code{mu.mean} and \code{mu.sd}
are suitable for effects mu on the log odds (ratio) scale.
%%  ~~ If necessary, more details than the description above ~~
}
\references{
Ott, M., Plummer, M., Roos, M. How vague is vague? How informative is
informative? Reference analysis for
Bayesian meta-analysis. Manuscript submitted to \emph{Statistics in Medicine}. 2020.
}
\section{Warning }{
This function may take several minutes to run if the parameter \code{m_J} and/or \code{M_J}
is not specified,
especially if the desired precision is \code{digits.J=2} or even \code{digits.J=3}.
}
\seealso{
\code{\link{fit_models_RA}}, \code{\link{plot_RA_fits}}
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