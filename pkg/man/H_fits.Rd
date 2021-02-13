\name{H_fits}
\alias{H_fits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hellinger distance between marignal posterior densities of two bayesmeta fits
%%  ~~function to do ... ~~
}
\description{
    Computes the Hellinger distance between the
    marginal posterior densities extracted from two bayesmeta fits,
    for the selected parameter of interest.
    Offers the option to select between integral-based and approximate moment-based computation of Hellinger distances.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
} 
\usage{
H_fits(fit1, fit2, parameter = "mu", individual = NA,
       method = "integral")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{fit1}{
  first model fit of class bayesmeta, computed with
  the \code{bayesmeta} function in the package bayesmeta.
%%     ~~Describe \code{x} here~~
}
\item{fit2}{
  second model fit of class bayesmeta, computed with
  the \code{bayesmeta} function in the package bayesmeta.
%%     ~~Describe \code{x} here~~
}
\item{parameter}{
  parameter for which the the marginal posterior densities should be compared. Possible values are \code{"mu", "tau", "theta_new"} and \code{"theta"} (\code{"theta"} refers to the random effects and is used in combination  with the individual argment).
%%     ~~Describe \code{x} here~~
}
\item{individual}{
  integer \eqn{<= k}, where \eqn{k} is the number of studies in the meta-analysis data set. Specifies which random effect parameter to consider (according to the order of the studies in the data set). Is used only if \code{parameter = "theta"}.
%%     ~~Describe \code{x} here~~
}
\item{method}{
    method for Hellinger distance computation. Either \code{"integral"} for integral-based computation or \code{"moment"} for approximate moment-based calculation using a normal approximation. Defaults to \code{"integral"}.
%%     ~~Describe \code{x} here~~
}
}
\details{
In the random effects model for meta-analysis, the parameter mu is the overall mean parameter, tau is the between-study (heterogeneity) standard deviation, theta_new is the predicted effect for a new study and theta_i is the random effect parameter for study i.

If integral-based computation (\code{method = "integral"}) of Hellinger distances is selected (the default), numerical integration is applied to obtain the Hellinger distance between the two marginal posterior densities (by using the function \code{H}). 
If moment-based computation (\code{method = "moment"}) is selected, the marginal densities are first approximated by normal densities with the same means and standard deviations and then the Hellinger distance between these normal densities can be obtained by an analytical formula (implemented in the function \code{H_normal}).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A real number in [0,1].
}
\section{Warning }{
If the integral-based method is used to compute Hellinger distances (\code{method = "integral"}),
numerical problems may occur in some cases, which may lead to implausible outputs.
Therefore, we generally recommend to double-check the results of the integral-based method using the moment-based method (\code{method = "moment"}) - especially if the former results are implausibe. If large differences between the two methods are observed, we recommend to rely on the moment-based method unless a normal approximation of the involved densities is inappropriate.
}
\seealso{
\code{\link{H}},
\code{\link{H_normal}}
}
\examples{
data(aa)
# fit under Jeffreys reference heterogeneity prior
fit1 <- bayesmeta(y = aa[ , "y"], sigma = aa[, "sigma"], 
                  mu.prior.mean = 0, mu.prior.sd = 4, 
                  tau.prior = "Jeffreys")
# fit under standard half-normal heterogeneity prior
fit2 <- bayesmeta(y = aa[ , "y"], sigma = aa[, "sigma"], 
                  mu.prior.mean = 0, mu.prior.sd = 4, 
                  tau.prior = function(t) dhalfnormal(t, scale = 1))
# Hellinger distance between the marginal posteriors for mu
H_fits(fit1, fit2)
# moment-based calculation
H_fits(fit1, fit2, method = "moment")
# Hellinger distance between the marginal posteriors for theta_2
# (random effect parameter)
H_fits(fit1, fit2, parameter = "theta", individual = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.