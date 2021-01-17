\name{fit_models_RA}
\alias{fit_models_RA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model fitting for reference analysis using 2 benchmarks:
Posterior inference for benchmark and actual heterogeneity priors
%%  ~~function to do ... ~~
}
\description{
Computes the posterior distribution of the
parameters in a random-effects meta-analysis (expressed
as a normal-normal hierarchical model) for 
two benchmark heterogeneity priors and
the actual heterogeneity prior(s) specified.
Applies the function \code{bayesmeta} from the package bayesmeta. 
%Also provides functions to evaluate marginal posterior distributions.
%for computations.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fit_models_RA(df, tau.prior = list(), scale.hn0 = 1/500,
              mu.mean = 0, mu.sd = 4, interval.type = "central")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{df}{
    data frame with one column "y" containing the (transformed) effect estimates for the individual studies and
    one column "sigma" containing the standard errors of these estimates.
%%     ~~Describe \code{x} here~~
}
\item{tau.prior}{
    list of prior specifications, which are \emph{either} functions returning the probability densities of the actual priors of interest for the heterogeneity parameter tau \emph{or} character strings specifying priors implemented in the \code{bayesmeta} function. See the documentation of the argument \code{tau.prior} of the \code{bayesmeta} function for details.
}
\item{scale.hn0}{
  scale parameter of the half-normal benchmark prior(usually small, so that the benchmark is conservative).
}
\item{mu.mean}{
  mean of the normal prior for the effect mu.
}
\item{mu.sd}{
  standard deviation of the normal prior for the effect mu.
}
\item{interval.type}{
  the type of (credible, prediction, shrinkage) interval to be returned by default.
  Either "central" for equal-tailed intervals or "shortest" for shortest intervals. Defaults to "central". See also the corresponding argument in the \code{bayesmeta} function.
}
}
\details{
The two heterogeneity benchmark priors used are
introduced in Ott et al. (2021, Section 3.4),
where they are denoted by 
HN0 and J.
Note that "J" refers to Jeffreys reference prior,
which is improper, but leads to a proper posterior if
there are at least two studies in the meta-analysis data set.
HN0 is a half-normal prior with scale parameter \code{scale.hn0}.

Decreasing the scale parameter \code{scale.hn0} of the half-normal benchmarks
leads to a more conservative (i.e. its mass is more concentrated near 0) HN0 benchmark prior.
However, \code{scale.hn0} cannot be chosen arbitrarily small
since too small values lead to numerical problems in the \code{bayesmeta} function
used to fit the models.
To verify how conservative the HN0 benchmark is,
one can compare the marginal posterior for the overall mean parameter mu
with the corresponding posterior for the fixed effects model,
e.g. by using the function \code{plot_RA}.
The better the match between these two marginal posteriors,
the more conservative the HN0 benchmark is.

The default values for \code{mu.mean} and \code{mu.sd}
are suitable for effects mu on the log odds (ratio) scale.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list containing the model fits, namely 
 a list of lists of class \code{bayesmeta}. 
This list has length \eqn{2} + length(\code{tau.prior})
and contains one element for each heterogeneity prior considered
(2 benchmark priors and the actual priors specified),
in the following order:
\item{fit.hn0 }{ for the half-normal HN0 benchmark prior with scale parameter \code{scale.hn0}}
\item{fit.j }{ for Jeffreys' (improper) reference prior}
\item{fit.actual_1 }{ for the first prior in the list \code{tau.prior} (if specified)}
\item{fit.actual_2 }{ for the second prior in the list \code{tau.prior} (if specified)}
\item{...}{ ...}
\item{fit.actual_n }{ for the nth prior in the list \code{tau.prior} (if specified)}
See \code{\link[bayesmeta]{bayesmeta}} in the package bayesmeta
for information on the structure of the lists of class \code{bayesmeta}.
}
\references{
Ott, M., Plummer, M., Roos, M. How vague is vague? How informative is
informative? Reference analysis for
Bayesian meta-analysis. Manuscript revised for \emph{Statistics in Medicine}. 2021.
}
\section{Warning }{
If \code{fit_models_RA} ends with an error or warning,
we recommend to increase the value of the parameter \code{scale.hn0}
for the conservative benchmark prior.

}
%\section{Remark }{
%}
\seealso{
\code{\link[bayesmeta]{bayesmeta}} in the package bayesmeta,
\code{\link{plot_RA}},
\code{\link{fit_models_RA_5bm}}
}
\examples{
# for aurigular acupuncture (AA) data set with one
# actual half-normal and the "DuMouchel" heterogeneity prior
data(aa)
fit_models_RA(df=aa, tau.prior=list(function(t)dhalfnormal(t, scale=1),
                                    "DuMouchel"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}