\name{fit_models_RA}
\alias{fit_models_RA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model fitting for reference analysis:
Posterior inference for benchmark and actual heterogeneity priors
%%  ~~function to do ... ~~
}
\description{
Computes the posterior distribution of the
parameters in a random-effects meta-analysis (expressed
as a normal-normal hierarchical model) for 
each of five benchmark heterogeneity priors and
the actual heterogeneity prior(s) specified.
Applies the function \code{bayesmeta}() from the package bayesmeta. 
%Also provides functions to evaluate marginal posterior distributions.
%for computations.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fit_models_RA(df, tau.prior=list(),
              m_J=NA, M_J=NA, upper.J=3, digits.J=2,
              m_inf=NA, M_inf=NA, rlmc0=0.0001, rlmc1=0.9999,
              mu.mean=0, mu.sd=4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{df}{
    data frame with one column "y" containing the (transformed) effect estimates for the individual studies and
    one column "sigma" containing the standard errors of these estimates.
%%     ~~Describe \code{x} here~~
}
\item{tau.prior}{
    list of prior specifications, which are \emph{either} functions returning the probability densities of the actual priors of interest for the heterogeneity parameter tau \emph{or} character strings specifying priors implemented in the \code{bayesmeta}() function. See the documentation of the argument \code{tau.prior} of the \code{bayesmeta}() function for details.
}
\item{m_J}{
    parameter value \eqn{m=m_J} of the SGC(\eqn{m}) prior,
    which induces a marignal posterior for the heterogeneity standard deviation tau   close to Jeffreys' reference posterior (wrt the Hellinger distance).
    If set to \code{NA} (the default), this parameter is computed 
    using the function \code{m_j_sgc} and the arguments
    \code{df}, \code{upper.J}, \code{digits.J}, \code{mu.mean} and \code{mu.sd}.
    %and some other parameters.
}
\item{M_J}{
    parameter value \eqn{M=M_J} of the SIGC(\eqn{M}) prior,
    which induces a marignal posterior for the heterogeneity standard deviation tau   close to Jeffreys' reference posterior (wrt the Hellinger distance).
    If set to \code{NA} (the default), this parameter is computed 
    using the function \code{M_j_sigc} and the arguments
    \code{df}, \code{upper.J}, \code{digits.J}, \code{mu.mean} and \code{mu.sd}.
    % some other parameters.
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
\item{mu.mean}{
  mean of the normal prior for the effect mu.
}
\item{mu.sd}{
  standard deviation of the normal prior for the effect mu.
}
}
\details{
The five heterogeneity benchmark priors used are
introduced in Ott et al. (2020, Sections 2.3.2, 2.5 and 2.7),
where they are denoted by 
SGC(\code{m_inf}), SIGC(\code{M_J}), SGC(\code{m_J}), SIGC(\code{M_inf})
and Jeffreys.
Note that "Jeffreys" refers to Jeffreys' reference prior,
which is improper, but leads to a proper posterior if
there are at least two studies in the meta-analysis data set.

The default values for \code{mu.mean} and \code{mu.sd}
are suitable for effects mu on the log odds (ratio) scale.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with two elements:
The first element named "fits" is a list containing the model fits and the
second element called "par" is a named vector giving the parameter values
of the benchmark priors.
The vector "par" has the following five components:
\code{m_inf}, \code{M_J}, \code{m_J}, \code{M_inf} and \code{C}.
 
The element "fits" is a list of lists of class \code{bayesmeta}. 
This list has length \eqn{5} + length(\code{tau.prior})
and contains one element for each heterogeneity prior considered
(5 benchmark priors and the actual priors specified),
in the following order:
\item{fit.SGC.m_inf }{ for SGC(\code{m_inf}) benchmark prior}
\item{fit.SIGC.M_J }{ for SIGC(\code{M_J}) benchmark prior}
\item{fit.SGC.m_J }{ for SGC(\code{m_J}) benchmark prior}
\item{fit.SIGC.M_inf }{ for SIGC(\code{M_inf}) benchmark prior}
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
Bayesian meta-analysis. Manuscript submitted to \emph{Statistics in Medicine}. 2020.
}
\section{Warning }{
This function may take several minutes to run if the parameter \code{m_J} and/or \code{M_J}
is not specified,
especially if the desired precision is \code{digits.J=2} or even \code{digits.J=3}.
}
\section{Remark }{
If \code{fit_models_RA} ends with an error or warning,
we recommend to check the values of the parameters \code{M_J}, \code{m_J} and \code{M_inf}
for the benchmark priors.
The choice of these parameters can be studied by using the functions 
\code{M_j_sigc}, \code{m_j_sgc} and \code{M_inf_sigc}.
}
\seealso{
\code{\link[bayesmeta]{bayesmeta}} in the package bayesmeta,
\code{\link{m_j_sgc}}, \code{\link{M_j_sigc}},
\code{\link{m_inf_sgc}}, \code{\link{M_inf_sigc}}
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