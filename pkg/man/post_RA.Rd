\name{post_RA}
\alias{post_RA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior reference analysis based on a data frame using 2 benchmarks
%%  ~~function to do ... ~~
}
\description{
% TODO: mention posterior estimates here
Computes a table of posterior estimates and informativeness values
for the marginal posterior distributions
of different parameters in the NNHM
induced by the actual heterogeneity priors specified in \code{tau.prior}.
Also provides the same estimates for the 
posterior benchmarks proposed in Ott et al. (2020).
%Hellinger distances between marginal posterior distributions
%for different parameters in the NNHM
%induced by the actual heterogeneity priors specified in \code{tau.prior} and
%posterior benchmarks proposed in Ott et al. (2020).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
post_RA(df, tau.prior = list(function(x) dhalfnormal(x, scale = 1)), 
        show.re = FALSE, estimate = "median", ci.method = "central", 
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
    list of prior specifications, which are \emph{either} functions returning the probability densities of the actual priors of interest for the heterogeneity parameter tau \emph{or} character strings specifying priors implemented in the \code{bayesmeta} function. See the documentation of the argument \code{tau.prior} of the \code{bayesmeta} function for details.
}
%\item{tau.prior}{
%    list of functions returning the probability densities of the actual priors of interest for the heterogeneity parameter tau.
%}
\item{show.re}{
    logical. If TRUE, the estimates for the random effect parameters \eqn{\theta_i} are also shown. Defaults to FALSE.
}
\item{estimate}{
    type of posterior point estimate. Either "median" for the posterior median or "mean" for the posterior mean. Defaults to "median".
}
\item{ci.method}{
    method for computation of 95\% credible intervals (CrIs). Either "central" for equal-tailed CrIs or "shortest" for shortest CrIs. Defaults to "central".
}
\item{scale.hn0}{
  scale parameter of the half-normal benchmark prior (usually a small number to 
  obtain a conservative benchmark which has most of its probability mass close to 0).
  Defaults to 0.002.
}
\item{mu.mean}{
  mean of the normal prior for the effect mu. Defaults to 0.
}
\item{mu.sd}{
  standard deviation of the normal prior for the effect mu. Defaults to 4.
}
}
\details{
The two posterior benchmarks used are
introduced in Ott et al. (2020, Section 3.4),
where they are denoted by 
\eqn{po_{J}(\Psi)} and \eqn{po_{HN0}(\Psi)}. 
Here, \eqn{\Psi \in \{ \mu, \tau, \theta_1, ..., \theta_k, \theta_{new}  \}}
denotes the parameter of interest in the NNHM,
where \eqn{\theta_{i}} is the random effect in the \eqn{i}th study and
\eqn{\theta_{new}} the predicted effect for a new study.
For the overall mean parameter \eqn{\mu},
we additionally consider the fixed-effects model benchmark
\eqn{po_{FE}(\mu)}.

Note that Jeffreys' reference posterior \eqn{po_{J}} is proper if
there are at least two studies in the meta-analysis data set.
It is based on the improper Jeffreys' reference prior,
which is minimally informative given the data.
The computation of the informativeness values is described in 
Ott et al. (2020, Section 3.6).

The absolute value of the signed informativeness quantifies how close the 
actual posterior and the reference posterior \eqn{po_{J}} are.
If the signed informativeness is negative, then the actual prior is conservative  with respect to Jeffreys reference prior J (that is puts more weight on smaller values for \eqn{\tau} than J).
If the signed informativeness is positive, then the actual prior is anti-conservative  with respect to Jeffreys reference prior J (that is puts more weight on larger values for \eqn{\tau} than J).

The default values for \code{mu.mean} and \code{mu.sd}
are suitable for effects mu on the log odds (ratio) scale.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix with 6 columns and a number of columns depending on the number
of actual heterogeneity priors specified and the parameters of interests in the NNHM
(if \code{show.re = FALSE}, then there are 3 parameters of interest (\eqn{\mu, \tau, \theta_{new}}) and the matrix has \eqn{3*(n+2)+1} rows, where
where \eqn{n}=length(\code{tau.prior}) is the number of actual heterogeneity priors specified;
if \code{show.re = TRUE}, then the matrix has \eqn{(k+3)*(n+2)+1} rows, where
\eqn{k} is the number of studies in the meta-analysis data set
(so that there are \eqn{k+3} parameters of interest).)

The row names specify the parameter in the NNHM for which the marginal posterior is considered,
followed by the heterogeneity prior used to compute that posterior.
HN0, J and FE denote the three benchmark priors introduced in Ott et al. (2020).
\code{pri_act_i} denotes the i-th prior in the \code{tau.prior} list.

The 6 columns provide the following estimates:
\item{\code{estimate} }{ posterior point estimate (median or mean)}
\item{\code{CrI_low} }{ lower limit of the 95\% credible interval}
\item{\code{CrI_up} }{ upper limit of the 95\% credible interval}
\item{\code{length_CrI} }{ length of the 95\% credible interval}
\item{\code{H(po_HN0, po_act)} }{ Hellinger distance between the posterior benchmark \eqn{po_{HN0}(\Psi)} and the marginal posterior induced by the heterogeneity prior listed on the left, for the parameter \eqn{\Psi} listed on the left}
\item{\code{signed_inf} }{ Signed informativeness, i.e. \eqn{sign(H(po_{HN0}, po_{act}) - H(po_{HN0}, po_{J})) H(po_{act}, po_{J})}, for the marginal posterior induced by the heterogeneity prior listed on the left, for the parameter \eqn{\Psi} listed on the left}
}
\references{
Ott, M., Plummer, M., Roos, M. How vague is vague? How informative is
informative? Reference analysis for
Bayesian meta-analysis. Manuscript revised for \emph{Statistics in Medicine}. 2020.
}
\seealso{
\code{\link{post_RA_fits}}, 
\code{\link{post_RA_3bm}}
}
\examples{
# for aurigular acupuncture (AA) data set with two
# actual half-normal heterogeneity priors
data(aa)
post_RA(df=aa, tau.prior=list(function(t)dhalfnormal(t, scale=0.5),
                              function(t)dhalfnormal(t, scale=1)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.