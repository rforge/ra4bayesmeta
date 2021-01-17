\name{post_RA_fits}
\alias{post_RA_fits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior reference analysis based on bayesmeta fits
%%  ~~function to do ... ~~
}
\description{
Computes a table of Hellinger distances between marginal posterior distributions
for different parameters in the NNHM
between the actual model fits in \code{fits.actual} and
the benchmark fits in \code{fits.bm}.
All fits should be based on the same data set.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
post_RA_fits(fits.actual, fits.bm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{fits.actual}{
    a list of model fits of class bayesmeta, computed with
    the \code{bayesmeta} function in the package bayesmeta.
}
\item{fits.bm}{
    a list of model fits of class bayesmeta, computed with
    the \code{bayesmeta} function in the package bayesmeta. 
    To be used as benchmarks.
}
}
\details{
Two alternative suggestions for posterior benchmarks are provided
in Ott et al. (2021, Section 3.4) and its Supplementary Material (Section 2.5) and they 
can be computed using the functions \code{fit_models_RA} and \code{fit_models_RA_5bm},
respectively.
}
\value{
A matrix of Hellinger distance estimates between marginal posteriors
with \eqn{n.bm} columns and \eqn{n.act*(k+3)} rows,
where \eqn{n.bm}=length(\code{fits.bm}) is the number of benchmark fits specified,
\eqn{n.act}=length(\code{fits.actual}) the number of actual fits specified
 and
\eqn{k} the number of studies in the meta-analysis data set
(so that there are \eqn{k+3} parameters \eqn{\Psi \in \{ \mu, \tau, \theta_1, ..., \theta_k, \theta_{new}  \}} of potential interest in the NNHM).
 
The columns of the matrix give the following Hellinger distance estimates
between two marginal posteriors (for the parameter of interest \eqn{\Psi} varying with rows)
induced by the following two heterogeneity priors
(from left to right):
\item{H(po_{bm_1}, po_act) }{ first benchmark prior bm_1 inducing the fit \code{fits.bm[[1]]} and actual prior}
\item{H(po_{bm_2}, po_act) }{ second benchmark prior bm_2 inducing the fit \code{fits.bm[[2]]} and actual prior}
\item{... }{ ... }
\item{H(po_{bm_{n.bm}}, po_act) }{ last benchmark prior bm_{n.bm} inducing the fit \code{fits.bm[[n.bm]]} and actual prior}

The actual heterogenity prior and the parameter of interest \eqn{\Psi} vary
with the rows in the following order:
\item{mu, pri_act_1 }{ \eqn{\Psi=\mu} and first actual prior in \code{tau.prior}}
\item{mu, pri_act_2 }{ \eqn{\Psi=\mu} and second actual prior in \code{tau.prior}}
\item{... }{ ...}
\item{mu, pri_act_n }{ \eqn{\Psi=\mu} and nth actual prior in \code{tau.prior}}
\item{tau, pri_act_1 }{ \eqn{\Psi=\tau} and first actual prior in \code{tau.prior}}
\item{... }{ ...}
\item{tau, pri_act_n }{ \eqn{\Psi=\tau} and nth actual prior}
\item{theta_1, pri_act_1 }{ \eqn{\Psi=\theta_1} and first actual prior}
\item{... }{ ...}
\item{theta_k, pri_act_n }{ \eqn{\Psi=\theta_k} and nth actual prior}
\item{theta_new, pri_act_1 }{ \eqn{\Psi=\theta_{new}} and first actual prior}
\item{... }{ ...}
\item{theta_new, pri_act_n }{ \eqn{\Psi=\theta_{new}} and nth actual prior}
}
\references{
Ott, M., Plummer, M., Roos, M. How vague is vague? How informative is
informative? Reference analysis for
Bayesian meta-analysis. Manuscript revised for \emph{Statistics in Medicine}. 2021.

Ott, M., Plummer, M., Roos, M. Supplementary Material:
How vague is vague? How informative is informative? Reference analysis for 
Bayesian meta-analysis. Revised for \emph{Statistics in Medicine}. 2021.
}
\seealso{
\code{\link[bayesmeta]{bayesmeta}} in the package bayesmeta,
\code{\link{fit_models_RA}},
\code{\link{post_RA}}, 
\code{\link{pri_RA_fits}}

}
\examples{
# for aurigular acupuncture (AA) data set
data(aa)
# compute the model fits
# actual standard half-normal and half-Cauchy heterogeneity priors
fits <- fit_models_RA(df=aa, tau.prior=
                             list(function(t)dhalfnormal(t, scale=1),
                                  function(t)dhalfcauchy(t, scale=1)))
# benchmark fits under HN0 and J (Jeffreys) priors
fits.bm.post <- fits[c(1,2)]
fits.actual <- fits[c(3,4)]
post_RA_fits(fits.actual=fits.actual, fits.bm=fits.bm.post)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.