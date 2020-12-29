\name{post_RA_3bm}
\alias{post_RA_3bm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior reference analysis based on a data frame using 3 benchmarks
%%  ~~function to do ... ~~
}
\description{
Computes a table of Hellinger distances between marginal posterior distributions
for different parameters in the NNHM
induced by the actual heterogeneity priors specified in \code{tau.prior} and
posterior benchmarks proposed in the Supplementary Material for Ott et al. (2020).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
post_RA_3bm(df, tau.prior=list(function(x) dhalfnormal(x, scale=1)),
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
    list of prior specifications, which are \emph{either} functions returning the probability densities of the actual priors of interest for the heterogeneity parameter tau \emph{or} character strings specifying priors implemented in the \code{bayesmeta} function. See the documentation of the argument \code{tau.prior} of the \code{bayesmeta} function for details.
}
%\item{tau.prior}{
%    list of functions returning the probability densities of the actual priors of interest for the heterogeneity parameter tau.
%}
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
The three posterior benchmarks used are
introduced in the Supplementary Material of Ott et al. (2020, Sections 2.2.1 and 2.5, see also Section 3.4 in Ott at al. (2020) for Jeffreys reference prior),
where they are denoted by 
\eqn{po_{m_{inf}}(\Psi)}, \eqn{po_{J}(\Psi)} and \eqn{po_{M_{inf}}(\Psi)}.
Here, \eqn{\Psi \in \{ \mu, \tau, \theta_1, ..., \theta_k, \theta_{new}  \}}
denotes the parameter of interest in the NNHM,
where \eqn{\theta_{i}} is the random effect in the \eqn{i}th study and
\eqn{\theta_{new}} the predicted effect for a new study.

Note that Jeffreys' reference posterior \eqn{po_{J}} is proper if
there are at least two studies in the meta-analysis data set.
It is based on the improper Jeffreys' reference prior,
which is minimally informative given the data.

The default values for \code{mu.mean} and \code{mu.sd}
are suitable for effects mu on the log odds (ratio) scale.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with two elements:
The first element named "table" is a matrix containing the Hellinger distance estimates and the
second element called "par" is a named vector giving the parameter values
of the benchmark priors.
The vector "par" has the following three components:
\code{m_inf}, \code{M_inf} and \code{C}.

The matrix "table" contains the Hellinger distance estimates between marginal posteriors
and has \eqn{3} columns and \eqn{n*(k+3)} rows,
where \eqn{n}=length(\code{tau.prior}) is the number of actual heterogeneity priors specified
and \eqn{k} is the number of studies in the meta-analysis data set
(so that there are \eqn{k+3} parameters of interest).
 
The columns of the matrix give the following Hellinger distance estimates
between two marginal posteriors (for the parameter of interest \eqn{\Psi} varying with rows)
induced by the following two heterogeneity priors
(from left to right):
\item{H(po_{m_inf}, po_act) }{ benchmark prior SGC(m_inf) and actual prior}
\item{H(po_J, po_act) }{ Jeffreys' reference prior \eqn{\pi_J} and actual prior}
\item{H(po_{M_inf}, po_act) }{ benchmark prior SIGC(M_inf) and actual prior}
The actual heterogenity prior and the parameter of interest \eqn{\Psi} vary
with the rows in the following order:
\item{mu, pri_act_1 }{ \eqn{\Psi=\mu} and first actual prior in \code{tau.prior}}
\item{mu, pri_act_2 }{ \eqn{\Psi=\mu} and second actual prior in \code{tau.prior}}
\item{... }{ }
\item{mu, pri_act_n }{ \eqn{\Psi=\mu} and nth actual prior in \code{tau.prior}}
\item{tau, pri_act_1 }{ \eqn{\Psi=\tau} and first actual prior in \code{tau.prior}}
\item{... }{ }
\item{tau, pri_act_n }{ \eqn{\Psi=\tau} and nth actual prior}
\item{theta_1, pri_act_1 }{ \eqn{\Psi=\theta_1} and first actual prior}
\item{... }{ }
\item{theta_k, pri_act_n }{ \eqn{\Psi=\theta_k} and nth actual prior}
\item{theta_new, pri_act_1 }{ \eqn{\Psi=\theta_{new}} and first actual prior}
\item{... }{ }
\item{theta_new, pri_act_n }{ \eqn{\Psi=\theta_{new}} and nth actual prior}
}
\references{
Ott, M., Plummer, M., Roos, M. Supplementary Material:
How vague is vague? How informative is informative? Reference analysis for 
Bayesian meta-analysis. Revised for \emph{Statistics in Medicine}. 2020.

Ott, M., Plummer, M., Roos, M. How vague is vague? How informative is
informative? Reference analysis for
Bayesian meta-analysis. Manuscript revised for \emph{Statistics in Medicine}. 2020.
}
\section{Warning }{
A warning message will be displayed if one of the parameters \code{m_inf} or \code{M_inf}
has a value larger than 5*10^6, since this may lead to numerical problems
in the \code{bayesmeta} function
used for computation of the marginal posteriors.
}
\seealso{
\code{\link{post_RA_fits}}, 
\code{\link{pri_RA_5bm}}
}
\examples{
# for aurigular acupuncture (AA) data set with two
# actual half-normal heterogeneity priors
data(aa)
post_RA_3bm(df=aa, tau.prior=list(function(t)dhalfnormal(t, scale=0.5),
                              function(t)dhalfnormal(t, scale=1)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.