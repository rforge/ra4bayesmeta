\name{m_inf_sgc}
\alias{m_inf_sgc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimization function for the SGC(m) prior: Adjust the prior to a target
relative latent model complexity (RLMC)
%%  ~~function to do ... ~~
}
\description{
  Computes the parameter value \eqn{m=m_{inf}} of the SGC(\eqn{m}) prior,
such that the relative latent model complexity (RLMC) with respect
to the reference threshold is approximately \code{rlmc}.
The reference threshold is chosen as the (1-\code{alpha})-quantile of the 
SGC(\eqn{m_{inf}}) prior.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
m_inf_sgc(rlmc, alpha=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{rlmc}{
    target RLMC value. Real number in \eqn{(0,1)}.
%%     ~~Describe \code{x} here~~
}
\item{alpha}{
   determines the (1-\code{alpha})-quantile of the SGC(\eqn{m}) prior,
   which is used as reference threshold. Defaults to 0.5 (i.e. the median).
%%     ~~Describe \code{x} here~~
}
}
\details{
See the Supplementary Material of Ott et al. (2020), Section 2.6,
for the formulas and explanations.
Note that the parameter value \eqn{m_{inf}} does not depend
on the data set considered.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Parameter value \eqn{m=m_{inf}} of the SGC(\eqn{m}) prior. Real number > 1.
}
\references{
Ott, M., Plummer, M., Roos, M. How vague is vague? How informative is
informative? Reference analysis for
Bayesian meta-analysis. Manuscript submitted to \emph{Statistics in Medicine}. 2020.
}
\seealso{
\code{\link{M_inf_sigc}}
}
\examples{
# extreme RLMC target value close to 0 used in Ott et al. (2020)
m_inf_sgc(rlmc=0.0001)

# 25\% quantile instead of the median as ref. threshold
m_inf_sgc(rlmc=0.0001, alpha=0.75)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.