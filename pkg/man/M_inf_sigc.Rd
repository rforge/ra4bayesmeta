\name{M_inf_sigc}
\alias{M_inf_sigc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimization function for the SIGC(M) prior: Adjust the prior to a target
relative latent model complexity (RLMC)
%%  ~~function to do ... ~~
}
\description{
Computes the parameter value \eqn{M=M_{inf}} of the SGC(\eqn{M}) prior,
such that the relative latent model complexity (RLMC) with respect
to the reference threshold is approximately \code{rlmc}.
The reference threshold is chosen as the (1-\code{alpha})-quantile of the 
SGC(\eqn{M_{inf}}) prior.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
M_inf_sigc(rlmc, df, alpha=0.5, truncation=5*10^6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{rlmc}{
    target RLMC value. Real number in \eqn{(0,1)}.
%%     ~~Describe \code{x} here~~
}
\item{df}{
  data frame with one column "sigma" containing the standard errors of the estimates for the individual studies.
}
\item{alpha}{
   determines the (1-\code{alpha})-quantile of the SIGC(\eqn{M}) prior,
   which is used as reference threshold. Defaults to 0.5 (i.e. the median).
%%     ~~Describe \code{x} here~~
}
\item{truncation}{
   upper bound for the parameter value \eqn{M}.
   Defaults to the empirically determined value 5*10^6.
%%     ~~Describe \code{x} here~~
}
}
\details{
See the Supplementary Material of Ott et al. (2020), Section 2.6,
for the formulas and explanations.
Note that the parameter value \eqn{M_{inf}} does depend
on the data set considered.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Parameter value \eqn{M=M_{inf}} of the SGC(\eqn{M}) prior. Real number > 1.
}
\references{
Ott, M., Plummer, M., Roos, M. How vague is vague? How informative is
informative? Reference analysis for
Bayesian meta-analysis. Manuscript submitted to \emph{Statistics in Medicine}. 2020.
}
\section{Warning }{
Occasionally, the formula for \eqn{M_{inf}} given in the Supplementary Material of 
Ott et al. (2020, Section 2.8)
yields values larger than 5*10^6. This can cause numerical problems in the
\code{bayesmeta} function.
Therefore, we truncate the parameter value at the
empirically determined threshold 5*10^6 by default.
}
\seealso{
\code{\link{m_inf_sgc}}
}
\examples{
# extreme RLMC target value close to 1 used in Ott et al. (2020)
# for the aurigular acupuncture (AA) data set
data(aa)
M_inf_sigc(df=aa, rlmc=0.9999)
# for the respiratory tract infections (RTI) data set
data(rti)
M_inf_sigc(df=rti, rlmc=0.9999)

# 75\% quantile instead of the median as ref. threshold
M_inf_sigc(df=rti, rlmc=0.9999, alpha=0.25)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.