\name{post_mu_fe}
\alias{post_mu_fe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normal posterior for the overall mean parameter in the fixed effects model
%%  ~~function to do ... ~~
}
\description{
This function computes the mean and the standard deviation of the
normal posterior distribution for the overall mean parameter mu in the fixed effects model.
%assuming a normal prior for mu.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
post_me_fe(df, mu.mean = 0, mu.sd = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
  data frame with one column "y" containing the (transformed) effect estimates for the individual studies and
    one column "sigma" containing the standard errors of these estimates. 
}
\item{mu.mean}{
   mean of the normal prior for the overall mean parameter mu. Defaults to 0.
}
\item{mu.sd}{
  standard deviation of the normal prior for the overall mean parameter mu. Defaults to 4.
}
}
\details{
For the fixed effects model, the normal posterior for the overall mean parameter mu can be computed analytically (under the common assumption of a normal prior on mu)
since this is a onjugate Bayesian normal-normal model.
See for example Ott et al. (2020), Equation (2) for the formula
(\code{mu.mean} corresponds to \eqn{nu} and \code{mu.sd} corresponds to \eqn{gamma}).
%This is a conjugate Bayesian normal-normal model.

The default values for \code{mu.mean} and \code{mu.sd}
are suitable for effects mu on the log odds (ratio) scale (Roever, 2020).
%For outcomes on the log odds ratio scale, Roever (2020) suggests to
%use the normal prior with nu=0 and gamma=4 if a "vague" prior for mu is desired.
%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with two elements: the first element "mean" and the second element "sd", which refer to the mean and the standard deviation of the normal posterior of mu.
}
\references{
Ott, M., Plummer, M., Roos, M. How vague is vague? How informative is
informative? Reference analysis for
Bayesian meta-analysis. Manuscript revised for \emph{Statistics in Medicine}. 2020.

Roever C. Bayesian random-effects meta-analysis using the bayesmeta R package (2020). \emph{Journal of Statistical Software} \bold{93(6)}, 1--51.
}

\examples{
# load the aurigular acupuncture (AA) data set 
data(aa)
# normal prior for log odds ratios suggested by Roever (2020)
post_mu_fe(df=aa, mu.mean=0, mu.sd=4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.