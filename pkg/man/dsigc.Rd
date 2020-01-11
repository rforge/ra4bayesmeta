\name{dsigc}
\alias{dsigc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Density function of the square-root inverse generalized conventional (SIGC) benchmark prior
%%  ~~function to do ... ~~
}
\description{
 Density function of the SIGC distribution described in Ott et al. (2020).   
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dsigc(x, M, C)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
   vector of quantiles.
%%     ~~Describe \code{x} here~~
}
\item{M}{
   real number in \eqn{(1,\infty)}.
%%     ~~Describe \code{x} here~~
}
\item{C}{
   non-negative real number. 
%%     ~~Describe \code{x} here~~
}
}
\details{
The density function with domain \eqn{[0, \infty)} is given by
\deqn{
   \pi(x) = 4(M-1)Cx^{-5}(1+Cx^{-4})^{-M}
}
for \eqn{x >= 0}.
This density is obtained
if the density function for
variance components given in equation (2.15) in Berger & Deely (1988)
is assigned to the precision (i.e. the inverse of the variance) and
then transformed to the standard deviation scale.
See Ott et al. (2020), Section 2.3.2 for more information.

For meta-analsis data sets, Ott et al. (2020) choose
\eqn{C=\sigma_{ref}^{-2}},
where \eqn{\sigma_{ref}} is the reference standard deviation (see function \code{\link{sigma_ref}}) of the
data set,
which is defined as the geometric mean of the standard deviations
of the individual studies.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Value of the density function at locations x, where \eqn{x >= 0}. Vector of non-negative real numbers.
}
\references{
Berger, J. O., Deely, J. (1988). A Bayesian approach to ranking and selection of
related means with alternatives to analysis-of-variance methodology. \emph{Journal of the
American Statistical Association} \bold{83(402)}, 364--373.

Ott, M., Plummer, M., Roos, M. How vague is vague? How informative is
informative? Reference analysis for
Bayesian meta-analysis. Manuscript submitted to \emph{Statistics in Medicine}. 2020.
}
\seealso{
\code{\link{dsgc}}
}
\examples{
dsigc(x=c(0.1,0.5,1), M=1.2, C=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keword{distribution}