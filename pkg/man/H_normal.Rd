\name{H_normal}
\alias{H_normal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Approximate moment-based Hellinger distance computation between two probability densities
%%  ~~function to do ... ~~
}
\description{
    Computes the Hellinger distance between the
    two probability densities with the given means and standard deviations
    using a normal approximation and moment-matching.
    Applies a moment-based analytical formula, which is exact for normal densities.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
} 
\usage{
H_normal(mean1, sd1, mean2, sd2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{mean1}{
  mean of the first density.
%%     ~~Describe \code{x} here~~
}
\item{sd1}{
  standard deviation of the first density.
%%     ~~Describe \code{x} here~~
}
\item{mean2}{
  mean of the second density.
%%     ~~Describe \code{x} here~~
}
\item{sd2}{
  standard deviation of the second density.
%%     ~~Describe \code{x} here~~
}
}
\details{
The Hellinger distance between two normal densities can be computed analytically and it only depends on the means and the standard deviations of these two densities.
Any density can be approximated to the first order by a normal distribution.
Thus, this moment-based formula can also be applied to non-normal densities to obtain an approximate estimate of the Hellinger distance.
This approximation relies on moment-matching, i.e. the non-normal densities are approximated by normal densities with the same means and standard deviations.
See Roos et al. (2020), Section 3.3, for the formula and for additional information.

Note that the Hellinger distance is invariant with respect to one-to-one transformations of both densities. For the between-study standard deviation \eqn{\tau}, the density of \eqn{log(\tau)} can be approximated by a normal density.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A real number in [0,1].
}
\references{
Roos, M., Hunanyan, S., Bakka, H., Rue, H. Sensitivity and identification quantification by a relative latent model complexity perturbation in Bayesian meta-analysis. Manuscript revised for \emph{Biometrical Journal}. 2020.
}
\seealso{
\code{\link{H}}
}
\examples{
H_normal(mean1 = 0, sd1 = 1, mean2 = 1, sd2 = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.