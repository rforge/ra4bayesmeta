\name{cal_h_dist}
\alias{cal_h_dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calibration of the Hellinger distance
%%  ~~function to do ... ~~
}
\description{
    By default (if \code{output="shift"}), this function returns
   the mean of a unit-variance normal distribution,
    such that the Hellinger distance between this distribution and
    the standard normal distribution equals the given value. 
    Offers the option to return the area of overlap (if \code{output="ao"}) between these two unit-variance
    normal distributions instead.
    Gives an intuitive interpretation of Hellinger distance values.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cal_h_dist(h, output="shift")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
  vector of Hellinger distances, consisting of real numbers in [0,1]
%%     ~~Describe \code{x} here~~
}
\item{output}{
  either \code{"shift"} or \code{"ao"}. Specifies if the output should be given as
  the shift between two unit-varaince normal distributions or as the area of overlap (AO) between
  these unit-varaince normal distributions
%%     ~~Describe \code{x} here~~
}
}
\details{
For a given Hellinger distance h, there is a mean \eqn{\mu(h)}, such that
\deqn{H(N(\mu(h), 1), N(0, 1))=h,}
where H denotes the Hellinger distance.
See Roos et al. (2015), Sect. 2.2 for details.

If \code{output="shift"}, the function returns the shift \eqn{\mu(h)} between
the two unit-variance normal distributions.
If \code{output="ao"}, the function returns
the area of overlap between the \eqn{N(\mu(h), 1)} and \eqn{N(0, 1)} distributions.
This area of overlap is given by
\deqn{AO(\mu(h)) = \Phi(\mu(h)/2 ;\mu(h), 1) + 1 - \Phi(\mu(h)/2 ;0, 1),}
where \eqn{\Phi(. ;\mu, \sigma^2)} denotes the cumulative distribution function of
the normal distribution with mean \eqn{\mu} and variance \eqn{\sigma^2}.
See Ott et al. (2021, Section 3.5) for more information on this area of overlap calibration.
%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of means (if \code{output="shift"}) or areas of overlap (if \code{output="ao"}), respectively.
}
\references{
Roos, M., Martins, T., Held, L., Rue, H. (2015). Sensitivity analysis for Bayesian hierarchical models. 
\emph{Bayesian Analysis} \bold{10(2)}, 321--349.
\url{https://projecteuclid.org/euclid.ba/1422884977}

Ott, M., Plummer, M., Roos, M. How vague is vague? How informative is
informative? Reference analysis for
Bayesian meta-analysis. Manuscript revised for \emph{Statistics in Medicine}. 2021.
}

\examples{
# calibration in terms of shifts 
cal_h_dist(h=c(0.1, 0.5, 0.9))
# calibration in terms of areas of overlap 
cal_h_dist(h=c(0.1, 0.5, 0.9), output="ao")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.