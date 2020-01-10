\name{cal_h_dist}
\alias{cal_h_dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calibration of the Hellinger distance
%%  ~~function to do ... ~~
}
\description{
    Returns the mean of a unit-variance normal distribution,
    such that the Hellinger distance between this distribution and
    the standard normal distribution equals the given value.
    Gives an intuitive interpretation of Hellinger distance values.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cal_h_dist(h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
  vector of Hellinger distances, consisting of real numbers in [0,1]
%%     ~~Describe \code{x} here~~
}
}
\details{
For a given Hellinger distance h, the function returns mu, such that
\deqn{H(N(\mu, 1), N(0, 1))=h,}{H(N(mu, 1), N(0, 1))=h,}
where H denotes the Hellinger distance.

See Roos et al. (2015), Sect. 2.2 for details.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of means.
}
\references{
Roos, M., Martins, T., Held, L., Rue, H. (2015). Sensitivity analysis for Bayesian hierarchical models. 
\emph{Bayesian Analysis} \bold{10(2)}, 321--349.
\url{https://projecteuclid.org/euclid.ba/1422884977}
}

\examples{
cal_h_dist(h=c(0.1, 0.5, 0.9))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.