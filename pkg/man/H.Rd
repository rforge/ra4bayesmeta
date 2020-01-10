\name{H}
\alias{H}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hellinger distance between two probability densities
%%  ~~function to do ... ~~
}
\description{
    Computes the Hellinger distance between the
    two probability density functions specfied.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
} 
\usage{
H(dens1, dens2, lower=-Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dens1}{
  function returning the first probability density.
%%     ~~Describe \code{x} here~~
}
\item{dens2}{
   function returning the second probability density.
%%     ~~Describe \code{x} here~~
}
\item{lower}{
  lower end of the common domain of two density functions. Usually equal to \eqn{-\infty} or 0.
%%     ~~Describe \code{x} here~~
}
}
\details{
The Hellinger distance is defined in Le Cam (1986) and also given in
Roos et al. (2015, Sect. 2.1).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A real number in [0,1].
}
\references{
Le Cam, L. (1986). Asymptotic Methods in Statistical Decision Theory. Springer.

Roos, M., Martins, T., Held, L., Rue, H. (2015). Sensitivity analysis for Bayesian hierarchical models. 
\emph{Bayesian Analysis} \bold{10(2)}, 321--349.
\url{https://projecteuclid.org/euclid.ba/1422884977}
%Ott, M., Plummer, M., Roos, M. How vague is vague? How informative is
%informative? Reference analysis for
%Bayesian meta-analysis. Manuscript in preparation. 2019.
}
\examples{
H(dens1=function(x) dhalfnormal(x, scale=0.5), 
  dens2=function(x) dhalfnormal(x, scale=1), lower=0)
H(dens1=function(x) dhalfnormal(x, scale=0.5), 
  dens2=function(x) dhalfcauchy(x, scale=1), lower=0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.