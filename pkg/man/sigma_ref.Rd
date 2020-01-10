\name{sigma_ref}
\alias{sigma_ref}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reference standard deviation
%%  ~~function to do ... ~~
}
\description{
Computes the reference standard deviation of the
given meta-analysis data set.
The reference standard deviation is defined as the geometric mean
of the standard deviations of the individual studies
(Sorbye & Rue 2014, equation (7)).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sigma_ref(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
  data frame with one column "sigma" containing the standard errors of the estimates for the individual studies.
}
}
\value{
The reference standard deviation of the data set. Non-negative real number.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Sorbye, S., Rue, H. (2014). Scaling intrinsic Gaussian Markov random field priors in
spatial modelling. \emph{Spatial Statistics} \bold{8}, 39--51.
\url{https://doi.org/10.1016/j.spasta.2013.06.004}
%% ~put references to the literature/web site here ~
}

\examples{
# reference sd for resporiratory tract infections (RTI) data
data(rti)
sigma_ref(df=rti)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")

