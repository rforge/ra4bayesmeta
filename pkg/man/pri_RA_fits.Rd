\name{pri_RA_fits}
\alias{pri_RA_fits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prior reference analysis based on bayesmeta fits
%%  ~~function to do ... ~~
}
\description{
Computes a table of Hellinger distance estimates between the
benchmark heterogeneity priors inducing the specified model fits in \code{fits.bm}
and the actual heterogenity priors inducing the specified actual model fits in \code{fits.actual}.
All fits should be based on the same data set.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pri_RA_fits(fits.actual, fits.bm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{fits.actual}{
    a list of model fits of class bayesmeta, computed with
    the \code{bayesmeta}() function in the package bayesmeta.
}
\item{fits.bm}{
    a list of model fits of class bayesmeta, computed with
    the \code{bayesmeta}() function in the package bayesmeta. 
    To be used as benchmarks.
}
}
\details{
Suggestions for prior benchmarks are provided
in Ott et al. (2020, Sections 2.5 and 2.7) and they 
can be computed using the function \code{fit_models_RA}.
}
\value{
A matrix of Hellinger distance estimates between actual and benchmark
heterogeneity priors
with \eqn{n.bm} columns and \eqn{n.act} rows,
where \eqn{n.bm}=length(\code{fit.bms}) is the number of benchmark fits specified and
\eqn{n.act}=length(\code{fit.actual}) the number of actual fits specified.
 
The columns of the matrix give the following Hellinger distance estimates
between two heterogeneity priors
(from left to right):
\item{H(pri_bm_1, pri_act) }{ first benchmark prior bm_1 inducing the fit \code{fits.bm[[1]]} and actual prior}
\item{H(pri_bm_2, pri_act) }{ second benchmark prior bm_2 inducing the fit \code{fits.bm[[2]]} and actual prior}
\item{... }{ ... }
\item{H(pri_bm_n.bm, pri_act) }{ last benchmark prior bm_n.bm inducing the fit \code{fits.bm[[n.bm]]} and actual prior}
Each row correspond to one actual heterogeneity prior
inducing one of the fits in \code{fits.actual}, in the same order
as in that list. Thus, the row names are:
\item{pri_act_1 }{ first actual prior in \code{tau.prior}}
\item{pri_act_2 }{ second actual prior in \code{tau.prior}}
\item{...}{ ...}
\item{pri_act_n }{ nth (last) actual prior in \code{tau.prior}}
}
\references{
Ott, M., Plummer, M., Roos, M. How vague is vague? How informative is
informative? Reference analysis for
Bayesian meta-analysis. Manuscript submitted to \emph{Statistics in Medicine}. 2020.
}
\seealso{
\code{\link[bayesmeta]{bayesmeta}} in the package bayesmeta,
\code{\link{fit_models_RA}},
\code{\link{post_RA_fits}},
\code{\link{pri_RA}} 
}
\examples{
# for aurigular acupuncture (AA) data set 
data(aa)
# compute the model fits
# two actual half-normal heterogeneity priors
fits <- fit_models_RA(df=aa, tau.prior=list(function(t)dhalfnormal(t, scale=0.5),
                                            function(t)dhalfnormal(t, scale=1)))[[1]]
# benchmark fits under SGC(m_inf), SIGC(M_J), SGC(m_J) and SIGC(M_inf) priors
fits.bm.pri <- fits[1:4]
# actual fits under HN priors
fits.actual <- fits[c(6,7)]
pri_RA_fits(fits.actual=fits.actual, fits.bm=fits.bm.pri)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.