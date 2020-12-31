\name{ra4bayesmeta-package}
\alias{ra4bayesmeta-package}
\docType{package}
\title{
\packageTitle{ra4bayesmeta}
}
\description{
%TODO: add notation for different priors/posteriors?
Functionality for performing a principled reference analysis in the Bayesian normal-normal hierarchical model (NNHM) used for Bayesian meta-analysis. 
Focuses on the choice of the prior for the between-study (heterogeneity) standard deviation \eqn{\tau}.
Implements the methodology proposed in Ott et al. (2020).

% TODO: decribe RA at posterior and prior level
The posterior reference analysis applies the Hellinger distance 
to quantify the informativeness
of an actual heterogeneity prior of interest by
comparing the corresponding marginal posteriors
(for different parameters in the NNHM)
with posterior benchmarks.
The prior reference analysis 
quantifies the informativeness
of an actual heterogeneity prior given the data by
computing its Hellinger distance to prior benchmarks.

The main posterior benchmark used is Jeffreys reference posterior,
induced by a minimally informative improper reference prior for the heterogeneity.
Two alternative proposals for additional posterior benchmarks are implemented.
The first, simpler proposal is described in Ott et al. (2020) and uses a benchmark induced by a conservative half-normal heterogeneity prior.
This benchmark allows us to discriminate between conservative and anti-conservative heterogeneity priors with respect to the reference prior.
The second, more involved proposal is introduced in the Supplementary Material of Ott et al. (2020) and uses up to four additional benchmarks induced by proper priors.
These benchmarks allow us to discriminate between strongly conservative, reference affine, and strongly anti-conservative heterogeneity priors.
All four benchmark priors belong to two closely related parametric families, the 
so-called "square-root generalized conventional" (SGC) and
"inverse square-root generalized conventional" (SIGC) distributions,
which are derived from prior distributions for variance components
described in Berger & Deely (1988).
Note that most of these prior benchmarks also depend on the data set considered.

Functions which implement this second proposal use 3 benchmarks for the posterior reference analysis (Jeffreys reference posterior and two additional benchmarks induced by a highly conservative/anti-conservative heterogeneity prior) and 5 benchmarks for the prior reference analysis (Jeffreys reference prior and 4 proper priors) and have a corresponding suffix \code{_3bm} or \code{_5bm}.
The corresponding functions without such a suffix implement the first proposal based on 2 heterogeneity benchmark priors (the half-normal and Jeffreys reference prior).


%Two additional extreme (in terms of the relative latent model complexity (RLMC)) posterior benchmarks are also considered, which allow to assess if a
%heterogeneity prior is more informative for the \eqn{\tau=0} or the \eqn{\tau=\infty} hypothesis.
%On the prior level, the impact of the improper Jeffreys reference prior
%on the marginal posteriors
%is bounded 
%by two proper prior benchmarks.
%In addition, the two extreme (in terms of RLMC) prior benchmarks are also included.
%All four benchmarks belong to two closely related parametric families, the 
%so-called "square-root generalized conventional" (SGC) and
%"inverse square-root generalized conventional" (SIGC) distributions,
%which are derived from prior distributions for variance components
%described in Berger & Deely (1988).
%Note that most of these prior benchmarks also depend on the data set considered.




%Computes a reference posterior, induced by a minimally informative improper reference prior for the between-study (heterogeneity) standard deviation. Determines two additional proper prior benchmarks bounding the impact of the improper reference prior. Includes functions for reference analyses at both the posterior and the prior level, which, given the data, quantify the informativeness of a heterogeneity prior of interest relative to the minimally informative reference prior and the two proper prior benchmarks.
The main functions \code{post_RA} and \code{post_RA_3bm}
perform the reference analysis at the posterior level based on a data frame as input.
The functions \code{fit_models_RA} and \code{fit_models_RA_5bm} compute the benchmark heterogeneity priors and 
corresponding (marginal) posteriors for a given data frame.
The more flexible functions \code{post_RA_fits} 
and \code{pri_RA_fits} 
perform
the reference analysis at the posterior level
and prior level
based on the provided benchmark and actual model fits (of class \code{bayesmeta}),
so that alternative benchmarks heterogeneity priors can also be used.
Functions for plotting the densities of the heterogeneity priors and
marginal posteriors for different parameters are also provided.
The functions operate on data frames which are compatible with the bayesmeta R package on CRAN.
}
\details{
Package: ra4bayesmeta

Type: Package

Title: Reference analysis for Bayesian meta-analysis

Version: 1.0-1

Date: 2020-12-31

Author: Manuela Ott [aut, cre], Malgorzata Roos [aut] 

Maintainer: Manuela Ott <manuela.c.ott@gmail.com>

Depends: bayesmeta

License: GPL (>=2)
%The DESCRIPTION file:
%\packageDESCRIPTION{pCalibrate}
%\packageIndices{pCalibrate}
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
}
\author{
Manuela Ott, Malgorzata Roos
%\packageAuthor{pCalibrate}
Maintainer: \packageMaintainer{ra4bayesmeta}
}
\references{
Ott, M., Plummer, M., Roos, M. How vague is vague? How informative is
informative? Reference analysis for
Bayesian meta-analysis. Manuscript revised for \emph{Statistics in Medicine}. 2020.

Ott, M., Plummer, M., Roos, M. Supplementary Material:
How vague is vague? How informative is informative? Reference analysis for 
Bayesian meta-analysis. Revised for \emph{Statistics in Medicine}. 2020.

Berger, J. O., Deely, J. (1988). A Bayesian approach to ranking and selection of
related means with alternatives to analysis-of-variance methodology. \emph{Journal of the
American Statistical Association} \bold{83(402)}, 364--373.
%submitted to \emph{Statisitcs in Medicine}
}
\keyword{ package
}
\keyword{htest}
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
# auricular acupuncture (AA) data set
data(aa)

# posterior reference analysis
# for standard HN and HC heterogeneity priors
# using the HN0 and J benchmarks
post_RA(df=aa, tau.prior=
               list(function(t)dhalfnormal(t, scale=1),
                    function(t)dhalfcauchy(t, scale=1)))

# posterior reference analysis
# for standard HN and HC heterogeneity priors
# using 3 benchmarks (SGC(m_inf), J and SIGC(M_inf) priors)
post_RA_3bm(df=aa, tau.prior=
                   list(function(t)dhalfnormal(t, scale=1),
                        function(t)dhalfcauchy(t, scale=1)))
          
# plot the marginal posteriors for the effect mu
# using the HN0 and J benchmarks
plot_RA(df=aa, tau.prior=
               list(function(t)dhalfnormal(t, scale=1),
                    function(t)dhalfcauchy(t, scale=1)),
        type="post.mu", xlim=c(-1,1.5), ylim=c(0,2.5),
        legend=TRUE, legend.tau.prior=c("HN(1)", "HC(1)"))   
}
