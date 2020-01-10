\name{aa}
\alias{aa}
\docType{data}
\title{
Auricular acupuncture data
%%   ~~ data name/kind ... ~~
}
\description{
Meta-analysis data set including 4 randomized, controlled trials comparing treatment completion among
cocaine dependents treated with auricular acupuncture versus sham acupuncture.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(aa)}
\format{
  A data frame with 4 observations on the following 7 variables:
  \describe{
    \item{\code{study}}{character string, label of the study}
    %TODO: specify the direction of the OR
    \item{\code{f.t}}{integer, number of persons who did not complete treatment among those treated with auricular acupuncture }
    \item{\code{n.t}}{integer, number of persons treated with auricular acupuncture}
    \item{\code{f.c}}{integer, number of persons who did not complete treatment among those treated with sham acupuncture}
    \item{\code{n.c}}{integer, number of persons treated with sham acupuncture (control group)}
    \item{\code{y}}{numeric, log odds ratio estimates for the individual trials}
    \item{\code{sigma}}{numeric, standard error of the log odds ratio estimate for the individual trials}
  }
}
\details{
This data set was originally analyzed by Gates et al. (2006) and reanalyzed by Bodnar et al. (2017).

The log odds ratio estimates and standard errors were computed using the \code{escalc} function in the
package \code{metafor}.
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Bodnar, O., Link, A., Arendacka, B., Possolo, A., Elster, C. (2017). Bayesian estimation in random effects meta-analysis using a non-informative prior. 
\emph{Statistics in Medicine} \bold{36}, 378--399.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Gates S, Smith LA, Foxcroft D. Auricular acupuncture for cocaine dependence. \emph{Cochrane Database of Systematic Reviews 2006},
Issue 1. Art. No.: CD005192. \url{https://doi.org/10.1002/14651858.CD005192.pub2}
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(aa)
str(aa)
# forest plot
forest(x=aa$y, sei=aa$sigma, xlab="log odds ratio")
}
\keyword{datasets}
